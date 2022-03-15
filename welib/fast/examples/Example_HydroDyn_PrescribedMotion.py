""" 
Compute hydrodynamics loads on a structure under rigid body motion using "Python HydroDyn"

NOTE:
   the code below can be replaced by a simple call to hydrodyn_driver.py:

       dfPH, dfOF, msy = hydroSimFromOpenFAST(fstFilename, tMax=None)


- Hydrodynamic parameters are read from an HydroDyn input file
- Motion (displacement, velocities, accelerations) are taken from an OpenFAST simulation
  NOTES:
      - requires latest dev branch of OpenFAST to get platform reference outputs for now
      - if ascii files are used, the output resolution "OutFmt" should have sufficient digits)
- The Motion is applied at the HydroDyn reference point 
- Loads are computed at the HydroDyn reference point using the python hydrodyn module.

"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
# Local 
import welib.weio as weio
from welib.weio.fast_input_deck import FASTInputDeck
from welib.weio.fast_input_file import FASTInputFile
from welib.fast.hydrodyn import HydroDyn
from welib.fast.fast_mesh import *
from welib.tools.clean_exceptions import *
from welib.tools.tictoc import Timer

MyDir=os.path.dirname(__file__)

def hydroSim(fstFilename, plot=True, json=True, tMax=None):
    """ 
    Compute hydrodynamics loads on a structure under rigid body motion based on an OpenFAST simulation
    fstFilename: OpenFAST or HydroDyn driver input file
    tMax: maximum time in time vector (to reduce the simulation length)
    """
    # --- Find out file and hydrodyn filenames from FST file
    fst = FASTInputFile(fstFilename)
    if 'HydroFile' in fst.keys():
        hdFilename = os.path.join(os.path.dirname(fstFilename), fst['HydroFile'].replace('"','') )
        outFilenames = [fstFilename.replace('.fst',ext) for ext in ['.outb','.out'] if os.path.exists(fstFilename.replace('.fst',ext))]
        outFilename = outFilenames[0]
    else:
        hdFilename = os.path.join(os.path.dirname(fstFilename), fst['HDInputFile'].replace('"','') )
        outFilename = fstFilename.replace('.dvr','.HD.out')
    base = os.path.basename(os.path.dirname(fstFilename)) + '_'+os.path.splitext(os.path.basename(fstFilename))[0]

    # --- Open OpenFAST output file
    dfOF = weio.read(outFilename).toDataFrame()
    if tMax is not None:
        dfOF=dfOF[dfOF['Time_[s]']<=tMax]

    # --- Initialize a python HydroDyn instance
    hd = HydroDyn(hdFilename)
    u, y = hd.init(Gravity = fst['Gravity'], WtrDens=fst['WtrDens'], WtrDpth=fst['WtrDpth'])
    hd.writeSummary(hdFilename.replace('.dat','.HD_python.sum'))
    umesh = u['Morison']['Mesh']
    ymesh = y['Morison']['Mesh']
    #print(hd.p)

    # --- Relevant columns from OpenFAST outputs
    time  = dfOF['Time_[s]'].values
    if 'PRPSurge_[m]' in dfOF.columns:
        qCol   = ['PRPSurge_[m]'    ,'PRPSway_[m]'    ,'PRPHeave_[m]'   ,'PRPRoll_[rad]'    ,'PRPPitch_[rad]'   ,'PRPYaw_[rad]']
        qdCol  = ['PRPTVxi_[m/s]'   ,'PRPTVyi_[m/s]'  ,'PRPTVzi_[m/s]'  ,'PRPRVxi_[rad/s]'  ,'PRPRVyi_[rad/s]'  ,'PRPRVzi_[rad/s]']
        qddCol = [ 'PRPTAxi_[m/s^2]','PRPTAyi_[m/s^2]','PRPTAzi_[m/s^2]','PRPRAxi_[rad/s^2]','PRPRAyi_[rad/s^2]','PRPRAzi_[rad/s^2]']
    else:
        raise NotImplementedError()
        qCol  = ['Q_Sg_[m]'   ,'Q_Sw_[m]'   ,'Q_Hv_[m]'   ,'Q_R_[rad]'   ,'Q_P_[rad]'   ,'Q_Y_[rad]']
        qdCol = ['QD_Sg_[m/s]','QD_Sw_[m/s]','QD_Hv_[m/s]','QD_R_[rad/s]','QD_P_[rad/s]','QD_Y_[rad/s]']
        qddCol = None

    # --- Prepare time stepping
    fh = np.zeros((len(time), 6))
    msy = MeshStorage(ymesh, time) # Store mesh at each time step

    # --- Time integration
    with Timer('Time integration'):
        for it, t in enumerate(time):
            # Reference point motion
            q  = dfOF[qCol].iloc[it].values
            qd = dfOF[qdCol].iloc[it].values
            omega = (qd[3],qd[4],qd[5]) # ...
            if qddCol is None:
                if it==0:
                    qdd = qd*0
                else:
                    qdd =(qd -  dfOF[qdCol].iloc[it-1].values) / (time[it]-time[it-1])
            else:
                qdd = dfOF[qddCol].iloc[it].values

            # Rigid body motion of the mesh
            umesh.rigidBodyMotion(u=(q[0],q[1],q[2]), theta=(q[3],q[4],q[5]), u_dot=(qd[0],qd[1],qd[2]), omega=omega, u_ddot=(qdd[0],qdd[1],qdd[2]), omega_dot = (qdd[3],qdd[4],qdd[5])  )
            # Calculate hydrodynamic loads at every nodes 
            y=hd.calcOutput(t, u=u, y=y)
            # Store mesh
            msy.store(ymesh, it)
            # Compute integral loads (force&moment) at the reference point (translated but not rotated)
            fh[it, :3], fh[it, 3:] = ymesh.mapLoadsToPoint((q[0],q[1],q[2]))
            print('f_Hydro {:16.3f}{:16.3f}{:16.3f}{:16.3f}{:16.3f}{:16.3f}'.format(*fh[it,:]))

    # --- Creating a dataframe for convenience
    lCols = ['HydroFxi_[N]','HydroFyi_[N]','HydroFzi_[N]','HydroMxi_[N-m]','HydroMyi_[N-m]','HydroMzi_[N-m]']
    dfPH = pd.DataFrame(data=np.column_stack((time,fh)), columns=['Time_[s]']+lCols)

    # --- Plot
    if plot:
        fig,axes = plt.subplots(6, 2, sharey=False, figsize=(12.8,8.5)) # (6.4,4.8)
        fig.subplots_adjust(left=0.11, right=0.95, top=0.95, bottom=0.07, hspace=0.40, wspace=0.22)
        # DOF
        for iCol, col in enumerate(qCol):
            axes[iCol,0].plot(time, dfOF[col].values  , '-', label='OpenFAST')
            axes[iCol,0].set_ylabel(col.replace('_',' '))
        # Forces
        for iCol, col in enumerate(lCols):
            axes[iCol,1].plot(time, dfPH[col].values  , label='Python non-linear')
            axes[iCol,1].plot(time, dfOF[col].values  , 'k:', label='OpenFAST')
            axes[iCol,1].set_ylabel(col.replace('_',' '))
        axes[0,1].legend()
        axes[5,0].set_xlabel('Time [s]')
        axes[5,1].set_xlabel('Time [s]')

    if json:
        msy.toJSON3D(ymesh, '_'+base+'_MeshMotion.json')

    return dfPH, dfOF

if __name__ == '__main__':
    fstFilename = os.path.join(MyDir, '../../../data/Spar/Main_Spar_ED_HydroExample.fst');
    dfPH, dfOF = hydroSim(fstFilename, tMax=None)
    plt.show()

if __name__ == '__test__':
    fstFilename = os.path.join(MyDir, '../../../data/Spar/Main_Spar_ED_HydroExample.fst');
    dfPH, dfOF = hydroSim(fstFilename, tMax=0.1, plot=False, json=False)
    np.testing.assert_almost_equal(dfPH['HydroFxi_[N]'].values  /1e6,dfOF['HydroFxi_[N]'].values  /1e6, 5)
    np.testing.assert_almost_equal(dfPH['HydroFyi_[N]'].values  /1e6,dfOF['HydroFyi_[N]'].values  /1e6, 5)
    np.testing.assert_almost_equal(dfPH['HydroFzi_[N]'].values  /1e6,dfOF['HydroFzi_[N]'].values  /1e6, 5)
    np.testing.assert_almost_equal(dfPH['HydroMxi_[N-m]'].values/1e6,dfOF['HydroMxi_[N-m]'].values/1e6, 3)
    np.testing.assert_almost_equal(dfPH['HydroMyi_[N-m]'].values/1e6,dfOF['HydroMyi_[N-m]'].values/1e6, 3)
    np.testing.assert_almost_equal(dfPH['HydroMzi_[N-m]'].values/1e6,dfOF['HydroMzi_[N-m]'].values/1e6, 3)
