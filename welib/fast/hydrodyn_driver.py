import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
# Local 
from welib.tools.tictoc import Timer

import welib.weio as weio
from welib.weio.fast_input_file import FASTInputFile
from welib.weio.fast_output_file import writeDataFrame
from welib.fast.hydrodyn import HydroDyn
from welib.fast.fast_mesh import MeshStorage
from welib.tools.clean_exceptions import *
from welib.tools.tictoc import Timer


def hydroSimLinFromOpenFAST(fstFilename, tMax=None, optsM=None, plot=True, json=True, out=True, png=True, verbose=False, base=None, MCKF=None, q0=None, fig=None):
    """ 
    """
    if base is None:
        base=fstFilename.replace('.fst','')
    if q0 is None:
        q0=np.zeros(6)

    # --- Initialize a python HydroDyn instance
    hd = HydroDyn(fstFilename)

    if MCKF is None:
        MCKF = hd.linearize_RigidMotion2Loads(q0)
    M,C,K,F0=MCKF

    # --- Open OpenFAST output file
    dfOF = weio.read(hd.outFilename).toDataFrame()
    if tMax is not None:
        dfOF=dfOF[dfOF['Time_[s]']<=tMax]
    if F0 is None:
        F0=np.zeros(6) # TODO mean value?

    # --- Relevant columns from OpenFAST outputs
    time  = dfOF['Time_[s]'].values
    if 'PRPSurge_[m]' in dfOF.columns:
        qCol   = ['PRPSurge_[m]'    ,'PRPSway_[m]'    ,'PRPHeave_[m]'   ,'PRPRoll_[rad]'    ,'PRPPitch_[rad]'   ,'PRPYaw_[rad]'     ]
        qdCol  = ['PRPTVxi_[m/s]'   ,'PRPTVyi_[m/s]'  ,'PRPTVzi_[m/s]'  ,'PRPRVxi_[rad/s]'  ,'PRPRVyi_[rad/s]'  ,'PRPRVzi_[rad/s]'  ]
        qddCol = ['PRPTAxi_[m/s^2]','PRPTAyi_[m/s^2]','PRPTAzi_[m/s^2]','PRPRAxi_[rad/s^2]','PRPRAyi_[rad/s^2]','PRPRAzi_[rad/s^2]']
    else:
        raise NotImplementedError()

    # --- Prepare time stepping
    fh = np.zeros((len(time), 6))

    # --- Time integration
    with Timer('Time integration'):
        for it, t in enumerate(time):
            # Reference point motion
            q   = dfOF[qCol].iloc[it].values
            qd  = dfOF[qdCol].iloc[it].values
            qdd = dfOF[qddCol].iloc[it].values

            fh[it,:] = -M.dot(qdd) - C.dot(qd) - K.dot(q)
            if verbose:
                print('f_Hydro {:12.4f} {:16.3f}{:16.3f}{:16.3f}{:16.3f}{:16.3f}{:16.3f}'.format(t, *fh[it,:]))

    # --- Creating a dataframe for convenience
    lCols = ['HydroFxi_[N]','HydroFyi_[N]','HydroFzi_[N]','HydroMxi_[N-m]','HydroMyi_[N-m]','HydroMzi_[N-m]']
    dfPH = pd.DataFrame(data=np.column_stack((time,fh)), columns=['Time_[s]']+lCols)

    # --- Adding operating point values if present
    for iCol, col in enumerate(lCols):
        dfPH[col]+=F0[iCol]


    # --- Plot
    if plot:
        if fig is None:
            fig,axes = plt.subplots(6, 2, sharey=False, figsize=(12.8,8.5)) # (6.4,4.8)
            fig.subplots_adjust(left=0.11, right=0.95, top=0.95, bottom=0.07, hspace=0.40, wspace=0.22)
            newFig = True
        else:
            axes= np.array(fig.axes).reshape(6,2)
            newFig = False

        # DOF
        if newFig:
            for iCol, col in enumerate(qCol):
                axes[iCol,0].plot(time, dfOF[col].values  , '-', label='OpenFAST')
                axes[iCol,0].set_ylabel(col.replace('_',' '))
        # Forces
        for iCol, col in enumerate(lCols):
            axes[iCol,1].plot(time, dfPH[col].values, '--', label='Python linear')
            if newFig:
                axes[iCol,1].plot(time, dfOF[col].values, 'k:', label='OpenFAST')
                mi=np.min(dfOF[col].values)
                mx=np.max(dfOF[col].values)
                mn=np.mean(dfOF[col].values)

        # Scale axes if tiny range
        for iCol, col in enumerate(lCols):
            mi, mx = axes[iCol,1].get_ylim()
            mn = (mx+mi)/2
            if np.abs(mx-mn)<1e-5:
                axes[iCol,1].set_ylim(mn-1, mn+1)
        axes[0,1].legend()
        axes[5,0].set_xlabel('Time [s]')
        axes[5,1].set_xlabel('Time [s]')
        if png:
            fig.savefig(base+'_hydroPyPrescrMotion_Lin.png')

    if out:
        writeDataFrame(dfPH, base + '_hydroPyPrescrMotion_Lin.outb')

    return dfPH, dfOF, fig


def hydroSimFromOpenFAST(fstFilename, tMax=None, optsM=None, plot=True, json=True, out=True, png=True, verbose=False, base=None, fig=None, motionRef='PRP', zRef=None):
    """ 
    Compute hydrodynamics loads on a structure under rigid body motion based on an OpenFAST simulation
    fstFilename: OpenFAST or HydroDyn driver input file
    tMax: maximum time in time vector (to reduce the simulation length)
    """
    if base is None:
        base=fstFilename.replace('.fst','')

    # --- Initialize a python HydroDyn instance
    hd = HydroDyn(fstFilename)
    #hd.writeSummary(hdFilename.replace('.dat','.HD_python.sum'))
    umesh = hd.u['Morison']['Mesh']
    ymesh = hd.y['Morison']['Mesh']
    #print(hd.p)

    # --- Open OpenFAST output file
    dfOF = weio.read(hd.outFilename).toDataFrame()
    if tMax is not None:
        dfOF=dfOF[dfOF['Time_[s]']<=tMax]


    # --- Relevant columns from OpenFAST outputs
    time  = dfOF['Time_[s]'].values
    if motionRef=='PRP':
        RefPoint=(0,0,0)
        if 'PRPSurge_[m]' in dfOF.columns:
            qCol   = ['PRPSurge_[m]'    ,'PRPSway_[m]'    ,'PRPHeave_[m]'   ,'PRPRoll_[rad]'    ,'PRPPitch_[rad]'   ,'PRPYaw_[rad]'     ]
            qdCol  = ['PRPTVxi_[m/s]'   ,'PRPTVyi_[m/s]'  ,'PRPTVzi_[m/s]'  ,'PRPRVxi_[rad/s]'  ,'PRPRVyi_[rad/s]'  ,'PRPRVzi_[rad/s]'  ]
            qddCol = [ 'PRPTAxi_[m/s^2]','PRPTAyi_[m/s^2]','PRPTAzi_[m/s^2]','PRPRAxi_[rad/s^2]','PRPRAyi_[rad/s^2]','PRPRAzi_[rad/s^2]']
        else:
            raise NotImplementedError()
    else:
        qCol   = ['Q_Sg_[m]'      ,'Q_Sw_[m]'      ,'Q_Hv_[m]'      ,'Q_R_[rad]'      ,'Q_P_[rad]'      ,'Q_Y_[rad]']
        qdCol  = ['QD_Sg_[m/s]'   ,'QD_Sw_[m/s]'   ,'QD_Hv_[m/s]'   ,'QD_R_[rad/s]'   ,'QD_P_[rad/s]'   ,'QD_Y_[rad/s]']
        qddCol = ['QD2_Sg_[m/s^2]','QD2_Sw_[m/s^2]','QD2_Hv_[m/s^2]','QD2_R_[rad/s^2]','QD2_P_[rad/s^2]','QD2_Y_[rad/s^2]']
        if zRef is None:
            raise Exception('When motionRef is not PRP, zRef needs to be provided')
        RefPoint=(0,0,zRef)

    # --- Prepare time stepping
    fh = np.zeros((len(time), 6))
    msy = MeshStorage(ymesh, time) # Store mesh at each time step

    # --- Time integration
    with Timer('Time integration'):
        for it, t in enumerate(time):
            # Reference point motion
            q  = dfOF[qCol].iloc[it].values
            qd = dfOF[qdCol].iloc[it].values
            qdd = dfOF[qddCol].iloc[it].values
            omega = (qd[3],qd[4],qd[5]) # ...

            # Rigid body motion of the mesh
            #umesh.rigidBodyMotion(u=(q[0],q[1],q[2]), theta=(q[3],q[4],q[5]), u_dot=(qd[0],qd[1],qd[2]), omega=omega, u_ddot=(qdd[0],qdd[1],qdd[2]), omega_dot = (qdd[3],qdd[4],qdd[5])  )
            umesh.rigidBodyMotion(q=q, qd=qd, qdd=qdd, RefPoint=RefPoint)
            # Calculate hydrodynamic loads at every nodes 
            hd.y=hd.calcOutput(t, u=hd.u, y=hd.y, optsM=optsM)
            # Store mesh
            msy.store(ymesh, it)
            # Compute integral loads (force&moment) at the reference point (translated but not rotated)
            fh[it, :3], fh[it, 3:] = ymesh.mapLoadsToPoint((q[0],q[1],q[2]))
            if verbose:
                print('f_Hydro {:12.4f} {:16.3f}{:16.3f}{:16.3f}{:16.3f}{:16.3f}{:16.3f}'.format(t, *fh[it,:]))

    # --- Creating a dataframe for convenience
    lCols = ['HydroFxi_[N]','HydroFyi_[N]','HydroFzi_[N]','HydroMxi_[N-m]','HydroMyi_[N-m]','HydroMzi_[N-m]']
    dfPH = pd.DataFrame(data=np.column_stack((time,fh)), columns=['Time_[s]']+lCols)

    # --- Plot
    if plot:
        if fig is None:
            fig,axes = plt.subplots(6, 2, sharey=False, figsize=(12.8,8.5)) # (6.4,4.8)
            fig.subplots_adjust(left=0.11, right=0.95, top=0.95, bottom=0.07, hspace=0.40, wspace=0.22)
            newFig=True
        else:
            newFig=False
            axes= np.array(fig.axes).reshape(6,2)
        # DOF
        if newFig:
            for iCol, col in enumerate(qCol):
                axes[iCol,0].plot(time, dfOF[col].values  , '-', label='OpenFAST')
                axes[iCol,0].set_ylabel(col.replace('_',' '))
        # Forces
        for iCol, col in enumerate(lCols):
            axes[iCol,1].plot(time, dfPH[col].values  , label='Python non-linear')
            if newFig:
                axes[iCol,1].plot(time, dfOF[col].values  , 'k:', label='OpenFAST')
                axes[iCol,1].set_ylabel(col.replace('_',' '))
        axes[0,1].legend()
        axes[5,0].set_xlabel('Time [s]')
        axes[5,1].set_xlabel('Time [s]')
        if png:
            fig.savefig(base+'_hydroPyPrescrMotion_NL.png')

    if json:
        msy.toJSON3D(ymesh, base + '_hydroPyPrescrMotion_NL_MeshMotion.json')

    if out:
        writeDataFrame(dfPH, base + '_hydroPyPrescrMotion_NL.outb')

    return dfPH, dfOF, msy, fig


if __name__ == '__main__':
    fstFilename = os.path.join(MyDir, '../../data/Spar/Main_Spar_ED_HydroExample.fst');
    dfPH, dfOF, msy = hydroSimFromOpenFAST(fstFilename, tMax=None)
    plt.show()
