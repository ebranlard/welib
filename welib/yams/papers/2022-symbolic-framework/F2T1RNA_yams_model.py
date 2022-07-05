""" 
1) Generate equations of motion for an onshore wind turbine with the folowing degrees of freedom:
  - Platform Surge and pitch 
  - 1 fore-aft shape function
  - shaft rotation

2) Linearize the equations of motion

3) Simulate the non linear and linear equations of motion based on a OpenFAST model, and compare results

"""
import unittest
import os
import numpy as np    
from welib.yams.rotations import *
from sympy import trigsimp
import matplotlib.pyplot as plt
import importlib
# yams
from welib.yams.models.FTNSB_sympy import *
from welib.yams.models.FTNSB_sympy_symbols import *

def main(runSim=True, runFAST=False):

    model = get_model('F2T1RNA', mergeFndTwr=False, linRot=False,
                      yaw='zero', tilt='fixed',tiltShaft=True,
                      rot_elastic_type='SmallRot', #rot_elastic_type='Body', 'Body' or 'SmallRot'
                      orderMM=1,
                      orderH=1,
                      twrDOFDir=['x','y','x','y'], # Order in which the flexible DOF of the tower are set
                     )
    extraSubs=model.shapeNormSubs # shape functions normalized to unity
    smallAngles  = [(model.twr.vcList, 2)]
    smallAngles += [([theta_tilt, phi_y]    , 1)]
    replaceDict={'theta_tilt':('tilt',None)}
    model.exportPackage(path='_F2T1RNA', extraSubs=extraSubs, smallAngles=smallAngles, replaceDict=replaceDict, pathtex='_F2T1RNA')

    # --- Run non linear and linear simulation using a FAST model as input
    if runSim:
        # TODO TODO all this can be replaced with a call to yams.models.simulator
        #sim = SimulatorFromOF(WT, modelName=modelName, packageDir='py')
        #time, dfFS, p = sim.setupSim(tMax=tMax)

        # --- Import the python module that was generated
        model_pkg = importlib.import_module('_F2T1RNA')

        # --- Load the wind turbine model, and extract relevant parameters "p"
        MyDir=os.path.dirname(__file__)
        #fstFilename = os.path.join(MyDir, '../../../data/NREL5MW/Main_Onshore.fst')
        fstFilename = os.path.join(MyDir, 'F2T1RNA_SmallAngle/Main_Spar_ED.fst')
        from welib.yams.windturbine import FASTWindTurbine
        WT = FASTWindTurbine(fstFilename, twrShapes=[0,2], nSpanTwr=50)
        p = WT.yams_parameters()

        # --- Perform time integration
        if os.path.exists(fstFilename.replace('.fst','.outb')):
            import welib.weio as weio
            dfFS = weio.read(fstFilename.replace('.fst','.outb')).toDataFrame()
            time =dfFS['Time_[s]'].values
        else:
            time = np.linspace(0,50,1000)
            dfFS = None
        resNL, sysNL, dfNL = WT.simulate_py    (model_pkg, p, time)
        resLI, sysLI, dfLI = WT.simulate_py_lin(model_pkg, p, time)

        # --- Simple Plot
        fig,axes = plt.subplots(3, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
        fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
        for i,ax in enumerate(axes):
            chan=dfNL.columns[i+1]
            ax.plot(dfNL['Time_[s]'], dfNL[chan], '-'  , label='non-linear')
            ax.plot(dfLI['Time_[s]'], dfLI[chan], '--' , label='linear')
            if dfFS is not None:
                ax.plot(dfFS['Time_[s]'], dfFS[chan], 'k:' , label='OpenFAST')
            ax.set_xlabel('Time [s]')
            ax.set_ylabel(chan)
            ax.tick_params(direction='in')
        ax.legend()
    else:
        dfNL=None
        dfLI=None
        dfFS=None
    return dfNL, dfLI, dfFS



if __name__=="__main__":
    dfNL, dfLI, dfFS = main(runSim=True)
    plt.show()
 
if __name__=="__test__":
    dfNL, dfLI, dfFS = main(runSim=True)

    from welib.tools.stats import mean_rel_err
    vb = False
    method='minmax'
    eps1= mean_rel_err(y1=dfNL['PtfmSurge_[m]'],   y2=dfFS['PtfmSurge_[m]']  , method=method, verbose=vb)
    eps2= mean_rel_err(y1=dfNL['PtfmPitch_[deg]'], y2=dfFS['PtfmPitch_[deg]'], method=method, verbose=vb)
    eps3= mean_rel_err(y1=dfNL['Q_TFA1_[m]'],      y2=dfFS['Q_TFA1_[m]']     , method=method, verbose=vb)
    np.testing.assert_array_less(eps1, 0.58)
    np.testing.assert_array_less(eps2, 0.59)
    np.testing.assert_array_less(eps3, 0.60)

    eps1= mean_rel_err(y1=dfLI['PtfmSurge_[m]'],   y2=dfFS['PtfmSurge_[m]']  , method=method, verbose=vb)
    eps2= mean_rel_err(y1=dfLI['PtfmPitch_[deg]'], y2=dfFS['PtfmPitch_[deg]'], method=method, verbose=vb)
    eps3= mean_rel_err(y1=dfLI['Q_TFA1_[m]'],      y2=dfFS['Q_TFA1_[m]']     , method=method, verbose=vb)
    np.testing.assert_array_less(eps1, 1.03)
    np.testing.assert_array_less(eps2, 1.02)
    np.testing.assert_array_less(eps3, 0.73)

    try:
        os.remove('./_F2T1RNA.py')
        os.remove('./_F2T1RNA.tex')
    except:
        pass
