""" 
1) Generate equations of motion for an onshore wind turbine with the folowing degrees of freedom:
  - 1 fore-aft shape function
  - 1 side-side shape function
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
from welib.yams.models.FTNSB_sympy import get_model
from welib.yams.models.FTNSB_sympy_symbols import *

def main(runSim=True, runFAST=False, create=True):

    if create:
        model = get_model('F0T2N0S1', mergeFndTwr=True, linRot=False,
                          yaw='zero', tilt='fixed',tiltShaft=True,
                          rot_elastic_type='SmallRot', #rot_elastic_type='Body', 'Body' or 'SmallRot'
                          orderMM=1,
                          orderH=1,
                          twrDOFDir=['x','y','x','y'], # Order in which the flexible DOF of the tower are set
                         )
        extraSubs=model.shapeNormSubs # shape functions normalized to unity
        smallAngles  = [(model.twr.vcList, 2)]
        smallAngles += [([theta_tilt]    , 1)]
        replaceDict={'theta_tilt':('tilt',None)}
        model.exportPackage(path='_F0T2N0S1', extraSubs=extraSubs, smallAngles=smallAngles, replaceDict=replaceDict, pathtex='_F0T2N0S1')



    # --- Run non linear and linear simulation using a FAST model as input
    if runSim:
        # --- Import the python module that was generated
        model_pkg = importlib.import_module('_F0T2N0S1')

        # --- Load the wind turbine model, and extract relevant parameters "p"
        MyDir=os.path.dirname(__file__)
        #fstFilename = os.path.join(MyDir, '../../../data/NREL5MW/Main_Onshore.fst')
        fstFilename = os.path.join(MyDir, 'F0T2N0S1/Main_Spar_ED.fst')
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
    dfNL, dfLI, dfFS = main(runSim=True, create=True)
    plt.show()

if __name__=="__test__":
    dfNL, dfLI, dfFS = main(runSim=True)

    from welib.tools.stats import mean_rel_err
    vb = False
    method='minmax'
    eps1= mean_rel_err(y1=dfNL['Azimuth_[deg]'],   y2=dfFS['Azimuth_[deg]']  , method=method, verbose=vb)
    eps2= mean_rel_err(y1=dfNL['Q_TFA1_[m]'],      y2=dfFS['Q_TFA1_[m]']     , method=method, verbose=vb)
    eps3= mean_rel_err(y1=dfNL['Q_TSS1_[m]'],      y2=dfFS['Q_TSS1_[m]']     , method=method, verbose=vb)
    np.testing.assert_array_less(eps1, 0.57)
    np.testing.assert_array_less(eps2, 0.58)
    np.testing.assert_array_less(eps3, 0.59)


    eps1= mean_rel_err(y1=dfLI['Azimuth_[deg]'],   y2=dfFS['Azimuth_[deg]']  , method=method, verbose=vb)
    eps2= mean_rel_err(y1=dfLI['Q_TFA1_[m]'],      y2=dfFS['Q_TFA1_[m]']     , method=method, verbose=vb)
    eps3= mean_rel_err(y1=dfLI['Q_TSS1_[m]'],      y2=dfFS['Q_TSS1_[m]']     , method=method, verbose=vb)
    np.testing.assert_array_less(eps1, 0.57)
    np.testing.assert_array_less(eps2, 0.58)
    np.testing.assert_array_less(eps3, 0.59)
    try:
        os.remove('./_F0T2N0S1.py')
        os.remove('./_F0T2N0S1.tex')
    except:
        pass
