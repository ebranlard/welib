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
from welib.yams.models.OneRigidBody_sympy import get_model_one_body
from welib.yams.models.FTNSB_sympy_symbols import *

def main(runSim=True, runFAST=False):

    model = get_model('F2T1RNA', mergeFndTwr=False, linRot=False,
                      yaw='zero', tilt='fixed',tiltShaft=True,
                      rot_elastic_type='SmallRot', #rot_elastic_type='Body', 'Body' or 'SmallRot'
                      orderMM=1,
                      orderH=1,
                      twrDOFDir=['x','y','x','y'], # Order in which the flexible DOF of the tower are set
                     )
    # --- Compute Kane's equations of motion
    model.kaneEquations(Mform='TaylorExpanded')
    EOM = model.to_EOM()

    # --- Additional substitution 
    extraSubs=model.shapeNormSubs # shape functions normalized to unity
    EOM.subs(extraSubs)

    #  --- Small angle approximations
    EOM.smallAngleApprox(model.twr.vcList, order = 2)
    EOM.smallAngleApprox([theta_tilt, phi_y], order = 1)
    EOM.simplify()

    # --- Separate EOM into mass matrix and forcing
    EOM.mass_forcing_form() # EOM.M and EOM.F

    # --- Linearize equation of motions 
    EOM.linearize(noAcc=True) # EOM.M0, EOM.K0, EOM.C0, EOM.B0

    # --- Export equations
    replaceDict={'theta_tilt':('tilt',None)}
    EOM.savePython(folder='./', prefix='_', replaceDict=replaceDict)
    EOM.saveTex   (folder='./', prefix='_', variables=['M','F','M0','K0','C0','B0'])

    # --- Run non linear and linear simulation using a FAST model as input
    if runSim:
        # --- Import the python module that was generated
        model_pkg = importlib.import_module('_F2T1RNA')

        # --- Load the wind turbine model, and extract relevant parameters "p"
        MyDir=os.path.dirname(__file__)
        #fstFilename = os.path.join(MyDir, '../../../data/NREL5MW/Main_Onshore_OF2.fst')
        fstFilename = os.path.join(MyDir, '../examples/_F2T1RNA_SmallAngle/Main_Spar_ED.fst')
        from welib.yams.windturbine import FASTWindTurbine
        WT = FASTWindTurbine(fstFilename, twrShapes=[0,2], nSpanTwr=50)
        p = WT.yams_parameters()

        # --- Perform time integration
        if os.path.exists(fstFilename.replace('.fst','.out')):
            import welib.weio as weio
            dfFS = weio.read(fstFilename.replace('.fst','.out')).toDataFrame()
            time =dfFS['Time_[s]'].values
        else:
            time = np.linspace(0,50,1000)
            dfFS = None
        resNL, sysNL = WT.simulate_py    (model_pkg, p, time)
        resLI, sysLI = WT.simulate_py_lin(model_pkg, p, time)
        dfNL = sysNL.toDataFrame(WT.channels, WT.FASTDOFScales)
        dfLI = sysLI.toDataFrame(WT.channels, WT.FASTDOFScales)

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

if __name__=="__main__":
    main(runSim=True)
    plt.show()

if __name__=="__test__":
    main(runSim=False)
    try:
        os.remove('./_F2T1RNA.py')
        os.remove('./_F2T1RNA.tex')
    except:
        pass
