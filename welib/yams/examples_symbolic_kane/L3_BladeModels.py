""" 

"""
import unittest
import os
import glob    
import numpy as np    
from welib.yams.rotations import *
from sympy import trigsimp
import matplotlib.pyplot as plt
import importlib
# yams
from welib.yams.models.FTNSB_sympy import get_model
from welib.yams.models.FTNSB_sympy_symbols import *
from welib.tools.clean_exceptions import *

def main(model, simFile=None, extraSubs=None, smallAngles=None, op_point=None,
        generate=True, runSim=True, latex=False):

    if extraSubs is None:
        extraSubs=[]
    if smallAngles is None:
        smallAngles=[]


    if generate:
        # --- Compute Kane's equations of motion
        model.kaneEquations(Mform='TaylorExpanded')
        EOM = model.to_EOM()

        # --- Additional substitution 
        EOM.subs(extraSubs)

        #  --- Small angle approximations
        for angles, order in smallAngles:
            EOM.smallAngleApprox(angle, order = order)
        EOM.simplify()

        # --- Separate EOM into mass matrix and forcing
        EOM.mass_forcing_form() # EOM.M and EOM.F

        # --- Linearize equation of motions 
        EOM.linearize(noAcc=True, op_point=op_point) # EOM.M0, EOM.K0, EOM.C0, EOM.B0

        # --- Export equations
        replaceDict={'theta_tilt':('tilt',None)}
        EOM.savePython(folder='./', prefix='_', replaceDict=replaceDict)
        if latex:
            EOM.saveTex   (folder='./', prefix='_', variables=['M','F','M0','K0','C0','B0'])

    # --- Run non linear and linear simulation using a FAST model as input
    if runSim:
        # --- Import the python module that was generated
        model_pkg = importlib.import_module('_'+model.name)

        # --- Load the wind turbine model, and extract relevant parameters "p"
        MyDir=os.path.dirname(__file__)
        fstFilename = os.path.join(MyDir, simFile)
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
        print(WT)
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
    else:
        dfNL=None
        dfLI=None
        dfFS=None
    return dfNL, dfLI, dfFS, WT

if __name__=="__main__":
    np.set_printoptions(linewidth=300, precision=4)
    smallAngles  = []
    op_point =[]
    extraSubs = []

#     model = get_model('F0T2N0S1', mergeFndTwr=True, linRot=False,
#                       yaw='zero', tilt='fixed',tiltShaft=True,
#                       rot_elastic_type='SmallRot', #rot_elastic_type='Body', 'Body' or 'SmallRot'
#                       orderMM=1,
#                       orderH=1,
#                       twrDOFDir=['x','y','x','y'], # Order in which the flexible DOF of the tower are set
#                      )
#     model = get_model('R1S0B100', mergeFndTwr=True, linRot=False,
#                   nB=1, yaw='zero', tilt='zero',azimuth_init='fixed', 
#                   pitch='fixed', cone='fixed', r_hub='fixed', coneAtRotorCenter=True,
#                   tiltShaft=True,
#                   rot_elastic_type='Body', #rot_elastic_type='Body', 'Body' or 'SmallRot'
#                   rot_elastic_smallAngle=False, # Keep False nu needs to be second order
#                   orderMM=1, orderH=1,
#                   twrDOFDir=['x','y','x','y'], # Order in which the flexible DOF of the tower are set
#                   verbose=True,
#                  )

    # --- Template
    #model = get_model('R3S0B100', azimuth_init='fixed', pitch='fixed', cone='fixed', r_hub='fixed', coneAtRotorCenter=True, orderMM=1, orderH=1, verbose=True)
    #simDir='_R1S0B100/'
    #extraSubs    = model.shapeNormSubs # shape functions normalized to unity
    #extraSubs   +=[(theta_pitch,0)]
    #extraSubs   +=[(theta_cone,0)]
    #extraSubs   +=[(psi_0,0)]
    #extraSubs   +=[(gravity,0)]
    #smallAngles  = [    ]
    #smallAngles += [ (model.twr.vcList, 2)]
    #smallAngles += [ ([theta_tilt]    , 1)]
    #op_point =[]
    #op_point +=[  (diff(q_psi, time), Symbol('Omega_op')) ,(q_psi, Symbol('psi_op'))]
    #for b in model.blds:
    #    if isinstance(b, YAMSFlexibleBody):
    #        op_point +=[(diff(b.q[0], time),0)]
    #dfNL, dfLI, dfFS = main(model, simDir, extraSubs, smallAngles, generate=True, runSim=True)

    # --- Flap only
    model = get_model('R3S0B100', azimuth_init='fixed', pitch='fixed', cone='fixed', r_hub='fixed', coneAtRotorCenter=True, orderMM=2, orderH=1, verbose=True)
    simFile='_R3S0BXXX/MainFlap.fst'
    extraSubs   +=[(theta_pitch,0)]
    extraSubs   +=[(theta_cone,0)]
    extraSubs   +=[(psi_0,0)]
    dfNL, dfLI, dfFS, WT = main(model, simFile, extraSubs, smallAngles, op_point, generate=True, runSim=True)


    # --- Edge only
#     model = get_model('R3S0B010', azimuth_init='fixed', pitch='fixed', cone='fixed', r_hub='fixed', coneAtRotorCenter=True, orderMM=1, orderH=1, verbose=True)
#     extraSubs   +=[(theta_pitch,0)]
#     extraSubs   +=[(theta_cone,0)]
#     extraSubs   +=[(psi_0,0)]
#     simFile='_R3S0BXXX/MainEdge.fst'
#     dfNL, dfLI, dfFS, WT = main(model, simFile, extraSubs, smallAngles, op_point, generate=True, runSim=True)


    # --- Flap with gravity
    #model = get_model('R3S0B100', azimuth_init='fixed', pitch='fixed', cone='fixed', r_hub='fixed', coneAtRotorCenter=True, orderMM=2, orderH=1, verbose=True)
    #simFile='_R3S0BXXX/MainFlapGravity.fst'
    #extraSubs   +=[(theta_pitch,0)]
    #extraSubs   +=[(theta_cone,0)]
    #extraSubs   +=[(psi_0,0)]
    #dfNL, dfLI, dfFS, WT = main(model, simFile, extraSubs, smallAngles, op_point, generate=True, runSim=True)

    # --- Edge gravity
#     model = get_model('R3S0B010', azimuth_init='fixed', pitch='fixed', cone='fixed', r_hub='fixed', coneAtRotorCenter=True, orderMM=1, orderH=1, verbose=True)
#     extraSubs   +=[(theta_pitch,0)]
#     extraSubs   +=[(theta_cone,0)]
#     extraSubs   +=[(psi_0,0)]
#     simFile='_R3S0BXXX/MainEdgeGravity.fst'
#     dfNL, dfLI, dfFS, WT = main(model, simFile, extraSubs, smallAngles, op_point, generate=True, runSim=True)

    plt.show()

