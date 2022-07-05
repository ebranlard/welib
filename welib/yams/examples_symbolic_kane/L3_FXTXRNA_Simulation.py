import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
# Local 
from welib.yams.flexibility import GeneralizedMCK_PolyBeam, GMBeam 
from welib.yams.windturbine import FASTWindTurbine
from welib.tools.clean_exceptions import *
from welib.system.mech_system import MechSystem
from welib.system.eva import *
import welib.weio as weio

MyDir=os.path.dirname(__file__)

def FAST2StructureInputs(FST_file, model_name=None):
    #WT = FASTWindTurbine(FST_file, twrShapes=[0,2], nSpanTwr=4, algo='ElastoDyn')
    WT = FASTWindTurbine(FST_file, twrShapes=[0,2], nSpanTwr=50)
    p = WT.yams_parameters()
    return p,WT

def readEDSummaryFile(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()
        i=0
        for l in lines:
            if l.find('Full system mass matrix')>=0:
                M=np.zeros((24,24))
                for j in np.arange(24):
                    v=np.asarray(lines[j+i+1].strip()[2:].split()).astype(float)
                    M[j,:]=v
                I = ~np.all(M==0, axis=0)
                M = M[I][:,I]
                return M
            i+=1
def readFSTLin(filename):
    try:
        lin = weio.FASTLinearizationFile(filename)
    except:
        return None
    A = lin['A']

    return lin


def simulate(fstFilename, model_name, sims, sim_name):
    # ---
    p, WT = FAST2StructureInputs(fstFilename, model_name)

    print('>>>>>',p['tilt'])
    import importlib
    model= importlib.import_module('_py.{}'.format(model_name))


    MM_ED = readEDSummaryFile(fstFilename.replace('.fst','.ED.sum'))
    lin   = readFSTLin(fstFilename.replace('.fst','.1.lin'))
#     freq_d2, zeta2, Q2, freq02  = eigA(lin['A'])

    # --- Reference simulation
    df=weio.read(fstFilename.replace('.fst','.out')).toDataFrame()
    #time = np.linspace(0,50,5000)
    time = df['Time_[s]'].values

    # --- Initial conditions
    nDOFExpected=np.sum([int(s) for s in model_name if s.isdigit()])
    print('DOFs:', WT.DOFname, 'Model:',model_name, 'nDOF:',nDOFExpected )
    if len(WT.DOFname)!=nDOFExpected:
        raise Exception('Inconsistency in number of DOFs')
    q0  = WT.q0
    qd0 = WT.qd0
    print('q0 :',q0)
    print('qd0:',qd0)
    q0l = WT.q0*0
    qd0l= WT.qd0*0
    DOFs = WT.activeDOFs
    print('q0l :',q0l)
    print('qd0l:',qd0l)

    # --- Evaluate linear structural model
    u0=dict() # Inputs at operating points
    u0['T_a']= 0 # thrust at operating point # TODO
    u0['M_y_a']= 0 # aerodynamic tilting moment at operating point
    u0['M_z_a']= 0 #+0*np.sin(0.1*t)  # aerodynamic "yawing" moment
    u0['F_B']= 0 # Buoyancy t(at operating point

    M_lin   = model.M_lin(q0l,p)
    C_lin   = model.C_lin(q0l,qd0l,p,u0)
    K_lin   = model.K_lin(q0l,qd0l,p,u0) 
    B_lin   = model.B_lin(q0l,qd0l,p,u0)

    # --- Print linearized mass damping 
#     print('--------------------')
#     print('Linear Mass Matrix: ')
#     print(M_lin)
#     print('--------------------')
#     print('Linear Damping Matrix: ')
#     print(C_lin)
#     print('--------------------')
#     print('Linear Stifness Matrix: ')
#     print(K_lin)
#     print('--------------------')
#     print('Linear RHS: ')
#     print(B_lin)

    # --- Non linear
    u=dict()
    u['T_a']= lambda t: 0 #+0*np.sin(0.1*t)  # Thrust as function of time # TODO
    u['M_y_a']= lambda t: 0 #+0*np.sin(0.1*t)  # aerodynamic tilting moment
    u['M_z_a']= lambda t: 0 #+0*np.sin(0.1*t)  # aerodynamic "yawing" moment
    u['F_B']= lambda t: 0 #+0*np.sin(0.1*t)  # Thrust as function of time # TODO
    #u['F_B']= lambda t: (p['M_F']  + p['M_RNA'] + p['MM_T'][0,0])*p['g']
    t=0
    MM      = model.mass_matrix(q0,p)
    forcing = model.forcing(t,q0,qd0,p,u)
    #print(WT.WT_rigid)
    print('--------------------')
    print('Mass Matrix: ')
    print(MM)
    print(MM_ED)
    M_Ref=MM_ED.copy()
    print(' Rel Error')
    MM   [np.abs(MM   )<1e-1]=1
    MM_ED[np.abs(MM_ED)<1e-1]=1
    M_Ref[np.abs(M_Ref)<1e-1]=1
    print(np.around(np.abs((MM_ED-MM))/M_Ref*100,1))
    print('--------------------')
    print('Forcing: ')
    print(forcing)

    if sims is False:
        return p, WT, None, None, None, None 


    # --- integrate non-linear system
    fM = lambda x: model.mass_matrix(x, p)
    fF = lambda t,x,xd: model.forcing(t, x, xd, p=p, u=u)
    sysNL = MechSystem(fM, F=fF, x0=q0, xdot0=qd0 )
    resNL=sysNL.integrate(time, method='RK45')

    # --- integrate linear system
    fF = lambda t,x,xd: np.array([0]*len(q0))
    sysLI = MechSystem(M=M_lin, K=K_lin, C=C_lin, F=fF, x0=q0, xdot0=qd0)
    resLI=sysLI.integrate(time, method='RK45') # **options):
    
    # --- Convert results to dataframe and save to file
    channels = WT.channels
    DOFscales= WT.FASTDOFScales
    dfNL = sysNL.toDataFrame(WT.channels, DOFscales)
    dfLI = sysLI.toDataFrame(WT.channels, DOFscales)
    sysNL.save(fstFilename.replace('.fst','_NonLinear.csv'), WT.channels, DOFscales)
    sysLI.save(fstFilename.replace('.fst','_Linear.csv'), WT.channels, DOFscales)

    try:
        print('>>>>> A')
        print(sysLI.A)
        print('>>>>> A')
        print(lin['A'])
        freq_d, zeta, Q, freq0  = eigA(sysLI.A)
        freq_d2, zeta2, Q2, freq02  = eigA(lin['A'])
        print(freq_d)
        print(freq_d2)
        print(zeta)
        print(zeta2)
    except:
        print('>> No lin')

    for c in dfNL.columns:
        if c.find('[deg]')>1: 
            if np.max(dfNL[c])>90:
                dfNL[c]= np.mod(dfNL[c],360)
                dfLI[c]= np.mod(dfLI[c],360)


    # --- Plot
    # sys.plot()
    legDone=False
    nDOF=sysNL.nDOF
    fig,axes = plt.subplots(nDOF, 2, sharey=False, sharex=True, figsize=(12.0,8.0)) # (6.4,4.8)
    axes = axes.reshape(nDOF,2)
    fig.subplots_adjust(left=0.08, right=0.99, top=0.98, bottom=0.06, hspace=0.07, wspace=0.18)
    for idof, dof in enumerate(WT.activeDOFs):
        # Positions
        chan=channels[idof]
        axes[idof,0].plot(dfNL['Time_[s]'], dfNL[chan], '-'  , label='non-linear')
        axes[idof,0].plot(dfLI['Time_[s]'], dfLI[chan], '--' , label='linear')
        if chan in df.columns:
            axes[idof,0].plot(df['Time_[s]'], df[chan], 'k:' , label='OpenFAST')
            if not legDone:
                legDone=True
                axes[idof,0].legend(loc='upper right')
        axes[idof,0].tick_params(direction='in')
        axes[idof,0].set_ylabel(chan.replace('_',' '))

        # Velocities
        vdof = idof+nDOF
        chan=channels[vdof]
        axes[idof,1].plot(dfNL['Time_[s]'], dfNL[chan], '-'  , label='non-linear')
        axes[idof,1].plot(dfLI['Time_[s]'], dfLI[chan], '--' , label='linear')
        if chan in df.columns:
            axes[idof,1].plot(df['Time_[s]'], df[chan], 'k:' , label='OpenFAST')
        axes[idof,1].tick_params(direction='in')
        axes[idof,1].set_ylabel(chan.replace('_',' '))

        if idof==nDOF-1:
            axes[idof,0].set_xlabel('Time [s]')
            axes[idof,1].set_xlabel('Time [s]')

    fig.savefig('_figs/{}.png'.format(sim_name))
    plt.show()

    return p, WT, sysNL, dfNL, sysLI, dfLI


def main():
    # --- Rigid "F2"
    #fstFilename = '_F2T0RNANoRefH/Main_Spar_ED.fst' ;model_name='F2T0RNA_fnd';sim_name='F2T0RNA_NoRefH'
    #fstFilename = '_F2T0_NoRNA_NoRefH/Main_Spar_ED.fst' ;model_name='F2T0RNA_fnd';sim_name='F2T0RNA_NoRNA_NoRefH'
    #fstFilename = '_F2T0RNA/Main_Spar_ED.fst'      ;model_name='F2T0RNA_fnd';sim_name='F2T0RNA'

    #fstFilename = '_F2T0N0S1/Main_Spar_ED.fst'; model_name='F2T0N0S1_fnd'; sim_name='F2T0N0S1'

    # Phi x phi z
    #fstFilename = '_F000101T0RNA/Main_Spar_ED.fst'; model_name='B000101'; sim_name=model_name
    #fstFilename = '_F000101T0RNA/Main_Spar_ED.fst'; model_name='F000101T0RNA_fnd'; sim_name=model_name
    #fstFilename = '_F000101T0N0S1/Main_Spar_ED.fst'; model_name='F000101T0N0S1_fnd'; sim_name=model_name

    # Phi y phi z
    #fstFilename = '_F000011T0RNA/Main_Spar_ED.fst'; model_name='B000011'; sim_name=model_name

    # Phi x phi y phi z 
    #fstFilename = '_F000111T0RNA/Main_Spar_ED.fst'; model_name='B000111'; sim_name=model_name
    #fstFilename = '_F000111T0RNA/Main_Spar_ED.fst'; model_name='F000111T0RNA_fnd'; sim_name=model_name

    # --- Rotor
    #fstFilename = '_F000111T0N0S1/Main_Spar_ED.fst'; model_name='F000111T0N0S1_fnd'; sim_name=model_name
    #fstFilename = '_F5T0N0S1/Main_Spar_ED.fst'; model_name='F5T0N0S1_fnd'; sim_name='F5T0N0S1'
    fstFilename = '_F5T1N0S1/Main_Spar_ED.fst'; model_name='F5T1N0S1_fnd'; sim_name='F5T1N0S1'

    # --- Flexibility "T1, T2"
    #fstFilename = '_F0T1RNA/Main_Spar_ED.fst'; model_name='F0T1RNA'; sim_name='F0T1RNA'

    #fstFilename = '_F0T2_NoRNA_sym/Main_Spar_ED.fst'; model_name='F0T2RNA'; sim_name='F0T2_NoRNA_sym'  # NOTE: Works fine large disp, symmetric shapes, with HubMass and NacMass, Twr2Shaft, detoriate slightly with overhang 
    #fstFilename = '_F0T2_NoRNA/Main_Spar_ED.fst'; model_name='F0T2RNA'; sim_name='F0T2_NoRNA'  # NOTE: with asymmetric shape functions, cannot achieve as good a result somehow. Wrong alpha???

    #fstFilename = '_F0T2RNA/Main_Spar_ED.fst'; model_name='F0T2RNA'; sim_name='F0T2RNA'
    #fstFilename = '_F0T2RNA_sym/Main_Spar_ED.fst'; model_name='F0T2RNA'; sim_name='F0T2RNA_sym'

    fstFilename = '_F0T2N0S1/Main_Spar_ED.fst'; model_name='_F0T2N0S1'; sim_name=model_name;

    # --- Floater + Flexibility "F2T1"
    #fstFilename = '_F2T1RNANoRefH/Main_Spar_ED.fst'; model_name='F2T1RNA_fnd'; sim_name='F2T1RNA_NoRefH'
    #fstFilename = '_F2T1RNA_SmallAngle/Main_Spar_ED.fst'; model_name='F2T1RNA_fnd'; sim_name='F2T1RNA_SmallAngle'
    #fstFilename = '_F2T1RNA/Main_Spar_ED.fst'; model_name='F2T1RNA_fnd'; sim_name='F2T1RNA_LargeAngle'

    # --- "F3 T1"
    #fstFilename = '_F3T1RNA/Main_Spar_ED.fst'; model_name='F3T1RNA_fnd'; sim_name='F3T1RNA'
    # --- "F5 T1"
    #fstFilename = '_F5T1RNA/Main_Spar_ED.fst'; model_name='F5T1RNA_fnd'; sim_name='F5T1RNA'


    sim = True
    p, WT, sysNL, resNL, sysL, resL  = simulate(fstFilename, model_name, sim, sim_name)

    # --- Print parameters
#     print('--------------------')
#     print('Strucural Parameters: ')
#     for k,v in p.items():
#         if hasattr(v,'__len__'):
#             print('{:10s}:\n{}'.format(k,v))
#         else:
#             print('{:10s}:{}'.format(k,v))

if __name__ == '__main__':
    np.set_printoptions(linewidth=300, precision=5)
    main()
