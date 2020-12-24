import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
# Local 
from welib.yams.flexibility import GeneralizedMCK_PolyBeam, GMBeam 
from welib.yams.windturbine import FASTWindTurbine
from welib.tools.clean_exceptions import *
from welib.system.mech_system import MechSystem
import welib.weio as weio

MyDir=os.path.dirname(__file__)

def FAST2StructureInputs(FST_file, model_name=None):
    
    #WT = FASTWindTurbine(FST_file, twrShapes=[0,2], nSpanTwr=4, algo='ElastoDyn')
    WT = FASTWindTurbine(FST_file, twrShapes=[0,2], nSpanTwr=50)
    print(WT.RNA)
    # --- Dict needed by structural script 
    p = dict()
    p['z_FG']     = WT.fnd.masscenter[2]
    p['M_F']      = WT.fnd.mass
    p['J_xx_F']   = WT.fnd.masscenter_inertia[0,0]
    p['J_yy_F']   = WT.fnd.masscenter_inertia[1,1]
    p['J_zz_F']   = WT.fnd.masscenter_inertia[2,2]
    p['g']        = WT.ED['Gravity']
    p['tilt']     =-WT.ED['ShftTilt']
    p['x_NR']     = WT.r_NR_inN[0]                    # x-coord from N to R in nac-coord
    p['z_NR']     = WT.r_NR_inN[2]                    # z-coord from N to R in nac-coord
    p['x_RNAG']   = WT.RNA.masscenter[0]            # x-coord from N to RNA_G in nac-coord
    p['z_RNAG']   = WT.RNA.masscenter[2]            # z-coord from N to RNA_G in nac-coord
    p['M_RNA']    = WT.RNA.mass                   # Total mass of RNA
    p['J_xx_RNA'] = WT.RNA.masscenter_inertia[0,0]           # Inertia of RNA at RNA_G in nac-coord
    p['J_yy_RNA'] = WT.RNA.masscenter_inertia[1,1]           # Inertia of RNA at RNA_G in nac-coord
    p['J_zz_RNA'] = WT.RNA.masscenter_inertia[2,2]           # Inertia of RNA at RNA_G in nac-coord
    p['J_zx_RNA'] = WT.RNA.masscenter_inertia[0,2]           # Inertia of RNA at RNA_G in nac-coord
    p['L_T']      = WT.twr.length
    p['z_OT']     = WT.twr.pos_global[2]         # distance from "Origin" (MSL) to tower base
    p['M_T']      = WT.twr.MM[0,0]
    p['z_TG']     = WT.twr.masscenter[2]
    p['J_xx_T']   = WT.twr.masscenter_inertia[0,0]
    p['J_yy_T']   = WT.twr.masscenter_inertia[1,1]
    p['J_zz_T']   = WT.twr.masscenter_inertia[2,2]
    p['Oe_T']     = WT.twr.Oe6
    p['Gr_T']     = WT.twr.Gr
    p['Ge_T']     = WT.twr.Ge
    p['MM_T']     = WT.twr.MM
    p['v_yT1c']   = WT.twr.Bhat_t_bc[1,0]  # Mode 1  3 x nShapes
    p['v_xT2c']   = WT.twr.Bhat_t_bc[0,1]  # Mode 2
    p['DD_T']     = WT.twr.DD
    p['KK_T']     = WT.twr.KK
    # Rotor (blades + hub)
    p['Jxx_R']    = WT.rotgen.masscenter_inertia[0,0]
    p['JO_R']     = WT.rotgen.masscenter_inertia[1,1]
    p['M_R']      = WT.rot.mass
    # Nacelle 
    p['J_xx_N']   = WT.nac.masscenter_inertia[0,0]
    p['J_yy_N']   = WT.nac.masscenter_inertia[1,1]
    p['J_zz_N']   = WT.nac.masscenter_inertia[2,2]
    p['M_N']      = WT.nac.mass
    p['x_NG']     = WT.nac.masscenter[0]            # x-coord from N to nac G in nac-coord
    p['z_NG']     = WT.nac.masscenter[2]            # z-coord from N to nac G in nac-coord
    print(WT.rot)

    # Mooring restoring
    p['z_TM']      = 0
    p['K_x_M']     = 0
    p['K_y_M']     = 0
    p['K_z_M']     = 0
    p['K_phi_x_M'] = 0
    p['K_phi_y_M'] = 0
    p['K_phi_z_M'] = 0
    # Buoyancy
    p['z_TB']      = 0

    #print('MT\n',WT.twr.MM[6:,6:])

#     print('>>>>>>>>> HACK')
#     p['Jxx_R']=4.370255E+07
#     p['Jxx_R']=3.863464E+07
    # 38677056
#  8.066928E+06   -7.943743E+08   0.000000E+00 
# -7.943743E+08    9.692683E+10   0.000000E+00
#  0.000000E+00    0.000000E+00   4.370255E+07  
#      x                y           tx            ty             tz                psi
#[x   8.066928E+06   0.000000E+00  0.000000E+00 -7.943743E+08   0.000000E+00    0.000000E+00
#[y   0.000000E+00   8.066928E+06  7.943743E+08  0.000000E+00 >-1.415862E+05    0.000000E+00
#[tx  0.000000E+00   7.943743E+08  9.694190E+10 -1.093750E-01 > 8.551839E+06  > 3.858149E+07
#[ty -7.943743E+08   0.000000E+00 -1.093750E-01  9.692683E+10   0.000000E+00    0.000000E+00
#[tz  0.000000E+00 >-1.415862E+05 >8.551839E+06  0.000000E+00 > 1.853499E+08  >-3.375442E+06
#[psi 0.000000E+00  0.000000E+00  >3.858149E+07  0.000000E+00 >-3.375442E+06    4.370255E+07


#         y           tx              tz                psi
#[y    8.066928E+06  7.943743E+08  >-1.415862E+05    0.000000E+00
#[tx   7.943743E+08  9.694190E+10  > 8.551839E+06  > 3.858149E+07
#[tz >-1.415862E+05 >8.551839E+06  > 1.853499E+08  >-3.375442E+06
#[psi 0.000000E+00  >3.858149E+07  >-3.375442E+06    4.370255E+07

    return p,WT

def simulate(fstFilename, model_name, sims, sim_name):

    # ---
    p, WT = FAST2StructureInputs(fstFilename, model_name)

    import importlib
    model= importlib.import_module('_py.{}'.format(model_name))

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
    u0['F_B']= 0 # Buoyancy at operating point

    M_lin   = model.M_lin(q0l,p)
    C_lin   = model.C_lin(q0l,qd0l,p,u0)
    K_lin   = model.K_lin(q0l,qd0l,p,u0) 
    B_lin   = model.B_lin(q0l,qd0l,p,u0)

    # --- Print linearized mass damping 
    print('--------------------')
    print('Linear Mass Matrix: ')
    print(M_lin)
    print('--------------------')
    print('Linear Damping Matrix: ')
    print(C_lin)
    print('--------------------')
    print('Linear Stifness Matrix: ')
    print(K_lin)
    print('--------------------')
    print('Linear RHS: ')
    print(B_lin)

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
    print('--------------------')
    print('Mass Matrix: ')
    print(MM)
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
    DOFscales=[]
    for s in channels:
        if s.find('[deg]')>0:
            DOFscales.append(180/np.pi)
        elif s.find('TSS1')>0:
            DOFscales.append(-1)
        elif s.find('[rpm]')>0:
            DOFscales.append(60/(2*np.pi))
        else:
            DOFscales.append(1)

    dfNL = sysNL.toDataFrame(WT.channels, DOFscales)
    dfLI = sysLI.toDataFrame(WT.channels, DOFscales)
    sysNL.save(fstFilename.replace('.fst','_NonLinear.csv'), WT.channels, DOFscales)
    sysLI.save(fstFilename.replace('.fst','_Linear.csv'), WT.channels, DOFscales)

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

    # --- Rigid "F5"
    fstFilename = '_F5T0N0S1/Main_Spar_ED.fst'; model_name='F5T0N0S1_fnd'; sim_name='F5T0N0S1'

    # --- Flexibility "T1, T2"
    #fstFilename = '_F0T1RNA/Main_Spar_ED.fst'; model_name='F0T1RNA'; sim_name='F0T1RNA'

    #fstFilename = '_F0T2_NoRNA_sym/Main_Spar_ED.fst'; model_name='F0T2RNA'; sim_name='F0T2_NoRNA_sym'  # NOTE: Works fine large disp, symmetric shapes, with HubMass and NacMass, Twr2Shaft, detoriate slightly with overhang 
    #fstFilename = '_F0T2_NoRNA/Main_Spar_ED.fst'; model_name='F0T2RNA'; sim_name='F0T2_NoRNA'  # NOTE: with asymmetric shape functions, cannot achieve as good a result somehow. Wrong alpha???

    #fstFilename = '_F0T2RNA/Main_Spar_ED.fst'; model_name='F0T2RNA'; sim_name='F0T2RNA'
    #fstFilename = '_F0T2RNA_sym/Main_Spar_ED.fst'; model_name='F0T2RNA'; sim_name='F0T2RNA_sym'

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
