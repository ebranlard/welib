import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
# Local 
from welib.yams.flexibility import GeneralizedMCK_PolyBeam, GMBeam 
from welib.yams.windturbine import FASTWindTurbine
from welib.tools.clean_exceptions import *

MyDir=os.path.dirname(__file__)


def FAST2StructureInputs(FST_file):

    WT = FASTWindTurbine(FST_file)

    # --------------------------------------------------------------------------------}
    # --- Setting Wisdem parameters
    # --------------------------------------------------------------------------------{
#     pW=dict()
#     # Twr
#     pW['TowerSpan']     = twrProp['HtFract_[-]'].values*(ED['TowerHt']-ED['TowerBsHt']) # from 0 to tower length
#     pW['TowerMassDens'] = twrProp['TMassDen_[kg/m]'].values
#     pW['TowerFAStiff']  = twrProp['TwFAStif_[Nm^2]'].values
#     pW['TowerSSStiff']  = twrProp['TwSSStif_[Nm^2]'].values
#     pW['TowerDamp']     = np.array([ twr['TwrFADmp(1)'], twr['TwrFADmp(2)'], twr['TwrSSDmp(1)'], twr['TwrSSDmp(2)']]) # structural damping ratio 
#     pW['TowerCoeffs']   = np.array([[ twr['TwFAM1Sh(2)'], twr['TwFAM2Sh(2)'], twr['TwSSM1Sh(2)'], twr['TwSSM2Sh(2)']],
#                                     [ twr['TwFAM1Sh(3)'], twr['TwFAM2Sh(3)'], twr['TwSSM1Sh(3)'], twr['TwSSM2Sh(3)']],
#                                     [ twr['TwFAM1Sh(4)'], twr['TwFAM2Sh(4)'], twr['TwSSM1Sh(4)'], twr['TwSSM2Sh(4)']],
#                                     [ twr['TwFAM1Sh(5)'], twr['TwFAM2Sh(5)'], twr['TwSSM1Sh(5)'], twr['TwSSM2Sh(5)']],
#                                     [ twr['TwFAM1Sh(6)'], twr['TwFAM2Sh(6)'], twr['TwSSM1Sh(6)'], twr['TwSSM2Sh(6)']]])
# 
#     pW['TowerExp']  = np.arange(2,7) # Exponents for shape functions polynomial, OpenFAST typically use 2,3,4,5,6
#     pW['TowerHt']   = ED['TowerHt']
#     pW['TowerBsHt'] = ED['TowerBsHt']
#     pW['Gravity']   = ED['Gravity']
# 
#     # Fnd
#     pW['PtfmCMzt']  = ED['PtfmCMzt']
#     pW['PtfmMass']  = ED['PtfmMass']
#     pW['PtfmRIner'] = ED['PtfmRIner']
#     pW['PtfmPIner'] = ED['PtfmPIner']
#     pW['PtfmYIner'] = ED['PtfmYIner']
# 
#     # RNA
#     pW['ShftTilt']  = ED['ShftTilt']  
#     pW['x_NR']      = r_NR_inN[0]     # x-coord from N to R in nac-coord
#     pW['z_NR']      = r_NR_inN[2]     # z-coord from N to R in nac-coord
#     pW['x_RNAG']    = r_SGhub_inS[0]  # x-coord from N to RNA_COG in nac-coord
#     pW['z_RNAG']    = r_SGhub_inS[2]  # z-coord from N to RNA_COG in nac-coord
#     pW['M_RNA']     = M_RNA           # Total mass of RNA
#     pW['J_xx_RNA']  = IG_RNA[0,0]     # Inertia of RNA at RNA_COG in nac-coord # NOTE diagonal for now..
#     pW['J_yy_RNA']  = IG_RNA[1,1]     # Inertia of RNA at RNA_COG in nac-coord
#     pW['J_zz_RNA']  = IG_RNA[2,2]     # Inertia of RNA at RNA_COG in nac-coord

    # --- Dict needed by structural script 
    p = dict()
#     p['z_FG']     = pW['PtfmCMzt']
    p['M_F']      = WT.fnd.mass
#     p['J_xx_F']   = pW['PtfmRIner']
#     p['J_yy_F']   = pW['PtfmPIner']
#     p['J_zz_F']   = pW['PtfmYIner']
#     p['g']        = pW['Gravity']
#     p['tilt']     =-pW['ShftTilt']
#     p['x_NR']     = pW['x_NR']                    # x-coord from N to R in nac-coord
#     p['z_NR']     = pW['z_NR']                    # z-coord from N to R in nac-coord
#     p['x_RNAG']   = pW['x_RNAG']                  # x-coord from N to RNA_G in nac-coord
#     p['z_RNAG']   = pW['z_RNAG']                  # z-coord from N to RNA_G in nac-coord
#     p['M_RNA']    = pW['M_RNA']                   # Total mass of RNA
#     p['J_xx_RNA'] = pW['J_xx_RNA']                # Inertia of RNA at RNA_G in nac-coord
#     p['J_yy_RNA'] = pW['J_xx_RNA']                # Inertia of RNA at RNA_G in nac-coord
#     p['J_zz_RNA'] = pW['J_xx_RNA']                # Inertia of RNA at RNA_G in nac-coord
#     p['L_T']      = pW['TowerHt']-pW['TowerBsHt']
#     p['z_OT']     = pW['TowerBsHt']               # distance from "Origin" (MSL) to tower base
#     p['M_T']      = pTwr['MM'][0,0]
#     p['z_TG']     = pTwr['s_OG'][2]
#     p['J_xx_T']   = pTwr['J_G'][0,0]
#     p['J_yy_T']   = pTwr['J_G'][1,1]
#     p['J_zz_T']   = pTwr['J_G'][2,2]
    p['MM_T']     = WT.twr.MM
    p['Oe_T']     = WT.twr.Oe6
    p['Gr_T']     = WT.twr.Gr
    p['Ge_T']     = WT.twr.Ge
    p['v_yT1c']   = WT.twr.Bhat_t_bc[1,0] 
    p['DD_T']     = WT.twr.DD
    p['KK_T']     = WT.twr.KK
    return p


def main(model_name='F2T0RNA_fnd'):

    # --- Generate necessary Wisdem inputs (pW) to the frequency domain component
    # NOTE: This will be handled by Wisdem
    MyDir=os.path.dirname(__file__)
    fstFilename = os.path.join(MyDir,'../../../data/_Spar2DOFNoHydroNoAero/Main_Spar_ED.fst')
    #pW = FAST2WisdemInputsWT(fstFilename)

    # --- Convert Wisdem inputs to inputs necessary for the structural part 
    # NOTE: Keep this interface, this is internal to frequency domain component
    #p  = WisdemInputs2StructureInputs(pW)
    p  = FAST2StructureInputs(fstFilename)

    # --- Print parameters
    print('--------------------')
#     print('Strucural Parameters: ')
#     for k,v in p.items():
#         if hasattr(v,'__len__'):
#             print('{:10s}:\n{}'.format(k,v))
#         else:
#             print('{:10s}:{}'.format(k,v))
    
    import importlib
    model= importlib.import_module('_py.{}'.format(model_name))

    # --- Evaluate linear structural model
    u0=dict() # Inputs at operating points
    u0['T_a']= 0 # thrust at operating point
    q0  = np.zeros(7) # [x, y, z, phi_x, phi_y, phi_z, q_t1] at operating point
    qd0 = np.zeros(7) # velocities at operating point
    M_lin   = model.M_lin(q0,p)
    C_lin   = model.C_lin(q0,qd0,p,u0)
    K_lin   = model.K_lin(q0,qd0,p,u0) 
    B_lin   = model.B_lin(q0,qd0,p,u0)

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
    #u['T_a']= lambda t: 1000+10*np.sin(0.1*t)  # Thrust as function of time
    #t=0
    #MM      = model.mass_matrix(q0,p)
    #forcing = model.forcing(t,q0,qd0,p,u)
    #print('--------------------')
    #print('Mass Matrix: ')
    #print(MM)
    #print('--------------------')
    #print('Forcing: ')
    #print(forcing)


if __name__ == '__main__':
    np.set_printoptions(linewidth=300, precision=2)
    main()
