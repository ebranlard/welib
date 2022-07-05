"""
Use a finite element method (FEM) formulation using frame elements to compute the mode shapes 
of a uniform cantilevered beam with and without top mass.


NOTE: Script might be unfinished

NOTE: an ElastoDyn wrapper should be added in the future in welib/fast/examples
"""
import unittest
import os
import numpy as np
import matplotlib.pyplot as plt

from welib.FEM.fem_beam import *
import welib.weio as weio
from welib.tools.clean_exceptions import *

MyDir=os.path.dirname(__file__)

def ModeShapesElastoDynTower():
    # --- Parameters
    nel      = 100            # Number of elements along the beam
    BC       = 'clamped-free' # Boundary condition: free-free or clamped-free
    element  = 'frame3d'      # Type of element used in FEM
    useFASTModel = True       # Use FAST model to find relevant quantities (RNA mass, tower length, etc.)

    if useFASTModel:
        # --- Option 1 obtain data from FAST model
        # fst file
        fstFile=os.path.join(MyDir,'./../../../data/NREL5MW/Main_Onshore.fst')
        from welib.yams.windturbine import WindTurbineStructure
        WT = WindTurbineStructure().fromFAST(fstFile)
        print(WT.twr)
        print(WT.RNA)
        TowerLen = WT.twr.length
        RNAMass  = WT.RNA.mass
        RNA_J_G  = WT.RNA.masscenter_inertia
        RNA_G    = WT.RNA.masscenter
    else:
        # --- Option 2 user inputs
        TowerLen = 87.6
        RNAMass  = 349606.5
        RNA_G    = [-0.407744, 0. ,      1.966427]
        RNA_J_G = np.array([[42620975.3 , -0.       ,  -924869.9],
                            [      -0.  ,24376664.3 ,   -0.     ],
                            [ -924869.9 , -0.       , 22652099.3]])

    # --- Convert for OpenFAST coord to FEM coordinate system
    R_OF2FEM = np.array(
            [[0,0,1],
             [1,0,0],
             [0,1,0]])
    RNA_J_G= R_OF2FEM.dot(RNA_J_G).dot(R_OF2FEM.T)
    RNA_G  = R_OF2FEM.dot(RNA_G)
    print('RNA_J_G\n',np.around(RNA_J_G,6))
    print('RNA_G\n',RNA_G)

    # --- Set RNA rigid body mass matrix at tower top
    #M_tip= None
    M_tip= rigidBodyMassMatrixAtP(m=RNAMass, J_G=RNA_J_G*0, Ref2COG=RNA_G)
    print('M_Tip\n',np.around(M_tip,2))

    TwrFile=os.path.join(MyDir,'./../../../data/NREL5MW/5MW_Baseline/NRELOffshrBsline5MW_Onshore_ElastoDyn_Tower.dat')
    twr = weio.FASTInputFile(TwrFile).toDataFrame()
    x   = twr['HtFract_[-]']*(TowerLen)
    m   = twr['TMassDen_[kg/m]']
    EIz = twr['TwFAStif_[Nm^2]']        # FA/Flap = IyOpenFAST = Iz(when x is axis)
    EIy = twr['TwSSStif_[Nm^2]']        # SS/Edge = IxOpenFAST = Iy(when x is axis)

    # --- Create Beam FEM model
    # Derived parameters
    A  = m*0+100                         # Area
    Kt = m*0+100                         # Saint Venant torsion
    E   = 211e9                           # Young modulus  [N/m^2]
    G   = 79.3e9                         # Shear modulus. Steel: 79.3  [Pa] [N/m^2]
    Ip  = 2*(EIz/E)                      # Polar second moment of area [m^4]
    EIx = E*Ip 

    # --- Compute FEM model and mode shapes
    FEM=cbeam(x,m=m,EIx=EIx,EIy=EIy,EIz=EIy,EA=E*A,A=A,E=E,G=G,Kt=Kt,
            element=element, nel=nel, BC=BC, M_tip=M_tip)

    return FEM

def compareFEMwithTheory():
    # --- Compare modes with theory
    fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
    fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
    COLORS=plt.rcParams['axes.prop_cycle'].by_key()['color']

    x =FEM['xNodes'][0,:]
    Q =FEM ['Q']
    QY1 = [Q [1::6, iMode] for iMode in range(10) if FEM ['modeNames'][iMode].startswith('uy')]

    for iModeY in range(3):
        ax.plot(x   , QY1[iModeY], 's', c=COLORS[iModeY], label='FEM (no top mass) Mode {}'.format(iModeY+1)) # FEM no top mass
    ax.set_xlabel('Beam span [m]')
    ax.set_ylabel('Deflection [m]')
    ax.tick_params(direction='in')
    ax.set_title('FEM - mode shapes of a beam')
    ax.legend()

    return FEM

def plotModes(FEM):
    nModesPlot=8
    # --- Show frequencies to screen
    print('Mode   Frequency  Label ')
    for i in np.arange(nModesPlot):
        print('{:4d} {:10.3f}   {:s}'.format(i+1,FEM['freq'][i],FEM['modeNames'][i]))

    # --- Plot mode components for first few modes
    x=FEM['xNodes'][0,:]
    Q=FEM['Q']

    fig,axes = plt.subplots(1, nModesPlot, sharey=False, figsize=(12.4,2.5))
    fig.subplots_adjust(left=0.04, right=0.98, top=0.91, bottom=0.11, hspace=0.40, wspace=0.30)
    for i in np.arange(nModesPlot):
        axes[i].plot(x, Q[0::6,i]  ,'-'  , label='ux')
        axes[i].plot(x, Q[1::6,i]  ,'-'  , label='uy')
        axes[i].plot(x, Q[2::6,i]  ,'-'  , label='uz')
        axes[i].plot(x, Q[3::6,i]  ,':'  , label='vx')
        axes[i].plot(x, Q[4::6,i]  ,':'  , label='vy')
        axes[i].plot(x, Q[5::6,i]  ,':'  , label='vz')
        axes[i].set_xlabel('')
        axes[i].set_ylabel('')
        axes[i].set_title(FEM['modeNames'][i])
        if i==0:
            axes[i].legend()

if __name__=='__main__':
    FEM  = ModeShapesElastoDynTower()

    plotModes(FEM)

    plt.show()

if __name__=='__test__':
    FEM  = ModeShapesElastoDynTower()
#     MM        = FEM['MM']
#     KK        = FEM['KK']
#     xNodes    = FEM['xNodes']
#     Q         = FEM['Q']
#     freq      = FEM['freq']
#     modeNames = FEM['modeNames']
# 
#     np.testing.assert_almost_equal(MM[0,0], 68400.0, 5)
#     np.testing.assert_almost_equal(MM[-1,-1], 97714.285714, 5)
#     np.testing.assert_almost_equal(MM[11,11], 195428.571428, 5)
#     np.testing.assert_almost_equal(MM[20,20], 76217.142857, 5)
#     np.testing.assert_almost_equal(KK[7,7]/1e10, 3.9696, 5)
#     np.testing.assert_almost_equal(KK[10,10]/1e11, 13.232, 5)
#     np.testing.assert_almost_equal(KK[11,11]/1e12, 1.3232, 5)
#     np.testing.assert_almost_equal(freq[0],     0.71050, 5)
#     np.testing.assert_almost_equal(freq[1],     0.71050, 5)
#     np.testing.assert_almost_equal(freq[2],     4.45278, 5)
#     np.testing.assert_almost_equal(freq[-1], 1209.98056, 5)
#     if modeNames[0] not in ['uz1','uy1']:
#         raise Exception('First mode not labelled properly')
#     if modeNames[2] not in ['uz2','uy2']:
#         raise Exception('Second mode not labelled properly')
#     np.testing.assert_equal(modeNames[7], 'vx1')

if __name__=='__export__':
    pass
#     from welib.tools.repo import export_figs_callback
#     FEM= compareFEMwithTheory()
#     export_figs_callback(__file__)
