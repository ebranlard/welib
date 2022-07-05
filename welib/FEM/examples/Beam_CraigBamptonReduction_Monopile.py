"""
Setup a FEM model of a cantilevered beam
A Craig-Bampton reduction is performed
The Guyan and Craig-Bampton modes are plotted

For the example, the beam has constant properties along its span, but it can readily be adapted to a generic beam


TODO: Work in progress

"""
import unittest
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from welib.FEM.fem_beam import *
from welib.FEM.reduction import CraigBampton
from welib.beams.theory import *

MyDir=os.path.dirname(__file__)

def CraigBamptonReduction(TopMass=False):
    # --- Parameters
    nel      = 10             # Number of elements along the beam
    BC       = 'clamped-free' # Boundary condition: free-free or clamped-free
    element  = 'frame3d'      # Type of element used in FEM
    if TopMass:
        # NOTE: you can use welib.yams.windturbine to compute RNA mass and inertia
        Mtop = 50000  # Top mass [kg]
        M_tip= rigidBodyMassMatrixAtP(m=Mtop, J_G=None, Ref2COG=None)
    else:
        M_tip=None

    # --- Structural data for uniform beam
    # TODO read SubDyn model of monopile instead
    E   = 210e9     # Young modulus [Pa] [N/m^2]
    G   = 79.3e9    # Shear modulus. Steel: 79.3  [Pa] [N/m^2]
    L   = 100       # Beam Length [m]
    rho = 7850      # Material density kg/m^3
    D  =  8.1966401  # Diameter
    t  =  0.04476779 # thickness
    A   = np.pi*( (D/2)**2 - (D/2-t)**2) # Area for annulus [m^2] 
    I   = np.pi/64*(D**4-(D-2*t)**4) # Second moment of area for annulus (m^4)
    Kt  = I                     # Torsion constant, same as I for annulus [m^4]
    Ip  = 2*I                   # Polar second moment of area [m^4]

    nNodes=30
    x   = np.linspace(-50,50,nNodes)
    EIy = E*I   * np.ones(nNodes)
    EIz = E*I   * np.ones(nNodes)
    m   = rho*A * np.ones(nNodes)
    EIx = E*Ip  * np.ones(nNodes)
    Kt =  Kt    * np.ones(nNodes)
    A  =  A     * np.ones(nNodes)

    # --- Compute FEM model and mode shapes
    FEM=cbeam(x,m=m,EIx=EIx,EIy=EIy,EIz=EIy,EA=E*A,A=A,E=E,G=G,Kt=Kt,
            element=element, BC=BC, M_tip=M_tip)


    # --- Craig-Bampton reduction
    MM     = FEM['MM']
    KK     = FEM['KK']
    xNodes = FEM['xNodes']
    IDOF_tip  = FEM['Nodes2DOF'][FEM['Elem2Nodes'][-1,:][1],:]
    Ileader=IDOF_tip-6 # TODO TODO TODO BC has removed DOF...
    print(Ileader)
    print(MM.shape)

    MMr, KKr, Phi_G, Phi_CB, f_G, f_CB,_,_ = CraigBampton(MM, KK, Ileader, nModesCB=4, Ifollow=None, F=None, DD=None)

    # TODO TODO TODO
    x=FEM['xNodes'][0,1:-1] 
    Q=Phi_G # NEED to Add x=0 and Tip
    print(Q.shape)
    print(xNodes.shape)
    GuyModes=np.column_stack((x,Q[0::6,0],Q[1::6,1], Q[2::6,2], Q[3::6,3], Q[2::6,4], Q[1::6,5]))
    nModesPlot=6
    print(GuyModes)
    df=pd.DataFrame(data=GuyModes, columns=['x','Phix','Phiy','Phiz','PhiVx','PhiVy','PhiVz'])
    df.to_csv('_Guyan.csv',index=False)



    # --- Show frequencies to screen
    print('Mode   Frequency  Label ')
    for i in np.arange(nModesPlot):
        print('{:4d} {:10.3f}   {:s}'.format(i+1,FEM['freq'][i],FEM['modeNames'][i]))

    # --- Plot mode components for first few modes
    #Q=FEM['Q']
    print(x.shape)
    print(Q[0::6,i].shape)

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
        #axes[i].set_title(FEM['modeNames'][i])
        if i==0:
            axes[i].legend()

    return FEM


if __name__=='__main__':
    FEM = CraigBamptonReduction() 

    plt.show()

if __name__=='__test__':
    pass
#     FEM, _ = UniformBeam()
#     MM        = FEM['MMr']
#     KK        = FEM['KKr']
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

#if __name__=='__export__':
#    from welib.tools.repo import export_figs_callback
#    FEM = CraigBamptonReduction() 
#    export_figs_callback(__file__)
