"""
Use a finite element method (FEM) formulation using frame elements to compute the mode shapes 
of a uniform beam cantilevered with and without top mass.
"""
import unittest
import os
import numpy as np
import matplotlib.pyplot as plt

from welib.FEM.fem_beam import *
from welib.beams.theory import *

MyDir=os.path.dirname(__file__)

def UniformBeam(TopMass=False):
    # --- Parameters
    nel      = 10             # Number of elements along the beam
    BC       = 'clamped-free' # Boundary condition: free-free or clamped-free
    element  = 'frame3d'      # Type of element used in FEM
    if TopMass:
        Mtop = 50000  # Top mass [kg]
        M_tip= rigidBodyMassMatrixAtP(m=Mtop, J_G=None, Ref2COG=None)
    else:
        M_tip=None

    # --- Structural data for uniform beam
    E   = 210e9     # Young modulus [Pa] [N/m^2]
    G   = 79.3e9    # Shear modulus. Steel: 79.3  [Pa] [N/m^2]
    L   = 100       # Beam Length [m]
    EIy0= 1.654e+12 # Planar second moment of area [m^4]
    m0  = 1.026e+04 # Mass per length [kg/m]
    EIx0= EIy0*2    # Polar second moment of area [m^4]
    A   = 1.00      # Area [m^2] 
    Kt  = EIy0/E*10 # Torsion constant [m^4]

    # --- Compute FEM model and mode shapes
    FEM=cbeam(L,m=m0,EIx=EIx0,EIy=EIy0,EIz=EIy0,EA=E*A,A=A,E=E,G=G,Kt=Kt,
            element=element, nel=nel, BC=BC, M_tip=M_tip)

    # --- Theory
    if TopMass:
        Theory = UniformBeamBendingModes('unloaded-topmass-clamped-free',EIy0,m0,A=1,L=L,nModes=3,Mtop=Mtop)
    else:
        Theory = UniformBeamBendingModes('unloaded-clamped-free',EIy0,m0,A=1,L=L,nModes=3)

    return FEM, Theory

def compareFEMwithTheory():
    FEM , Theory  = UniformBeam()
    FEM2, Theory2 = UniformBeam(TopMass=True)

    # --- Compare modes with theory
    fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
    fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
    COLORS=plt.rcParams['axes.prop_cycle'].by_key()['color']

    x =FEM['xNodes'][0,:]
    Q =FEM ['Q']
    Q2=FEM2['Q']
    QY1 = [Q [1::6, iMode] for iMode in range(10) if FEM ['modeNames'][iMode].startswith('uy')]
    QY2 = [Q2[1::6, iMode] for iMode in range(10) if FEM2['modeNames'][iMode].startswith('uy')]

    for iModeY in range(3):
        if iModeY==0:
            lblTM='FEM (top mass)'; lblTH='Theory (no top mass)'; lblTH2='Theory (top mass)'
        else:
            lblTM=None; lblTH=None; lblTH2=None
        x_th = Theory[1]
        U_th = Theory[2]
        x_th2= Theory2[1]
        U_th2= Theory2[2]
        ax.plot(x_th , U_th[iModeY] , 'k-', label=lblTH) # Theory No top mass 
        ax.plot(x_th2, U_th2[iModeY], 'k--', label=lblTH2) # Theory Top mass 
        ax.plot(x   , QY2[iModeY], 'o', c=COLORS[iModeY], label=lblTM)    # FEM top mass
        ax.plot(x   , QY1[iModeY], 's', c=COLORS[iModeY], label='FEM (no top mass) Mode {}'.format(iModeY+1)) # FEM no top mass
    ax.set_xlabel('Beam span [m]')
    ax.set_ylabel('Deflection [m]')
    ax.tick_params(direction='in')
    ax.set_title('FEM - mode shapes of a beam')
    ax.legend()

    return FEM, FEM2

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
    FEM, FEM2= compareFEMwithTheory()

    plotModes(FEM)

    plt.show()

if __name__=='__test__':
    FEM, _ = UniformBeam()
    MM        = FEM['MM']
    KK        = FEM['KK']
    xNodes    = FEM['xNodes']
    Q         = FEM['Q']
    freq      = FEM['freq']
    modeNames = FEM['modeNames']

    np.testing.assert_almost_equal(MM[0,0], 68400.0, 5)
    np.testing.assert_almost_equal(MM[-1,-1], 97714.285714, 5)
    np.testing.assert_almost_equal(MM[11,11], 195428.571428, 5)
    np.testing.assert_almost_equal(MM[20,20], 76217.142857, 5)
    np.testing.assert_almost_equal(KK[7,7]/1e10, 3.9696, 5)
    np.testing.assert_almost_equal(KK[10,10]/1e11, 13.232, 5)
    np.testing.assert_almost_equal(KK[11,11]/1e12, 1.3232, 5)
    np.testing.assert_almost_equal(freq[0],     0.71050, 5)
    np.testing.assert_almost_equal(freq[1],     0.71050, 5)
    np.testing.assert_almost_equal(freq[2],     4.45278, 5)
    np.testing.assert_almost_equal(freq[-1], 1209.98056, 5)
    if modeNames[0] not in ['uz1','uy1']:
        raise Exception('First mode not labelled properly')
    if modeNames[2] not in ['uz2','uy2']:
        raise Exception('Second mode not labelled properly')
    np.testing.assert_equal(modeNames[7], 'vx1')

if __name__=='__export__':
    from welib.tools.repo import export_figs_callback
    FEM= compareFEMwithTheory()
    export_figs_callback(__file__)
