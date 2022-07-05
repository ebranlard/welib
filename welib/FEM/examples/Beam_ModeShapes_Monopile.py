"""
Use a finite element method (FEM) formulation using frame elements to compute the mode shapes 
of a uniform beam cantilevered with and without top mass.

NOTE: all of this can be done in few lines using the class: welib.fast.subdyn.SubDyn

NOTE: Script might be unfinished

NOTE: a SubDyn wrapper with a higher level interface should be available in the future in welib/fast/examples

"""
import unittest
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from welib.FEM.fem_beam import *
from welib.FEM.reduction import CraigBampton, augmentModes

import welib.weio as weio

MyDir=os.path.dirname(__file__)


def MonopileFEM(TopMass=False):
    # --- Parameters
    BC       = 'clamped-free' # Boundary condition: free-free or clamped-free
    element  = 'frame3d'      # Type of element used in FEM
    UseSubDynModel=True
    TopMass = False

    # Add an optional top mass and ineria
    if TopMass:
        # NOTE: you can use welib.yams.windturbine to compute RNA mass and inertia
        Mtop = 50000  # Top mass [kg]
        M_tip= rigidBodyMassMatrixAtP(m=Mtop, J_G=None, Ref2COG=None)
    else:
        M_tip=None

    if UseSubDynModel:
        # --- Option 1 - Read data from SubDyn
        # Read SubDyn file
        sdFilename=os.path.join(MyDir,'../../../data/Monopile/MT100_SD.dat')
        sd = weio.read(sdFilename)
        # Convert to "welib.fem.Graph" class to easily handle the model (overkill for a monopile)
        graph = sd.toGraph()
        graph.divideElements(sd['NDiv'])
        graph.sortNodesBy('z')
        df = graph.nodalDataFrame()
        x   = df['z'] # NOTE: FEM uses "x" as main axis
        D   = df['D'] # Diameter [m]
        t   = df['t'] # thickness [m]
        E   = df['E'] # Young modules [N/m^2]
        G   = df['G'] # Shear modules [N/m^2]
        rho = df['rho'] # material density [kg/m^3]
    else:
        # --- Option 2 - User specify data
        x    = np.linspace(-50,50,100) # Vertical distance along Monopile [m]
        ONES = np.ones(x.shape)
        E   = ONES*210e9     # Young modulus [Pa] [N/m^2]
        G   = ONES*79.3e9    # Shear modulus. Steel: 79.3  [Pa] [N/m^2]
        rho = ONES*7850      # Material density [kg/m^3]
        D  =  ONES*8.1966401  # Diameter [m]
        t  =  ONES*0.04476779 # thickness [m]

    # Derive section properties for a hollow cylinder based on diameter and thickness
    A   = np.pi*( (D/2)**2 - (D/2-t)**2) # Area for annulus [m^2] 
    I   = np.pi/64*(D**4-(D-2*t)**4) # Second moment of area for annulus (m^4)
    Kt  = I                     # Torsion constant, same as I for annulus [m^4]
    Ip  = 2*I                   # Polar second moment of area [m^4]
    L = np.max(x)-np.min(x) # Monopile length

    # --- Compute FEM model and mode shapes
    FEM=cbeam(x,m=rho*A,EIx=E*Ip,EIy=E*I,EIz=E*I,EA=E*A,A=A,E=E,G=G,Kt=Kt,
                element=element, BC=BC, M_tip=M_tip)

    # --- Perform Craig-Bampton reduction, fixing the top node of the beam
    Q_G,_Q_CB, df_G, df_CB, Modes_G, Modes_CB, CB = CB_topNode(FEM, nCB=8, element='frame3d', main_axis='x')

    df_CB.to_csv('_CB.csv',index=False)
    df_G.to_csv('_Guyan.csv',index=False)

    nModesPlot=8

    # --- Show frequencies to screen
    print('Mode   Frequency  Label ')
    for i in np.arange(8):
        print('{:4d} {:10.3f}   {:s}'.format(i+1,FEM['freq'][i],FEM['modeNames'][i]))

    # --- Plot mode components for first few modes
    print(x.shape)
    #Q=FEM['Q'] ; modeNames = FEM['modeNames']
    #Q=Q_CB ;modeNames = names_CB
    Modes=Modes_CB
    nModesPlot=min(len(Modes),nModesPlot)

    fig,axes = plt.subplots(1, nModesPlot, sharey=False, figsize=(12.4,2.5))
    fig.subplots_adjust(left=0.04, right=0.98, top=0.91, bottom=0.11, hspace=0.40, wspace=0.30)
    for i in np.arange(nModesPlot):
        key= list(Modes.keys())[i]

        axes[i].plot(x, Modes[key]['comp'][:,0]  ,'-'  , label='ux')
        axes[i].plot(x, Modes[key]['comp'][:,1]  ,'-'  , label='uy')
        axes[i].plot(x, Modes[key]['comp'][:,2]  ,'-'  , label='uz')
        axes[i].plot(x, Modes[key]['comp'][:,3]  ,':'  , label='vx')
        axes[i].plot(x, Modes[key]['comp'][:,4]  ,':'  , label='vy')
        axes[i].plot(x, Modes[key]['comp'][:,5]  ,':'  , label='vz')
        axes[i].set_xlabel('')
        axes[i].set_ylabel('')
        axes[i].set_title(Modes[key]['label'])
        if i==0:
            axes[i].legend()

if __name__=='__main__':
    FEM = MonopileFEM() 

    plt.show()


