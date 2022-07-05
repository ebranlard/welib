""" 
TODO

NOTE: a SubDyn wrapper with a higher level interface should be available in the future in welib/fast/examples
"""

import unittest
import os
import scipy
import numpy as np
from numpy.linalg import inv

from welib.FEM.fem_beam import applyBC, generalizedMassMatrix, shapeIntegrals
from welib.FEM.fem_beam import geometricalStiffening
from welib.FEM.fem_beam import orthogonalizeModePair, normalize_to_last
from welib.FEM.fem_beam import cbeam_assembly_frame3dlin, cbeam_frame3dlin_Kg
from welib.FEM.fem_beam import cbeam_assembly_frame3d
from welib.FEM.fem_beam import cbeam_assembly
from welib.FEM.frame3dlin import frame3dlin_Mcross, frame3dlin_Kg

from welib.tools.clean_exceptions import *

import welib.weio as weio
from welib.system.eva import eig
from welib.yams.sid import FEMBeam2SID
import matplotlib.pyplot as plt


MyDir=os.path.dirname(__file__)

np.set_printoptions(linewidth=300, precision=4)

def OpenFASTIsolatedTower():
    # --- Read data from NREL5MW tower
    TowerHt=87.6;
    TowerBs=0;
    TwrFile=os.path.join(MyDir,'./../../../data/NREL5MW/5MW_Baseline/NRELOffshrBsline5MW_Onshore_ElastoDyn_Tower.dat')
    twr = weio.FASTInputFile(TwrFile).toDataFrame()
    z = twr['HtFract_[-]']*(TowerHt-TowerBs)
    m = twr['TMassDen_[kg/m]']  # mu
    EIy = twr['TwFAStif_[Nm^2]'] 
    EIz = twr['TwSSStif_[Nm^2]']  # TODO actually EIx

    # --- Create Beam FEM model
    # Derived parameters
    A  = m*0+100                         # Area
    Kv = m*0+100                         # Saint Venant torsion
    E  = 214e9                           # Young modulus  [N/m^2]
    Iy = EIy/E                           # Area moment [m^4]
    Iz = EIz/E                           # Area moment [m^4]
    nNodes   = len(z)
    nElem    = nNodes-1

    # Nodes positions
    xNodes = np.zeros((3,nNodes))
    xNodes[2,:]=z

    # Assembly
    MM, KK, xNodes, DCM, Elem2Nodes, Nodes2DOF, Elem2DOF = cbeam_assembly_frame3dlin(xNodes, m, Iy, Iz=Iz, A=A, Kv=Kv, E=E)

    # --- Constraints/ BC
    MMr, KKr, Tr,_,_ = applyBC(MM, KK, Elem2Nodes, Nodes2DOF, BC_root=[0,0,0,0,0,0], BC_tip=[1,1,1,1,1,1])
    iStart= 0; 

    # --- Eigenvalues/vectors
    [Q, freq]= eig(KKr, MMr, freq_out=True)

    # --- Orthogonalization/ normalization of modes
    Imodes=[0,1]
    #Q[:,0],Q[:,1] = orthogonalizeModePair(Q[:,0],Q[:,1], iStart)
    #Q= normalize_to_last(Q, Imodes, iStart);


    # --- Export Modes
    U1 = np.concatenate(([0],Q[0::6,0] )) # Deflection mode 1, along x
    V1 = np.concatenate(([0],Q[4::6,0] )) # Slope mode 1 , theta y
    U2 = np.concatenate(([0],Q[1::6,1] )) # Deflection mode 2, along y
    V2 = np.concatenate(([0],Q[3::6,1] )) # Slope mode 2, theta x
    #print(U1)
    #print(U2)
    #print(Q[:,0])
    #print(Q[:,1])
    #print(Q[:,2])
    #print(Q[:,3])
    #M=np.column_stack([z,U1,V2,U2,V2])
    # np.savetxt('out.csv',M)

    # --- Generalized mass matrix
    # Selecting modes
    Imodes=[0,1]
    if len(Imodes)>0:
        Se= Tr.dot(Q[:, Imodes]) # nDOF_tot x nShapes
    else:
        Se= Tr # All
    Mtt, J0, Mrt, Mgt, Mgr, Mgg, St, Sr= generalizedMassMatrix(xNodes, MM, Se)


    #fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
    #fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
#   #  ax.plot(z, U1    , label='Mode 1')
    #ax.plot(z, U2    , label='Mode 2')
    #ax.set_xlabel('')
    #ax.set_ylabel('')
    #ax.legend()

    return MM, KK, Q, freq

if __name__=='__main__':
    MM, KK, Q, freq = OpenFASTIsolatedTower()
    plt.show()

if __name__=='__test__':
    MM, KK, Q, freq = OpenFASTIsolatedTower()

    np.testing.assert_almost_equal(MM[0,0], 17921.9563543, 5)
    np.testing.assert_almost_equal(MM[-1,-1], 7590.2188  , 5)
    np.testing.assert_almost_equal(MM[11,11], 30565.98330, 5)
    np.testing.assert_almost_equal(MM[20,20], 26585.67290, 5)
    np.testing.assert_almost_equal(KK[7,7]/1e10, 1.91655, 5)
    np.testing.assert_almost_equal(KK[10,10]/1e11, 4.893305, 5)
    np.testing.assert_almost_equal(KK[11,11]/1e12, 1.87917, 5)
    np.testing.assert_almost_equal(freq[0],     0.891449, 5)
    np.testing.assert_almost_equal(freq[1],     0.891449, 5)
    np.testing.assert_almost_equal(freq[-1], 5250.756553, 5)


