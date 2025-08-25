""" Unsteady lifting surface"""

import numpy as np
from scipy.linalg import lu_factor, lu_solve
import os
from welib.tools.tictoc import Timer
import welib.weio as weio
from vtk import VTK_Misc, WrVTK_Lattice

# Constants
PI = 3.141592654
RCUT = 1.0e-10

# Parameters
PARAMS = {
    'NTMAX': 640,  # Max time steps
    'nChord': 4,    # Max chordwise panels
    'NSMAX': 13,   # Max spanwise panels
}

# Global data
class GLOBAL_DATA:
    def __init__(self):
        self.SNO = None
        self.CSO = None
        self.GAMA1J = None
        self.QF = None  # Wing nodes
        self.CP = None
        self.A = None
        self.DW = None
        self.GAMA1 = None
        self.DLY = None
        self.GAMA = None # Wing Gammas
        self.DL = None
        self.DP = None
        self.DS = None
        self.DLT = None
        self.DD = None
        self.A1 = None
        self.QW = None     # Wake nodes
        self.VORTIC = None # Wake vorticity Gammas
        self.UVW = None
        self.QW1 = None
        self.VORT1 = None
        self.US = None
        self.WW = None
        self.WTS = None
        self.IP = None
        self.SX = 0.0
        self.SZ = 0.0
        self.CS1 = 0.0
        self.SN1 = 0.0
        self.DXW = 0.0
        self.nChord = 0
        self.NS = 0
        self.nChord = 0
        self.nSpan_half = 0

    def initialize_arrays(self, nChord, nSpan_half, NTMAX):
        # Spanwise arrays
        self.DLY = np.zeros(nSpan_half)
        self.US  = np.zeros(nSpan_half)
        # Chord arrays
        self.SNO = np.zeros(nChord + 1)
        self.CSO = np.zeros(nChord + 1)
        self.GAMA1J = np.zeros(nChord + 1)
        # 
        self.QF = np.zeros((nChord + 1, nSpan_half + 1, 3))
        self.CP = np.zeros((nChord, nSpan_half, 3))
        self.GAMA = np.zeros((nChord, nSpan_half))
        self.DL = np.zeros((nChord, nSpan_half))
        self.DP = np.zeros((nChord, nSpan_half))
        self.DS = np.zeros((nChord, nSpan_half))
        self.DLT = np.zeros((nChord, nSpan_half))
        self.DD = np.zeros((nChord, nSpan_half))
        self.A1 = np.zeros((nChord, nSpan_half))
        self.QW = np.zeros((NTMAX, nSpan_half + 1, 3))
        self.VORTIC = np.zeros((NTMAX, nSpan_half))
        self.UVW = np.zeros((NTMAX, nSpan_half + 1, 3))
        self.QW1 = np.zeros((NTMAX, nSpan_half + 1, 3))
        self.VORT1 = np.zeros((NTMAX, nSpan_half))
        self.WTS = np.zeros((nChord, nSpan_half))
        # System arrays
        self.WW    = np.zeros(nChord * nSpan_half)
        self.A     = np.zeros((nChord * nSpan_half, nChord * nSpan_half))
        self.DW    = np.zeros(nChord * nSpan_half)
        self.GAMA1 = np.zeros(nChord * nSpan_half)
        self.IP    = np.zeros(nChord * nSpan_half, dtype=np.int32)

# Instantiate global data

def vortex(X, Y, Z, X1, Y1, Z1, X2, Y2, Z2, GAMA):
    R1R2X = (Y - Y1) * (Z - Z2) - (Z - Z1) * (Y - Y2)
    R1R2Y = -((X - X1) * (Z - Z2) - (Z - Z1) * (X - X2))
    R1R2Z = (X - X1) * (Y - Y2) - (Y - Y1) * (X - X2)
    SQUARE = R1R2X * R1R2X + R1R2Y * R1R2Y + R1R2Z * R1R2Z
    R1 = np.sqrt((X - X1)**2 + (Y - Y1)**2 + (Z - Z1)**2)
    R2 = np.sqrt((X - X2)**2 + (Y - Y2)**2 + (Z - Z2)**2)
    if R1 < RCUT or R2 < RCUT or SQUARE < RCUT:
        return 0.0, 0.0, 0.0
    R0R1 = (X2 - X1) * (X - X1) + (Y2 - Y1) * (Y - Y1) + (Z2 - Z1) * (Z - Z1)
    R0R2 = (X2 - X1) * (X - X2) + (Y2 - Y1) * (Y - Y2) + (Z2 - Z1) * (Z - Z2)
    COEF = GAMA / (4.0 * PI * SQUARE) * (R0R1 / R1 - R0R2 / R2)
    U = R1R2X * COEF
    V = R1R2Y * COEF
    W = R1R2Z * COEF
    return U, V, W

def wake(X, Y, Z, IT, Gammas, QW, NS):
    U = 0.0
    V = 0.0
    W = 0.0
    I1 = IT - 1
    for I in range(I1):
        for J in range(NS):
            Gamma = Gammas[I, J]
            U1, V1, W1 = vortex(X, Y, Z, QW[I, J, 0], QW[I, J, 1], QW[I, J, 2], QW[I+1, J, 0], QW[I+1, J, 1], QW[I+1, J, 2], Gamma)
            U2, V2, W2 = vortex(X, Y, Z, QW[I+1, J, 0], QW[I+1, J, 1], QW[I+1, J, 2], QW[I+1, J+1, 0], QW[I+1, J+1, 1], QW[I+1, J+1, 2], Gamma)
            U3, V3, W3 = vortex(X, Y, Z, QW[I+1, J+1, 0], QW[I+1, J+1, 1], QW[I+1, J+1, 2], QW[I, J+1, 0], QW[I, J+1, 1], QW[I, J+1, 2], Gamma)
            U4, V4, W4 = vortex(X, Y, Z, QW[I, J+1, 0], QW[I, J+1, 1], QW[I, J+1, 2], QW[I, J, 0], QW[I, J, 1], QW[I, J, 2], Gamma)
            U += U1 + U2 + U3 + U4
            V += V1 + V2 + V3 + V4
            W += W1 + W2 + W3 + W4
    return U, V, W

def veloce(X, Y, Z, IT, GD):
    X1 = (X - GD.SX) * GD.CS1 + (Z - GD.SZ) * GD.SN1
    Y1 = Y
    Z1 = -(X - GD.SX) * GD.SN1 + (Z - GD.SZ) * GD.CS1
    U1, V1, W1 = wake(X,  Y, Z, IT, GD.VORTIC, GD.QW, GD.nSpan_half)
    U2, V2, W2 = wake(X, -Y, Z, IT, GD.VORTIC, GD.QW, GD.nSpan_half)
    U3, V3, W3 = wing(X1,  Y1, Z1, GD.GAMA, GD.QF)
    U4, V4, W4 = wing(X1, -Y1, Z1, GD.GAMA, GD.QF)
    U = U1 + U2 + GD.CS1 * (U3 + U4) - GD.SN1 * (W3 + W4)
    V = V1 - V2 + V3 - V4
    W = W1 + W2 + GD.SN1 * (U3 + U4) + GD.CS1 * (W3 + W4)
    return U, V, W

def wing(X, Y, Z, GAMA, QF, A1=None, SNO=None, CSO=None):
    U = 0.0
    V = 0.0
    W = 0.0
    for I in range(QF.shape[0]-1):
        for J in range(QF.shape[1]-1):
            U1, V1, W1 = vortex(X, Y, Z, QF[I, J, 0], QF[I, J, 1], QF[I, J, 2], QF[I, J+1, 0], QF[I, J+1, 1], QF[I, J+1, 2], GAMA[I, J])
            U2, V2, W2 = vortex(X, Y, Z, QF[I, J+1, 0], QF[I, J+1, 1], QF[I, J+1, 2], QF[I+1, J+1, 0], QF[I+1, J+1, 1], QF[I+1, J+1, 2], GAMA[I, J])
            U3, V3, W3 = vortex(X, Y, Z, QF[I+1, J+1, 0], QF[I+1, J+1, 1], QF[I+1, J+1, 2], QF[I+1, J, 0], QF[I+1, J, 1], QF[I+1, J, 2], GAMA[I, J])
            U4, V4, W4 = vortex(X, Y, Z, QF[I+1, J, 0], QF[I+1, J, 1], QF[I+1, J, 2], QF[I, J, 0], QF[I, J, 1], QF[I, J, 2], GAMA[I, J])
            U0 = U1 + U2 + U3 + U4
            V0 = V1 + V2 + V3 + V4
            W0 = W1 + W2 + W3 + W4
            if A1 is not None:
                A1[I, J] = U0 * SNO[I] + W0 * CSO[I]
            U += U0
            V += V0
            W += W0
    return U, V, W

def wingl(X, Y, Z, GAMA, QF):
    U = 0.0
    V = 0.0
    W = 0.0
    for I in range(QF.shape[0]-1):
        for J in range(QF.shape[1]-1):
            U2, V2, W2 = vortex(X, Y, Z, QF[I, J+1, 0], QF[I, J+1, 1], QF[I, J+1, 2], QF[I+1, J+1, 0], QF[I+1, J+1, 1], QF[I+1, J+1, 2], GAMA[I, J])
            U4, V4, W4 = vortex(X, Y, Z, QF[I+1, J, 0], QF[I+1, J, 1], QF[I+1, J, 2], QF[I, J, 0], QF[I, J, 1], QF[I, J, 2], GAMA[I, J])
            U += U2 + U4
            V += V2 + V4
            W += W2 + W4
    I = QF.shape[0]-2
    for J in range(QF.shape[1]-1):
        U3, V3, W3 = vortex(X, Y, Z, QF[I+1, J+1, 0], QF[I+1, J+1, 1], QF[I+1, J+1, 2], QF[I+1, J, 0], QF[I+1, J, 1], QF[I+1, J, 2], GAMA[I, J])
        U += U3
        V += V3
        W += W3
    return U, V, W

def rectangularWingPanelling(nChord, nSpan,  bSpan, chord, vSpan=None, vChord=None, NW_length=0.5):
    if vSpan is not None:
        chord = vChord[1] - vChord[0]
        bSpan = vSpan[1] - vSpan[0]
    else:
        vSpan = np.linspace(0, bSpan, nSpan+1)
        vChord = np.linspace(0, chord, nChord+1)

    S = bSpan * chord
    AR = 2.0 * bSpan**2 / S

    vChord_QP = vChord+0.25*chord/nChord
    vSpan_CP  = (vSpan[:-1] + vSpan[1:]) / 2
    vChord_CP = (vChord_QP[:-1] + vChord_QP[1:]) / 2

    QF = np.zeros((nChord + 1, nSpan + 1, 3))
    CP = np.zeros((nChord, nSpan, 3))
    DS = np.zeros((nChord, nSpan))
    for J in range(nSpan+1):
        for I in range(nChord+1):
            QF[I, J, 0] = vChord_QP[I]
            QF[I, J, 1] = vSpan[J]
            QF[I, J, 2] = 0
    # Last chord point governed by NW_length
    QF[nChord, :, 0] = chord + NW_length
    QF[nChord, :, 1] = vSpan
    QF[nChord, :, 2] = 0 

    for J in range(nSpan):
        for I in range(nChord):
            CP[I, J, 0] = vChord_CP[I]
            CP[I, J, 1] = vSpan_CP[J]
            CP[I, J, 2] = 0
    # Compute differential surface areas
    for J in range(nSpan):
        for I in range(nChord):
            DS[I, J] = (vChord[I+1]- vChord[I]) * (vSpan[J+1]- vSpan[J])
    return QF, CP, DS, vSpan, vChord, S, AR

def rotateWing(QF, CP, ALFA):
    nChord, nSpan = CP.shape[0], CP.shape[1]
    # Rotate coordinates
    SIN, COS = np.sin(-ALFA), np.cos(-ALFA)
    for I in range(nChord+1):
        for J in range(nSpan+1):
            QF1 = QF[I, J, 0]
            QF[I, J, 0] = QF1 * COS - QF[I, J, 2] * SIN
            QF[I, J, 2] = QF1 * SIN + QF[I, J, 2] * COS
    for I in range(nChord):
        for J in range(nSpan):
            CP1 = CP[I, J, 0]
            CP[I, J, 0] = CP1 * COS - CP[I, J, 2] * SIN
            CP[I, J, 2] = CP1 * SIN + CP[I, J, 2] * COS
    return QF, CP


def computeA(ALFA, GD):
    # --- Compute A
    K = 0
    for I in range(GD.nChord):
        for J in range(GD.nSpan_half):
            U, V, W = wing(GD.CP[I, J, 0], GD.CP[I, J, 1], GD.CP[I, J, 2], GD.GAMA, GD.QF, GD.A1, GD.SNO, GD.CSO)
            L = 0
            for I1 in range(GD.nChord):
                for J1 in range(GD.nSpan_half):
                    GD.A[K, L] = GD.A1[I1, J1]
                    L += 1
            U, V, W = wing(GD.CP[I, J, 0], -GD.CP[I, J, 1], GD.CP[I, J, 2], GD.GAMA, GD.QF, GD.A1, GD.SNO, GD.CSO)
            L = 0
            for I1 in range(GD.nChord):
                for J1 in range(GD.nSpan_half):
                    GD.A[K, L] += GD.A1[I1, J1]
                    L += 1
            K += 1


def main(nChord=4, nSpan_half=13, nStep=10, chord=1, span=8.0, outputDir=''):
    GD = GLOBAL_DATA()
    GD.nChord     = nChord
    GD.nSpan_half = nSpan_half 
    # Initialize arrays
    GD.initialize_arrays(GD.nChord, GD.nSpan_half, PARAMS['NTMAX'])

    # Input data
    RO = 1.0
    BH = 0.0
    OM = 0.0
    VINF = 10.0

    chord = 1.0
    B = span/2.0 # TODO
    DX = chord / GD.nChord
    ALFA1 = 5.0
    ALFAO = 0.0
    ALFA = (ALFA1 + ALFAO) * PI / 180.0
    DT = DX / VINF / 4.0
    T = -DT
    NW_length = 0.3 * VINF * DT

    # Initialize constants
    K = 0
    for I in range(GD.nChord):
        for J in range(GD.nSpan_half):
            GD.WW[K] = 0.0
            GD.DLT[I, J] = 0.0
            GD.VORTIC[I, J] = 0.0
            GD.VORT1[I, J] = 0.0
            GD.GAMA[I, J] = 1.0  # For influence matrix calculations
            K += 1

    # --- Wing geometry
    # Mesh wing and calculate collocation points
    GD.QF, GD.CP, GD.DS, vSpan, vChord, S, AR = rectangularWingPanelling(GD.nChord, GD.nSpan_half, bSpan = B, chord=chord, NW_length=NW_length)
    GD.QF, GD.CP = rotateWing(GD.QF, GD.CP, ALFA)
    GD.dl_LL = np.diff(vSpan)

    if nChord==4 and nSpan_half==13 and span==8:
        np.testing.assert_almost_equal(S, 4.0)
        np.testing.assert_almost_equal(GD.QF[0,:,0], [0.06226217]*14)
        np.testing.assert_almost_equal(GD.QF[1,:,0], [0.3113108]*14)
        np.testing.assert_almost_equal(GD.QF[-1,:,0], [1.0148733]*14) 
        np.testing.assert_almost_equal(GD.QF[:,0,1], [0]*5 )
        np.testing.assert_almost_equal(GD.QF[:,1,1], [0.3076923]*5 )
        np.testing.assert_almost_equal(GD.QF[:,-1,1], [4]*5 )
        np.testing.assert_almost_equal(np.unique(GD.QF[:,:,2]), [-0.08878991, -0.07081404, -0.04902511, -0.02723617, -0.00544723])
        np.testing.assert_almost_equal(GD.CP[0,:,0], [0.1867865]*13)
        np.testing.assert_almost_equal(GD.CP[1,:,0], [0.4358352]*13)
        np.testing.assert_almost_equal(GD.CP[-1,:,0], [0.9339325]*13) 
        np.testing.assert_almost_equal(GD.CP[:,0,1], [0.1538462]*4)
        np.testing.assert_almost_equal(np.unique(GD.CP[:,:,2]), [-0.0817085, -0.0599196, -0.0381306, -0.0163417])
        np.testing.assert_almost_equal(GD.DS.flatten(), [0.0769231]*13*4)


    # --- Compute A
    # TODO, could that change? NOTE: some local inclination was added in original code
    GD.SNO[:GD.nChord] = np.sin(ALFA)
    GD.CSO[:GD.nChord] = np.cos(ALFA)
    computeA(ALFA, GD)

    # Output initial geometry
    print("-" * 56)
    print(f" ALFA: {ALFA1:10.2f}  B : {B:10.2f}  C : {chord:13.2f}")
    print(f" S : {S:10.2f}  AR : {AR:13.2f}")
    print(f" nChord : {GD.nChord:10d}  NS : {GD.nSpan_half:10d} \n")
    NC1 = GD.nChord + 1
    NS1 = GD.nSpan_half + 1

    # Ensure output directory exists
    os.makedirs("_outputs", exist_ok=True)
    with open("_outputs/uLS.csv", "w") as f:
        f.write("#SX,T,SZ,VINF,TETA,OMEGA,CL,L,CM,CD,CL_rel,Gamma_rel\n")




    GD.rWingh= np.zeros_like(GD.QF)
    # Main time-stepping loop
    for IT in range(1, nStep + 1):
        T += DT

        # Path information
        GD.SX = -VINF * T
        DSX = -VINF
        SZ = BH * np.sin(OM * T)
        DSZ = BH * OM * np.cos(OM * T)
        TETA = 0.0
        OMEGA = 0.0
        VINF = -np.cos(TETA) * DSX - np.sin(TETA) * DSZ
        GD.SN1 = np.sin(TETA)
        GD.CS1 = np.cos(TETA)
        WT = GD.SN1 * DSX - GD.CS1 * DSZ
        GD.SNO[:GD.nChord] = np.sin(ALFA) # NOTE: add local inclination
        GD.CSO[:GD.nChord] = np.cos(ALFA)


        GD.rWingh[:,:, 0] = GD.QF[:, :, 0] * GD.CS1 - GD.QF[:, :, 2] * GD.SN1 + GD.SX
        GD.rWingh[:,:, 1] = GD.QF[:, :, 1]
        GD.rWingh[:,:, 2] = GD.QF[:, :, 0] * GD.SN1 + GD.QF[:, :, 2] * GD.CS1 + SZ






        # Vectorized wake shedding points
        GD.QW[IT-1, :, 0] = GD.QF[NC1-1, :, 0] * GD.CS1 - GD.QF[NC1-1, :, 2] * GD.SN1 + GD.SX
        GD.QW[IT-1, :, 1] = GD.QF[NC1-1, :, 1]
        GD.QW[IT-1, :, 2] = GD.QF[NC1-1, :, 0] * GD.SN1 + GD.QF[NC1-1, :, 2] * GD.CS1 + SZ

        # Aerodynamic calculations
        K = 0
        if IT > 1:
            for I in range(GD.nChord):
                for J in range(GD.nSpan_half):
                    W11 = 0.0
                    XX1 = GD.CP[I, J, 0] * GD.CS1 - GD.CP[I, J, 2] * GD.SN1 + GD.SX
                    ZZ1 = GD.CP[I, J, 0] * GD.SN1 + GD.CP[I, J, 2] * GD.CS1 + SZ
                    U, V, W =    wake(XX1,  GD.CP[I, J, 1], ZZ1, IT, GD.VORTIC, GD.QW, GD.nSpan_half)
                    U1, V1, W1 = wake(XX1, -GD.CP[I, J, 1], ZZ1, IT, GD.VORTIC, GD.QW, GD.nSpan_half)
                    U = U + U1
                    W = W + W1
                    U11 = U * GD.CS1 + W * GD.SN1
                    W11 = -U * GD.SN1 + W * GD.CS1
                    GD.WW[K] = U11 * GD.SNO[I] + W11 * GD.CSO[I]
                    GD.DW[K] = -VINF * GD.SNO[I] + GD.CP[I, J, 0] * OMEGA - WT
                    GD.WTS[I, J] = W11
                    K += 1
        else:
            for I in range(GD.nChord):
                for J in range(GD.nSpan_half):
                    W11 = 0.0
                    GD.DW[K] = -VINF * GD.SNO[I] + GD.CP[I, J, 0] * OMEGA - WT
                    GD.WTS[I, J] = W11
                    K += 1

        # Solve the linear system
        K1 = GD.nChord * GD.nSpan_half
        GD.GAMA1[:K1] = GD.DW[:K1] - GD.WW[:K1]
        if IT == 1:
            LU, PIV = lu_factor(GD.A[:K1, :K1])
        GD.GAMA1[:K1] = lu_solve((LU, PIV), GD.GAMA1[:K1])

        # Wing vortex lattice listing
        K = 0
        for I in range(GD.nChord):
            for J in range(GD.nSpan_half):
                GD.GAMA[I, J] = GD.GAMA1[K]
                K += 1

        # Wake shedding
        GD.VORTIC[IT-1, :GD.nSpan_half] = GD.GAMA[GD.nChord-1, :]
        GD.VORTIC[IT, :GD.nSpan_half] = 0.0

        # Wake rollup calculation
        IW = 1
        NWMAX = 5    # Max wake roll up elements
        if IT >= 2:
            if IT >= NWMAX:
                IW = IT - NWMAX + 1
            for I in range(IW-1, IT-1):
                for J in range(NS1):
                    U, V, W = veloce(GD.QW[I, J, 0], GD.QW[I, J, 1], GD.QW[I, J, 2], IT, GD)
                    GD.UVW[I, J, :] = [U * DT, V * DT, W * DT]
            GD.QW[IW-1:IT-1, :NS1, :] += GD.UVW[IW-1:IT-1, :NS1, :]

        # Force calculations
        FL = FD = FM = FG = 0.0
        QUE = 0.5 * RO * VINF * VINF
        for J in range(GD.nSpan_half):
            SIGMA = 0.0
            SIGMA1 = 0.0
            GD.DLY[J] = 0.0
            for I in range(GD.nChord):
                if I == 0:
                    GAMAIJ = GD.GAMA[I, J]
                else:
                    GAMAIJ = GD.GAMA[I, J] - GD.GAMA[I-1, J]
                DXM = (GD.QF[I, J, 0] + GD.QF[I, J+1, 0]) / 2.0
                SIGMA1 = (0.5 * GAMAIJ + SIGMA) * DX
                SIGMA = GD.GAMA[I, J]
                DFDT = (SIGMA1 - GD.DLT[I, J]) / DT
                GD.DLT[I, J] = SIGMA1
                GD.DL[I, J] = RO * (VINF * GAMAIJ + DFDT) * GD.dl_LL[J] * GD.CSO[I]
                U1, V1, W1 = wingl(GD.CP[I, J, 0],  GD.CP[I, J, 1], GD.CP[I, J, 2], GD.GAMA, GD.QF)
                U2, V2, W2 = wingl(GD.CP[I, J, 0], -GD.CP[I, J, 1], GD.CP[I, J, 2], GD.GAMA, GD.QF)
                W8 = W1 + W2
                CTS = -(GD.WTS[I, J] + W8) / VINF
                DD1 = RO * GD.dl_LL[J] * DFDT * GD.SNO[I]
                DD2 = RO * GD.dl_LL[J] * VINF * GAMAIJ * CTS
                GD.DD[I, J] = DD1 + DD2
                GD.DP[I, J] = GD.DL[I, J] / GD.DS[I, J] / QUE
                GD.DLY[J] += GD.DL[I, J]
                FL += GD.DL[I, J]
                FD += GD.DD[I, J]
                FM += GD.DL[I, J] * DXM
                FG += GAMAIJ * GD.dl_LL[J]
        CL = FL / (QUE * S)
        CD = FD / (QUE * S)
        CM = FM / (QUE * S * chord)
        CLOO = 2.0 * PI * ALFA / (1.0 + 2.0 / AR)
        if abs(CLOO) < 1.0e-20:
            CLOO = CL
        CLT = CL / CLOO
        CFG = FG / (0.5 * VINF * S) / CLOO

        # Output results
        #print(f"{IT}/{nStep}")
        #print(f" T={T:10.2f}  SX={GD.SX:10.2f}  SZ={SZ:10.2f}  VINF={VINF:10.2f}  TETA={TETA:10.2f}  OMEGA={OMEGA:10.2f}")
        print(f",CL={CL:10.4f}  L={FL:10.4f}  CM={CM:10.4f}  CD={CD:10.4f}  L/L(INF)={CLT:10.4f}  GAMA/GAMA(INF)={CFG:10.4f}")
        with open("_outputs/uLS.csv", "a") as f:
            f.write(f"{-GD.SX},{T},{SZ},{VINF},{TETA},{OMEGA},{CL},{FL},{CM},{CD},{CLT},{CFG}\n")

        # After updating GD.QW and GD.VORTIC for the current time step IT:
        #write_wing_vtk(IT, GD, outputDir=outputDir)
        #write_wake_vtk(IT, GD, outputDir=outputDir)
    print("[ OK ] Program 16")

def write_wake_vtk(it, GD, outputDir=''):
    """
    Write the current wake lattice and vorticity to a VTK file for visualization.
    """
    mvtk = VTK_Misc()
    filename = os.path.join(outputDir, f"_vtk/_wake_{it:05d}.vtk")
    # GD.QW shape: (NTMAX, NSMAX+1, 3)
    # GD.VORTIC shape: (NTMAX, NSMAX)
    # Reshape to match expected input (nSpan, nDepth, 3)
    QW    = GD.QW[:it, :, :]
    Gamma = GD.VORTIC[:it-1, :]
    QWt   = np.transpose(QW, (1, 0, 2))
    Gammat = np.transpose(Gamma, (1, 0))
    WrVTK_Lattice(filename, mvtk, QWt, Gammat)

def write_wing_vtk(it, GD, outputDir=''):
    """
    Write the current lifting surface lattice and circulation to a VTK file for visualization.
    """
    mvtk     = VTK_Misc()
    filename = os.path.join(outputDir, f"_vtk/_wing_{it:05d}.vtk")

    # GD.QF shape: (nChord+1, NSMAX+1, 3)
    # GD.GAMA shape: (nChord, NSMAX)
    # For time step it, use the current geometry and gamma
    # Reshape QF to (nChord+1, NSMAX+1, 3) -> (nChord+1, NSMAX+1, 3) if not already
    QF    = GD.rWingh[:GD.nChord+1, :GD.nSpan_half+1, :]  # (nChord, nSpan, 3)
    Gamma  = GD.GAMA[:GD.nChord, :GD.nSpan_half]       # (nChord-1, nSpan-1)
    QFt   = np.transpose(QF, (1, 0, 2))    # (nSpan, nDepth, 3)
    Gammat = np.transpose(Gamma, (1, 0))
    WrVTK_Lattice(filename, mvtk, QFt, Gammat)

if __name__ == "__main__":
    scriptDir = os.path.dirname(os.path.abspath(__file__))
    if True:
        with Timer():
            main(nChord=4, nSpan_half=13, nStep=10, chord=1, span=8.0, outputDir=scriptDir)
        # ---
        nt = 10
        df_ref = weio.read(os.path.join(scriptDir, './_outputs/uLS.csv')).toDataFrame()
        df = weio.read(os.path.join(scriptDir, './_outputs/uLS_ref.csv')).toDataFrame()
        np.testing.assert_almost_equal(df['CL'].values[:nt], df_ref['CL'].values[:nt], 6)
        np.testing.assert_almost_equal(df['CD'].values[:nt], df_ref['CD'].values[:nt], 6)
        np.testing.assert_almost_equal(df['CM'].values[:nt], df_ref['CM'].values[:nt], 6)
        np.testing.assert_almost_equal(df['Gamma_rel'].values[:nt], df_ref['Gamma_rel'].values[:nt], 6)
    if False:
        with Timer():
            main(nChord=1, nSpan_half=3, nStep=50, chord=1, span=8.0, outputDir=scriptDir)
        # ---
        nt = 10
        df_ref = weio.read(os.path.join(scriptDir, './_outputs/uLS.csv')).toDataFrame()
        df = weio.read(os.path.join(scriptDir, './_outputs/uLS_ref.csv')).toDataFrame()
        np.testing.assert_almost_equal(df['CL'].values[:nt], df_ref['CL'].values[:nt], 6)
        np.testing.assert_almost_equal(df['CD'].values[:nt], df_ref['CD'].values[:nt], 6)
        np.testing.assert_almost_equal(df['CM'].values[:nt], df_ref['CM'].values[:nt], 6)
        np.testing.assert_almost_equal(df['Gamma_rel'].values[:nt], df_ref['Gamma_rel'].values[:nt], 6)
    print('[ OK ] test pass')
