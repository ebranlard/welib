""" Unsteady lifting surface"""

import numpy as np
from scipy.linalg import lu_factor, lu_solve
import os
from welib.tools.tictoc import Timer
import welib.weio as weio

# Constants
PI = 3.141592654
RCUT = 1.0e-10

# Parameters
PARAMS = {
    'NTMAX': 640,  # Max time steps
    'NCMAX': 4,    # Max chordwise panels
    'NSMAX': 13,   # Max spanwise panels
    'NWMAX': 5     # Max wake elements
}

# Global data
class GLOBAL_DATA:
    def __init__(self):
        self.ALF = None
        self.SNO = None
        self.CSO = None
        self.GAMA1J = None
        self.QF = None  # Wing nodes
        self.CP = None
        self.A = None
        self.DW = None
        self.GAMA1 = None
        self.BB = None
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
        self.NC = 0
        self.NS = 0
        self.IW = 0
        self.CH = 0.0
        self.LU = None
        self.PIV = None

    def initialize_arrays(self):
        self.ALF = np.zeros(PARAMS['NCMAX'] + 1)
        self.SNO = np.zeros(PARAMS['NCMAX'] + 1)
        self.CSO = np.zeros(PARAMS['NCMAX'] + 1)
        self.GAMA1J = np.zeros(PARAMS['NCMAX'] + 1)
        self.QF = np.zeros((PARAMS['NCMAX'] + 1, PARAMS['NSMAX'] + 1, 3))
        self.CP = np.zeros((PARAMS['NCMAX'], PARAMS['NSMAX'], 3))
        self.A = np.zeros((PARAMS['NCMAX'] * PARAMS['NSMAX'], PARAMS['NCMAX'] * PARAMS['NSMAX']))
        self.DW = np.zeros(PARAMS['NCMAX'] * PARAMS['NSMAX'])
        self.GAMA1 = np.zeros(PARAMS['NCMAX'] * PARAMS['NSMAX'])
        self.BB = np.zeros(PARAMS['NSMAX'])
        self.DLY = np.zeros(PARAMS['NSMAX'])
        self.GAMA = np.zeros((PARAMS['NCMAX'], PARAMS['NSMAX']))
        self.DL = np.zeros((PARAMS['NCMAX'], PARAMS['NSMAX']))
        self.DP = np.zeros((PARAMS['NCMAX'], PARAMS['NSMAX']))
        self.DS = np.zeros((PARAMS['NCMAX'], PARAMS['NSMAX']))
        self.DLT = np.zeros((PARAMS['NCMAX'], PARAMS['NSMAX']))
        self.DD = np.zeros((PARAMS['NCMAX'], PARAMS['NSMAX']))
        self.A1 = np.zeros((PARAMS['NCMAX'], PARAMS['NSMAX']))
        self.QW = np.zeros((PARAMS['NTMAX'], PARAMS['NSMAX'] + 1, 3))
        self.VORTIC = np.zeros((PARAMS['NTMAX'], PARAMS['NSMAX']))
        self.UVW = np.zeros((PARAMS['NTMAX'], PARAMS['NSMAX'] + 1, 3))
        self.QW1 = np.zeros((PARAMS['NTMAX'], PARAMS['NSMAX'] + 1, 3))
        self.VORT1 = np.zeros((PARAMS['NTMAX'], PARAMS['NSMAX']))
        self.US = np.zeros(PARAMS['NSMAX'])
        self.WW = np.zeros(PARAMS['NCMAX'] * PARAMS['NSMAX'])
        self.WTS = np.zeros((PARAMS['NCMAX'], PARAMS['NSMAX']))
        self.IP = np.zeros(PARAMS['NCMAX'] * PARAMS['NSMAX'], dtype=np.int32)

# Instantiate global data
GD = GLOBAL_DATA()

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

def veloce(X, Y, Z, IT, JS1, JS2):
    X1 = (X - GD.SX) * GD.CS1 + (Z - GD.SZ) * GD.SN1
    Y1 = Y
    Z1 = -(X - GD.SX) * GD.SN1 + (Z - GD.SZ) * GD.CS1
    U1, V1, W1 = wake(X,  Y, Z, IT, GD.VORTIC, GD.QW, GD.NS)
    U2, V2, W2 = wake(X, -Y, Z, IT, GD.VORTIC, GD.QW, GD.NS)
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

def geo(B, C, NC, NS, DX, DY, DGAP, ALFA):
    NC1 = NC + 1
    NS1 = NS + 1
    SN = np.sin(GD.ALF[:NC1])
    CS = np.cos(GD.ALF[:NC1])
    CTG1 = np.tan(PI / 2.0 - 90.0 * PI / 180.0)
    CTG2 = np.tan(PI / 2.0 - 90.0 * PI / 180.0)
    CTIP = C + B * (CTG2 - CTG1)
    S = B * (C + CTIP) / 2.0
    AR = 2.0 * B * B / S

    BJ = 0.0
    for J in range(NS1):
        if J > 0:
            BJ += GD.BB[J-1]
        Z1 = 0.0
        DC1 = BJ * CTG1
        DC2 = BJ * CTG2
        DX1 = (C + DC2 - DC1) / NC
        for I in range(NC):
            GD.QF[I, J, 0] = DC1 + DX1 * (I+1 - 0.75)
            GD.QF[I, J, 1] = BJ
            GD.QF[I, J, 2] = Z1 - 0.25 * DX1 * SN[I]
            Z1 -= DX1 * SN[I]
        GD.QF[NC, J, 0] = C + DC2 + GD.DXW
        GD.QF[NC, J, 1] = GD.QF[NC-1, J, 1]
        GD.QF[NC, J, 2] = Z1 - GD.DXW * SN[NC-1]

    for J in range(NS):
        Z1 = 0.0
        BJ = GD.QF[0, J, 1] + GD.BB[J] / 2.0
        DC1 = BJ * CTG1
        DC2 = BJ * CTG2
        DX1 = (C + DC2 - DC1) / NC
        for I in range(NC):
            GD.CP[I, J, 0] = DC1 + DX1 * (I+1 - 0.25)
            GD.CP[I, J, 1] = BJ
            GD.CP[I, J, 2] = Z1 - 0.75 * DX1 * SN[I]
            Z1 -= DX1 * SN[I]
            GD.DS[I, J] = DX1 * GD.BB[J]

    # Rotate coordinates
    SN1 = np.sin(-ALFA)
    CS1 = np.cos(-ALFA)
    for I in range(NC1):
        for J in range(NS1):
            QF1 = GD.QF[I, J, 0]
            GD.QF[I, J, 0] = QF1 * CS1 - GD.QF[I, J, 2] * SN1
            GD.QF[I, J, 2] = QF1 * SN1 + GD.QF[I, J, 2] * CS1
            if I == NC or J >= NS:
                continue
            CP1 = GD.CP[I, J, 0]
            GD.CP[I, J, 0] = CP1 * CS1 - GD.CP[I, J, 2] * SN1
            GD.CP[I, J, 2] = CP1 * SN1 + GD.CP[I, J, 2] * CS1
    return S, AR


def computeA(ALFA):
    # --- Compute A
    K = 0
    for I in range(GD.NC):
        for J in range(GD.NS):
            U, V, W = wing(GD.CP[I, J, 0], GD.CP[I, J, 1], GD.CP[I, J, 2], GD.GAMA, GD.QF, GD.A1, GD.SNO, GD.CSO)
            L = 0
            for I1 in range(GD.NC):
                for J1 in range(GD.NS):
                    GD.A[K, L] = GD.A1[I1, J1]
                    L += 1
            U, V, W = wing(GD.CP[I, J, 0], -GD.CP[I, J, 1], GD.CP[I, J, 2], GD.GAMA, GD.QF, GD.A1, GD.SNO, GD.CSO)
            L = 0
            for I1 in range(GD.NC):
                for J1 in range(GD.NS):
                    GD.A[K, L] += GD.A1[I1, J1]
                    L += 1
            K += 1


def main():
    # Initialize arrays
    GD.initialize_arrays()

    # Input data
    GD.NC = PARAMS['NCMAX']
    GD.NS = PARAMS['NSMAX']
    NSTEPS = 40  # For DT / 4
    RO = 1.0
    BH = 0.0
    OM = 0.0
    VINF = 10.0
    C = 1.0
    B = 4.0  # AR = 8
    DX = C / GD.NC
    DY = B / GD.NS
    GD.CH = 10000.0 * C
    ALFA1 = 5.0
    ALFAO = 0.0
    ALFA = (ALFA1 + ALFAO) * PI / 180.0
    GD.ALF[:GD.NC] = 0.0
    GD.ALF[GD.NC] = GD.ALF[GD.NC-1]
    DT = DX / VINF / 4.0
    print(f"DT={DT}")
    T = -DT
    GD.DXW = 0.3 * VINF * DT
    GD.BB[:] = DY

    # Initialize constants
    K = 0
    for I in range(GD.NC):
        for J in range(GD.NS):
            GD.WW[K] = 0.0
            GD.DLT[I, J] = 0.0
            GD.VORTIC[I, J] = 0.0
            GD.VORT1[I, J] = 0.0
            GD.GAMA[I, J] = 1.0  # For influence matrix calculations
            K += 1

    # Calculate collocation points
    S, AR = geo(B, C, GD.NC, GD.NS, DX, DY, 0.0, ALFA)
    # --- Compute A
    # TODO, could that change?
    GD.SNO[:GD.NC] = np.sin(ALFA + GD.ALF[:GD.NC])
    GD.CSO[:GD.NC] = np.cos(ALFA + GD.ALF[:GD.NC])
    computeA(ALFA)

    # Output initial geometry
    print("\nWING LIFT DISTRIBUTION CALCULATION (WITH GROUND EFFECT)")
    print("-" * 56)
    print(f" ALFA: {ALFA1:10.2f}  B : {B:10.2f}  C : {C:13.2f}")
    print(f" S : {S:10.2f}  AR : {AR:13.2f}")
    print(f" NC : {GD.NC:10d}  NS : {GD.NS:10d}  L.E. HEIGHT: {GD.CH:6.2f}\n")
    for I in range(GD.NC):
        print(f" ALF({I+1:2d})={GD.ALF[I] * 180.0 / PI:10.4f}")
    for I in range(0, GD.NS, 2):
        print(f" BB({I+1:3d})={GD.BB[I]:10.4f}")
    NC1 = GD.NC + 1
    NS1 = GD.NS + 1

    # Ensure output directory exists
    os.makedirs("_outputs", exist_ok=True)
    with open("_outputs/uLS.csv", "w") as f:
        f.write("#SX,T,SZ,VINF,TETA,OMEGA,CL,L,CM,CD,CL_rel,Gamma_rel\n")




    # Main time-stepping loop
    for IT in range(1, NSTEPS + 1):
        print(f"{IT}/{NSTEPS}")
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
        GD.SNO[:GD.NC] = np.sin(ALFA + GD.ALF[:GD.NC])
        GD.CSO[:GD.NC] = np.cos(ALFA + GD.ALF[:GD.NC])

        # Vectorized wake shedding points
        GD.QW[IT-1, :, 0] = GD.QF[NC1-1, :, 0] * GD.CS1 - GD.QF[NC1-1, :, 2] * GD.SN1 + GD.SX
        GD.QW[IT-1, :, 1] = GD.QF[NC1-1, :, 1]
        GD.QW[IT-1, :, 2] = GD.QF[NC1-1, :, 0] * GD.SN1 + GD.QF[NC1-1, :, 2] * GD.CS1 + SZ

        # Aerodynamic calculations
        K = 0
        if IT > 1:
            for I in range(GD.NC):
                for J in range(GD.NS):
                    W11 = 0.0
                    XX1 = GD.CP[I, J, 0] * GD.CS1 - GD.CP[I, J, 2] * GD.SN1 + GD.SX
                    ZZ1 = GD.CP[I, J, 0] * GD.SN1 + GD.CP[I, J, 2] * GD.CS1 + SZ
                    U, V, W =    wake(XX1,  GD.CP[I, J, 1], ZZ1, IT, GD.VORTIC, GD.QW, GD.NS)
                    U1, V1, W1 = wake(XX1, -GD.CP[I, J, 1], ZZ1, IT, GD.VORTIC, GD.QW, GD.NS)
                    U = U + U1
                    W = W + W1
                    U11 = U * GD.CS1 + W * GD.SN1
                    W11 = -U * GD.SN1 + W * GD.CS1
                    GD.WW[K] = U11 * GD.SNO[I] + W11 * GD.CSO[I]
                    GD.DW[K] = -VINF * GD.SNO[I] + GD.CP[I, J, 0] * OMEGA - WT
                    GD.WTS[I, J] = W11
                    K += 1
        else:
            for I in range(GD.NC):
                for J in range(GD.NS):
                    W11 = 0.0
                    GD.DW[K] = -VINF * GD.SNO[I] + GD.CP[I, J, 0] * OMEGA - WT
                    GD.WTS[I, J] = W11
                    K += 1

        # Solve the linear system
        K1 = GD.NC * GD.NS
        GD.GAMA1[:K1] = GD.DW[:K1] - GD.WW[:K1]
        if IT == 1:
            GD.LU, GD.PIV = lu_factor(GD.A[:K1, :K1])
        GD.GAMA1[:K1] = lu_solve((GD.LU, GD.PIV), GD.GAMA1[:K1])

        # Wing vortex lattice listing
        K = 0
        for I in range(GD.NC):
            for J in range(GD.NS):
                GD.GAMA[I, J] = GD.GAMA1[K]
                K += 1

        # Wake shedding
        GD.VORTIC[IT-1, :GD.NS] = GD.GAMA[GD.NC-1, :]
        GD.VORTIC[IT, :GD.NS] = 0.0

        # Wake rollup calculation
        GD.IW = 1
        if IT >= 2:
            if IT >= PARAMS['NWMAX']:
                GD.IW = IT - PARAMS['NWMAX'] + 1
            for I in range(GD.IW-1, IT-1):
                for J in range(NS1):
                    U, V, W = veloce(GD.QW[I, J, 0], GD.QW[I, J, 1], GD.QW[I, J, 2], IT, 0, 0)
                    GD.UVW[I, J, :] = [U * DT, V * DT, W * DT]
            GD.QW[GD.IW-1:IT-1, :NS1, :] += GD.UVW[GD.IW-1:IT-1, :NS1, :]

        # Force calculations
        FL = FD = FM = FG = 0.0
        QUE = 0.5 * RO * VINF * VINF
        for J in range(GD.NS):
            SIGMA = 0.0
            SIGMA1 = 0.0
            GD.DLY[J] = 0.0
            for I in range(GD.NC):
                if I == 0:
                    GAMAIJ = GD.GAMA[I, J]
                else:
                    GAMAIJ = GD.GAMA[I, J] - GD.GAMA[I-1, J]
                DXM = (GD.QF[I, J, 0] + GD.QF[I, J+1, 0]) / 2.0
                SIGMA1 = (0.5 * GAMAIJ + SIGMA) * DX
                SIGMA = GD.GAMA[I, J]
                DFDT = (SIGMA1 - GD.DLT[I, J]) / DT
                GD.DLT[I, J] = SIGMA1
                GD.DL[I, J] = RO * (VINF * GAMAIJ + DFDT) * GD.BB[J] * GD.CSO[I]
                U1, V1, W1 = wingl(GD.CP[I, J, 0],  GD.CP[I, J, 1], GD.CP[I, J, 2], GD.GAMA, GD.QF)
                U2, V2, W2 = wingl(GD.CP[I, J, 0], -GD.CP[I, J, 1], GD.CP[I, J, 2], GD.GAMA, GD.QF)
                W8 = W1 + W2
                CTS = -(GD.WTS[I, J] + W8) / VINF
                DD1 = RO * GD.BB[J] * DFDT * GD.SNO[I]
                DD2 = RO * GD.BB[J] * VINF * GAMAIJ * CTS
                GD.DD[I, J] = DD1 + DD2
                GD.DP[I, J] = GD.DL[I, J] / GD.DS[I, J] / QUE
                GD.DLY[J] += GD.DL[I, J]
                FL += GD.DL[I, J]
                FD += GD.DD[I, J]
                FM += GD.DL[I, J] * DXM
                FG += GAMAIJ * GD.BB[J]
        CL = FL / (QUE * S)
        CD = FD / (QUE * S)
        CM = FM / (QUE * S * C)
        CLOO = 2.0 * PI * ALFA / (1.0 + 2.0 / AR)
        if abs(CLOO) < 1.0e-20:
            CLOO = CL
        CLT = CL / CLOO
        CFG = FG / (0.5 * VINF * S) / CLOO

        # Output results
        print(f" T={T:10.2f}  SX={GD.SX:10.2f}  SZ={SZ:10.2f}  VINF={VINF:10.2f}  TETA={TETA:10.2f}  OMEGA={OMEGA:10.2f}")
        print(f",CL={CL:10.4f}  L={FL:10.4f}  CM={CM:10.4f}  CD={CD:10.4f}  L/L(INF)={CLT:10.4f}  GAMA/GAMA(INF)={CFG:10.4f}")
        with open("_outputs/uLS.csv", "a") as f:
            f.write(f"{-GD.SX},{T},{SZ},{VINF},{TETA},{OMEGA},{CL},{FL},{CM},{CD},{CLT},{CFG}\n")

        I2 = 5
        if IT == I2:
            print("=" * 118)
            print(" I     DL    II          DCP          I I          GAMA")
            print(" I            I= 1     2     3     4     I I     1     2     3     4")
            print("=" * 118)
            for J in range(GD.NS):
                for I in range(1, GD.NC):
                    GD.GAMA1J[I+1] = GD.GAMA[I, J] - GD.GAMA[I-1, J]
                DLYJ = GD.DLY[J] / GD.BB[J]
                print(f"{J+1:3d} I{DLYJ:9.3f} II{GD.DP[0,J]:9.3f} I{GD.DP[1,J]:9.3f} I{GD.DP[2,J]:9.3f} I{GD.DP[3,J]:9.3f} II{GD.GAMA[0,J]:9.3f} I{GD.GAMA1J[2]:9.3f} I{GD.GAMA1J[3]:9.3f} I{GD.GAMA1J[4]:9.3f} I")
            print(" WAKE ELEMENTS")
            for I in range(IT):
                print(f" VORTIC(IT={I+1:3d})={GD.VORTIC[I,:PARAMS['NSMAX']]}")
                for J in range(3):
                    print(f" QW({J+1:2d})={GD.QW[I,:PARAMS['NSMAX']+1,J]}")
    print("[ OK ] Program 16")

if __name__ == "__main__":
    with Timer():
        main()
    # ---
    df_ref = weio.read('./_outputs/uLS.csv').toDataFrame()
    df = weio.read('./_outputs/uLS_ref.csv').toDataFrame()
    np.testing.assert_almost_equal(df['CL'].values, df_ref['CL'].values, 6)
    np.testing.assert_almost_equal(df['CD'].values, df_ref['CD'].values, 6)
    np.testing.assert_almost_equal(df['CM'].values, df_ref['CM'].values, 6)
    np.testing.assert_almost_equal(df['Gamma_rel'].values, df_ref['Gamma_rel'].values, 6)
    print('[ OK ] test pass')
