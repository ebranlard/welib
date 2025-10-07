""" Unsteady lifting surface"""

import numpy as np
from scipy.linalg import lu_factor, lu_solve
import os
from welib.tools.tictoc import Timer
import welib.weio as weio
from vtk import VTK_Misc, WrVTK_Lattice
#from welib.vortilib.elements.VortexSegment import vs_u_raw

# Constants
PI = 3.141592654
RCUT = 1.0e-10

# Global data
class GLOBAL_DATA:
    def __init__(self):
        self.SNO = None
        self.CSO = None
         # Wing nodes
        self.CP = None # Control points
        self.AA = None
        self.DW = None
        self.Gamma_LS = None # Wing Gammas
        self.Gamma_LS_prev = None # Wing Gammas
        self.Area = None
        self.dLT = None
        self.Uind_wake_on_wake = None
        self.WW = None
        self.Uind_wake_on_LS = None
        self.nChord = 0
        self.nSpan = 0
        # --- Wake
        self.Gamma_NW = None # Wake vorticity Gammas
        self.r_NW = None     # Wake nodes (was QW)

    def initialize_arrays(self, nChord, nSpan, NTMAX):
        self.nChord = nChord
        self.nSpan  = nSpan 
        # Spanwise arrays
        # Chord arrays
        self.SNO = np.zeros(nChord + 1)
        self.CSO = np.zeros(nChord + 1)
        # Lifting surface
        self.QF              = np.zeros((nChord + 1, nSpan + 1, 3))
        self.CP              = np.zeros((nChord, nSpan, 3))
        self.Gamma_LS        = np.zeros((nChord, nSpan))
        self.Gamma_LS_prev   = np.zeros((nChord, nSpan))
        self.Area            = np.zeros((nChord, nSpan))
        self.dLT             = np.zeros((nChord, nSpan))
        self.Uind_wake_on_LS = np.zeros((nChord, nSpan, 3))
        # Wake
        self.r_NW     = np.zeros((NTMAX, nSpan + 1, 3))
        self.Gamma_NW = np.zeros((NTMAX, nSpan))
        self.Uind_wake_on_wake = np.zeros((NTMAX, nSpan + 1, 3))
        # System arrays
        self.WW    = np.zeros(nChord * nSpan)
        self.AA    = np.zeros((nChord * nSpan, nChord * nSpan))
        self.DW    = np.zeros(nChord * nSpan)
        
def vs_u_raw(CP, P1, P2, GAMA, **kwargs):
    X, Y, Z = CP
    X1, Y1, Z1 = P1
    X2, Y2, Z2 = P2
    R1R2X = (Y - Y1) * (Z - Z2) - (Z - Z1) * (Y - Y2)
    R1R2Y = -((X - X1) * (Z - Z2) - (Z - Z1) * (X - X2))
    R1R2Z = (X - X1) * (Y - Y2) - (Y - Y1) * (X - X2)
    SQUARE = R1R2X * R1R2X + R1R2Y * R1R2Y + R1R2Z * R1R2Z
    R1 = np.sqrt((X - X1)**2 + (Y - Y1)**2 + (Z - Z1)**2)
    R2 = np.sqrt((X - X2)**2 + (Y - Y2)**2 + (Z - Z2)**2)
    if R1 < RCUT or R2 < RCUT or SQUARE < RCUT:
        return np.array([0.0, 0.0, 0.0])
    R0R1 = (X2 - X1) * (X - X1) + (Y2 - Y1) * (Y - Y1) + (Z2 - Z1) * (Z - Z1)
    R0R2 = (X2 - X1) * (X - X2) + (Y2 - Y1) * (Y - Y2) + (Z2 - Z1) * (Z - Z2)
    COEF = GAMA / (4.0 * PI * SQUARE) * (R0R1 / R1 - R0R2 / R2)
    U = R1R2X * COEF
    V = R1R2Y * COEF
    W = R1R2Z * COEF
    return np.array([U, V, W])

def UI_from_wake(CP, IT, r_NW, Gammas_NW, flip=False):
    if flip:
        CP = CP.copy()
        CP[1] *= -1
    U = np.zeros(3)
    for I in range(IT):
        for J in range(Gammas_NW.shape[1]):
            Gamma = Gammas_NW[I, J]
            U1 = vs_u_raw(CP, r_NW[I  , J   , :], r_NW[I+1, J  , :], Gamma, rcut_den=RCUT, rcut_norm=RCUT)
            U2 = vs_u_raw(CP, r_NW[I+1, J   , :], r_NW[I+1, J+1, :], Gamma, rcut_den=RCUT, rcut_norm=RCUT)
            U3 = vs_u_raw(CP, r_NW[I+1, J+1 , :], r_NW[I  , J+1, :], Gamma, rcut_den=RCUT, rcut_norm=RCUT)
            U4 = vs_u_raw(CP, r_NW[I  , J+1 , :], r_NW[I  , J  , :], Gamma, rcut_den=RCUT, rcut_norm=RCUT)
            U += U1 + U2 + U3 + U4
    if flip:
        U[1] *= -1
    return U

def UI_from_all(CP, IT, r_NW, Gamma_NW, r_LS, Gamma_LS):
    U1 = UI_from_wake(CP , IT, r_NW, Gamma_NW)
    U2 = UI_from_wake(CP, IT, r_NW, Gamma_NW, flip=True)
    U3 = UI_from_wing(CP, Gamma_LS, r_LS)
    U4 = UI_from_wing(CP, Gamma_LS, r_LS, flip=True)
    U = U1 + U2 + U3 + U4
    return U

def UI_mat_from_wing(CP, Gamma_LS, r_LS, n_hat, flip=False):
    if flip:
        CP = CP.copy()
        CP[1] *= -1

    nChord, nSpan = Gamma_LS.shape[0:2]
    Un = np.zeros((nChord, nSpan))

    for I in range(r_LS.shape[0]-1):
        for J in range(r_LS.shape[1]-1):
            U1 = vs_u_raw(CP, r_LS[I  , J  , :], r_LS[I  , J+1, :], Gamma_LS[I, J], rcut_den=RCUT, rcut_norm=RCUT)
            U2 = vs_u_raw(CP, r_LS[I  , J+1, :], r_LS[I+1, J+1, :], Gamma_LS[I, J], rcut_den=RCUT, rcut_norm=RCUT)
            U3 = vs_u_raw(CP, r_LS[I+1, J+1, :], r_LS[I+1, J  , :], Gamma_LS[I, J], rcut_den=RCUT, rcut_norm=RCUT)
            U4 = vs_u_raw(CP, r_LS[I+1, J  , :], r_LS[I  , J  , :], Gamma_LS[I, J], rcut_den=RCUT, rcut_norm=RCUT)
            U = U1 + U2 + U3 + U4
            Un[I, J] = U[0] * n_hat[I, J, 0] + U[2] * n_hat[I, J, 2]
    return Un


def UI_from_wing(CP, Gamma_LS, r_LS, flip=False):
    if flip:
        CP = CP.copy()
        CP[1] *= -1
    U = np.zeros(3)
    for I in range(r_LS.shape[0]-1):
        for J in range(r_LS.shape[1]-1):
            U1 = vs_u_raw(CP, r_LS[I  , J  , :], r_LS[I  , J+1, :], Gamma_LS[I, J], rcut_den=RCUT, rcut_norm=RCUT)
            U2 = vs_u_raw(CP, r_LS[I  , J+1, :], r_LS[I+1, J+1, :], Gamma_LS[I, J], rcut_den=RCUT, rcut_norm=RCUT)
            U3 = vs_u_raw(CP, r_LS[I+1, J+1, :], r_LS[I+1, J  , :], Gamma_LS[I, J], rcut_den=RCUT, rcut_norm=RCUT)
            U4 = vs_u_raw(CP, r_LS[I+1, J  , :], r_LS[I  , J  , :], Gamma_LS[I, J], rcut_den=RCUT, rcut_norm=RCUT)
            U0 = U1 + U2 + U3 + U4
            U += U0
    if flip:
        U[1] *= -1
    return U

def UI_from_wing_trailed(CP, Gamma, r_LS, flip=False):
    if flip:
        CP = CP.copy()
        CP[1] *= -1
    U = np.zeros(3)
    for I in range(r_LS.shape[0]-1):
        for J in range(r_LS.shape[1]-1):
            U2 = vs_u_raw(CP, r_LS[I  , J+1, :], r_LS[I+1, J+1, :], Gamma[I, J], rcut_den=RCUT, rcut_norm=RCUT)
            U4 = vs_u_raw(CP, r_LS[I+1, J  , :], r_LS[I  , J  , :], Gamma[I, J], rcut_den=RCUT, rcut_norm=RCUT)
            #print(f"WG{I+1:5d}{J+1:5d}{Gamma[I, J]:15.6f}")
            #print(f"WQ{I+1:5d}{J+1:5d}{r_LS[I, J+1, 0]:15.6f}{r_LS[I, J+1, 1]:15.6f}{r_LS[I, J+1, 2]:15.6f}")
            #print(f"WQ{I+1:5d}{J+1:5d}{r_LS[I+1, J, 0]:15.6f}{r_LS[I+1, J, 1]:15.6f}{r_LS[I+1, J, 2]:15.6f}")
            #print(f"WQ{I+1:5d}{J+1:5d}{r_LS[I+1, J+1, 0]:15.6f}{r_LS[I+1, J+1, 1]:15.6f}{r_LS[I+1, J+1, 2]:15.6f}")
            #print(f"WX{I+1:5d}{J+1:5d}{CP[0]:15.6f}{CP[1]:15.6f}{CP[2]:15.6f}")
            #print(f"WU{I+1:5d}{J+1:5d}{U2[0]:15.6f} {U2[1]:15.6f} {U2[2]:15.6f}")
            U += U2 + U4
    I = r_LS.shape[0]-2
    for J in range(r_LS.shape[1]-1):
        U3 = vs_u_raw(CP, r_LS[I+1, J+1, :], r_LS[I+1, J, :], Gamma[I, J], rcut_den=RCUT, rcut_norm=RCUT)
        U += U3
    if flip:
        U[1] *= -1
#   subroutine WINGL(X, Y, Z, GAMA, U, V, W)
#     ! Calculate induced velocity for induced drag calculation
#     U = 0.0_MK
#     V = 0.0_MK
#     W = 0.0_MK
#     do I = 1, NC
#       do J = 1, NS
#         call VORTEX(X, Y, Z, QF(I,J+1,1), QF(I,J+1,2), QF(I,J+1,3), QF(I+1,J+1,1), QF(I+1,J+1,2), QF(I+1,J+1,3), GAMA(I,J), U2, V2, W2)
#         call VORTEX(X, Y, Z, QF(I+1,J,1), QF(I+1,J,2), QF(I+1,J,3), QF(I,J,1), QF(I,J,2), QF(I,J,3), GAMA(I,J), U4, V4, W4)
#         U = U + U2 + U4
#         V = V + V2 + V4
#         W = W + W2 + W4
#       end do
#     end do
#     ! Add influence of latest unsteady wake element
#     I = NC
#     do J = 1, NS
#       call VORTEX(X, Y, Z, QF(I+1,J+1,1), QF(I+1,J+1,2), QF(I+1,J+1,3), QF(I+1,J,1), QF(I+1,J,2), QF(I+1,J,3), GAMA(I,J), U3, V3, W3)
#       U = U + U3
#       V = V + V3
#       W = W + W3
#     end do
#   end subroutine WINGL


    return np.array(U)

def rectangularWingPanelling(nChord, nSpan,  bSpan, chord, vSpan=None, vChord=None, NW_length=0.5):
    if vSpan is not None:
        chord = vChord[1] - vChord[0]
        bSpan = vSpan[1] - vSpan[0]
    else:
        vSpan = np.linspace(0, bSpan, nSpan+1)
        vChord = np.linspace(0, chord, nChord+1)

    vChord_QP = vChord+0.25*chord/nChord
    vSpan_CP  = (vSpan[:-1] + vSpan[1:]) / 2
    vChord_CP = (vChord_QP[:-1] + vChord_QP[1:]) / 2

    r_LS = np.zeros((nChord + 1, nSpan + 1, 3))
    CP = np.zeros((nChord, nSpan, 3))
    Area = np.zeros((nChord, nSpan))
    for J in range(nSpan+1):
        for I in range(nChord+1):
            r_LS[I, J, 0] = vChord_QP[I]
            r_LS[I, J, 1] = vSpan[J]
            r_LS[I, J, 2] = 0
    # Last chord point governed by NW_length
    r_LS[nChord, :, 0] = chord + NW_length
    r_LS[nChord, :, 1] = vSpan
    r_LS[nChord, :, 2] = 0 

    for J in range(nSpan):
        for I in range(nChord):
            CP[I, J, 0] = vChord_CP[I]
            CP[I, J, 1] = vSpan_CP[J]
            CP[I, J, 2] = 0
    # Compute differential surface areas
    for J in range(nSpan):
        for I in range(nChord):
            Area[I, J] = (vChord[I+1]- vChord[I]) * (vSpan[J+1]- vSpan[J])
    return r_LS, CP, Area, vSpan, vChord

def rotAndTrans_abs(P0, theta, translation=(0,0,0)):
    """
    Positive Rotation about y
    ------------------------- 
    [X_rot] =  [ C  S ] [X]
    [Z_rot] =  [-S  C ] [Z]
    Change of coordinate from local to inertial:
    [X]  =  [ C  S ] [x] 
    [Z]i =  [-S  C ] [z]b
    """
    P = np.zeros_like(P0)
    CS = np.cos(theta)
    SN = np.sin(theta)
    P[:,:, 0] =  P0[:, :, 0] * CS + P0[:, :, 2] * SN + translation[0]
    P[:,:, 1] =  P0[:, :, 1]                         + translation[1]
    P[:,:, 2] = -P0[:, :, 0] * SN + P0[:, :, 2] * CS + translation[2]
    return P

def rotAndTrans_abs_opp(P0, theta, translation=(0,0,0)):
    # NOTE: I disaggree with this sign convention
    P = np.zeros_like(P0)
    CS = np.cos(theta)
    SN = np.sin(theta)
    P[:,:, 0] = P0[:, :, 0] * CS - P0[:, :, 2] * SN + translation[0]
    P[:,:, 1] = P0[:, :, 1]                         + translation[1]
    P[:,:, 2] = P0[:, :, 0] * SN + P0[:, :, 2] * CS + translation[2]
    return P

def inertial_to_body_theta(Vi, theta):
    Vb = np.zeros_like(Vi)
    CS = np.cos(theta)
    SN = np.sin(theta)
    Vb[0] = Vi[0] * CS - Vi[2] * SN
    Vb[1] = Vi[1]
    Vb[2] = Vi[0] * SN + Vi[2] * CS
    return Vb

def computeA(CP, r_LS, GD, n_hat):
    # --- Compute A
    K = 0
    nChord, nSpan = CP.shape[0:2]
    Gamma_LS = np.ones((nChord, nSpan)) # Unit circulation
    A  = np.zeros((nChord * nSpan, nChord * nSpan))
    for I in range(nChord):
        for J in range(nSpan):
            A1 = UI_mat_from_wing(CP[I, J, :], Gamma_LS, r_LS, n_hat)
            L = 0
            for I1 in range(nChord):
                for J1 in range(nSpan):
                    A[K, L] = A1[I1, J1]
                    L += 1
            A1 = UI_mat_from_wing(CP[I,J, :], Gamma_LS, r_LS, n_hat, flip=True)
            L = 0
            for I1 in range(nChord):
                for J1 in range(nSpan):
                    A[K, L] += A1[I1, J1]
                    L += 1
            K += 1
    return A


def panlInfo(r, spanFirst=True):
    if not spanFirst:
        r   = np.transpose(r, (1, 0, 2))  # (nSpan, nDepth, 3)
    nSpanP1, nChordP1, _ = r.shape
    nSpan = nSpanP1- 1
    nChord = nChordP1- 1
    Tang = np.zeros((nSpan, nChord, 3))
    Norm = np.zeros((nSpan, nChord, 3))
    Orth = np.zeros((nSpan, nChord, 3))
    dl   = np.zeros((nSpan, nChord, 3))
    Area = np.zeros((nSpan, nChord))
    for iChord in range(nChord):
        for iSpan in range(nSpan):
            P1     = r[iSpan  , iChord  , :]
            P2     = r[iSpan  , iChord+1, :]
            P3     = r[iSpan+1, iChord+1, :]
            P4     = r[iSpan+1, iChord  , :]
            P8     = (P1+P4)/2
            P6     = (P2+P3)/2
            P5     = (P1+P2)/2
            P7     = (P4+P3)/2
            P9     = 0.75*P1+0.25*P2
            P10    = 0.75*P4+0.25*P3
            DP1    = P6-P8
            DP2    = P10-P9
            DP3    = P7-P5
            Tang[iSpan, iChord] = (DP1)/np.linalg.norm(DP1) # tangential unit vector, along chord
            Norm[iSpan, iChord] = np.cross(DP1,DP2)
            Norm[iSpan, iChord] /= np.linalg.norm(Norm[iSpan, iChord])
            dl  [iSpan, iChord]  = DP2
            Orth[iSpan, iChord]  = np.cross(Norm[iSpan, iChord], Tang[iSpan, iChord]) # orthogonal vector to N and T
            Area[iSpan] = np.linalg.norm(np.cross(DP1,DP3))
    if not spanFirst:
        Norm = np.transpose(Norm, (1, 0, 2))
        Tang = np.transpose(Tang, (1, 0, 2))
        Orth = np.transpose(Orth, (1, 0, 2))
        dl   = np.transpose(dl, (1, 0, 2))
        Area = np.transpose(Area, (1, 0))
    return Norm, Tang, Orth, Area, dl


def main(nChord=4, nSpan=13, nStep=10, chord=1, span=8.0, alpha=5, 
         omega_pitch=0, A_pitch=0,
         omega_heave=0, A_heave=0, 
         U0=10.0, U0_wind=0, U0_body=0, rho=1.0, outputDir='', simName='default', motionType='body',
         debug_print=False
         ):
    

    V_wind_i = np.zeros(3) # Wind velocity in inertial frame
    U0_body_b = 0

    if motionType=='body':
        U0_body_b = U0
    elif motionType=='wind':
        U0_body_b = 0
        V_wind_i[0]  = U0
    elif motionType=='mixed':
        V_wind_i[0]  = U0_wind
        U0_body_b = U0_body
        U0 = U0_wind # We chose wind as a ref velocity
    else:
        raise NotImplementedError('Unknown motionType '+motionType)


    GD = GLOBAL_DATA()
    # Initialize arrays
    GD.initialize_arrays(nChord, nSpan, nStep)

    # --- Derived parameters
    Vref = U0
    VINF = U0 # TODO
    span_half = span/2.0 # TODO
    DX = chord / GD.nChord
    alpha = alpha * PI / 180.0
    DT = DX / Vref / 4.0
    T = -DT
    NW_length = 0.3 * Vref * DT

    S_half = span_half * chord
    AR = 2.0 * span_half**2 / S_half

    # --- Reference data for wing
    CL_ref = 2.0 * PI * alpha / (1.0 + 2.0 / AR)
    if abs(CL_ref) < 1.0e-20:
        CL_ref = CL
    FG_ref = (0.5 * Vref * S_half) * CL_ref

    # --- Wing geometry
    # Mesh wing and calculate collocation points
    GD.r_LS, GD.CP, GD.Area, vSpan, vChord = rectangularWingPanelling(GD.nChord, GD.nSpan, bSpan = span_half, chord=chord, NW_length=NW_length)
    # GD.r_LS2, GD.CP2 = rotAndTrans(GD.r_LS, -alpha, )
    GD.r_LS = rotAndTrans_abs(GD.r_LS, alpha) # Positive rotation about y
    GD.CP    = rotAndTrans_abs(GD.CP , alpha)
    GD.dl_LL = np.diff(vSpan)
    Norm0, Tang0, Orth0, Area0, dl0 = panlInfo(GD.r_LS, spanFirst=False)

    # Backup Initial position for rigid body motion
    GD.CP0 = GD.CP.copy()
    GD.r_LS0 = GD.r_LS.copy()

    if nChord==4 and nSpan==13 and span==8:
        test_geometry(GD, S_half)

    # --- Output initial geometry to screen
    print("-" * 56)
    print(f" alpha: {alpha*180/np.pi:10.2f}  span_half : {span_half:10.2f}  C : {chord:13.2f}")
    print(f" Heave  : {omega_heave:10.2f}rad/s   A: {A_heave:13.2f} k={omega_heave*chord/(2*U0)}")
    print(f" S_half : {S_half:10.2f}  AR : {AR:13.2f}")
    print(f" nChord : {GD.nChord:10d}  NS : {GD.nSpan:10d} \n")
    # Ensure output directory exists
    os.makedirs("_outputs", exist_ok=True)
    os.makedirs("_vtk", exist_ok=True)
    with open("_outputs/uLS_{}.csv".format(simName), "w") as f:
        f.write("#SX,T,SZ,VINF,TETA,OMEGA,CL,L,CM,CD,CL_rel,Gamma_rel,omT\n")

    # --- Compute A
    # TODO, could that change? NOTE: some local inclination was added in original code
    GD.SNO[:GD.nChord] = np.sin(alpha)
    GD.CSO[:GD.nChord] = np.cos(alpha)
    GD.AA = computeA(GD.CP0, GD.r_LS0, GD, Norm0)
    if debug_print:
        for K in range(GD.AA.shape[0]):
            for L in range(GD.AA.shape[1]):
                print(f"AA {K+1:5d} {L+1:5d} {GD.AA[K, L]:15.6f}")

    # --- Main time-stepping loop
    for IT in range(nStep):

        # --- Wake rollup calculation between T and T+DT  - Update Continuous States from T to T+DT
        # Freestream rollup contribution
        GD.Uind_wake_on_wake[:, :, 0] = V_wind_i[0]
        GD.Uind_wake_on_wake[:, :, 1] = V_wind_i[1]
        GD.Uind_wake_on_wake[:, :, 2] = V_wind_i[2]
        # Wake self-induced contribution
        IT_loc = IT-1
        if IT_loc >= 1:
            IW = 0
            NWMAX = 5    # Max wake roll up elements
            if IT_loc >= NWMAX-1:
                IW = IT_loc - NWMAX + 1
            for iAge in range(IW, IT_loc):
                for iSpan in range(GD.nSpan+1):
                    U = UI_from_all(GD.r_NW[iAge, iSpan, :], IT_loc, GD.r_NW, GD.Gamma_NW, GD.r_LS, GD.Gamma_LS)
                    GD.Uind_wake_on_wake[iAge, iSpan, :] += U
        GD.r_NW += GD.Uind_wake_on_wake*DT #GD.r_NW[IW-1:IT-1, :nSpanP1, :] += GD.Uind_wake_on_wake[IW-1:IT-1, :nSpanP1, :]*DT

        T += DT

        # --- Rigid body motion - Positions at t+DT - Inputs at T+DT
        # Absolute structural displacement of body origin in inertial frame
        rO_str_i = np.zeros(3) 
        rO_str_i[0] = -U0_body_b * T # We are going along negative x
        rO_str_i[2] = A_heave * np.sin(omega_heave * T)
        # V0_str_i = dr0/dt - Absolute structural velocity of body origin in inertial frame
        # TODO there is a choice as to what U0_body_b means. It sounds like KP defines it ias always along x_body
        VO_str_i = np.zeros(3) 
        VO_str_i[0] = - U0_body_b  # TODO or V0_str_b[0] = -U0_body_b
        VO_str_i[2] = A_heave * omega_heave * np.cos(omega_heave * T)

        # Rotate
        TETA = A_pitch * np.sin(omega_pitch * T) # Pitch angle, NOTE: I THINK IT'S USING A NEGATIVE CONVENTION
        SN1 = np.sin(TETA)
        CS1 = np.cos(TETA)
        # ###  VINF = -np.cos(TETA) * DSX - np.sin(TETA) * DSZ # TODO COMMENTED OUT
        VO_str_b = inertial_to_body_theta(VO_str_i, TETA)
        VO_str_b[0] =- U0_body_b # TODO TODO HACK


        # NOTE: I DISAGREE WITH THIS, LIKELY WRONG TETA CONVENTION
        WT = - ( - SN1 * VO_str_i[0] + CS1 * VO_str_i[2])

        # DSX = -VINF
        # DSZ = BH * OM * cos(OM * T)
        # WT =  SN1 * DSX - CS1 * DSZ

        GD.SNO[:GD.nChord] = np.sin(alpha) # NOTE: add local inclination
        GD.CSO[:GD.nChord] = np.cos(alpha)

        GD.r_LS = rotAndTrans_abs_opp(GD.r_LS0, TETA, rO_str_i) # TODO I think convention is wrong for TETA
        GD.CP   = rotAndTrans_abs_opp(GD.CP0  , TETA, rO_str_i)

        # Vectorized wake shedding points
        GD.r_NW[IT, :, :] = GD.r_LS[GD.nChord, :, :]  # Last wing point is first NW point

        # --- Solving circulation on LS and NW at T+DT - Solve Constraint at T+DT
        # Induced velocities on LS
        K = 0
        if IT > 0:
            for I in range(GD.nChord):
                for J in range(GD.nSpan):
                    U =  UI_from_wake(GD.CP[I, J, :], IT,  GD.r_NW, GD.Gamma_NW)
                    U1 = UI_from_wake(GD.CP[I, J, :], IT,  GD.r_NW, GD.Gamma_NW, flip=True)
                    U = U + U1
                    # From inertial to non-pitched body frame
                    U11, _, W11 = inertial_to_body_theta(U, -TETA) # TODO WRONG CONVENTION
                    #U11 =  U[0] * CS1 + U[2] * SN1 # Theta projection if any pitching
                    #W11 = -U[0] * SN1 + U[2] * CS1
                    # From body-frame to normal of each panel
                    # OK: Un = u|b * sin("alpha") + w|b * cos("alpha")
                    GD.WW[K] = U11 * GD.SNO[I] + W11 * GD.CSO[I] # np.dot( Norm0[I,J], (U11, 0, W11))
                    GD.Uind_wake_on_LS[I, J, 2] = W11 # NOTE: projected by theta
                    K += 1
        # RHS - 
        K = 0
        for I in range(GD.nChord):
            for J in range(GD.nSpan):
                # NOTE: this is about -n it seems
                # NOTE: MANU: I DISAGREE WITH MOST OF THE TERMS BELOW
                #   - NOTE: omega r below is for rotation about origin (LE), with r taken at t=0 where the angle of attack is already included
                GD.DW[K] = 0
                GD.DW[K] += -VINF * GD.SNO[I]    # -V0_str_i . n_hat, OK only if V0_str_b[0] = VINF
                # GD.DW[K] += VO_str_b[0] * GD.SNO[I]  # V0_str_i . n_hat, OK only if V0_str_b[0] = VINF
                # GD.DW[K] += - V_wind_b[0] * GD.SNO[I] - V_wind_b[2] * GD.CSO[I] # - Vwind . n_hat

                GD.DW[K] += GD.CP0[I, J, 0] * omega_pitch  # I DISAGREE, it should be CP0 before rotation by alpha
                GD.DW[K] += - WT # TODO I THINK IT NEEDS further projection using SNO and CSO
                K += 1

        # Solve the linear system
        RHS = GD.DW - GD.WW # NOTE: I'd like the opposite, but that's likely because of the Gamma Convention of KP
        if IT == 0:
            LU, PIV = lu_factor(GD.AA)
        Gamma_LS_flat = lu_solve((LU, PIV), RHS)
        # Reshape the solution on the lifting surface lattice
        GD.Gamma_LS = Gamma_LS_flat.reshape(GD.nChord, GD.nSpan)
        # Assign Circulation of first Near wake panel to match Kutta condition
        GD.Gamma_NW[IT, :] = GD.Gamma_LS[GD.nChord-1, :]

        if debug_print:
            for K in range(nSpan * nChord):
                print(f"K {IT+1:4d} {K+1:4d} {GD.DW[K]:10.4f} {GD.WW[K]:10.4f} {Gamma_LS_flat[K]:10.4f}")
            if IT == 9:
                raise ValueError("Test failed")



        # --- Force calculations - CalcOutput at T+DT
        FL, FD, FM, CL, CD, CM, Fl, Fd, Fm, dL_IJ, dD_IJ, dM_IJ, FG = LS_calcForce(rho, chord, Vref, S_half, DT, DX, GD.Gamma_LS, GD.r_LS, GD.CP, GD.Gamma_LS_prev, GD, TETA)
        # Write force to screen and output file
        CLT = CL / CL_ref
        CFG = FG / FG_ref
        # Output results
        #print(f" T={T:10.2f}  SX={rO_str_i[0]:10.2f}  SZ={rO_str_i[2]:10.2f}  VINF={VINF:10.2f}  TETA={TETA:10.2f}  OMEGA={omega_pitch:10.2f}")
        print(f"CL={CL:10.4f}  L={FL:10.4f}  CM={CM:10.4f}  CD={CD:10.4f}  L/L(INF)={CLT:10.4f}  GAMA/GAMA(INF)={CFG:10.4f}")
        with open("_outputs/uLS_{}.csv".format(simName) , "a") as f:
            f.write(f"{-rO_str_i[0]},{T},{rO_str_i[2]},{Vref},{TETA},{omega_pitch},{CL},{FL},{CM},{CD},{CLT},{CFG},{np.mod(omega_heave*T, 2*np.pi)}\n")

        # After updating GD.QW and GD.Gamma_NW for the current time step IT:
        write_wing_vtk(IT+1, GD, outputDir=outputDir, simName=simName)
        write_wake_vtk(IT+1, GD, outputDir=outputDir, simName=simName)

        # --- Prepare for next time step
        GD.Gamma_LS_prev = GD.Gamma_LS


def LS_calcForce(rho, chord, Vref, S_ref, dt, DX, Gamma_LS, r_LS, CP, Gamma_LS_prev, GD, TETA):
    nChord, nSpan = Gamma_LS.shape
    dL_IJ    = np.zeros((nChord, nSpan))
    dD_IJ    = np.zeros((nChord, nSpan))
    dM_IJ    = np.zeros((nChord, nSpan))
    FG = 0.0
    qdyn = 0.5 * rho * Vref * Vref # dynamic pressure
    for J in range(nSpan):
        for I in range(nChord):
            # --- Lift
            if I == 0:
                dGamma_shed = Gamma_LS[I, J] # shed segment between two chordwise panels on LS
                SIGMA1 = (0.5 * dGamma_shed ) * DX
            else:
                dGamma_shed = Gamma_LS[I, J] - Gamma_LS[I-1, J] # shed segment between two chordwise panels on LS
                SIGMA1 = (0.5 * dGamma_shed + Gamma_LS[I-1, J]) * DX
            DFDT = (SIGMA1 - GD.dLT[I, J]) / dt # some kind of DGamma/Dt
            GD.dLT[I, J] = SIGMA1
            dL_IJ[I, J] = rho * (Vref * dGamma_shed + DFDT) * GD.dl_LL[J] * GD.CSO[I]
            # --- Drag
            U = UI_from_wing_trailed(CP[I, J, :], Gamma_LS, r_LS)
            U2 = UI_from_wing_trailed(CP[I, J, :], Gamma_LS, r_LS, flip=True)
            _, _, W1 = inertial_to_body_theta(U, -TETA) # TODO WRONG CONVENTION
            _, _, W2 = inertial_to_body_theta(U2, -TETA) # TODO WRONG CONVENTION
            W8 = W1 + W2
            #W8 = U[2] + U2[2]
            Wtot = GD.Uind_wake_on_LS[I, J, 2] + W8 # Total downwash from wake and LS (trailed)
            CTS = -Wtot / Vref
            DD1 = rho * GD.dl_LL[J] * DFDT * GD.SNO[I]         # Drag from varying lift
            DD2 = rho * GD.dl_LL[J] * Vref * dGamma_shed * CTS # Drag from downwash
            dD_IJ[I, J] = DD1 + DD2
            # Moment
            DXM = (GD.r_LS0[I, J, 0] + GD.r_LS0[I, J+1, 0]) / 2.0 
            dM_IJ[I,J] = dL_IJ[I, J] * DXM
            # KJ force?
            FG += dGamma_shed * GD.dl_LL[J]
            # debug_print
            #print(f"Drag{I+1:4}{J+1:4}{DFDT:10.3f}{DD1:10.3f}{GD.Uind_wake_on_LS[I, J, 2]:10.3f}{Gamma_LS[I,J]:10.3f}{W1:10.3f}{W2:10.3f}{W8:10.3f}{CTS:10.3f}{DD2:10.3f}")
    # Spanwise forces
    Fl = np.sum(dL_IJ, axis=0) # TODO check axis
    Fd = np.sum(dD_IJ, axis=0) # TODO check axis
    Fm = np.sum(dM_IJ, axis=0) # TODO check axis
    # Total forces
    FL = np.sum(dL_IJ.flatten())
    FD = np.sum(dD_IJ.flatten())
    FM = np.sum(dM_IJ.flatten())
    # Force coefficients
    CL = FL / (qdyn * S_ref)
    CD = FD / (qdyn * S_ref)
    CM = FM / (qdyn * S_ref * chord)
    return FL, FD, FM, CL, CD, CM, Fl, Fd, Fm, dL_IJ, dD_IJ, dM_IJ, FG


def write_wake_vtk(it, GD, outputDir='', simName='default'):
    """
    Write the current wake lattice and vorticity to a VTK file for visualization.
    """
    mvtk = VTK_Misc()
    filename = os.path.join(outputDir, f"_vtk/{simName}_wake_{it:05d}.vtk")
    # r_NW shape: (NTMAX, NSMAX+1, 3)
    # Gamma_NW shape: (NTMAX, NSMAX)
    QW    = GD.r_NW[:it, :, :]
    Gamma = GD.Gamma_NW[:it-1, :]
    QWt   = np.transpose(QW, (1, 0, 2))
    Gammat = np.transpose(Gamma, (1, 0))
    WrVTK_Lattice(filename, mvtk, QWt, Gammat)

def write_wing_vtk(it, GD, outputDir='', simName='default'):
    """
    Write the current lifting surface lattice and circulation to a VTK file for visualization.
    """
    mvtk     = VTK_Misc()
    filename = os.path.join(outputDir, f"_vtk/{simName}_wing_{it:05d}.vtk")

    rWing  = GD.r_LS[:GD.nChord+1, :GD.nSpan+1, :]  # (nChord, nSpan, 3)
    Gamma = GD.Gamma_LS[:GD.nChord, :GD.nSpan]       # (nChord-1, nSpan-1)
    rLSt   = np.transpose(rWing, (1, 0, 2))    # (nSpan, nDepth, 3)
    Gammat = np.transpose(Gamma, (1, 0))
    WrVTK_Lattice(filename, mvtk, rLSt, Gammat)



def test_geometry(GD, S_half):
    np.testing.assert_almost_equal(S_half, 4.0)
    np.testing.assert_almost_equal(GD.r_LS0[0,:,0], [0.06226217]*14)
    np.testing.assert_almost_equal(GD.r_LS0[1,:,0], [0.3113108]*14)
    np.testing.assert_almost_equal(GD.r_LS0[-1,:,0], [1.0148733]*14) 
    np.testing.assert_almost_equal(GD.r_LS0[:,0,1], [0]*5 )
    np.testing.assert_almost_equal(GD.r_LS0[:,1,1], [0.3076923]*5 )
    np.testing.assert_almost_equal(GD.r_LS0[:,-1,1], [4]*5 )
    np.testing.assert_almost_equal(np.unique(GD.r_LS0[:,:,2]), [-0.08878991, -0.07081404, -0.04902511, -0.02723617, -0.00544723])
    np.testing.assert_almost_equal(GD.CP0[0,:,0], [0.1867865]*13)
    np.testing.assert_almost_equal(GD.CP0[1,:,0], [0.4358352]*13)
    np.testing.assert_almost_equal(GD.CP0[-1,:,0], [0.9339325]*13) 
    np.testing.assert_almost_equal(GD.CP0[:,0,1], [0.1538462]*4)
    np.testing.assert_almost_equal(np.unique(GD.CP0[:,:,2]), [-0.0817085, -0.0599196, -0.0381306, -0.0163417])
    np.testing.assert_almost_equal(GD.Area.flatten(), [0.0769231]*13*4)



if __name__ == "__main__":

    scriptDir = os.path.dirname(os.path.abspath(__file__))

    if True:
        nt = 10
        with Timer():
            # omega = 2 , k = 0.1
            # omega =10 , k = 0.5
            main(nChord=4, nSpan=13, nStep=nt, chord=1, span=4.0, alpha=-5, omega_pitch =10, A_pitch=0.1745, outputDir=scriptDir, simName='pitch_AR=4_k=0.5', motionType='body')
        df = weio.read(os.path.join(scriptDir, './_outputs/uLS_pitch_AR=4_k=0.5.csv')).toDataFrame()
        df_ref = weio.read(os.path.join(scriptDir, './_outputs/uLS_pitch_AR=4_k=0.5_ref.csv')).toDataFrame()
        np.testing.assert_almost_equal(df['CL'].values[:nt], df_ref['CL'].values[:nt], 6)
        np.testing.assert_almost_equal(df['CD'].values[:nt], df_ref['CD'].values[:nt], 6)
        np.testing.assert_almost_equal(df['CM'].values[:nt], df_ref['CM'].values[:nt], 6)
        np.testing.assert_almost_equal(df['Gamma_rel'].values[:nt], df_ref['Gamma_rel'].values[:nt], 6)
        print('[ OK ] test pass')

    if True:
        nt = 10
        with Timer():
            main(nChord=4, nSpan=13, nStep=nt, chord=1, span=8.0, outputDir=scriptDir, simName='acc_body', motionType='body')
        # ---
        df = weio.read(os.path.join(scriptDir, './_outputs/uLS_acc_body.csv')).toDataFrame()
        df_ref = weio.read(os.path.join(scriptDir, './_outputs/uLS_ref.csv')).toDataFrame()
        np.testing.assert_almost_equal(df['CL'].values[:nt], df_ref['CL'].values[:nt], 6)
        np.testing.assert_almost_equal(df['CD'].values[:nt], df_ref['CD'].values[:nt], 6)
        np.testing.assert_almost_equal(df['CM'].values[:nt], df_ref['CM'].values[:nt], 6)
        np.testing.assert_almost_equal(df['Gamma_rel'].values[:nt], df_ref['Gamma_rel'].values[:nt], 6)
        print('[ OK ] test pass')

        with Timer():
            main(nChord=4, nSpan=13, nStep=nt, chord=1, span=8.0, outputDir=scriptDir, simName='acc_wind', motionType='wind')
        # ---
        df = weio.read(os.path.join(scriptDir, './_outputs/uLS_acc_wind.csv')).toDataFrame()
        df_ref = weio.read(os.path.join(scriptDir, './_outputs/uLS_ref.csv')).toDataFrame()
        np.testing.assert_almost_equal(df['CL'].values[:nt], df_ref['CL'].values[:nt], 6)
        np.testing.assert_almost_equal(df['CD'].values[:nt], df_ref['CD'].values[:nt], 6)
        np.testing.assert_almost_equal(df['CM'].values[:nt], df_ref['CM'].values[:nt], 6)
        np.testing.assert_almost_equal(df['Gamma_rel'].values[:nt], df_ref['Gamma_rel'].values[:nt], 6)
        print('[ OK ] test pass')

    if True:
        nt = 10
        with Timer():
            # omega = 2 , k = 0.1
            # omega =10 , k = 0.5
            main(nChord=4, nSpan=13, nStep=nt, chord=1, span=4.0, alpha=-5, omega_heave =10, A_heave=0.1, outputDir=scriptDir, simName='heave_AR=4_k=0.5_wind', motionType='wind')
        df = weio.read(os.path.join(scriptDir, './_outputs/uLS_heave_AR=4_k=0.5_wind.csv')).toDataFrame()
        df_ref = weio.read(os.path.join(scriptDir, './_outputs/uLS_heave_AR=4_k=0.5_ref.csv')).toDataFrame()
        np.testing.assert_almost_equal(df['CL'].values[:nt], df_ref['CL'].values[:nt], 6)
        np.testing.assert_almost_equal(df['CD'].values[:nt], df_ref['CD'].values[:nt], 6)
        np.testing.assert_almost_equal(df['CM'].values[:nt], df_ref['CM'].values[:nt], 6)
        np.testing.assert_almost_equal(df['Gamma_rel'].values[:nt], df_ref['Gamma_rel'].values[:nt], 6)
        print('[ OK ] test pass')

        with Timer():
            main(nChord=4, nSpan=13, nStep=nt, chord=1, span=4.0, alpha=-5, omega_heave =10, A_heave=0.1, outputDir=scriptDir, simName='heave_AR=4_k=0.5_body', motionType='body')
        df = weio.read(os.path.join(scriptDir, './_outputs/uLS_heave_AR=4_k=0.5_body.csv')).toDataFrame()
        df_ref = weio.read(os.path.join(scriptDir, './_outputs/uLS_heave_AR=4_k=0.5_ref.csv')).toDataFrame()
        np.testing.assert_almost_equal(df['CL'].values[:nt], df_ref['CL'].values[:nt], 6)
        np.testing.assert_almost_equal(df['CD'].values[:nt], df_ref['CD'].values[:nt], 6)
        np.testing.assert_almost_equal(df['CM'].values[:nt], df_ref['CM'].values[:nt], 6)
        np.testing.assert_almost_equal(df['Gamma_rel'].values[:nt], df_ref['Gamma_rel'].values[:nt], 6)
        print('[ OK ] test pass')

    if True:
        nt = 10
        with Timer():
            # omega = 2 , k = 0.1
            # omega =10 , k = 0.5
            main(nChord=4, nSpan=13, nStep=nt, chord=1, span=4.0, alpha=-5, omega_heave =2, A_heave=0.1, outputDir=scriptDir, simName='heave_AR=4_k=0.1')
        df = weio.read(os.path.join(scriptDir, './_outputs/uLS_heave_AR=4_k=0.1.csv')).toDataFrame()
        df_ref = weio.read(os.path.join(scriptDir, './_outputs/uLS_heave_AR=4_k=0.1_ref.csv')).toDataFrame()
        np.testing.assert_almost_equal(df['CL'].values[:nt], df_ref['CL'].values[:nt], 6)
        np.testing.assert_almost_equal(df['CD'].values[:nt], df_ref['CD'].values[:nt], 6)
        np.testing.assert_almost_equal(df['CM'].values[:nt], df_ref['CM'].values[:nt], 6)
        np.testing.assert_almost_equal(df['Gamma_rel'].values[:nt], df_ref['Gamma_rel'].values[:nt], 6)
        print('[ OK ] test pass')

    #if True:
    #    with Timer():
    #        main(nChord=1, nSpan=3, nStep=50, chord=1, span=8.0, outputDir=scriptDir)
    #    # ---
    #    nt = 10
    #    df = weio.read(os.path.join(scriptDir, './_outputs/uLS.csv')).toDataFrame()
    #    df_ref = weio.read(os.path.join(scriptDir, './_outputs/uLS_ref.csv')).toDataFrame()
    #    np.testing.assert_almost_equal(df['CL'].values[:nt], df_ref['CL'].values[:nt], 6)
    #    np.testing.assert_almost_equal(df['CD'].values[:nt], df_ref['CD'].values[:nt], 6)
    #    np.testing.assert_almost_equal(df['CM'].values[:nt], df_ref['CM'].values[:nt], 6)
    #    np.testing.assert_almost_equal(df['Gamma_rel'].values[:nt], df_ref['Gamma_rel'].values[:nt], 6)
