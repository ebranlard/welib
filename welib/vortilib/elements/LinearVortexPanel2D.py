import numpy as np

def lvp_u11(rCP, rS1, rS2, gamma1=1, gamma2=0, tol=1e-12):
    # TODO BUGGY DEBUG ME
    """
    Velocity induced at rCP by a linear-varying vortex panel from rS1 to rS2.
    gamma1: vortex strength at rS1
    gamma2: vortex strength at rS2
    """
    import numpy as np
    xCP, yCP = rCP
    x1, y1   = rS1
    x2, y2   = rS2
    dx       = x2 - x1
    dy       = y2 - y1
    L        = np.hypot(dx, dy)
    if L < tol: return 0.0, 0.0
    # Panel frame
    C        = dx / L
    S        = dy / L
    # Transform control point to panel frame
    xC =  (xCP - x1) * C + (yCP - y1) * S
    yC = -(xCP - x1) * S + (yCP - y1) * C
    # Endpoints in panel frame
    xA, xB = 0.0, L
    yA, yB = 0.0, 0.0
    # Influence integrals
    def I1(xa, xb, y):
        return np.arctan2(xb - xC, y - yC) - np.arctan2(xa - xC, y - yC)
    def I2(xa, xb, y):
        ra2 = (xa - xC)**2 + (y - yC)**2
        rb2 = (xb - xC)**2 + (y - yC)**2
        return 0.5 * np.log(rb2 / ra2)
    # Linear variation: gamma(s) = gamma1 + (gamma2-gamma1)*s/L
    # Induced velocity (Katz & Plotkin, XFoil, etc.)
    I1v = I1(xA, xB, yA)
    I2v = I2(xA, xB, yA)
    u_panel = (gamma1 * I1v + (gamma2 - gamma1) * ( (xB - xC) * I1v - L * I2v ) / L ) / (2 * np.pi)
    v_panel = (gamma1 * I2v + (gamma2 - gamma1) * ( (xB - xC) * I2v + L * I1v ) / L ) / (2 * np.pi)
    # Rotate back to global frame
    u =  u_panel * C - v_panel * S
    v =  u_panel * S + v_panel * C
    return u, v

def lvp_u11_kp(rCP, rS1, rS2, gamma1=1, gamma2=0, tol=1e-12):

    M      = len(PT1)
    N      = M + 1

    # TODO convert this into a single panel function
    HOLDA = 0.0
    HOLDB = 0.0

    for j in [0, M-1]:
        XT  = CO[i,0] - PT1[j,0]
        ZT  = CO[i,1] - PT1[j,1]
        X2T = PT2[j,0] - PT1[j,0]
        Z2T = PT2[j,1] - PT1[j,1]
        cth = np.cos(TH[j])
        sth = np.sin(TH[j])
        X   = XT*cth + ZT*sth
        Z   = -XT*sth + ZT*cth
        X2  = X2T*cth + Z2T*sth
        # Z2 = 0 always
        R1  = np.sqrt(X**2 + Z**2)
        R2  = np.sqrt((X-X2)**2 + Z**2)
        TH1 = np.arctan2(Z, X)
        TH2 = np.arctan2(Z, X-X2)
        if i == j:
            U1L = -0.5*(X-X2)/(X2)
            U2L = 0.5*(X)/(X2)
            W1L = -1/(2*np.pi)
            W2L = 1/(2*np.pi)
        else:
            U1L = -(Z*np.log(R2/R1) + X*(TH2-TH1) - X2*(TH2-TH1)) / (2*np.pi*X2)
            U2L = (Z*np.log(R2/R1) + X*(TH2-TH1)) / (2*np.pi*X2)
            W1L = -((X2-Z*(TH2-TH1)) - X*np.log(R1/R2) + X2*np.log(R1/R2)) / (2*np.pi*X2)
            W2L = ((X2-Z*(TH2-TH1)) - X*np.log(R1/R2)) / (2*np.pi*X2)

        cti = np.cos(-TH[j])
        sti = np.sin(-TH[j])
        U1  = U1L*cti + W1L*sti
        U2  = U2L*cti + W2L*sti
        W1  = -U1L*sti + W1L*cti
        W2  = -U2L*sti + W2L*cti

        if j == 0:
            A[i,0]   = -U1*np.sin(TH[i]) + W1*np.cos(TH[i])
            HOLDA    = -U2*np.sin(TH[i]) + W2*np.cos(TH[i])
            B[i,0]   =  U1*np.cos(TH[i]) + W1*np.sin(TH[i])
            HOLDB    =  U2*np.cos(TH[i]) + W2*np.sin(TH[i])
        elif j == M-1:
            A[i,M-1] = -U1*np.sin(TH[i]) + W1*np.cos(TH[i]) + HOLDA
            A[i,N-1] = -U2*np.sin(TH[i]) + W2*np.cos(TH[i])
            B[i,M-1] =  U1*np.cos(TH[i]) + W1*np.sin(TH[i]) + HOLDB
            B[i,N-1] =  U2*np.cos(TH[i]) + W2*np.sin(TH[i])
        else:
            A[i,j]   = -U1*np.sin(TH[i]) + W1*np.cos(TH[i]) + HOLDA
            HOLDA    = -U2*np.sin(TH[i]) + W2*np.cos(TH[i])
            B[i,j]   =  U1*np.cos(TH[i]) + W1*np.sin(TH[i]) + HOLDB
            HOLDB    =  U2*np.cos(TH[i]) + W2*np.sin(TH[i])

# SUBROUTINE VLINV(XP,YP, X1,Y1, X2,Y2, GAM1, GAM2, U, V)
# ! Computes velocity (U,V) at (XP,YP) due to a linear-varying vortex panel
# ! from (X1,Y1) to (X2,Y2) with strengths GAM1 and GAM2 at each end

# REAL XP, YP, X1, Y1, X2, Y2, GAM1, GAM2, U, V
# REAL DX, DY, S, C, SINE, COSINE
# REAL X, Y, XA, XB, YA, YB
# REAL A, B, C1, D, E, F, G, H, I, J

# DX = X2 - X1
# DY = Y2 - Y1
# S  = SQRT(DX*DX + DY*DY)
# IF (S .EQ. 0.0) RETURN

# C = DX / S
# SINE = DY / S

# ! Transform control point to panel frame
# X =  ( (XP-X1)*C + (YP-Y1)*SINE )
# Y = (-(XP-X1)*SINE + (YP-Y1)*C )

# XA = 0.0
# XB = S
# YA = 0.0
# YB = 0.0

# ! Influence integrals
# A = ATAN2(XB-X, Y-YA) - ATAN2(XA-X, Y-YA)
# B = 0.5 * LOG( ((XB-X)**2 + (Y-YA)**2) / ((XA-X)**2 + (Y-YA)**2) )

# ! Induced velocities in panel frame
# U1 = (GAM1 * A + (GAM2-GAM1) * ( (XB-X)*A - S*B ) / S ) / (2.0*PI)
# V1 = (GAM1 * B + (GAM2-GAM1) * ( (XB-X)*B + S*A ) / S ) / (2.0*PI)

# ! Rotate back to global frame
# U = U1*C - V1*SINE
# V = U1*SINE + V1*C

# RETURN
# END
