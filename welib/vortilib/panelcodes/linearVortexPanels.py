import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
from welib.vortilib.panelcodes.panel_tools import panel_geometry
from welib.vortilib.elements.LinearVortexPanel2D import lvp_u11 # TODO

def LVP_solve(XP, YP, Uxy=None, fU=None, closed=True, verbose=False, alpha=None, x_Cm_ref=0.25):
    """
    Solve panel method for a given geometry and external velocity

        N linear vortex panels

    OUTPUTS:
     - out: storage for multiple variables, like Cp
    """
    out = {}
    # --- Velocity function
    if Uxy is not None:
        fU =lambda X, Y : (X*0+Uxy[0], X*0+Uxy[1])

    # --- Geometry
    PP, mids, dP, ds, t_hat, n_hat, phi, ns = panel_geometry(XP, YP, closed_expected=closed, force_clockwise=True)

    # --- Panels / Control points    
    CP = mids # Positions of control points
    nPanels = len(PP) - 1    
    
    # --- Build influence matrix

    # --- Manu's method, unfinished    
    A_UnV, A_UtV = build_influence_matrix(CP, PP, n_hat, t_hat)
    M = np.zeros((nPanels+1, nPanels+1))
    M[:-1, :] = A_UnV
    M2        = M # Temporary backup for comparison

    # --- Katz Plotkin
    TH     = np.arctan2(dP[:,1], dP[:,0])
    #TH[TH<0] = TH[TH<0] + 2*np.pi
    M, B  = build_influence_matrix_KP(CP, PP[:-1], PP[1:], TH)

    # --- Kutta condition
    M[-1,0]  = 1
    M[-1,-1] = 1

    # --- Right hand side
    Ux, Uy = fU(CP[:,0], CP[:,1])
    U0_n = (Ux*n_hat[:,0] + Uy*n_hat[:,1])
    U0_t = (Ux*t_hat[:,0] + Uy*t_hat[:,1])
    rhs = np.zeros(nPanels+1)
    rhs[:-1] = -U0_n

    # --- Solve
    gamma_nodes = np.linalg.solve(M, rhs)

    # Tangential velocities on panel
    UV_t = B @ gamma_nodes # Tangential velocities on panels due to vortices
    Utot_t = U0_t + UV_t
    # Normal velocities on panel (residual check)
    UV_n =  M[:-1,:] @ gamma_nodes # Normal velocities on panels due to vortices
    Utot_n = U0_n + UV_n

    Vwall = Utot_n[:,np.newaxis]*n_hat + Utot_t[:,np.newaxis]*t_hat

    # Pressure coefficient
    Vinf2 = Ux**2 + Uy**2 # NOTE: only fine for constant velocity
    Cp = 1-(Utot_t)**2/Vinf2

    # --- Outputs
    # Output: Geometry
    out['x']        = XP
    out['y']        = YP
    out['theta']    = np.arctan2(YP, XP)
    out['ds']       = ds        # Panel lengths
    out['n']        = n_hat
    out['t']        = t_hat
    out['CP']       = CP
    out['theta_CP'] = np.arctan2(CP[:,1], CP[:,0])
    # Output: System
    out['rhs']      = rhs       # Right hand side
    out['M']        = M         # System matrix
    # Output: Solution
    out['gammas']   = gamma_nodes
    # Output: Velocity at wall
    out['Vwall'] = Vwall
    out['Un'] = Utot_n
    out['Ut'] = Utot_t
    # Output: Cp
    out['Cp'] = Cp #1 - (Vwall[:,0]**2 + Vwall[:,1]**2)/(Ux**2+Uy**2)
    # Output: Loads
    if alpha is not None:
        # TODO TODO TODO
        # angle of panel normal w.r.t. horizontal and include AoA
        delta                = phi + (np.pi/2) # Angle from x-axis to normal vector
        beta                 = delta - alpha   # Angle between freestream and normal vector
        beta[beta > 2*np.pi] = beta[beta > 2*np.pi] - 2*np.pi
        out['Cn'] = -Cp*ds*np.sin(beta) # Normal to chord
        out['Cc'] = -Cp*ds*np.cos(beta) # Along chord
        out['Cl'] = sum(out['Cn']*np.cos(alpha)) - sum(out['Cc']*np.sin(alpha))
    out['Cx'] = sum(Cp*ds*n_hat[:,0])
    out['Cy'] = sum(Cp*ds*n_hat[:,1])
    out['Cm'] = sum(Cp*(CP[:,0]-x_Cm_ref)*ds*np.cos(phi))
    #out['Cl_KJ'] = 2*sum(gamma*ds)
    return out






def build_influence_matrix(CP, PP, n_hat, t_hat):
    # TODO lvp_u11 unvalidated yet!
    nPanels = len(PP) - 1
    A_UnV   = np.zeros((nPanels, nPanels+1))
    A_UtV   = np.zeros((nPanels, nPanels+1))
    for i in range(nPanels):
        for j in range(nPanels+1):
            if j == 0:
                gamma1, gamma2 = 1, 0
            elif j == nPanels:
                gamma1, gamma2 = 0, 1
            else:
                gamma1, gamma2 = 0, 0
            if j < nPanels:
                u, v = lvp_u11(CP[i], PP[j], PP[j+1], gamma1=1, gamma2=0)
                A_UnV[i, j]   = u*n_hat[i,0] + v*n_hat[i,1]
                A_UtV[i, j]   = u*t_hat[i,0] + v*t_hat[i,1]
            if j > 0:
                u, v = lvp_u11(CP[i], PP[j-1], PP[j], gamma1=0, gamma2=1)
                A_UnV[i, j]  += u*n_hat[i,0] + v*n_hat[i,1]
                A_UtV[i, j]  += u*t_hat[i,0] + v*t_hat[i,1]
    return A_UnV, A_UtV

def build_influence_matrix_KP(CO, PT1, PT2, TH):
    """
    Katz Plotkin influence matrix for linear vortex panels

    Build the influence matrix A and tangential matrix B for the linear vortex panel method,
    following the logic and formulas of KatzProg6_2DLinearVortex.f90.

    CO: (M,2) array of collocation points
    PT1, PT2: (M,2) arrays of panel start/end points
    TH: (M,) array of panel angles
    Returns:
        A: (M+1, M+1) influence matrix (normal velocity)
        B: (M, M+1) tangential velocity matrix
    """
    M      = len(PT1)
    N      = M + 1
    A      = np.zeros((N, N))
    B      = np.zeros((M, N))
    for i in range(M):
        # TODO convert this into a single panel function
        HOLDA = 0.0
        HOLDB = 0.0
        for j in range(M):
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
    return A, B


if __name__ == "__main__":
    from welib.vortilib.panelcodes.panel_tools import plot_Cp, plot_pressure_force_bars
    scriptDir = os.path.dirname(__file__)

    testcase=2 # 2: Van de Vooren (Katz-Plotkin), 3: XFoil

    if testcase==2:
        # Test case 2 - Van de Vooren (Katz Plotkin)  # Sharp
        # NUMBER OF AIRFOIL PANELS, M    :   90
        # THE ANGLE OF ATTACK IN DEGREES :   5
        # THICKNESS COEFF. Eps (<1)      :   0.075
        # T.E. ANGLE COEFF. K (1-2)      :   1.90555555555
        airfoil_file = os.path.join(scriptDir, 'data/VonDeVooren_esp0.075_k1.906_AFOIL2.csv') # KatzPlotkin example - VanDeVooren
        df = pd.read_csv(airfoil_file)
        XP, YP = df['x'].values, df['y'].values
        df_ref = pd.read_csv(os.path.join(scriptDir, 'data/VonDeVooren_Cp_lvortex.csv')) # KatzPlotkin example - VanDeVooren
        Cp_ref = df_ref['Cp_[-]'].values
    elif testcase==3:
        # Test case 3 - XFoil NACA 2412 PPAR N 170 P 4 T 1 R 1 # Slightly blunt
        airfoil_file = os.path.join(scriptDir,'data/NACA2412.txt')
        df = pd.read_csv(airfoil_file)
        XP, YP = df['x'].values, df['y'].values
        import pickle
        with open(os.path.join(scriptDir, 'tests/LSN_LV1_NACA2412.pkl'), 'rb') as f:
            data = pickle.load(f) 
        Cp_ref = data['Cp']


    alpha = 5 * np.pi / 180  # Angle of attack in radians
    Vinf  = 1.0              # Freestream velocity (m/s)
    Uxy   = np.array([Vinf*np.cos(alpha), Vinf*np.sin(alpha)])
    fU =lambda X, Y : (X*0+Uxy[0], X*0+Uxy[1]) # External velocity function

    # --- Panel method
    out = LVP_panel_solve(XP, YP, fU=fU, alpha=alpha)
    x = out['CP'][:,0]
    Cp = out['Cp']
    CP = out['CP']

    # --- Misc outputs and plotting
    #df_num = pd.DataFrame(data=np.column_stack([x, Cp]), columns=['x','Cp'])
    #df_num.to_csv(os.path.join(scriptDir,'CP_output_new.csv'))

    ax = plot_Cp(CP[:,0], Cp, Cp_ref=Cp_ref, simple=False, label='linear vortex')

    plt.show()
