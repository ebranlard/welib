"""
Setup:
 - n panels consisting of both consant source and constant vortex panels
 - n control points (CP) at the midpoints of the panels

Unknowns (n+1):
 - n unknowns sigma_i for the source panels
 - 1 unknown gamma for the vortex panel (they all have the same gamma)

Equation (n+1)
 - n no flow through conditions at the control points
 - 1 Kutta condition 

This is also referred to as Hess/Smith panel method, or Source Panel Vortex Panel (SPVP) method.
"""
import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
# Local 
import pickle

import numpy as np
from welib.essentials import *
from welib.vortilib.panelcodes.panel_tools import panel_geometry
from welib.vortilib.elements.SourcePanel2D import *
from welib.vortilib.elements.VortexPanel2D import *

def CSPN_CVP1_solve(XP, YP, Uxy=None, fU=None, verbose=False, alpha=None, x_Cm_ref=0.25):
    """
    Solve panel method for a given geometry and external velocity
     N linear source panels
     1 linear vorticity panel

    OUTPUTS:
     - out: storage for multiple variables, like Cp
    """
    out = {}
    # --- Velocity function
    if Uxy is not None:
        fU =lambda X, Y : (X*0+Uxy[0], X*0+Uxy[1])

    # --- Geometry
    PP, mids, dP, ds, t_hat, n_hat, phi, ns = panel_geometry(XP, YP, closed_expected=True, force_clockwise=True)

    # --- Build influence matrix
    CP = mids # Positions of control points
    nPanels = len(PP)-1
    A_UnS = np.zeros((len(CP), nPanels))
    A_UtS = np.zeros((len(CP), nPanels))
    for i in range(len(CP)):
        for j in range(nPanels):
            if (i == j): 
                #u, v = csp_u11(CP[i,:], PP[j,:], PP[j+1,:], principal=True)
                A_UnS[i,j] = 0.5  #  Principal value
                A_UtS[i,j] = 0 
            else:
                u, v = csp_u11(CP[i,:], PP[j,:], PP[j+1,:])
                A_UnS[i,j] = u*n_hat[i,0] + v*n_hat[i,1]
                A_UtS[i,j] = u*t_hat[i,0] + v*t_hat[i,1]


    # --- Build Influence matrix for vortex panels
    A_UnV = np.zeros((len(CP), nPanels))
    A_UtV = np.zeros((len(CP), nPanels))
    for i in range(len(CP)):
        for j in range(nPanels):
            if (i == j): 
                #u, v = cvp_u11(CP[i,:], PP[j,:], PP[j+1,:], principal=True)
                #print('PVV',u*n_hat[i,0] + v*n_hat[i,1], u*t_hat[i,0] + v*t_hat[i,1])
                A_UnV[i,j] = 0.0  #  Principal value
                A_UtV[i,j] = 0.5  #  Principal value # TODO TODO TODO TODO PRINCIPAL VALUE HERE IS 0.5
            else:
                u, v = cvp_u11(CP[i,:], PP[j,:], PP[j+1,:])
                A_UnV[i,j] = u*n_hat[i,0] + v*n_hat[i,1]
                A_UtV[i,j] = u*t_hat[i,0] + v*t_hat[i,1]
    if A_UtV[0,0] != 0.5:
        print('[WARN] Principal value of A_UtV[0,0] is not 0.5, it is {}'.format(A_UtV[0,0]))


    # --- Build system matrix.
    # We use sources for most of it, and a sum of gammas for the last column
    M = np.zeros((len(CP)+1, nPanels+1))
    M[:-1,:-1] = A_UnS
    # Column for main Gamma (superposition of all gammas)
    for i in range(len(CP)):
        M[i, nPanels] = - np.sum(A_UnV[i,:] )
    # Last row for Kutta condition, tangential velocity at 0 and last panels 
    for j in range(nPanels):
        M[nPanels, j] = A_UtS[0,j] + A_UtS[nPanels-1,j]
    if A_UtV[0,0] != 0.5:
        M[nPanels, nPanels] = -(sum(A_UtV[0,:] + A_UtV[nPanels-1,:])) + 1 
    else:
        M[nPanels, nPanels] = -(sum(A_UtV[0,:] + A_UtV[nPanels-1,:])) + 2

    # --- Right hand side
    Ux, Uy = fU(CP[:,0], CP[:,1])
    U0_n = (Ux*n_hat[:,0] + Uy*n_hat[:,1])
    U0_t = (Ux*t_hat[:,0] + Uy*t_hat[:,1])
    rhs = np.zeros(len(CP)+1)
    rhs[:-1] = - U0_n
    # Kutta condition, tangential freestream velocity on both panels
    Ux0, Uy0 = fU(CP[0,0], CP[0,1])
    Uxn, Uyn = fU(CP[-1,0], CP[-1,1])
    rhs[-1] = -(Ux0*t_hat[0,0] + Uy0*t_hat[0,1])  -(Uxn*t_hat[-1,0] + Uyn*t_hat[-1,1])
    rhs[-1] = -(U0_t[0] + U0_t[-1])  # Kutta condition, tangential velocity at 0 and last panels

    # --- Solve
    sigmas_gamma = np.linalg.solve(M, rhs)
    sigmas = sigmas_gamma[:len(CP)]
    gamma = sigmas_gamma[-1]

    # Tangential velocities on panel
    US_t = A_UtS @ sigmas
    UV_t = gamma - A_UtV @ np.array([gamma]*nPanels) # TODO problem in convention here
    Utot_t = U0_t + US_t + UV_t
    # Normal velocities on panel (residual check)
    US_n = A_UnS @ sigmas
    UV_n = - A_UnV @ np.array([gamma]*nPanels) # TODO problem in convention here
    Utot_n = U0_n + US_n + UV_n

    Vwall = Utot_n[:,np.newaxis]*n_hat + Utot_t[:,np.newaxis]*t_hat

    # Pressure coefficient
    Vinf2 = Ux**2 + Uy**2 # NOTE: only fine for constant velocity
    Cp = 1-(Utot_t)**2/Vinf2

    # --- Outputs
    # Output: Geometry
    out['x']        = XP
    out['y']        = YP
    out['theta']    = np.arctan2(YP, XP)
    out['ds']       = ds
    out['n']        = n_hat
    out['t']        = t_hat
    out['CP']       = CP
    out['theta_CP'] = np.arctan2(CP[:,1], CP[:,0])
    # Output: System
    out['rhs']      = rhs
    out['M']        = M
    # Output: Solution
    out['sigmas']   = sigmas
    out['gamma']    = gamma
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
    out['Cl_KJ'] = 2*sum(gamma*ds)
    return out


def CSPN_CVP1_u(X, Y, PP, sigmas, gamma, debug=False):
    nS = len(sigmas)
    U = np.zeros_like(X)
    V = np.zeros_like(Y)
    us, vs = csp_u(X, Y, PP, sigmas)
    gammas = np.array([gamma]*nS)
    uv, vv = cvp_u(X, Y, PP, gammas)
    U = us + uv
    V = vs + vv
    return U, V


if __name__ == '__main__':
    from welib.CFD.flows2D import flowfield2D, flowfield2D_plot
    from welib.vortilib.panelcodes.panel_tools import plot_Cp, plot_pressure_force_bars
    # Set numpy print options for nicer array display
    np.set_printoptions(precision=2, linewidth=140) #, suppress=True)

    scriptDir = os.path.dirname(__file__)
    # --- Parameters
    airfoil_file = os.path.join(scriptDir,'tests/NACA2412.txt')
    closed=True
    verbose=True

    Vinf = 1                                                                        # Freestream velocity [] (just leave this at 1)
    AoA  = 5                                                                        # Angle of attack [deg]
    alpha = AoA*(np.pi/180)                                                          # Angle of attack [rad]
    Uxy = np.array([Vinf*np.cos(alpha), Vinf*np.sin(alpha)])  # Freestream velocity vector [m/s]

    # --- External velocity function
    fU =lambda X, Y : (X*0+Uxy[0], X*0+Uxy[1])

    df = pd.read_csv(airfoil_file)
    XP, YP = df['x'].values, df['y'].values

    # --- Panel method
    out = CSPN_CVP1_solve(XP, YP, fU=fU, closed=True, verbose=True, alpha=alpha)
    CP = out['CP']
    Cp = out['Cp']
    n_hat = out['n']

    plot_pressure_force_bars(CP, Cp, n_hat, scale=0.1, ax=None)
    plot_Cp(CP, Cp)

    # --- Flow field
    PP  = np.column_stack((XP, YP))
    vel = lambda X, Y : CSPN_CVP1_u(X, Y, PP, out['sigmas'], -out['gamma']) # TODO gamma sign issue
    X, Y, U, V =  flowfield2D(vel, xmax=1.5, xmin=-0.5, ymin=-0.3, ymax=0.3, nx=140, ny=120, U0x=Uxy[0], U0y=Uxy[1], rel=False)
    #bIn = points_inside(XP, YP, X, Y)
    #U[bIn] = 0
    #V[bIn] = 0

    ax =  flowfield2D_plot(X, Y, U, V, ax=None, minVal=0, maxVal=2, bounded=False, rel=False)
    ax.fill(XP, YP, 'k')
    ax.set_aspect('equal')

    # --- Cp field
    Speed = np.sqrt(U**2 + V**2) 
    CpXY = 1 - Speed**2/Vinf**2 # (Vxy/Vinf)**2 
    ax =  flowfield2D_plot(X, Y, U, V, Speed=CpXY, ax=None, nLevels=31, minVal=-2, maxVal=1, bounded=False, rel=False, cmap='jet')
    ax.fill(XP,YP,'k')
    ax.set_aspect('equal')


    plt.show()
