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

def CVP_solve(XP, YP, Uxy=None, fU=None, verbose=False, alpha=None, x_Cm_ref=0.25):
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

    CP = mids # Positions of control points
    nPanels = len(PP)-1

    # --- Build Influence matrix for vortex panels
    A_UnV = np.zeros((len(CP), nPanels))
    A_UtV = np.zeros((len(CP), nPanels))
    for i in range(len(CP)):
        for j in range(nPanels):
            #if (i == j): 
            #    #u, v = cvp_u11(CP[i,:], PP[j,:], PP[j+1,:], principal=i==j)
            #    #print('PVV',u*n_hat[i,0] + v*n_hat[i,1], u*t_hat[i,0] + v*t_hat[i,1])
            #    A_UnV[i,j] = 0.0  #  Principal value
            #    A_UtV[i,j] =-0.5  #  Principal value 
            #else:
            u, v = cvp_u11(CP[i,:], PP[j,:], PP[j+1,:], principal=i==j)
            A_UnV[i,j] = u*n_hat[i,0] + v*n_hat[i,1]
            A_UtV[i,j] = u*t_hat[i,0] + v*t_hat[i,1]


    # --- Build system matrix.
    # We use sources for most of it, and a sum of gammas for the last column
    M = A_UnV

    # --- Kutta condition
    # Katz method, replace a random row with a Kutta condition
    # Kutta condition: gamma0 and gammaN must have oppostite signs, e.g.:
    #      gamma0 + gammaN = 0   or  -gamma0 - gammaN = 0
    iKutta = nPanels//4-1
    M[iKutta,:] = 0
    M[iKutta, 0] = +1
    M[iKutta, -1] = +1

    # --- Right hand side
    Ux, Uy = fU(CP[:,0], CP[:,1])
    U0_n = (Ux*n_hat[:,0] + Uy*n_hat[:,1])
    U0_t = (Ux*t_hat[:,0] + Uy*t_hat[:,1])
    rhs = -U0_n
    rhs[iKutta] = 0

    # --- Solve
    gammas = np.linalg.solve(M, rhs)

    # Tangential velocities on panel
    UV_t = A_UtV @ gammas # TODO TODO TODO TODO TODO TODO TODO TODO sign
    # Smoothing 
    UV_t_smooth = UV_t.copy()
    UV_t_smooth[1:] = (UV_t[0:-1]+UV_t[1:])/2
    gammas_smooth = gammas.copy()
    gammas_smooth[1:] = (gammas[0:-1]+gammas[1:])/2
    Utot_t = U0_t + UV_t_smooth
    # Normal velocities on panel (residual check)
    UV_n =  A_UnV @ gammas
    Utot_n = U0_n + UV_n

    Vwall = Utot_n[:,np.newaxis]*n_hat + Utot_t[:,np.newaxis]*t_hat
    Vwall_smooth = Vwall.copy()
    Vwall_smooth[1:] = (Vwall[0:-1]+Vwall[1:])/2

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
    out['gammas']    = gammas
    # Output: Velocity at wall
    out['Vwall'] = Vwall_smooth
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
    out['Cl_KJ'] = -2*sum(gammas_smooth*ds)
    return out


def CVP_u(X, Y, PP, gammas, debug=False):
    U, V = cvp_u(X, Y, PP, gammas)
    return U, V


if __name__ == '__main__':
    from welib.CFD.flows2D import flowfield2D, flowfield2D_plot
    from welib.vortilib.panelcodes.panel_tools import plot_Cp, plot_pressure_force_bars, plot_airfoil
    from welib.vortilib.panelcodes.panel_examples import getCase
    # Set numpy print options for nicer array display
    np.set_printoptions(precision=4, linewidth=140) #, suppress=True)
    scriptDir = os.path.dirname(__file__)


    # --- Parameters
    cas = getCase('VonDeVooren_lift', mid_out=True, solver=None) #, m=90, U0=1, alpha=5, solver='LVP', mid_out=True)
    XP, YP = cas['XP'], cas['YP']
    Uxy = cas['Uxy']

    XP=XP[::-1]
    YP=YP[::-1]
    #plot_airfoil(XP, YP)
    #plt.show()

    # --- Panel method
    out = CVP_solve(XP, YP, Uxy=Uxy, alpha=cas['alpha'])
    CP = out['CP']
    Cp = out['Cp']
    n_hat = out['n']

    print('Cl   :',  out['Cl'])
    print('Cl_KJ:',  out['Cl_KJ'])

    #plot_pressure_force_bars(CP, Cp, n_hat, scale=0.1, ax=None)
    ax = plot_Cp(CP[:,0], Cp, Cp_ref=cas['Cp'])
    #ax.set_ylim([-1.8, 1])

    # --- Flow field
    PP  = np.column_stack((XP, YP))
    vel = lambda X, Y : CVP_u(X, Y, PP, out['gammas'])
    X, Y, U, V =  flowfield2D(vel, xmax=2.0, xmin=-1.5, ymin=-0.6, ymax=0.6, nx=81, ny=50, U0x=Uxy[0], U0y=Uxy[1], rel=False)
    #bIn = points_inside(XP, YP, X, Y)
    #U[bIn] = 0
    #V[bIn] = 0

    ax =  flowfield2D_plot(X, Y, U, V, ax=None, minVal=0, maxVal=2, bounded=False, rel=False)
    ax.fill(XP, YP, 'k')
    ax.set_aspect('equal')

    # --- Cp field
    Speed = np.sqrt(U**2 + V**2) 
    CpXY = 1 - Speed**2/cas['U0']**2 # (Vxy/Vinf)**2 
    ax =  flowfield2D_plot(X, Y, U, V, Speed=CpXY, ax=None, nLevels=31, minVal=-2, maxVal=1, bounded=False, rel=False, cmap='jet')
    ax.fill(XP,YP,'k')
    ax.set_aspect('equal')


    if False:
        # Compare with Katz
        A2   = np.loadtxt("VortexCode2D_Katz/cvortex_A.DAT")
        B2   = np.loadtxt("VortexCode2D_Katz/cvortex_B.DAT")
        rhs2 = np.loadtxt("VortexCode2D_Katz/cvortex_RHS.DAT")
        g2   = np.loadtxt("VortexCode2D_Katz/cvortex_solution.DAT")
        v2   = np.loadtxt("VortexCode2D_Katz/cvortex_VEL.DAT")
        UVt2 = B2 @ g2
        B    = A_UtV
        print(A2.shape, M.shape)
        print(B2.shape, B.shape)
        print(rhs2.shape, rhs.shape)
        print(g2.shape, gammas.shape)

        fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8))
        ax.plot(rhs , '-', label='rhs')
        ax.plot(rhs2 , '--', label='rhs2')
        ax.legend()
        fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8))
        ax.plot(gammas , '-', label='gammas')
        ax.plot(- g2 , '--', label='minus g2')
        ax.plot(gammas_smooth , 'k-', label='gammas_smooth')
        ax.legend()

        fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8))
        ax.plot(Utot_t , '-', label='Utot_t')
        ax.plot(v2, '--', label='v2')
        ax.legend()

        fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8))
        ax.plot(UV_t , '-', label='UV_t')
        ax.plot(UVt2, '--', label='UVt2')
        ax.legend()





    plt.show()
