"""
Setup:
 - n panels consisting of constant doublet
 - 1 wake panel extending to infinity

Unknowns (n+1):
  -  +n unknowns mu for the doublet panels
  -  +1 unknown  mu for the wake panel

Equations (n+1)
 - n no flow through conditions at the control points accoutning for the relationship between mu_wake 
 - 1 Kutta condition links mu_1 mu_n and mu_wake

NOTE: 
  - mu_wake can be elliminated if needed by direct manipulation of the system matrices

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
from welib.vortilib.elements.DoubletPanel2D import *

def CDP_solve(XP, YP, Uxy=None, fU=None, verbose=False, alpha=None, x_Cm_ref=0.25):
    """
    Solve panel method for a given geometry and external velocity
     N constant doublet panels
     1 doublet panel in the wake

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
    nPanels = len(PP)-1 # NOTE: we will have one more wake panel

    # --- Wake panel
    WP1 = (PP[0,:] + PP[-1,:])/2
    WP2 = [100, 0] # Far away
    WP = np.vstack((WP1, WP2)) # Wake panel geometry

    # TEMP
    theta    = np.arctan2(dP[:,1], dP[:,0])

    # --- Build Influence matrix for vortex panels
    # NOTE: KatzPlotkin use a different sign convention. Also, they use 1/2pi=0.15916 which makes the system less ill-conditioned
    A_UnV = np.zeros((len(CP), nPanels+1))
    A_UtV = np.zeros((len(CP), nPanels+1))
    A_UnKP = np.zeros((len(CP), nPanels+1))
    A_UtKP = np.zeros((len(CP), nPanels+1))
    for i in range(len(CP)):
        for j in range(nPanels):
            if (i == j): 
                u, v = cdp_u11(CP[i,:], PP[j,:], PP[j+1,:], principal=i==j)
                #print('PVV',u*n_hat[i,0] + v*n_hat[i,1], u*t_hat[i,0] + v*t_hat[i,1], 2/(np.pi*ds[j]))
                A_UnV[i,j] = 2/(np.pi*ds[j]) #  Principal value
                A_UtV[i,j] = 0.0  #  Principal value TODO, need gradient
            else:
                u, v = cdp_u11(CP[i,:], PP[j,:], PP[j+1,:], principal=i==j)
            A_UnV[i,j] = u*n_hat[i,0] + v*n_hat[i,1]
            A_UtV[i,j] = u*t_hat[i,0] + v*t_hat[i,1]

            u_kp, w_kp = cdp_u11_kp_raw(CP[i,:], PP[j,:], PP[j+1,:], principal=i==j)
            A_UnKP[i,j] = -u_kp*np.sin(theta[i]) + w_kp*np.cos(theta[i])
            A_UtKP[i,j] =  u_kp*np.cos(theta[i]) + w_kp*np.sin(theta[i])

        # --- Wake panel
        j = nPanels
        u, v = cdp_u11(CP[i,:], WP1, WP2)
        A_UnV[i,j] = u*n_hat[i,0] + v*n_hat[i,1]
        A_UtV[i,j] = u*t_hat[i,0] + v*t_hat[i,1]

        # --- Wake panel KP (vortex point)
        R = np.sqrt((CP[i,0] - WP1[0])**2 + (CP[i,1] - WP1[1])**2)
        U = - 1/(2*np.pi)  * (CP[i,1] / (R**2))
        W =   1/(2*np.pi)  * (CP[i,0] - WP1[0]) / (R**2)
        A_UnKP[i,j]  = -U * np.sin(theta[i]) + W * np.cos(theta[i])
        A_UtKP[i,j]  =  U * np.cos(theta[i]) + W * np.sin(theta[i])
        #A_UnV[i,j]   = -U * np.sin(theta[i]) + W * np.cos(theta[i]) 
        #A_UtV[i,j]   =  U * np.cos(theta[i]) + W * np.sin(theta[i])

    # --- Build system matrix.
    # We use sources for most of it, and a sum of gammas for the last column
    M = np.zeros((len(CP)+1, nPanels+1))
    #M[:-1,:] = A_UnKP
    M[:-1,:] = A_UnV

    # --- Kutta condition
    # Kutta condition: mu1 - muN + mu_wake = 0
    iKutta = nPanels
    M[iKutta,:] = 0
    M[iKutta, 0] = +1
    M[iKutta, -2] = -1
    M[iKutta, -1] = +1

    # --- Right hand side
    Ux, Uy = fU(CP[:,0], CP[:,1])
    U0_n = (Ux*n_hat[:,0] + Uy*n_hat[:,1])
    U0_t = (Ux*t_hat[:,0] + Uy*t_hat[:,1])
    rhs = np.zeros(nPanels+1)
    rhs[:-1] = -U0_n
    rhs[iKutta] = 0

    # --- Solve
    # NOTE: With constant doublet panels + Neumann (no‐through) BC we have a rank-1 singular system.
    #  - The no-through equations depend only on jumps of μ between adjacent panels (since a constant μ on a panel has zero tangential derivative).
    #   Hence adding a constant offset to all μ leaves every no-through equation unchanged.
    # - The Kutta relation, μ1 − μN + μW = 0, also involves only differences, so it doesn't remove that nullspace.
    # Potential fix:
    # - Set one of the mu to zero (removing one column and one row)
    # - Set sum of mus to zero (adding one row, overdetermined system)
    # - Use Dirichlet formulation of the potential (phi=cst on body)
    # - Use Sources

    #mus = np.linalg.solve(M, rhs)
    # Add a row to the system matrix so that sum of mus is zero
    #M = np.vstack([M, np.ones(M.shape[1])])
    #M[-1,-1] = 0
    #rhs = np.append(rhs, 0)  # Add a zero to the rhs
    mus = np.linalg.lstsq(M, rhs, rcond=None)[0]
    #print('Mean mus: ', np.mean(mus[:-1]))
    #expected_mean = -3.78860454883324
    #mus = mus - np.mean(mus[:-1]) + expected_mean
    #A2   = np.loadtxt("VortexCode2D_Katz/cdoublets_A.DAT")
    #B2   = np.loadtxt("VortexCode2D_Katz/cdoublets_B.DAT")
    #rhs2 = np.loadtxt("VortexCode2D_Katz/cdoublets_RHS.DAT")
    #mus2 = np.loadtxt("VortexCode2D_Katz/cdoublets_solution.DAT")
    #Cp2  = np.loadtxt("VortexCode2D_Katz/Cp_cdoublets.csv", delimiter=',', skiprows=1)
    # mus = mus2




    # --- Principal value of tangential velocity due to doublet distribution
    mup = mus[:-1] # mu on panels
    R         = np.zeros(len(CP))
    UV_t_PV   = np.zeros(len(CP))
    # Interior points - Central difference around control point j
    R[1:-1]       = np.sqrt((CP[2:,0] - CP[:-2,0])**2 + (CP[2:,1] - CP[:-2,1])**2)
    UV_t_PV[1:-1] = (mup[2:] - mup[:-2]) / R[1:-1]
    # First point - Forward difference
    R[0]       = np.sqrt((CP[1,0] - CP[0,0])**2 + (CP[1,1] - CP[0,1])**2)
    UV_t_PV[0] = (mup[1] - mup[0]) / R[0]
    # Last point - Backward difference
    R[-1]       = np.sqrt((CP[-1,0] - CP[-2,0])**2 + (CP[-1,1] - CP[-2,1])**2)
    UV_t_PV[-1] = (mup[-1] - mup[-2]) / R[-1]
    # Delta Gamma
    dGamma = (mup[1:] - mup[:-1])
    Gamma_TE = mus[0]-mus[-2]+ mus[-1] # Kutta condition
    Gamma_tot = np.sum(dGamma)
    #print('Gamma TE: ', Gamma_TE, 'Gamma tot: ', Gamma_tot, mus[-1])
    #print('dGamma: ', dGamma)

    # Tangential velocities on panel
    UV_t = A_UtV @ mus
    Utot_t = U0_t + UV_t - UV_t_PV/2
    # Normal velocities on panel (residual check)
    UV_n =  A_UnV @ mus
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
    out['ds']       = ds
    out['n']        = n_hat
    out['t']        = t_hat
    out['CP']       = CP
    out['theta_CP'] = np.arctan2(CP[:,1], CP[:,0])
    out['WP']       = WP # wake panel
    # Output: System
    out['rhs']      = rhs
    out['M']        = M
    # Output: Solution
    out['mus']    = mus
    # Output: Velocity at wall
    out['Vwall'] = Vwall
    out['Un'] = Utot_n
    out['Ut'] = Utot_t
    # Output: Cp
    out['Cp'] = Cp #1 - (Vwall[:,0]**2 + Vwall[:,1]**2)/(Ux**2+Uy**2)
    # Output: Loads
    chord = max(XP)-min(XP)
    Vref = np.mean(np.sqrt(Vinf2)) # TODO
    # TODO TODO TODO Missing some chord and Uref - Do Dimension loads first
    if alpha is not None:
        # angle of panel normal w.r.t. horizontal and include AoA
        delta                = phi + (np.pi/2) # Angle from x-axis to normal vector
        beta                 = delta - alpha   # Angle between freestream and normal vector
        beta[beta > 2*np.pi] = beta[beta > 2*np.pi] - 2*np.pi
        out['Cn'] = -Cp*ds*np.sin(beta)/chord # Normal to chord
        out['Cc'] = -Cp*ds*np.cos(beta)/chord # Along chord
        out['Cl'] = sum(out['Cn']*np.cos(alpha)) - sum(out['Cc']*np.sin(alpha))
    out['Cx'] = sum(Cp*ds*n_hat[:,0])/chord
    out['Cy'] = sum(Cp*ds*n_hat[:,1])/chord
    out['Cm'] = sum(Cp*(CP[:,0]-x_Cm_ref)*ds*np.cos(phi))
    out['Cl_KJ'] = -2*Gamma_tot / (Vref * chord)
    return out


def CDP_u(X, Y, PP, WP, mus, debug=False):

    U, V = cdp_u(X, Y, PP, mus[:-1])

    Uw, Vw = cdp_u(X, Y, WP, [mus[-1]])
    U += Uw
    V += Vw 
    return U, V


if __name__ == '__main__':
    from welib.CFD.flows2D import flowfield2D, flowfield2D_plot
    from welib.vortilib.panelcodes.panel_tools import plot_Cp, plot_pressure_force_bars, plot_airfoil
    from welib.vortilib.panelcodes.panel_examples import getCase
    # Set numpy print options for nicer array display
    np.set_printoptions(precision=4, linewidth=140) #, suppress=True)
    scriptDir = os.path.dirname(__file__)

    # --- Parameters
    #cas = getCase('NACA2412_lift', mid_out=True, solver=None); xref = cas['x_ref'] #, m=90, U0=1, alpha=5, solver='LVP', mid_out=True)
    #XP, YP = cas['XP'], cas['YP']
    #XP = np.concatenate(([1.04], XP, [1.04]))
    #YP = np.concatenate(([0.00], YP, [0.00]))

    cas = getCase('VonDeVooren_lift', mid_out=True, solver=None); xref=None; #, m=90, U0=1, alpha=5, solver='LVP', mid_out=True)
    XP, YP = cas['XP'], cas['YP']
    XP=XP[::-1]
    YP=YP[::-1]

    Uxy = cas['Uxy']
    #plot_airfoil(XP, YP)

    # --- Panel method
    out = CDP_solve(XP, YP, Uxy=Uxy, alpha=cas['alpha'])
    CP = out['CP']
    Cp = out['Cp']
    n_hat = out['n']

    print('Cl   :',  out['Cl'])
    print('Cl_KJ:',  out['Cl_KJ'])

    #plot_pressure_force_bars(CP, Cp, n_hat, scale=0.1, ax=None)
    ax = plot_Cp(CP[:,0], Cp, Cp_ref=cas['Cp'], x_ref=xref)
    #ax.set_ylim([-1.8, 1])

    # --- Flow field
    PP  = np.column_stack((XP, YP))
    vel = lambda X, Y : CDP_u(X, Y, PP, out['WP'],out['mus'])
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
        #A2   = np.loadtxt("VortexCode2D_Katz/cdoublets_A.DAT")
        #B2   = np.loadtxt("VortexCode2D_Katz/cdoublets_B.DAT")
        #rhs2 = np.loadtxt("VortexCode2D_Katz/cdoublets_RHS.DAT")
        #mus2 = np.loadtxt("VortexCode2D_Katz/cdoublets_solution.DAT")
        #Cp2  = np.loadtxt("VortexCode2D_Katz/Cp_cdoublets.csv", delimiter=',', skiprows=1)
        # Compare with Katz
        fig,axes = plt.subplots(3, 1, sharey=False, figsize=(6.4,4.8))
        ax=axes[0]
        ax.plot(mus , '-', label='mus')
        ax.plot(mus2, '--', label='mus2')
        ax.legend()
        ax=axes[1]
        Vt2 = B2 @ mus2
        ax.plot(UV_t , '-', label='VT')
        ax.plot(Vt2, '--', label='VT2')
        ax.legend()
        ax=axes[2]
        ax.plot(Cp , '-', label='Cp')
        ax.plot(Cp2[:,1], '--', label='Cp2')

        #import pdb; pdb.set_trace() 


    plt.show()
