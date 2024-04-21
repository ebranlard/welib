""" 
Implement a vortex point panel method.

"""
import numpy as np
import pandas as pd
from welib.tools.clean_exceptions import *
from welib.vortilib.panelcodes.panel_tools import airfoil_params, plot_airfoil
# TODO look at welib.airfoils as well.. 


def vp_u(DX, DY, Gamma=1, regParam=0, regMethod=None): 
    """ 
    Induced velocity by one 2D vortex point on one Control Point (CP)
    CP: position of control point
    Pv: position of vortex
    """
    r2 = DX ** 2 + DY ** 2
    tX = -DY
    tY =  DX
    if regMethod is None:
        U = Gamma/(2*np.pi) * tX/r2   # [1] Eq 32.7
        V = Gamma/(2*np.pi) * tY/r2 
    else:                  
        U = Gamma/(2*np.pi) * tX/r2 * (1 - np.exp(- r2 / regParam ** 2))
        V = Gamma/(2*np.pi) * tY/r2 * (1 - np.exp(- r2 / regParam ** 2))
    return U,V


def kutta(M, rhs, iTE=None, verbose=False):
    M2 = M.copy()
    rhs2 = rhs.copy()
    if iTE==0:
        if verbose:
            print('Kutta with iTE=0')
        M2[:,-1] -= M2[:,0]
        M2[:,0:-1] = M2[:,1:]
        M2[:,-1]=-9

        M2[-1,:] -= M2[0,:]
        M2[0:-1,:] = M2[1:,:]
        M2[-1,:]=-9

        rhs2[-1]  -= rhs2[0]  # Subtracting row iTE+1 to col iTE
        rhs2[0:-1] = rhs2[1:] # Shifting indices
        rhs2[-1]   = -9
    elif iTE==-1:
        if verbose:
            print('Kutta with iTE=-1')
        M2[:,0] -= M2[:,-1]
        M2[:,-1]=-9

        M2[-1,:] -= M2[0,:]
        M2[-1,:]=-9

        rhs2[-1]  -= rhs2[0]  # Subtracting row iTE+1 to col iTE
        rhs2[-1]   = -9
    else:
        if verbose:
            print(f'Kutta with iTE={iTE}')
        M2[:,iTE]    -=M2[:,iTE+1]  # Subtracting col iTE+1 to col iTE
        M2[:,iTE+1:-1]=M2[:,iTE+2:] # Shifting indices 
        M2[:,-1]=-9

        M2[iTE,:]     -= M2[iTE+1,:]  # Subtracting row iTE+1 to col iTE
        M2[iTE+1:-1,:] = M2[iTE+2:,:] # Shifting indices
        M2[-1,:]       = -9

        rhs2[iTE]     -= rhs2[iTE+1]  # Subtracting row iTE+1 to col iTE
        rhs2[iTE+1:-1] = rhs2[iTE+2:] # Shifting indices
        rhs2[-1]     = -9

    M2   = M2[:-1,:-1]
    rhs2 = rhs2[:-1]
    return M2, rhs2

def kutta_revert(sol, iTE=None, start='LE'):
    """" revert the Kutta condition """
    sol2 = np.zeros(len(sol)+1)
    if iTE==0:
        sol2[-1] = -sol[0]
        sol2[:-1]= sol
    elif iTE==-1:
        sol2[-1]= -sol[0]
        sol2[:-1] = sol
    else:
        sol2[:iTE+1] =  sol[:iTE+1]
        sol2[iTE+1]  = -sol[iTE]
        sol2[iTE+2:] = sol[iTE+1:]
    return sol2

def backDiagonalCorrection(M, ds):
    """ Apply backdiagonal correction to enforce zero internal circulation"""
    m = M.shape[0]
    for i in range(m):
        vsum=0
        for j in range(m):
            if j!=m-i-1:
                vsum+= M[j,i]*ds[j]
        M[m-i-1, i] = -vsum/ ds[m-i-1]
    return M


def panel_solve_vps(XP, YP, Ux, Uy, lift=True, iTE=0, curv_method='Menger', backDiagCorr=True, verbose=False):
    """ 
    """

    # --- Geometry
    normals, tangents, mids, ds, curvature, ax =  airfoil_params(XP, YP, plot=False, ntScale=0.3, curv_method=curv_method)

    if verbose:
        print('ds  ', ds)
        print('crv ', curvature)
        #print('cds2', curvature*ds*2)

    # --- Right hand side
    # We implement the "Dirichlet" condition, no flow tangential
    rhs = -(Ux*tangents[:,0] + Uy*tangents[:,1])

    # --- Build matrix
    m = len(XP)-1
    M = np.zeros((m,m))
    for i in range(m):
        for j in range(m): # 
            if i==j:
                # Self weight and curvature
                # See Lewis
                M[i,i] = -0.5 - curvature[i]*ds[i]/(4*np.pi)
            else:
               ri=mids[i,:] # 
               rj=mids[j,:] # 
               dr = ri-rj
               u,v = vp_u(dr[0], dr[1], Gamma=1, regParam=0, regMethod=None)
               M[j,i] = (u*tangents[j,0] + v*tangents[j,1]) * ds[i]

    # --- Back diagonal correction
    # See Lewis
    if lift:
        if backDiagCorr:
            M = backDiagonalCorrection(M, ds)

    # --- KUTTA condition
    if lift: 
        M, rhs = kutta(M, rhs, iTE)

    # --- SOLVE
    gammas_r = np.linalg.solve(M, rhs)

    if lift: 
        gammas = kutta_revert(gammas_r, iTE=iTE)
        rhs    = kutta_revert(rhs, iTE=iTE)
    else:
        gammas = gammas_r


    # --- Outputs
    out = {}
    # Geometry
    out['x']    = XP
    out['y']    = YP
    out['theta'] = np.arctan2(YP, XP)
    out['ds']   = ds
    out['n']    = normals
    out['t']    = tangents
    out['CP']   = mids
    out['VP']   = mids
    out['curv'] = curvature
    out['theta_CP'] = np.arctan2(mids[:,1], mids[:,0])
    out['rhs']  = rhs

    # --- Output: Velocity at wall
    # TODO
    #Vtheta = gammas
    Vwall = np.zeros_like(mids)
    # NOTE: using the fact that tangent condition should be satisfied
    Vwall[:,0]=gammas*tangents[:,0] + 0*Ux
    Vwall[:,1]=gammas*tangents[:,1] + 0*Uy
    out['Vwall'] = Vwall

    # --- Output: Cp
    out['Cp'] = 1-(Vwall[:,0]**2+Vwall[:,1]**2)/(Ux**2+Uy**2)


    # --- Output: Lift coefficient
    # TODO
    # LEwis Eq. 2.30

    # --- Output: Stream function
    #psi = Ux*y - Uy*x + 1/(2*np.pi) * sum(gamma_i ds_i ln(r_ii))

    return gammas, out



if __name__ == '__main__':
    case = 'file'
#     case = 'cylinder'
#     case = 'ellipse_lift'
    alpha=0
    ge=None
    curv_method='zero'
#     curv_method='Menger'
    backDiagCorr=False
#     backDiagCorr=True
    coords_filename =None
    if case=='file':
        U0    = 1
        alpha = 9*np.pi/180
#         coords_filename = 'FFA-W3-301-coords.csv'
        coords_filename = 'geom-KarmanTrefftz_less.csv'
#         coords_filename = 'geom-KarmanTrefftz.csv'
#         coords_filename = 'Diamond-coords.csv'
        coords_filename = 'NACA0012-n399-coords.csv'
        # --- Read data
        df = pd.read_csv(coords_filename)
        XP = df['x'].values
        YP = df['y'].values
        lift=True
    elif case=='cylinder':
        m=21
        U0    = 1
        R=2
        theta   =-np.linspace(0,2*np.pi, m+1)
        XP = R*np.cos(theta)
        YP = R*np.sin(theta)
        dtheta=theta[1]-theta[0]
        theta_mid   = theta[:-1]+dtheta/2
        ge =  2.0*U0*np.sin(theta_mid) # Gamma_exact = Utheta at CP
        lift=False
    elif case=='ellipse_lift':
        m=300
        U0    = 1
        major=1
        ratio=0.05
        lift=True
        alpha = 30*np.pi/180
        theta   =-np.linspace(0,2*np.pi, m+1)
        abyr= np.sqrt((1.-ratio)/(1.+ratio))# ! see Lewis p 50
        XP = 0.5*major       *(np.cos(theta))
        YP = 0.5*major*ratio*  np.sin(theta)
        dtheta=theta[1]-theta[0]
        theta_mid   = theta[:-1]+dtheta/2
        f = np.sqrt(1.0+ (abyr**4) - 2.*(abyr**2)*np.cos(2.*theta_mid) )
        ge =  (2.0*U0*np.sin(theta_mid-alpha)+2*np.sin(alpha))/f          # See Lewis p50


    Ux    = U0*np.cos(alpha)
    Uy    = U0*np.sin(alpha)

    gammas, out = panel_solve_vps(XP, YP, Ux, Uy, lift=lift, iTE=0, curv_method=curv_method, verbose=True, backDiagCorr=backDiagCorr)
    print(gammas)
    VP = out['VP']

    #plot_airfoil(XP, YP, Uwall=out['Vwall'], ntScale=0.1, UScale=0.8)
    plot_airfoil(XP, YP, nt=True, ntScale=0.1, UScale=0.8)
    
#     print('MaxError ',np.max(np.abs(ge-gn)))
#     print('')
#     print(f'Case: {iCase} - Curv Method: {curv_method} - iCoord: {iCoord} - backCorr:{backCorr} - m:{m}')
    import matplotlib.pyplot as plt
    fig,axes = plt.subplots(1, 2, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
    fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
    ax=axes[0]
    if ge is not None:
        ax.plot(VP[:,0], ge    , 'k-', label='Theory')
    ax.plot(out['VP'][:,0], gammas, '--', label='gammas')
    ax=axes[1]
    if ge is not None:
        ax.plot(out['VP'][:,0], 1-ge**2/U0**2, '--', label='Theory')
    ax.plot(out['VP'][:,0], out['Cp'], '--', label='Cp')
    ax.legend()

    if coords_filename is not None:
        if 'Karman' in coords_filename:
            M = np.column_stack((out['VP'][:,0], out['Cp']))
            df = pd.DataFrame(data=M,columns=['x','Cp'])
            df.to_csv('CP-KarmanTrefftz_5_num_py.csv',sep=',', index=False)


#     ax.plot(theta_mid, gammas, '--', label='gammas')
#     ax.plot(theta_mid, ge    , '-', label='gammas')
#     ax.set_xlabel('phi')
#     ax=axes[1]
#     ax.plot(theta_mid*180/np.pi,rhs, label='rhs')
#     ax.plot(theta_mid*180/np.pi,ge , label='ge')
#     ax.plot(theta_mid*180/np.pi,gn ,'--', label='gn')
#     ax.set_xlabel('theta')
#     ax.set_ylabel('')
#     ax.legend()
#     # plt.show()


    plt.show()
