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
    bOK=r2>1e-8
    U = np.zeros_like(DX)
    V = np.zeros_like(DX)
    if regMethod is None:
        U[bOK] = Gamma/(2*np.pi) * tX[bOK]/r2[bOK]   # [1] Eq 32.7
        V[bOK] = Gamma/(2*np.pi) * tY[bOK]/r2[bOK] 
    else:                  
        U[bOK] = Gamma/(2*np.pi) * tX[bOK]/r2[bOK] * (1 - np.exp(- r2[bOK] / regParam ** 2))
        V[bOK] = Gamma/(2*np.pi) * tY[bOK]/r2[bOK] * (1 - np.exp(- r2[bOK] / regParam ** 2))
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


def panel_solve_vps(XP, YP, Ux, Uy, lift=True, iTE=0, curv_method='Menger', backDiagCorr=True, verbose=False, GammaConvention='z'):
    """ 
    """

    # --- Geometry
    n_hat, t_hat, mids, ds, curvature, ax =  airfoil_params(XP, YP, plot=False, ntScale=0.3, curv_method=curv_method)

    if verbose:
        print('ds  ', ds)
        print('crv ', curvature)
        print('crv ', np.max(curvature), np.min(curvature))
#         curvature[0]  =curvature[1]
#         curvature[-1] =-100
        #print('cds2', curvature*ds*2)
    if GammaConvention=='z':
        sgn = 1
    else:
        sgn =-1

    # --- Right hand side
    # We implement the "Dirichlet" condition, no flow tangential
    rhs = -(Ux*t_hat[:,0] + Uy*t_hat[:,1])

    # --- Build matrix
    nCP = len(XP)-1
    nV  = len(XP)-1
    VP = mids # Position of vortices
    CP = mids
    M = np.zeros((nCP,nV))
    for j in range(nCP): # 
        for i in range(nV):
            if i==j:
                # Self weight and curvature
                # See Lewis
                M[i,i] = sgn*(0.5 + curvature[i]*ds[i]/(4*np.pi))
            else:
               rj = CP[j,:] # 
               ri = VP[i,:] # 
               dr = sgn * (rj-ri) 
               u,v = vp_u(dr[0], dr[1], Gamma=1, regParam=0, regMethod=None)
               M[j,i] = (u*t_hat[j,0] + v*t_hat[j,1]) * ds[i]

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

    Gammas = gammas*ds

    # --- Outputs
    out = {}
    # Geometry
    out['x']    = XP
    out['y']    = YP
    out['theta'] = np.arctan2(YP, XP)
    out['ds']   = ds
    out['n']    = n_hat
    out['t']    = t_hat
    out['CP']   = CP
    out['VP']   = mids
    out['curv'] = curvature
    out['theta_CP'] = np.arctan2(mids[:,1], mids[:,0])

    # ---
    out['Gammas']  = Gammas
    out['rhs']  = rhs

    # --- Output: Velocity at wall
    # TODO
    #Vtheta = gammas
    Vwall = np.zeros_like(mids)
    # NOTE: using the fact that tangent condition should be satisfied
    Vwall[:,0]= -sgn * gammas*t_hat[:,0] + 0*Ux
    Vwall[:,1]= -sgn * gammas*t_hat[:,1] + 0*Uy
    out['Vwall'] = Vwall
    out['Un'] = Vwall[:,0]*n_hat[:,0] + Vwall[:,1]*n_hat[:,1]
    out['Ut'] = Vwall[:,0]*t_hat[:,0] + Vwall[:,1]*t_hat[:,1]
    theta = out['theta_CP'] 
    r_hat = np.zeros_like(n_hat)
    th_hat = np.zeros_like(n_hat)
    r_hat[:,0] = np.cos(theta)
    r_hat[:,1] = np.sin(theta)
    th_hat[:,0] =-np.sin(theta)
    th_hat[:,1] = np.cos(theta)
    out['Ur']  = Vwall[:,0]*r_hat[:,0] + Vwall[:,1]*r_hat[:,1]
    out['Uth'] = Vwall[:,0]*th_hat[:,0] + Vwall[:,1]*th_hat[:,1]

    # --- Output: Cp
    out['Cp'] = 1-(Vwall[:,0]**2+Vwall[:,1]**2)/(Ux**2+Uy**2)


    # --- Output: Lift coefficient
    # TODO
    # LEwis Eq. 2.30

    # --- Output: Stream function
    #psi = Ux*y - Uy*x + 1/(2*np.pi) * sum(gamma_i ds_i ln(r_ii))

    return gammas, out

def VP_velocity(X, Y, Ux, Uy, VP, Gammas, regMethod=None, regParams=None):
    nV = len(Gammas)
    if regParams is None:
        regParams = [0]*nV
    U = np.zeros_like(X)
    V = np.zeros_like(Y)
    for i in range(nV):
        p = VP[i]
        u, v = vp_u(X-p[0], Y-p[1], Gamma=Gammas[i], regMethod=regMethod, regParam=regParams[i])
        U+=u
        V+=v
    U+=Ux
    V+=Uy
    return U, V


if __name__ == '__main__':
    import matplotlib.pyplot as plt
    case = 'file'
#     case = 'KT'
    case = 'cylinder'
#     case = 'ellipse_lift'
    alpha=0
    ge=None
    Cp_theory = None
    Uth_theory=None
    theta_mid=None
    curv_method='zero'
    curv_method='Menger'
#     curv_method='Lewis'
    backDiagCorr=False
    backDiagCorr=True
    coords_filename =None
    from welib.airfoils.karman_trefftz import KT_wall, KT_shape
    from welib.airfoils.shapes import AirfoilShape
    if case=='file':
        U0    = 1
        alpha = 5*np.pi/180
#         coords_filename = 'geom-KarmanTrefftz_less.csv'
#         coords_filename = 'geom-KarmanTrefftz.csv'
        coords_filename = 'KarmanTrefftz-300-NoScale.csv'
#         coords_filename = 'FFA-W3-301-coords.csv'
#         coords_filename = 'Diamond-coords.csv'
#         coords_filename = 'NACA0012-n399-coords.csv'
        df = pd.read_csv(coords_filename)
        XP = df['x'].values
        YP = df['y'].values
        m=len(XP)
        iTE=0
        XC      = -0.2# X-coord of cylinder center in Z-plane
        YC      = 0.1 # Y-coord of cylinder center in Z-plane
        tau_deg = 10
        l = 2 - tau_deg / 180 # l = 2 - tau/np.pi
        #XP2, YP2 = KT_shape(XC, YC, l=l, n=len(XP)) # ,Xg,Yg)
        XP2, YP2 = KT_shape(XC, YC, l=l, n=300) # ,Xg,Yg)
        arf = AirfoilShape(XP2, YP2)
        #import pdb; pdb.set_trace()
        XCp_th, YCp, Cp_theory, U_circ, V_circ = KT_wall(XC, YC, l=l, n=m-1, alpha=alpha) 

        lift=abs(alpha)>0
    elif case=='KT':
        m=300
        U0    = 1
        alpha = 5*np.pi/180
        XC      = -0.2# X-coord of cylinder center in Z-plane
        YC      = 0.1 # Y-coord of cylinder center in Z-plane
        tau_deg = 10
        l = 2 - tau_deg / 180 # l = 2 - tau/np.pi
        XP, YP = KT_shape(XC, YC, l=l, n=m) # ,Xg,Yg)
        XCp_th, YCp, Cp_theory, U_circ, V_circ = KT_wall(XC, YC, l=l, n=m-1, alpha=alpha) 
        iTE=0
        lift=abs(alpha)>0
    elif case=='cylinder':
        m=151
        U0    = 1
        R=2
        Gamma =-4*np.pi*U0*R*0.5
        theta   =-np.linspace(0,2*np.pi, m+1)
        theta_TE = np.arcsin(Gamma/(4*np.pi*U0*R))
        theta   += theta_TE
        iTE=0
        XP = R*np.cos(theta)
        YP = R*np.sin(theta)
        dtheta=theta[1]-theta[0]
        theta_mid   = theta[:-1]+dtheta/2
        Uth_theory = -2.0*U0*np.sin(theta_mid) + Gamma/(2*np.pi*R)# Utheta at CP
        ge         = -Uth_theory
        Cp_theory  = 1-(Uth_theory)**2/U0**2
        XCp_th  = (XP[0:-1]+XP[1:])/2
        lift=abs(Gamma)>0
    elif case=='ellipse_lift':
        m=300
        U0    = 1
        major=1
        ratio=0.05
        alpha = 30*np.pi/180
        theta   =-np.linspace(0,2*np.pi, m+1)
        abyr= np.sqrt((1.-ratio)/(1.+ratio))# ! see Lewis p 50
        XP = 0.5*major       *(np.cos(theta))
        YP = 0.5*major*ratio*  np.sin(theta)
        dtheta=theta[1]-theta[0]
        theta_mid   = theta[:-1]+dtheta/2
        f = np.sqrt(1.0+ (abyr**4) - 2.*(abyr**2)*np.cos(2.*theta_mid) )
        Uth_theory =-(-2.0*U0*np.sin(theta_mid-alpha)-2*np.sin(alpha))/f          # See Lewis p50
        ge         =-(-2.0*U0*np.sin(theta_mid-alpha)-2*np.sin(alpha))/f          # See Lewis p50
        iTE=0
        Cp_theory  = 1-(Uth_theory)**2/U0**2
        XCp_th  = (XP[0:-1]+XP[1:])/2
        lift=abs(alpha)>0

    Ux    = U0*np.cos(alpha)
    Uy    = U0*np.sin(alpha)

    gammas, out = panel_solve_vps(XP, YP, Ux, Uy, lift=lift, iTE=iTE, curv_method=curv_method, verbose=True, backDiagCorr=backDiagCorr)
    print('gammas:', gammas)
    print('>>> n', len(XP), case, 'alpha:',alpha*180/np.pi)
    ds_mean = np.mean(out['ds'])
    print('>>> ds mean', ds_mean)
    VP = out['VP']

    plot_airfoil(XP, YP, Uwall=out['Vwall'], ntScale=0.1, UScale=0.8)
    plot_airfoil(XP, YP, nt=True, ntScale=0.1, UScale=0.8)
    
#     print('MaxError ',np.max(np.abs(ge-gn)))
#     print('')
#     print(f'Case: {iCase} - Curv Method: {curv_method} - iCoord: {iCoord} - backCorr:{backCorr} - m:{m}')
    fig,axes = plt.subplots(1, 3, sharey=False, figsize=(12.8,4.8)) # (6.4,4.8)
    fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
    ax=axes[0]
    if ge is not None:
        ax.plot(VP[:,0], ge    , 'k-', label='Theory')
    ax.plot(out['VP'][:,0], gammas, '--', label='gammas')
    ax=axes[1]
    if theta_mid is not None:
        ax.plot(theta_mid, Uth_theory, 'k-', label='Utheta')
        ax.plot(theta_mid, out['Un'], '--', label='Un')
        ax.plot(theta_mid, out['Ur'], ':', label='Ur')
        ax.plot(theta_mid, out['Ut'],  '.', label='Ut')
        ax.plot(theta_mid, out['Uth'], '--', label='Utheta')

        ax.legend()

    ax=axes[2]
    if Cp_theory is not None:
        ax.plot(XCp_th, Cp_theory, 'k-', label='Theory', lw=2)
    ax.plot(out['VP'][:,0], out['Cp'], '--', label='Cp')
    ax.legend()

    if coords_filename is not None:
        if 'Karman' in coords_filename:
            M = np.column_stack((out['VP'][:,0], out['Cp']))
            df = pd.DataFrame(data=M,columns=['x','Cp'])
            df.to_csv('CP-KarmanTrefftz_5_num_py.csv',sep=',', index=False)
# 

    regParams = out['ds']*0.50
    regMethod = None
    regMethod = 1


    # --- CP num
    #Uw, Vw = VP_velocity(XP, YP, Ux, Uy, out['VP'], out['Gammas'], regMethod=regMethod, regParams=regParams); Xw=XP

    VP = out['VP']
    n_hat =out['n']
    ds  =out['ds']
    VP2 = VP + 1*(n_hat.T*ds).T 
    Uw, Vw = VP_velocity(VP2[:,0], VP2[:,1], Ux, Uy, out['VP'], out['Gammas'], regMethod=regMethod, regParams=regParams); Xw=VP2[:,0]; Xw=VP[:,0]

    Q = np.sqrt(Uw**2+Vw**2)
    CP2 = 1-(Q**2/U0**2)
    fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
    fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
    ax.plot(Xw, CP2, label='From Num eval')
    if Cp_theory is not None:
        ax.plot(XCp_th, Cp_theory, 'k-', label='Theory', lw=2)
    ax.plot(out['VP'][:,0], out['Cp'], '--', label='Cp from Gamma')
#     ax.plot(XP, YP, 'k-')
#     ax.plot(VP2[:,0], VP2[:,1], '+')
    ax.set_xlabel('')
    ax.set_ylabel('Cp')
    ax.legend()

    # --- Flow field
    xmax     = 3.5
    minSpeed = 0
    maxSpeed = 2.0
    nLevels  = 11
    rlim     = 2

    # --- Derived Parameters
    # Grid
    vg = np.linspace(-xmax, xmax, 100)
    vx = vg
    vy = vg
    X, Y = np.meshgrid(vx, vy)

    U, V   = VP_velocity(X, Y, Ux, Uy, out['VP'], out['Gammas'],regMethod=regMethod, regParams=regParams)

    # --- Plot velocity and streamlines from velocity field
    Speed = np.sqrt((U**2+V**2))/U0
    # Streamlines
    ys = np.linspace(-xmax, xmax, 15)
    start = np.array([ys*0-xmax, ys])
    # Plot
    fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
    im = ax.contourf(X, Y, Speed, levels=np.linspace(minSpeed, maxSpeed, nLevels), vmin=minSpeed, vmax=maxSpeed)
    cb=fig.colorbar(im)
    cb.set_label(r'$\sqrt{u^2+v^2}/U_\text{ref}$')
    sp = ax.streamplot(vx, vy, U, V, color='k',start_points=start.T, linewidth=0.7,density=30,arrowstyle='-')
    #sp = ax.streamplot(vx, vy, U, V, color='k', linewidth=0.7,density=30,arrowstyle='-')
    #qv = streamQuiver(ax, sp, spacing=1, offset=1.0, scale=40, angles='xy')
    #ax.plot(0, 0, 'ko',lw=3)
    ax.plot(XP, YP, 'k-',lw=3)
    ax.set_xlabel(r'$x$ [m]')
    ax.set_ylabel(r'$y$ [m]')
    ax.set_aspect('equal','box')
    ax.tick_params(direction='in', top=True, right=True)
    ax.set_title('Flow about a 2D point source - from velocity field')





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
