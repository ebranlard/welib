""" 
Implement a source point panel method.

"""
import numpy as np
import pandas as pd
from welib.tools.clean_exceptions import *
from welib.vortilib.panelcodes.panel_tools import airfoil_params, plot_airfoil

print('>>>> SOURCE POINTS IS UNFINISHED')

# --------------------------------------------------------------------------------}
# --- 2d point source 
# --------------------------------------------------------------------------------{
def source_flow(DX, DY, Sigma=1): 
    """ 
    Induced velocity by one 2D source point on multiple Control Points (CP)
    INPUTS:
    - CP: position of control point
    - Ps: position of source

    EXAMPLE:
        # [Ui]= sp2d_u(1, 0, 2*np.pi)
        # [Ui]= sp2d_u(0, 1, 2*np.pi)
    """
    r2 = DX ** 2 + DY ** 2
    rX = DX
    rY = DY
    bOK=r2>1e-5
    tU = np.zeros_like(DX)
    V = np.zeros_like(DX)
    U = np.zeros_like(DX)
    V = np.zeros_like(DX)
    U[bOK] = Sigma/(2*np.pi) * rX[bOK] / r2[bOK]   # [1] Eq. 32.3
    V[bOK] = Sigma/(2*np.pi) * rY[bOK] / r2[bOK] 
#     U[~bOK] = 0.5*Sigma* rX[~bOK] / r2[~bOK]
#     V[~bOK] = 0.5*Sigma* rY[~bOK] / r2[~bOK]
    return U, V

def SP_panel_solve(XP, YP, Ux, Uy, lift=True, verbose=False):
    """ 
    """
    # --- Geometry
    n_hat, t_hat, mids, ds, _, ax =  airfoil_params(XP, YP, plot=False, ntScale=0.3, curv_method='zero')

    if verbose:
        print('ds  ', ds)
        #print('cds2', curvature*ds*2)

    # --- Right hand side
    # We implement the "Dirichlet" condition, no flow tangential
    #rhs = -(Ux*t_hat[:,0] + Uy*t_hat[:,1])
    rhs = -(Ux*n_hat[:,0] + Uy*n_hat[:,1])

    # --- Build matrix
    nCP = len(XP)-1
    nS  = len(XP)-1
    SP = mids # Position of sources
    SP = mids - 0.1*(n_hat.T*ds).T #-0.05*t_hat
    CP = mids
#     CP = np.column_stack( (XP[:-1], YP[:-1]) ) # Positions of control points
    M = np.zeros((nCP,nS))
    for i in range(nS):
        for j in range(nCP): # 
#             if i==j:
#                 M[i,i] = 0.5 *ds[i]
#             else:
# #               if True:
               ri = SP[i,:] # 
               rj = CP[j,:] # 
               dr = rj-ri
               u,v = source_flow(dr[0], dr[1], Sigma=1 ) # ds[i])
               M[j,i] = (u*n_hat[j,0] + v*n_hat[j,1]) # *ds[i]

    # --- SOLVE
    Sigmas = np.linalg.solve(M, rhs)
#     Sigmas = Sigmas*ds

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
    out['SP']   = SP
#     out['curv'] = curvature
    out['theta_CP'] = np.arctan2(CP[:,1], CP[:,0])
    out['rhs']  = rhs

#     # --- Output: Velocity at wall
#     # TODO
#     #Vtheta = gammas
    Vwall = np.zeros_like(mids)
    Vwall[:,0], Vwall[:,1] = SP_velocity(CP[:,0], CP[:,1], Ux, Uy, SP, Sigmas)
#     # NOTE: using the fact that tangent condition should be satisfied
#     Vwall[:,0]=gammas*tangents[:,0] + 0*Ux
#     Vwall[:,1]=gammas*tangents[:,1] + 0*Uy
    out['Vwall'] = Vwall
    out['Un'] = Vwall[:,0]*n_hat[:,0] + Vwall[:,1]*n_hat[:,1]
    out['Ut'] = Vwall[:,0]*t_hat[:,0] + Vwall[:,1]*t_hat[:,1]
    # --- Output: Cp
    out['Cp'] = 1-(Vwall[:,0]**2+Vwall[:,1]**2)/(Ux**2+Uy**2)
#     # --- Output: Lift coefficient
    # --- Output: Stream function
    #psi = Ux*y - Uy*x + 1/(2*np.pi) * sum(gamma_i ds_i ln(r_ii))

    return Sigmas, out

def SP_velocity(X, Y, Ux, Uy, SP, Sigmas):
    nS = len(Sigmas)
    U = np.zeros_like(X)
    V = np.zeros_like(Y)
    for i in range(nS):
        sp = SP[i]
        us, vs = source_flow(X-sp[0], Y-sp[1], Sigma=Sigmas[i])
        U+=us
        V+=vs
    U+=Ux
    V+=Uy
    return U, V



if __name__ == '__main__':
    import matplotlib.pyplot as plt
    case = 'cylinder'
#     case = 'ellipse_lift'
    alpha=0
    if case=='cylinder':
        m = 30
        U0    = 1
        R=2
        theta   =-np.linspace(0,2*np.pi, m+1)
        XP = R*np.cos(theta)
        YP = R*np.sin(theta)
        dtheta=theta[1]-theta[0]
        theta_mid   = theta[:-1]+dtheta/2
        ge =  2.0*U0*np.sin(theta_mid) # Gamma_exact = Utheta at CP
    elif case=='ellipse_lift':
        m=30
        U0    = 1
        major=1
        ratio=0.05
        lift=True
        alpha = 0*np.pi/180
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

    Sigmas, out = SP_panel_solve(XP, YP, Ux, Uy, verbose=True)
    print(Sigmas)
    SP = out['SP']
# 
    ax = plot_airfoil(XP, YP, Uwall=out['Vwall'], ntScale=0.1, UScale=0.8)
    ax.plot(out['SP'][:,0], out['SP'][:,1], '+')
#     plot_airfoil(XP, YP, nt=True, ntScale=0.1, UScale=0.8)
#     
# #     print('MaxError ',np.max(np.abs(ge-gn)))
# #     print('')
# #     print(f'Case: {iCase} - Curv Method: {curv_method} - iCoord: {iCoord} - backCorr:{backCorr} - m:{m}')
#     import matplotlib.pyplot as plt
    fig,axes = plt.subplots(1, 2, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
    fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
    ax=axes[0]
#     if ge is not None:
#         ax.plot(SP[:,0], ge    , 'k-', label='Theory')
#     ax.plot(out['SP'][:,0], Sigmas, '--', label='Sigmas')
    ax.plot(theta_mid, out['Un'], '--', label='Un')
    ax.plot(theta_mid, out['Ut'], '--', label='Ut')
    ax.plot(theta_mid, ge, 'k+', label='Ut')
    ax.legend()
    ax=axes[1]
    if ge is not None:
        ax.plot(out['CP'][:,0], 1-ge**2/U0**2, 'k+', label='Theory')
    ax.plot(out['CP'][:,0], out['Cp'], '--', label='Cp')
    ax.legend()
# 
#     if coords_filename is not None:
#         if 'Karman' in coords_filename:
#             M = np.column_stack((out['VP'][:,0], out['Cp']))
#             df = pd.DataFrame(data=M,columns=['x','Cp'])
#             df.to_csv('CP-KarmanTrefftz_5_num_py.csv',sep=',', index=False)


    plt.show()
