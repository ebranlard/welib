import numpy as np
from welib.essentials import *
# from welib.airfoils.karman_trefftz import *
    
def vc_u(DX, DY, U0=1, R=1, Gamma=0, alpha=0, internal_flow=True):
    """
    2D Cylinder with possible free stream and circulation
    """
    if alpha !=0:
        print('[WARN] cylinder not ready for students')
    r     = np.sqrt(DX**2 + DY**2) # Radius at mesh locations
    theta = np.arctan2(DY, DX)          # 
    U = U0*np.cos(alpha) - (U0*((R/r)**2)*np.cos(2*theta - alpha))  - Gamma*np.sin(theta)/(2*np.pi*r)
    V = U0*np.sin(alpha) - (U0*((R/r)**2)*np.sin(2*theta - alpha))  + Gamma*np.cos(theta)/(2*np.pi*r)
    if not internal_flow:
        b=r<R-1e-3
        U[b]=np.nan
        V[b]=np.nan
    return U, V


def KT_comf_map(z, l=2, a=1):
    # Karman Trefftz conformal map
    Z = l * a * ((z + a)**l + (z - a)**l) / ((z + a)** l - (z - a)** l)
    dZdz = 4*(l*a)**2 * ( (z-a)**(l-1)*(z+a)**(l-1) )/( ((z + a)**l - (z-a)**l) ** 2)
    return Z, dZdz

def KT_comf_map_inv(Z, l=2, a=1):
    # Inverse Karman Trefftz conformal map
    z = - a *(((Z-l)/(Z+l))**(1/l) + 1) / ( ( ( (Z-l)/(Z+l) )**(1/l) )-1)
    dzdZ = 1/ (4*(l*a)**2 * ( (z-a)**(l-1)*(z+a)**(l-1) )/( ((z + a)**l - (z-a)**l) ** 2))
    return z, dzdZ

def cyl_params(xc, yc, a=1, U0=0, alpha=0):
    rc = np.sqrt((a - xc) ** 2 + yc ** 2)
    beta = np.arcsin(- yc / (rc))
    Gamma = 4 * np.pi * rc * U0 * np.sin(beta - alpha)
    return rc, beta, Gamma

def KT_shape(xc, yc, l=2, a=1, n=100):
    """
    returns Karman-Trefftz profile coordinates
    INPUT:
    -  xc,yc   # Circle Center Location
    -  n       # Number of points along mapped foil surface
    """
    ## Main parameters
    rc, beta, _ =  cyl_params(xc, yc, a=a)
    print('>>> beta', beta)
    # --- Circle
    print('>>> TODO decide on which theta')
    #theta = np.linspace(beta, -2*np.pi+beta, n+1)[:-1]
    theta = np.linspace(0, 2*np.pi, n+1)[:-1]
    zc = (xc + 1j * yc)
    z_circ = zc + rc * np.exp(1j * theta)
    # --- Karman-Trefftz Conformal map - Profile Shape
    Zp, dZdz = KT_comf_map(z_circ, l=l, a=a)
    Xp = np.real(Zp)
    Yp = np.imag(Zp)
    return Xp, Yp

def KT_wall(xc, yc, l=2, a=1, U0=1, n=10, alpha=0):
    """ 
    # OUTPUT:
    #   X and Y coordinates of airfoil
    #   Cp : pressure coeff at foil surface
    #   Xg,Yg,Ug,Vg,CP : grid points, velocity field and  pressure coeff
    # EXAMPLES:
    #   [X_p,Y_p]    = fProfileKarmanTrefftz(xc,yc,l,n);
    #   [X_p,Y_p,Cp] = fProfileKarmanTrefftz(xc,yc,l,n,U0,alpha_deg);
    #   [X_p,Y_p,Cp,Ug,Vg,CP] = fProfileKarmanTrefftz(xc,yc,l,n,U0,alpha_deg,Xg,Yg);

    """
    rc, beta, Gamma =  cyl_params(xc, yc, a=a, U0=U0, alpha=alpha)
    # --- Pressure distribution on the airfoil
    print('Gamma', Gamma)
    # Velocity at circle surface
    #theta = np.linspace(beta, -2*np.pi+beta, n+1)[:-1]
    theta = np.linspace(0, 2*np.pi, n+1)[:-1]
    zc = (xc + 1j * yc)
    z_circ = zc + rc * np.exp(1j * theta)
    Zp, dZdz = KT_comf_map(z_circ, l=l, a=a)
    Xp = np.real(Zp)
    Yp = np.imag(Zp)
    # ---- 
    u_circ, v_circ = vc_u(np.real(z_circ)-xc, np.imag(z_circ)-yc, R=rc, U0=U0, alpha=alpha, Gamma=Gamma)
    # Velocities, -Cp on surface
    Zp, dZdz = KT_comf_map(z_circ, l=l, a=a)
    W_circ = (u_circ - 1j * v_circ) / dZdz
    U_circ =   np.real(W_circ)
    V_circ = - np.imag(W_circ)
    Q = np.sqrt(U_circ ** 2 + V_circ ** 2)
    Cp = 1 - (Q / U0) ** 2
    return Xp, theta, Cp, U_circ, V_circ


def KT_flow(Xg, Yg, xc, yc, l=2, a=1, U0=1, alpha=0):

    rc, beta, Gamma=  cyl_params(xc, yc, a=a, U0=U0, alpha=alpha)

    Zg = Xg + 1j * Yg
    z, dzdZ = KT_comf_map_inv(Zg, l=l, a=a)
    x = np.real(z)
    y = np.imag(z)
    u, v = vc_u(x-xc, y-yc, R=rc, U0=U0, alpha=alpha, Gamma=Gamma, internal_flow=False)
    W = (u - 1j*v) * dzdZ
    Ug =   np.real(W)
    Vg = - np.imag(W)
    Q  = np.sqrt(Ug ** 2 + Vg ** 2)
    CP = 1 - (Q / U0) ** 2
    return Ug, Vg, CP 
    
if __name__ == '__main__':
    import matplotlib.pyplot as plt
    ## Parameters for Karman Trefftz airfoil
    U0 = 1
    alpha = 5*np.pi/180;
    xc = - 0.2
    yc = 0.1
    tau_deg = 10
#     xc = - 0.3
#     yc = 0.1
#     tau_deg = 20
#     # ---  Display parameters
    n = 100


    # --- 
    l = 2 - tau_deg / 180 # l = 2 - tau/np.pi
    rc, beta =  cyl_params(xc, yc, a=1)
    nx = 200
    ny = 200
#     nx = 400
#     ny = 400
    nStreamlines = 21
    XLIM = np.array([- 3.5,3.5])
    YLIM = np.array([- 3,3])
#     ## Main function call (with grid for contour plots, but not necessary)
    Xg,Yg = np.meshgrid(np.linspace(XLIM[0],XLIM[1],nx), np.linspace(YLIM[0], YLIM[1], ny))

    Xp, Yp = KT_shape(xc, yc, l=l, n=n) # ,Xg,Yg)
    XCp, theta, Cp, U_circ, V_circ = KT_wall(xc, yc, l=l, n=n, alpha=alpha) # ,Xg,Yg)
    Ug, Vg, CP = KT_flow(Xg, Yg, l=l, U0=U0, alpha=alpha)

    # 
    U, V, X, Y = Ug, Vg, Xg, Yg
    xa, ya = Xp, Yp
    chord=np.max(Xp)-np.min(Xp)


    CP[np.isnan(CP)]=1
    print('Cp min', np.min(Cp))
    print('Cp min/max', np.nanmin(CP.flatten()), np.nanmax(CP.flatten()))


    # --- Plot
#     ## Plotting pressure distribution about the airfoil
#     fig,axes = plt.subplots(1, 3, sharey=False, figsize=(12.8,4.8)) # (6.4,4.8)
#     fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
#     ax = axes[0]
#     ax.plot(Xp, Yp    ,'k.-')
#     ax.plot(Xp[0], Yp[0] , 'bo')
#     ax.plot(Xp[1], Yp[1] , 'rd')
#     ax.set_xlabel('')
#     ax.set_ylabel('')
#     ax.legend()
#     ax.set_aspect('equal','box')
# 
#     ax = axes[1]
#     ax.plot((XCp - np.amin(XCp)) / (np.amax(XCp) - np.amin(XCp)), -Cp    , label='')
#     ax.set_xlabel('')
#     ax.set_ylabel('-Cp')
# #     plt.title(sprintf('Karman-Trefftz C_p \alpha = %.1f deg.',alpha_deg))
# #     plt.xlim(np.array([0,1]))
# #     plt.axis('ij')
# 
#     ax=axes[2]
#     im =ax.contourf(Xg, Yg, CP, 15)
#     ax.fill(Xp, Yp  ) #  ,'k.-')
#     ax.set_aspect('equal','box')
#     fig.colorbar(im)
# #     fill(X_p,Y_p,'w')
# #     plt.plot(P(:,1) * chord + np.amin(X_p),P(:,2) * chord,'b.')
# #     plt.plot(P(1,1) * chord + np.amin(X_p),P(1,2) * chord,'bd')
# #     plt.plot(P(2,1) * chord + np.amin(X_p),P(2,2) * chord,'ro')
# #     Y0 = np.linspace(plt.ylim(1),plt.ylim(2),nStreamlines)
# #     X0 = Y0 * 0 + plt.xlim(1)
# #     hlines = streamline(stream2(Xg,Yg,Ug,Vg,X0,Y0))
# #     set(hlines,'Color','k')
# #     plt.ylim(YLIM)
# #     plt.xlim(XLIM)
# #     plt.title('Karman-Trefftz Cp and streamlines')


    from flows2D import vorticity2D, circulation2D, flow_interp2D
    # --- Circulation contour
    print('>>> Circulation')
    theta=np.linspace(0,2*np.pi,2000)
    r_circ = 2.6
    xc = r_circ*np.cos(theta)+0
    yc = r_circ*np.sin(theta)+0
    with Timer('Interp'):
        #uc, vc = flow_interp2D(xc, yc, U, V, X, Y, method='nearest', algo='TriInterpolator')
        uc, vc = flow_interp2D(xc, yc, U, V, X, Y, method='nearest')
    Gamma  =  circulation2D(xc, yc, uc, vc, verbose=True)
    Cl     = -2*Gamma/(U0*chord)
    Cl_lin = 2*np.pi*np.sin(alpha -beta)
    Gamma = 4 * np.pi * rc * U0 * np.sin(beta - alpha)
    Cl2 = 8*np.pi*rc/chord*np.sin(alpha-beta)

    print('> Gamma : {:8.3f}'.format(Gamma))
    print('> rc    : {:8.3f}'.format(rc))
    print('> chord : {:8.3f}'.format(chord))
    print('> Cl    : {:8.3f}'.format(Cl))
    print('> Cl2   : {:8.3f}'.format(Cl2))
    print('> Cl_lin: {:8.3f}'.format(Cl_lin))
    print('> Cl/Cl: {:8.3f}'.format(Cl_lin/Cl))
    print('> c/4rc : {:8.3f}'.format(chord/(4*rc)))



    # --- Vorticity
    omax = 70 ; nLevels=8;
    omax = 90 ; nLevels=10;
    omax = 55 ; nLevels=12;
    omax = 0.5 ; nLevels=10;
    # omax = 75 ; nLevels=16;
    # omax = 150 ; nLevels=10;
    omz = vorticity2D(U, V, X, Y)
    omz = np.clip(omz, -omax, omax)

    fig,ax = plt.subplots(1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
    fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
    # im =ax.contourf(X, Y, omz, np.linspace(-omax,omax, nLevels), cmap="bwr")
    im =ax.contourf(X, Y, omz, 100, cmap="bwr")
    # im =ax.contourf(X, Y, omz, 100)
    ax.plot(xa, ya, 'k',lw=3)
    ax.plot(xc, yc, 'k',lw=2)
    cb=fig.colorbar(im)
    ax.set_xlabel(r'$x/c$ [-]')
    ax.set_ylabel(r'$y/c$ [-]')
    ax.set_aspect('equal','box')





    plt.show()

