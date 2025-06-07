import numpy as np
from welib.essentials import *
from welib.airfoils.shapes import *
    
def vc_u(DX, DY, U0=1, R=1, Gamma=0, alpha=0, internal_flow=True):
    """
    2D Cylinder with possible free stream and circulation
    """
    r     = np.sqrt(DX**2 + DY**2) # Radius at mesh locations
    theta = np.arctan2(DY, DX)          # 
    U = U0*np.cos(alpha) - (U0*((R/r)**2)*np.cos(2*theta - alpha))  - Gamma*np.sin(theta)/(2*np.pi*r)
    V = U0*np.sin(alpha) - (U0*((R/r)**2)*np.sin(2*theta - alpha))  + Gamma*np.cos(theta)/(2*np.pi*r)
    if not internal_flow:
        b=r<R-1e-3
        U[b]=np.nan
        V[b]=np.nan
    return U, V

def KT_comf_map(Z, l=2, A=1):
    """ Karman Trefftz conformal map. Z: cylinder plane, z: airfoil plane"""
    z    = l * A * ((Z + A)**l + (Z - A)**l) / ((Z + A)** l - (Z - A)** l)
    dzdZ = 4*(l*A)**2 * ( (Z-A)**(l-1)*(Z+A)**(l-1) )/( ((Z + A)**l - (Z-A)**l) ** 2)
    return z, dzdZ

def KT_comf_map_inv(z, l=2, A=1):
    """ Inverse Karman Trefftz conformal map: Z: cylinder plane, z:airfoil plane"""
    Z    = - A *(((z-l)/(z+l))**(1/l) + 1) / ( ( ( (z-l)/(z+l) )**(1/l) )-1)
    dZdz = 1/ (4*(l*A)**2 * ( (Z-A)**(l-1)*(Z+A)**(l-1) )/( ((Z + A)**l - (Z-A)**l) ** 2))
    return Z, dZdz

def cyl_params(XC, YC, A=1, U0=0, alpha=0):
    """ Geometrical parameters of the Z-plane cylinder """
    R     = np.sqrt((A - XC) ** 2 + YC ** 2)
    beta  = np.arcsin(- YC /R)
    Gamma = 4 * np.pi * R * U0 * np.sin(beta - alpha)
    return R, beta, Gamma

def KT_shape(XC, YC, l=2, A=1, n=100, thetaStart=None):
    """
    Returns Karman-Trefftz profile coordinates
    INPUT:
     - XC,YC: Center of Cylinder (in Z-plane)
     - l, A : Karman-Treffz parameters
     - n    : Number of points along mapped foil surface
    OUTPUTS:
     - xa,ya: Coordinates of airfoil (in z-plane)
    """
    # Cylinder parameters
    R, beta, Gamma =  cyl_params(XC, YC, A=A, U0=0, alpha=0)
    # Cylinder coordinates in the Z-plane
    if thetaStart is None:
        theta = np.linspace(beta, -2*np.pi+beta, n+1)[:-1]
    else:
        theta = np.linspace(0, 2*np.pi, n+1)[:-1]
    ZC = (XC + 1j * YC)
    Z_cyl = ZC + R * np.exp(1j * theta)
    # Transform back to the z-plane
    za, dzdZ = KT_comf_map(Z_cyl, l=l, A=A)
    xa = np.real(za)
    ya = np.imag(za)
    return xa, ya

def KT_wall(XC, YC, l=2, A=1, n=10, U0=1, alpha=0, thetaStart=None):
    """ 
    Returns the Karman-Trefftz profile coordinates and pressure coefficient
    on the wall surface.

    INPUTS:
     - XC,YC: Center of Cylinder [m]
     - l, A : Karman-Treffz parameters
     - n    : Number of points along mapped surface
     - U0   : Freestream velocity [m/s]
     - alpha: angle of attack [rad]
    OUTPUTS:
     - xa, ya:  coordinates of airfoil (in z-plane) [m]
     - Cp_w : pressure coefficient at airfoil surface [-]
     - u_w  : x-velocity at arifoil surface [m/s]
     - v_w  : y-velocity at arifoil surface [m/s]
    """
    # Cylinder parameters
    R, beta, Gamma =  cyl_params(XC, YC, A=A, U0=U0, alpha=alpha)
    # Cylinder coordinates in the Z-plane
    if thetaStart is None:
        theta = np.linspace(beta, -2*np.pi+beta, n+1)[:-1]
    else:
        theta = np.linspace(0, 2*np.pi, n+1)[:-1]
    ZC    = (XC + 1j * YC)
    Z_cyl = ZC + R * np.exp(1j * theta)
    # Transform back to the z-plane
    za, dzdZ = KT_comf_map(Z_cyl, l=l, A=A)
    xa = np.real(za)
    ya = np.imag(za)
    # Velocities on cylinder surface in Z plane
    U_cyl, V_cyl = vc_u(np.real(Z_cyl)-XC, np.imag(Z_cyl)-YC, R=R, U0=U0, alpha=alpha, Gamma=Gamma)
    # Back to z-plane
    w_cyl = (U_cyl - 1j * V_cyl) / dzdZ
    u_cyl =   np.real(w_cyl)
    v_cyl = - np.imag(w_cyl)
    Q = np.sqrt(u_cyl ** 2 + v_cyl ** 2)
    Cp = 1 - (Q / U0) ** 2
    return xa, ya, Cp, u_cyl, v_cyl

def KT_flow(x, y, XC, YC, l=2, A=1, U0=1, alpha=0):
    # Cylinder parameters
    R, _, Gamma =  cyl_params(XC, YC, A=A, U0=U0, alpha=alpha)
    # Transform z to Z plane
    z = x + 1j * y
    Z, dZdz = KT_comf_map_inv(z, l=l, A=A)
    X = np.real(Z)
    Y = np.imag(Z)
    # --- Compute cylinder velocity in Z-plane
    U, V = vc_u(X-XC, Y-YC, R=R, U0=U0, alpha=alpha, Gamma=Gamma, internal_flow=False)
    # Back to z-plane
    w = (U - 1j*V) * dZdz
    u =   np.real(w)
    v = - np.imag(w)
    Q  = np.sqrt(u ** 2 + v ** 2)
    CP = 1 - (Q / U0) ** 2
    return u, v, CP 
    
if __name__ == '__main__':
    import matplotlib.pyplot as plt
    # --- Parameters
    U0 = 1
    alpha = 5*np.pi/180;
    XC      = -0.2              # X-coord of cylinder center in Z-plane
    YC      = 0.1               # Y-coord of cylinder center in Z-plane
    tau_deg = 10
#     xc = - 0.3
#     yc = 0.1
#     tau_deg = 20
#     # ---  Display parameters
    n       = 300               # Number of points for airfoil


    # --- 
    l = 2 - tau_deg / 180 # l = 2 - tau/np.pi
    rc, beta, _ =  cyl_params(XC, YC, A=1)
    nx = 200
    ny = 200
#     nx = 400
#     ny = 400
    nStreamlines = 21
    XLIM = np.array([- 3.5,3.5])
    YLIM = np.array([- 3.5,3.5])
#     ## Main function call (with grid for contour plots, but not necessary)
    vx = np.linspace(XLIM[0],XLIM[1],nx)
    vy = np.linspace(YLIM[0], YLIM[1], ny)
    Xg,Yg = np.meshgrid(vx, vy)

    Xp, Yp = KT_shape(XC, YC, l=l, n=n) # ,Xg,Yg)
    XCp, theta, Cp, U_circ, V_circ = KT_wall(XC, YC, l=l, n=n, alpha=alpha) # ,Xg,Yg)
    Ug, Vg, CP = KT_flow(Xg, Yg, XC, YC, l=l, U0=U0, alpha=alpha)

    # 
    U, V, X, Y = Ug, Vg, Xg, Yg
    xa, ya = Xp, Yp
    xmax=XLIM[1]
    chord=np.max(Xp)-np.min(Xp)


    CP[np.isnan(CP)]=1
    print('Cp min', np.min(Cp))
    print('Cp min/max', np.nanmin(CP.flatten()), np.nanmax(CP.flatten()))

    # --- AirfoilShpae
    arf = AirfoilShape(xa, ya)
    arf.write('../vortilib/panelcodes/KarmanTrefftz-{}-NoScale.csv'.format(n), comment_char='',delim=',')
    #arf.plot()


    # --- Plot
    ## Plotting pressure distribution about the airfoil
    fig,axes = plt.subplots(1, 3, sharey=False, figsize=(12.8,4.8)) # (6.4,4.8)
    fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
    ax = axes[0]
    ax.plot(Xp, Yp    ,'k.-')
    ax.plot(Xp[0], Yp[0] , 'bo')
    ax.plot(Xp[1], Yp[1] , 'rd')
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.legend()
    ax.set_aspect('equal','box')

    ax = axes[1]
    ax.plot((XCp - np.amin(XCp)) / (np.amax(XCp) - np.amin(XCp)), Cp    , label='')
    ax.set_xlabel('')
    ax.set_ylabel('-Cp')
    #     plt.title(sprintf('Karman-Trefftz C_p \alpha = %.1f deg.',alpha_deg))
    #     plt.xlim(np.array([0,1]))
    #     plt.axis('ij')

    ax=axes[2]
    im =ax.contourf(Xg, Yg, CP, 15)
    ax.fill(Xp, Yp  ) #  ,'k.-')
    ax.set_aspect('equal','box')
    fig.colorbar(im)


    # Streamlines
    ys = np.linspace(-xmax, xmax, 15)
    start = np.array([ys*0-xmax, ys])
    nLevels = 11
    minSpeed = 0
    maxSpeed = 2.0

    Q = np.sqrt(U**2+V**2)
    Speed = np.sqrt((U**2+V**2))/U0

    # Plot
    fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
    im = ax.contourf(X, Y, Speed, levels=np.linspace(minSpeed, maxSpeed, nLevels), vmin=minSpeed, vmax=maxSpeed)
    cb=fig.colorbar(im)
    cb.set_label(r'$\sqrt{u^2+v^2}/U_\text{ref}$')
    sp = ax.streamplot(vx, vy, U, V, color='k',start_points=start.T, linewidth=0.7,density=30,arrowstyle='-')
    ax.set_aspect('equal','box')
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.legend()

#     fill(X_p,Y_p,'w')
#     plt.plot(P(:,1) * chord + np.amin(X_p),P(:,2) * chord,'b.')
#     plt.plot(P(1,1) * chord + np.amin(X_p),P(1,2) * chord,'bd')
#     plt.plot(P(2,1) * chord + np.amin(X_p),P(2,2) * chord,'ro')
#     Y0 = np.linspace(plt.ylim(1),plt.ylim(2),nStreamlines)
#     X0 = Y0 * 0 + plt.xlim(1)
#     hlines = streamline(stream2(Xg,Yg,Ug,Vg,X0,Y0))
#     set(hlines,'Color','k')
#     plt.ylim(YLIM)
#     plt.xlim(XLIM)
#     plt.title('Karman-Trefftz Cp and streamlines')


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

