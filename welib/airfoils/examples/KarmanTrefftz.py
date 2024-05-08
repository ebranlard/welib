import numpy as np
import matplotlib.pyplot as plt

from welib.airfoils.karman_trefftz import *  

def main(IPlot):

    # --- Parameters
    alpha = 10 *np.pi/180
    U0    = 1
    # Main parameters affecting airfoil geometry
    XC      = -0.2              # X-coord of cylinder center in Z-plane [m]
    YC      = 0.1               # Y-coord of cylinder center in Z-plane [m]
    tau_deg = 10                # Trailing edge [deg]
    A       = 1                 # Cylinder X-intersect in Z-plane [m]

    # Plotting parameters
    n       = 100               # Number of points for airfoil

    # --- Derived parameters
    l = 2 - tau_deg / 180 # Karman-Trefftz "Lambda" parameter    

    # --- Airfoil shape
    [xa, ya] = KT_shape(XC, YC, l, A, n);
    XCp, theta, Cp, U_circ, V_circ = KT_wall(XC, YC, l=l, n=n, alpha=alpha) # ,Xg,Yg)

    # --- Transform a couple of points for fun
    R, _, _ = cyl_params(XC, YC, A=A, U0=U0, alpha=alpha)
    print('R',R)
    Z_C = XC + 1j*YC
    Z_O = 0 + 1j*0
    Z_A = A + 1j*0
    z_C, _ = KT_comf_map(Z_C, l=l, A=A)
    z_O, _ = KT_comf_map(Z_O, l=l, A=A)
    z_A, _ = KT_comf_map(Z_A, l=l, A=A)
    #print('>>> C', z_C)
    #print('>>> O', z_O)
    #print('>>> A', z_A)

    theta = np.linspace(0, 2*np.pi, 100)
    Xcirc = R* np.cos(theta) + XC
    Ycirc = R* np.sin(theta) + YC


    # --- Flow
    if 4 in IPlot:
        XLIM = np.array([-3.5, 3.5])
        YLIM = np.array([-3.5, 3.5])
        #     ## Main function call (with grid for contour plots, but not necessary)
        vx = np.linspace(XLIM[0], XLIM[1], 200 )
        vy = np.linspace(YLIM[0], YLIM[1], 200 )
        X,Y = np.meshgrid(vx, vy)
        U, V, CP = KT_flow(X, Y, XC, YC, l=l, U0=U0, alpha=alpha)
        Speed = np.sqrt((U**2+V**2))/U0

        # --- Plot CP / Flow field
        chord=np.max(xa)-np.min(xa)

        CP[np.isnan(CP)]=1
        print('Cp min', np.min(Cp))
        print('Cp min/max', np.nanmin(CP.flatten()), np.nanmax(CP.flatten()))

        # Streamlines
        ys = np.linspace(-XLIM[0], XLIM[0], 15)
        start = np.array([ys*0-XLIM[0], ys])
        nLevels = 11
        minSpeed = 0
        maxSpeed = 2.0

        # --- Plot Velocity norm
        fig,ax = plt.subplots(1, 1, sharey=False, figsize=(8.4,5.8)) # (6.4,4.8)
        fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
        im = ax.contourf(X, Y, Speed, levels=np.linspace(minSpeed, maxSpeed, nLevels), vmin=minSpeed, vmax=maxSpeed)
        cb=fig.colorbar(im)
        cb.set_label(r'$\sqrt{u^2+v^2}/U_\text{ref}$')
        sp = ax.streamplot(vx, vy, U, V, color='k',start_points=start.T, linewidth=0.7,density=30,arrowstyle='-')
        ax.fill(xa, ya , 'k' ) 
        ax.plot(xa, ya, 'k-')
        ax.set_aspect('equal','box')
        ax.set_xlabel('x [m]')
        ax.set_ylabel('y [m]')
        ax.set_title('Airfoil - Karman-Trefftz Velocity Field')


    # --- Plot
    if 0 in IPlot:
        fig, axes = plt.subplots(1, 2, sharex=True, figsize=(12.8,4.8))
        axes[0].plot(xa, ya)
        axes[0].plot(np.real(z_A), np.imag(z_A), '.')
    #     axes[0].plot(np.real(z_O), np.imag(z_O), 'o')
    #     axes[0].plot(np.real(z_C), np.imag(z_C), 'd')
        axes[0].grid(ls=':', lw=0.5)
        axes[0].set_xlim([-2.5,2.5])
        axes[0].set_aspect('equal','box')
        axes[0].set_xlim([-2.5,2.5])

        axes[1].plot(Xcirc, Ycirc)
        axes[1].plot(np.real(Z_A), np.imag(Z_A), '.')
    #     axes[1].plot(np.real(Z_O), np.imag(Z_O), 'o')
        axes[1].plot(np.real(Z_C), np.imag(Z_C), 'd')
        axes[1].set_xlim([-2.5,2.5])
        axes[1].set_aspect('equal','box')
        axes[1].grid(ls=':', lw=0.5)
        axes[1].set_xlim([-2.5,2.5])

    # --- Flow
    if 1 in IPlot:
        XLIM = np.array([-3.5, 3.5])
        YLIM = np.array([-3.5, 3.5])
        #     ## Main function call (with grid for contour plots, but not necessary)
        vx = np.linspace(XLIM[0], XLIM[1], 200 )
        vy = np.linspace(YLIM[0], YLIM[1], 200 )
        X,Y = np.meshgrid(vx, vy)
        U, V, CP = KT_flow(X, Y, XC, YC, l=l, U0=U0, alpha=alpha)
        Speed = np.sqrt((U**2+V**2))/U0



        # --- Plot CP / Flow field
        chord=np.max(xa)-np.min(xa)

        CP[np.isnan(CP)]=1
        print('Cp min', np.min(Cp))
        print('Cp min/max', np.nanmin(CP.flatten()), np.nanmax(CP.flatten()))


        # Streamlines
        ys = np.linspace(-XLIM[0], XLIM[0], 15)
        start = np.array([ys*0-XLIM[0], ys])
        nLevels = 11
        minSpeed = 0
        maxSpeed = 2.0



        # --- Plot
        ## Plotting pressure distribution about the airfoil
        fig,axes = plt.subplots(1, 3, sharey=False, figsize=(12.8,4.8)) # (6.4,4.8)
        fig.subplots_adjust(left=0.06, right=0.96, top=0.93, bottom=0.11, hspace=0.30, wspace=0.50)

        ax = axes[0]
        ax.plot((XCp - np.amin(XCp))/chord , Cp    , label='')
        ax.set_ylim([-3,1])
        ax.set_xlabel('x/c [-]')
        ax.set_ylabel('Cp [-]')
        ax.invert_yaxis()
        ax.set_title('Pressure coefficient on airfoil surface')
        #     plt.title(sprintf('Karman-Trefftz C_p \alpha = %.1f deg.',alpha_deg))
        #     plt.xlim(np.array([0,1]))
        #     plt.axis('ij')

        ax=axes[1]
        im =ax.contourf(X, Y, CP, np.linspace(-3, 1, 41))
        ax.fill(xa, ya , 'k' ) 
        ax.plot(xa, ya, 'k-')
        ax.set_aspect('equal','box')
        ax.set_xlabel('x [m]')
        ax.set_ylabel('y [m]')
        cb = fig.colorbar(im)
        cb.set_label(r'$C_p$')
        ax.set_title('Pressure coefficient')

        ax=axes[2]
        im = ax.contourf(X, Y, Speed, levels=np.linspace(minSpeed, maxSpeed, nLevels), vmin=minSpeed, vmax=maxSpeed)
        cb=fig.colorbar(im)
        cb.set_label(r'$\sqrt{u^2+v^2}/U_\text{ref}$')
        sp = ax.streamplot(vx, vy, U, V, color='k',start_points=start.T, linewidth=0.7,density=30,arrowstyle='-')
        ax.fill(xa, ya , 'k' ) 
        ax.plot(xa, ya, 'k-')
        ax.set_aspect('equal','box')
        ax.set_xlabel('x [m]')
        ax.set_ylabel('y [m]')
        ax.set_title('Velocity norm')

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

    if 3 in IPlot:
        from welib.CFD.flows2D import vorticity2D, circulation2D, flow_interp2D
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


if __name__ == '__main__':
    main([4])
    plt.show()
    
if __name__ == '__test__':
    main([0,1,2,3])
if __name__=="__export__":
    main([4])
    from welib.tools.repo import export_figs_callback
    export_figs_callback(__file__)


