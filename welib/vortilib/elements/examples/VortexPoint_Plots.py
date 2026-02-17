""" 
Velocity field about a 2d vortex point

Reproduce Figure 32.1 from [1]

Reference: 
    [1] Branlard (2017) Wind turbine aerodynamics and vorticity-based methods, Springer

"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from welib.tools.curves import streamQuiver
from welib.vortilib.elements.VortexPoint import *

def main(whichPlots=[0]):

    vpsi  = np.arange(-np.pi, np.pi+0.1, np.pi/7) # desired streamlines values
    vphi  = np.arange(-np.pi, np.pi+0.1, np.pi/7) # desired streamlines values
    rlim =2
    vg = np.arange(-2.5, 2.5+0.05, 0.05)       # grid
    vx = vg
    vy = vg
    X, Y = np.meshgrid(vx, vy)
    psi  = vp_psi(X, Y, Pv = [0,0], Gamma = 1)
    phi  = vp_phi(X, Y, Pv = [0,0], Gamma = 1)
    U, V = vp_u  (X, Y, Pv = [0,0], Gamma = 1, regMethod = 1, regParam = 0.10)
    Speed = np.sqrt((U**2+V**2))
    R = np.sqrt(X**2 + Y**2)
    phi[np.abs(Y)<0.0001] = np.nan # To avoid repetition of streamline on y axis
    psi[R>rlim]=np.nan
    phi[R>rlim]=np.nan
    
    # Streamlines
    minSpeed=0
    maxSpeed=1.20
    nLevels =25
    yseed = np.linspace(0.25 , np.max(vy)*0.95, 7) # TODO figure something out based on levels
    start = np.array([yseed*0, yseed])

    if 0 in whichPlots:
        # --- Figure 32.1 from [1]
        # Plot
        fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
        fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
        ax.contour(X, Y, phi, levels=vphi, colors='k', linestyles=':', linewidths=1.2)
        sp = ax.contour(X, Y, psi, levels=vpsi, colors='k', linestyles='-', linewidths=1.2)
        #qv = streamQuiver(ax, sp, n=7, scale=40, angles='xy')
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_aspect('equal','box')
        ax.set_axis_off()
        ax.set_title('Vortilib - Flow about a 2D point vortex') # [1] VortexStreamlines.pdf
        ax.plot([-rlim, rlim], [0,0], 'k:', lw=1.2)
        # Legend
        linePsi = Line2D([0], [0], label='Streamlines', color='k', ls='-')
        linePhi = Line2D([0], [0], label='Potential'  , color='k', ls=':')
        ax.legend(handles = [linePhi, linePsi])


    if 1 in whichPlots:
        # --- 
#         bSing        = (X**2**2+Y**2**2)<0.001
#         Speed[bSing] = np.nan
        fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
        im = ax.contourf(X, Y, Speed, levels=np.linspace(minSpeed, maxSpeed, nLevels), vmin=minSpeed, vmax=maxSpeed)
        cb=fig.colorbar(im)
        sp = ax.streamplot(vx, vy, U, V, color='k',start_points=start.T,linewidth=0.7,density=30,arrowstyle='-')
        qv = streamQuiver(ax, sp, n=7, scale=40, angles='xy')
        ax.plot(0, 0, 'ko',lw=3)
        ax.set_xlabel(r'$x$ [m]')
        ax.set_ylabel(r'$y$ [m]')
        ax.set_aspect('equal','box')
        ax.tick_params(direction='in', top=True, right=True)
        ax.set_title('Vortilib - Flow about a 2D point vortex - alt')


if __name__ == '__main__':
    main()
    plt.show()
if __name__=="__test__":
    main()
    pass
if __name__=="__export__":
    main(whichPlots=[0,1])
    from welib.tools.repo import export_figs_callback
    export_figs_callback(__file__)
