""" 
Velocity field about a source point

Reproduce Figure 32.1 from [1]

Reference: 
    [1] Branlard (2017) Wind turbine aerodynamics and vorticity-based methods, Springer

"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from welib.tools.curves import streamQuiver
from welib.vortilib.elements.SourcePoint import *

def main(whichPlots=[0,1]):
    dth = np.pi/12
    vpsi  = np.arange(-np.pi, np.pi, dth) # desired streamlines values
    vphi  = np.linspace(-1, 1, 10) # desired streamlines values
    rlim =2
    vx = np.linspace(-2.5, 2.5, 301)           # grid
    vy = np.linspace(-2.5, 2.5, 503)           # grid
    X, Y = np.meshgrid(vx, vy)
    U0=1

    psi = sp2d_psi(X, Y, Ps=[0,0], Sigma=2*np.pi)
    phi = sp2d_phi(X, Y, Ps=[0,0], Sigma=2*np.pi)
    U, V = sp2d_u(X, Y, Ps=[0,0], Sigma=2*np.pi, regMethod=1, regParam=0.30)
    R = np.sqrt(X**2 + Y**2)
    psi[np.abs(Y)<0.0001] = np.nan # To avoid repetition of streamline on y axis
    psi[R>rlim]=np.nan
    phi[R>rlim]=np.nan

    if 0 in whichPlots:
        # --- Figure 32.1 from [1]
        fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
        fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
        ax.contour(X, Y, phi, levels=vphi, colors='k', linestyles=':', linewidths=1.2)
        ax.contour(X, Y, psi, levels=vpsi, colors='k', linestyles='-', linewidths=1.2)
        ax.plot([-rlim, rlim], [0,0], 'k-', lw=1.2)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_aspect('equal','box')
        ax.set_axis_off()
        ax.set_title('Vortilib - Flow about a 2D point source') # [1] SourceStreamlines.pdf
        # Legend
        linePsi = Line2D([0], [0], label='Streamlines', color='k', ls='-')
        linePhi = Line2D([0], [0], label='Potential'  , color='k', ls=':')
        ax.legend(handles = [linePhi, linePsi])



    if 1 in whichPlots:
        # --- Plot velocity and streamlines
        minSpeed=0
        maxSpeed=3.00
        nLevels =25
        theta = np.linspace(0, 2*np.pi, 25)
        start = np.array([0.1*np.cos(theta), 0.1*np.sin(theta)])
        Speed        = np.sqrt((U**2+V**2))/U0

        fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
        im = ax.contourf(X, Y, Speed, levels=np.linspace(minSpeed, maxSpeed, nLevels), vmin=minSpeed, vmax=maxSpeed)
        cb=fig.colorbar(im)
        sp = ax.streamplot(vx, vy, U, V, color='k',start_points=start.T, linewidth=0.7,density=30,arrowstyle='-')
        qv = streamQuiver(ax, sp, spacing=1, offset=1.0, scale=40, angles='xy')
        ax.plot(0, 0, 'ko',lw=3)
        ax.set_xlabel(r'$x$ [m]')
        ax.set_ylabel(r'$y$ [m]')
        ax.set_aspect('equal','box')
        ax.tick_params(direction='in', top=True, right=True)
        ax.set_title('Vortilib - Flow about a 2D point source - alt')


if __name__ == '__main__':
    main(whichPlots=[0,1])
    plt.show()
if __name__=="__test__":
    main()
    pass
if __name__=="__export__":
    main(whichPlots=[0,1])
    from welib.tools.repo import export_figs_callback
    export_figs_callback(__file__)
