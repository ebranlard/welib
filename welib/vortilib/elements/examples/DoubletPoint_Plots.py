""" 
Velocity field about a Double point

Reproduce Figure 32.3 from [1]

Reference: 
    [1] Branlard (2017) Wind turbine aerodynamics and vorticity-based methods, Springer.

"""
import numpy as np
import matplotlib.pyplot as plt
from welib.tools.curves import streamQuiver
from welib.vortilib.elements.DoubletPoint import *

def main(whichPlots=[0,1]):
    v  = np.arange(-4, 4.1, 0.5)         # desired streamlines values
    vg = np.arange(-2.5, 2.5+0.05, 0.05) # grid
    vx = vg
    vy = vg
    X, Y = np.meshgrid(vx, vy)
    U0=1
    alpha = 30*np.pi/180
    psi = dp2d_psi(X, Y, Pd=[0,0], Mu=1, alpha=alpha)
    U, V = dp2d_u (X, Y, Pd=[0,0], Mu=1, alpha=alpha) #, regMethod=1, regParam=0.10)

    if 0 in whichPlots:
        # --- Figure 32.3 from [1]
        fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
        fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
        ax.contour(X, Y, psi, levels=v, colors='k', linestyles='-', linewidths=1.0)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_aspect('equal','box')
        ax.set_axis_off()
        ax.set_title('Vortilib - Flow about a 2D point doublet') # [1] DoubletStreamlines.pdf

    if 1 in whichPlots:
        # --- Plot velocity and streamlines
        minSpeed=0
        maxSpeed=10.20
        nLevels =25
        yseed = np.linspace(0.2 , np.max(vy)*0.95, 7) # TODO figure something out based on levels/mass flow
        yseed = np.concatenate((yseed,-yseed))
        start = np.array([yseed*np.cos(alpha+np.pi/2), yseed*np.sin(alpha+np.pi/2)])

        Speed = np.sqrt((U**2+V**2))/1

        fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
        im = ax.contourf(X, Y, Speed, levels=np.linspace(minSpeed, maxSpeed, nLevels), vmin=minSpeed, vmax=maxSpeed)
        cb=fig.colorbar(im)
        sp = ax.streamplot(vx, vy, U, V, color='k',start_points=start.T, linewidth=0.7,density=30,arrowstyle='-')
        #qv = streamQuiver(ax, sp, spacing=1, offset=1.0, scale=40, angles='xy')
        qv = streamQuiver(ax, sp, n=5, scale=40, angles='xy')
        ax.plot(0, 0, 'ko',lw=3)
        #ax.plot(start[0], start[1], 'ko',lw=3)
        ax.set_xlabel(r'$x$ [m]')
        ax.set_ylabel(r'$y$ [m]')
        ax.set_aspect('equal','box')
        ax.tick_params(direction='in', top=True, right=True)
        ax.set_title('Vortilib - Flow about a 2D point doublet - alt')


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
