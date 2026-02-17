""" 
Velocity field about the rankine nose (point source + freestream)
Note: methods are stored in SourcePoint package

Reference: 
    [1] Branlard (2017) Wind turbine aerodynamics and vorticity-based methods, Springer

"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from welib.tools.curves import streamQuiver
from welib.vortilib.elements.SourcePoint import *
from welib.vortilib.elements.RankineBody import *

def main(whichPlots=[0]):

    if 0 in whichPlots:
        rlim =2
        Sigma = 1 #2*np.pi
        vx = np.linspace(-2.5, 2.5, 301)           # grid
        vy = np.linspace(-2.5, 2.5, 302)           # grid
        X, Y = np.meshgrid(vx, vy)
        U0=1
        y_max = Sigma/(4*U0 )
        x0    = Sigma/(2*U0 )

        # --- Adding a freestream
        xc,yc = rn_stag(Sigma=Sigma, U0=U0)
        xs,ys = rn_coord(Sigma=Sigma, U0=U0, x_max = np.max(vx))
        U, V  = rn_u(X, Y, Ps=[xc,0], Sigma=Sigma, U0=U0)
        minSpeed=0
        maxSpeed=1.50
        nLevels =25

        yseed = np.concatenate( (np.linspace(np.min(vy)*0.85 , np.max(vy)*0.85, 7), [0.001]))
        start = np.array([yseed*0+np.min(vx)*0.85, yseed])

        yseed = np.linspace(-y_max*.8, y_max*.8, 3)
        start_i = np.array([yseed*0 +np.max(vx)*0.8,  yseed])

        Speed = np.sqrt((U**2+V**2))/U0
        Speed[Speed>maxSpeed] = maxSpeed

        fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
        im = ax.contourf(X, Y, Speed, levels=np.linspace(minSpeed, maxSpeed, nLevels), vmin=minSpeed, vmax=maxSpeed)
        cb=fig.colorbar(im)
        # Outer streamlines
        sp = ax.streamplot(vx, vy, U, V, color='k',start_points=start.T, linewidth=0.7,density=30,arrowstyle='-')
        qv = streamQuiver(ax, sp, n=5, scale=40, angles='xy')
        # Inner streamlines
        sp = ax.streamplot(vx, vy, U, V, color='k',start_points=start_i.T, linewidth=0.7,density=30,arrowstyle='-')
        qv = streamQuiver(ax, sp, n=3, scale=40, angles='xy')
        ax.plot(0, 0, 'ko',lw=3)
        ax.plot(xc, 0, 'kd',lw=3)
        ax.plot(xs+xc, ys, 'k',lw=2)
        ax.plot(xs+xc,-ys, 'k',lw=2)
        #ax.plot([-2,2], 2*[ Sigma/(2*U0)], 'r',lw=1)
        #ax.plot([-2,2], 2*[-Sigma/(2*U0)], 'r',lw=1)
        #ax.plot([xc,xc], [-Sigma/(4*U0 ),Sigma/(4*U0)], 'r',lw=1)
        ax.set_ylim([np.min(vy), np.max(vy)])
        ax.set_xlim([np.min(vx), np.max(vx)])
        ax.set_xlabel(r'$x$ [m]')
        ax.set_ylabel(r'$y$ [m]')
        ax.set_aspect('equal','box')
        ax.tick_params(direction='in', top=True, right=True)
        ax.set_title('Vortilib - Rankine nose')


if __name__ == '__main__':
    main(whichPlots=[0])
    plt.show()
if __name__=="__test__":
    main()
    pass
if __name__=="__export__":
    main(whichPlots=[0])
    from welib.tools.repo import export_figs_callback
    export_figs_callback(__file__)
