""" 
Velocity field about the rankine oval (point source + pointsource + freestream)

Reference: 
    [1] Branlard (2017) Wind turbine aerodynamics and vorticity-based methods, Springer

"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from welib.tools.curves import streamQuiver
from welib.vortilib.elements.RankineBody import *

def main(whichPlots=[0,1]):

    vx = np.linspace(-2.5, 2.5, 401)           # grid
    vy = np.linspace(-2.5, 2.5, 402)           # grid
    X, Y = np.meshgrid(vx, vy)
    U0 = 1
    R  = 1
    t  = 0.4

    # --- 
    # Find oval parameters
    Sigma, b = ro_params(R=R, t=t, U0=U0) #, method='ajkj')
    psi   = ro_psi(X, Y, Sigma=Sigma, b=b, U0=U0)
    U, V  = ro_u  (X, Y, Sigma=Sigma, b=b, U0=U0)
    x0, y0= ro_stag(b=b, Sigma=Sigma, U0=U0)
    xs, ys = ro_coord(b=b, Sigma=Sigma, U0=U0)

    minSpeed=0
    maxSpeed=1.50
    nLevels =25

    yseed = np.concatenate( (np.linspace(np.min(vy)*0.85 , np.max(vy)*0.85, 7), [0.001]))
    start = np.array([yseed*0+np.min(vx)*0.85, yseed])
    #yseed = np.linspace(-t*0.8, t*0.8, 3)
    yseed = np.linspace(-t*0.6, t*0.6, 3)
    start_i = np.array([yseed*0,  yseed])


    if 0 in whichPlots:
        # --- Figure 32.1 from [1]
        fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
        fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
        #ax.contour(X, Y, phi, levels=vphi, colors='k', linestyles=':', linewidths=1.2)
        ax.contour(X, Y, psi, levels=8, colors='k', linestyles='-', linewidths=1.2)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_aspect('equal','box')
        ax.set_axis_off()
        ax.set_title('Vortilib - Flow about a 2D Rankine oval') # [1] SourceStreamlines.pdf
        # Legend
        #linePsi = Line2D([0], [0], label='Streamlines', color='k', ls='-')
        #linePhi = Line2D([0], [0], label='Potential'  , color='k', ls=':')
        #ax.legend(handles = [linePhi, linePsi])


    if 1 in whichPlots:
        #U, V = sp2d_u(X, Y, Ps=[xc,0], Sigma=Sigma, U0=U0)
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
        ax.plot(x0, y0, 'ko',lw=3)
        ax.plot([-b,b], [0,0], 'kd',lw=1)
        ax.plot(xs, ys, 'k-',lw=2)
        ax.plot(xs,-ys, 'k-',lw=2)
        #ax.plot([-R,-R, R, R, -R], [-t, t, t, -t, -t], 'k--',lw=1)
        ax.set_ylim([np.min(vy), np.max(vy)])
        ax.set_xlim([np.min(vx), np.max(vx)])
        ax.set_xlabel(r'$x$ [m]')
        ax.set_ylabel(r'$y$ [m]')
        ax.set_aspect('equal','box')
        ax.tick_params(direction='in', top=True, right=True)
        ax.set_title('Vortilib - Flow about a 2D Rankine oval - alt')


if __name__ == '__main__':
    main(whichPlots=[1])
    plt.show()
if __name__=="__test__":
    main(whichPlots=[1])
    pass
if __name__=="__export__":
    main(whichPlots=[1])
    from welib.tools.repo import export_figs_callback
    export_figs_callback(__file__)
