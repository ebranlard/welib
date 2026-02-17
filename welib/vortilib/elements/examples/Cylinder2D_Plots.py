""" 
Velocity field about a 2D acyclic and cyclic cylinder point

Reference: 
    [1] Branlard (2017) Wind turbine aerodynamics and vorticity-based methods, Springer

"""
import numpy as np
import matplotlib.pyplot as plt
from welib.tools.curves import streamQuiver
from welib.vortilib.elements.VortexCylinder2D import *

def main(gammaFact=[0, 0.5, 1.1]):
    vg = np.arange(-2.5, 2.5+0.005, 0.005)           # grid
    vx = vg
    vy = vg
    X, Y = np.meshgrid(vx, vy)
    U0 = 1
    R =1 
    Gamma_ref =  -4*np.pi*U0*R
    alpha =0*np.pi/180

    for fact in gammaFact:
        Gamma = Gamma_ref*fact

        U, V   = vc_u(X, Y, P=[0,0], U0=U0, R=R, Gamma=Gamma, alpha=alpha) 
        xs, ys = vc_stag(U0=U0, R=R, Gamma=Gamma, alpha=alpha)
        xc, yc = vc_coords(R=R, ne=100)

        # --- Plot velocity and streamlines
        nLevels =25
        # Internal streamlines 
        ys2 = np.linspace(-0.8, 0.8,5)
        xs2 = ys2*0

        yseed = np.linspace(np.min(vy)*0.85 , np.max(vy)*0.85, 7)
        xs1   = yseed*0+np.min(vx)*0.85
        ys1   = yseed

        minSpeed=0
        maxSpeed=3.00

        if fact==0:
            pass
        elif fact<1:
            xs1 = np.concatenate( (xs+[+0.02,-0.02], xs1))
            ys1 = np.concatenate( (ys, ys1))
        else:
#             xs1 = np.concatenate( (xs+[+0.01,-0.01], xs1))
#             ys1 = np.concatenate( (ys, ys1))
            xs1 = np.concatenate( ([xs[1]     ,xs[1]     , xs[0]], xs1))
            ys1 = np.concatenate( ([ys[1]+0.01,ys[1]-0.01, ys[0]], ys1))

        start  = np.array( [xs1, ys1])
        starti = np.array([xs2, ys2])

        Speed        = np.sqrt((U**2+V**2))/U0
        Speed[Speed>maxSpeed] = maxSpeed

        fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
        im = ax.contourf(X, Y, Speed, levels=np.linspace(minSpeed, maxSpeed, nLevels), vmin=minSpeed, vmax=maxSpeed)
        cb=fig.colorbar(im)
        # External streamlines
        sp = ax.streamplot(vx, vy, U, V, color='k',start_points=start.T, linewidth=0.7,density=30,arrowstyle='-')
        qv = streamQuiver(ax, sp, n=5, scale=40, angles='xy')

        # Internal streamlines
        sp = ax.streamplot(vx, vy, U, V, color='k',start_points=starti.T, linewidth=0.7,density=30,arrowstyle='-')
        qv = streamQuiver(ax, sp, n=3, scale=40, angles='xy')


        ax.plot(xs, ys, 'ko',lw=3)

        ax.plot(xc, yc, 'k',lw=3)
        ax.set_xlabel(r'$x/R$ [-]')
        ax.set_ylabel(r'$y/R$ [-]')
        ax.set_aspect('equal','box')
        ax.tick_params(direction='in', top=True, right=True)
        ax.set_title('Vortilib - Flow about a 2D vortex cylinder - Gamma={:.1f}'.format(fact))

if __name__ == '__main__':
    main(gammaFact=[0, 0.5, 1.1])
    plt.show()
if __name__=="__test__":
    main()
    pass
if __name__=="__export__":
    main(gammaFact=[0, 0.5, 1.1])
    from welib.tools.repo import export_figs_callback
    export_figs_callback(__file__)
