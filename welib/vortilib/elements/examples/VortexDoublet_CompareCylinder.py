"""
Compare the velocity field from a vortex doublet line obtained using numerical or analytical integration

"""
# --- General
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable
# --- Local
from welib.vortilib.elements.VortexDoublet import *
from welib.vortilib.elements.VortexCylinder import *

from welib.tools.tictoc import Timer
from welib.tools.curves import streamQuiver

def main():
    # --- Parameters
    R      = 10
    gamma_t=-10
    XLIM  = np.asarray([0.01*R,10*R ]) # 
    ZLIM  = np.asarray([-10*R,0.01*R]) # 
    nx    = 210                   # Number of points for velocity evaluation
    nz    = nx
    bLog  = False # take logarithm of speed
    # bLog  = True # take logarithm of speed
    IPlot=[1,2]

    # --- Derived parameters
    dmz_dz = gamma_t * R**2 * np.pi # doublet intensity per length

    # --- Flow field and speed
    #zs = np.linspace(ZLIM[0],ZLIM[1],nx).astype(np.float32)
    #xs = np.linspace(XLIM[0]*1.1,XLIM[1]*1.1,nx).astype(np.float32)
    xbeg = XLIM[0]*1.1
    zbeg = ZLIM[0]
    xend = XLIM[1]*1.1
    zend = ZLIM[1]
    zs = np.mgrid[zbeg:zend+((zend-zbeg)/(nz-1)):(zend-zbeg)/(nz-1)]
    xs = np.mgrid[xbeg:xend+((xend-xbeg)/(nx-1)):(xend-xbeg)/(nx-1)]

    # remove singularities
    #bBad=np.abs(np.abs(xs)-R)<0.2*R
    #xs=xs[~bBad]
    [Z,X]=np.meshgrid(zs,xs)
    Y=X*0
    with Timer('Cylinder'):
        urt,uzt = vc_tang_u(X,Y,Z,gamma_t=gamma_t,R=R,polar_out=True)
    with Timer('Cylinder Doublet'):
        urt,uzt = vc_tang_u_doublet(X,X*0,Z,gamma_t=gamma_t,R=R)
    with Timer('Doublet'):
        urn,uzn = doublet_line_polar_u(X,Z,dmz_dz)


    if 1 in IPlot:
        if bLog:
            Speed_t=np.log(np.sqrt(uzt**2+urt**2))
            Speed_n=np.log(np.sqrt(uzn**2+urn**2))
        else:
            Speed_t=np.sqrt(uzt**2+urt**2)
            Speed_n=np.sqrt(uzn**2+urn**2)
        supermin=min(np.min(Speed_t),np.min(Speed_n))
        supermax=max(np.max(Speed_t),np.max(Speed_n))/2
        if supermax>10:
            supermax=-gamma_t
        print(supermin)
        print(supermax)

        # --- Plotting
        fig,ax = plt.subplots(1,2, sharex=True, sharey=True)
        levels=np.linspace(supermin,supermax,10)
        im1=ax[0].contourf(Z/R,X/R,Speed_t, levels=levels) # the easiest way to get contourf to comply with cmap is to give levels
        im2=ax[1].contourf(Z/R,X/R,Speed_n, levels=levels) # the easiest way to get contourf to comply with cmap is to give levels
        # ax[0].plot([0,ZLIM[1]],[0,0],'k-',lw=2)
        # ax[1].plot([0,ZLIM[1]],[0,0],'k-',lw=2)
        # ax[0].plot([0,ZLIM[1]],[0,0],'k-',lw=2)
        # ax[1].plot([0,ZLIM[1]],[0,0],'k-',lw=2)

        # Streamlines and quiver
        yseed=np.linspace(XLIM[0]*1.1/R,XLIM[1]*0.9/R,7)
        start=np.array([yseed*0+(ZLIM[0]+ZLIM[1])/(2*R),yseed])
        sp=ax[0].streamplot(zs/R,xs/R,uzt,urt,color='k',start_points=start.T,linewidth=0.7,density=30,arrowstyle='-')
        qv=streamQuiver(ax[0],sp,n=[5,5,5,5,5,5,5],scale=40,angles='xy')

        yseed=np.linspace(XLIM[0]*1.1/R,XLIM[1]*0.9/R,7)
        start=np.array([yseed*0+(ZLIM[0]+ZLIM[1])/(2*R),yseed])
        sp=ax[1].streamplot(zs/R,xs/R,uzn,urn,color='k',start_points=start.T,linewidth=0.7,density=30,arrowstyle='-')
        qv=streamQuiver(ax[1],sp,n=[5,5,5,5,5,5,5],scale=40,angles='xy')
        ax[0].set_xlabel('z/R [-]')
        ax[1].set_xlabel('z/R [-]')
        ax[0].set_ylabel('r/R [-]')

        # fig.subplots_adjust(right=0.8)
        # cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
        # fig.colorbar(im1, cax=cbar_ax)

        ax[0].set_xlim(ZLIM/R)
        ax[0].set_ylim(XLIM/R)
        ax[1].set_xlim(ZLIM/R)
        ax[1].set_ylim(XLIM/R)



    # --- Relative error
    if 2 in IPlot:
        relmin=0
        relmax=10

        RelErrUz = np.abs((uzt-uzn)/uzt)*100
        RelErrUr = np.abs((urt-urn)/urt)*100

        fig,ax = plt.subplots(1,2, sharex=True, sharey=True,figsize=(6.4,3))
        fig.subplots_adjust(left=0.08, right=0.85, top=0.98, bottom=0.15, hspace=0.0, wspace=0.08)
    #     levels=[0,0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 10.0]
        levels=[0,0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0]
    #     levels=[0,1.5,5.0,10.0,15.0]
        im1=ax[0].contourf(Z/R, X/R, RelErrUz, levels=levels) 
        im2=ax[1].contourf(Z/R, X/R, RelErrUr, levels=levels)
    #     ax[1].axis('equal')

        theta=np.linspace(np.pi/2,np.pi,30)
        for r0 in [1,5,10]:
            ax[0].plot(r0*np.cos(theta), r0*np.sin(theta),'k--',lw=2)
            ax[1].plot(r0*np.cos(theta), r0*np.sin(theta),'k--',lw=2)
        cbar_ax = fig.add_axes([0.88, 0.14, 0.03, 0.84])
        fig.colorbar(im1, cax=cbar_ax)

        ax[0].set_xlabel('z/R [-]')
        ax[1].set_xlabel('z/R [-]')
        ax[0].set_ylabel('r/R [-]')

        ax[0].set_xlim(ZLIM/R)
        ax[0].set_ylim(XLIM/R)
        ax[1].set_xlim(ZLIM/R)
        ax[1].set_ylim(XLIM/R)








if __name__ == '__main__':
    main()
    plt.show()

