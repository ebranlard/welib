"""
Compare the velocity field from a vortex doublet line obtained using numerical or analytical integration

"""
# --- General
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable
# --- Local
from welib.vortilib.elements.VortexDoublet import *
from welib.vortilib.elements.VortexCylinderSkewed import *
from welib.tools.curves import streamQuiver
from welib.tools.tictoc import Timer

def main():
    # --- Parameters
    bPolar = False
    m      = np.tan(-30*np.pi/180) 
    R      = 1
    U0     = 1
    a      = 0.3
    gamma_t=-2*a*U0
    XLIM  = np.asarray([-5*R,5*R ]) # 
    ZLIM  = np.asarray([-5*R,4*R]) # 
    ZMax  = 100*R  # Extent of line for numerical integration
    nQuad = 1000   # Numbero of quadrature points for numerical integration
    nx    = 200                    # Number of points for velocity evaluation
    nz    = nx

    # --- Derived parameters
    dmz_dz = gamma_t * R**2 * np.pi # doublet intensity per length

    # --- Flow field and speed
    zs = np.linspace(ZLIM[0],ZLIM[1],nx).astype(np.float32)
    xs = np.linspace(XLIM[0]*1.1,XLIM[1]*1.1,nx).astype(np.float32)
    [Z,X]=np.meshgrid(zs,xs)
    Y=X*0

    ux    = [0,0,0]
    uz    = [0,0,0]
    Speed = [0,0,0]

    labs=['Vortex doublet theory','Vortex doublet numerical','Vortex cylinder']
    with Timer(labs[0]):
        ux[0],_,uz[0] = doublet_line_u    (X,Y,Z,dmz_dz, m=m)
    with Timer(labs[1]):
        ux[1],_,uz[1] = doublet_line_u_num(X,Y,Z,dmz_dz, m=m, zmax=ZMax, nQuad=nQuad)
    with Timer(labs[2]):
        ux[2],_,uz[2] = svc_tang_u(X,Y,Z,gamma_t=gamma_t,R=R,m=m,ntheta=180,polar_out=False)

    for i in np.arange(3):
        Speed[i]=(uz[i] + U0)/U0
        supermin=np.min(Speed[i])
        supermax=np.max(Speed[i])
        print(supermin)
        print(supermax)
    supermin =  0
    supermax =  1.1



    # --- Plotting
    fig,ax = plt.subplots(1, 3, sharey=True, figsize=(14.4,4.8)) # (6.4,4.8)
    # fig.subplots_adjust(left=0.12, right=0.80, top=0.88, bottom=0.11, hspace=0.20, wspace=0.20)
    levels=np.linspace(supermin,supermax,25)
    for i in np.arange(3):
        im=ax[i].contourf(Z/R,X/R,Speed[i], levels=levels) # the easiest way to get contourf to comply with cmap is to give levels
    for i in np.arange(2):
        ax[i].plot([0,ZLIM[1]],[0,m*ZLIM[1]],'k-',lw=2)
    ax[2].plot([0,ZLIM[1]],[1,1+m*ZLIM[1]],'k-',lw=2)
    ax[2].plot([0,ZLIM[1]],[-1,-1+m*ZLIM[1]],'k-',lw=2)

    # Streamlines and quiver
    for i in np.arange(3):
        yseed=np.linspace(XLIM[0]*0.9/R,XLIM[1]*0.9/R,11)
        start=np.array([yseed*0+ZLIM[0]*0.8/(R),yseed])
        sp=ax[i].streamplot(zs/R,xs/R,uz[i]+U0,ux[i],color='k',start_points=start.T,linewidth=0.7,density=30,arrowstyle='-')
        qv=streamQuiver(ax[i],sp,n=[5]*len(yseed),scale=40,angles='xy')

    #     yseed=np.linspace(XLIM[0]*0.9/R,XLIM[1]*0.9/R,11)
    #     start=np.array([yseed*0+ZLIM[1]*0.5/(R),yseed])
    #     sp=ax[i].streamplot(zs/R,xs/R,uz[i]+U0,ux[i],color='k',start_points=start.T,linewidth=0.7,density=30,arrowstyle='-')
    #     qv=streamQuiver(ax[i],sp,n=[5]*len(yseed),scale=40,angles='xy')

        ax[i].set_xlabel('z/R [-]')
        ax[i].set_xlim(ZLIM/R)
        ax[i].set_ylim(XLIM/R)
        ax[i].set_title(labs[i])
    ax[0].set_ylabel('x/R [-]')

    fig.subplots_adjust(right=0.92,left=0.06)
    cbar_ax = fig.add_axes([0.93, 0.11, 0.01, 0.8]) #
    fig.colorbar(im, cax=cbar_ax)

if __name__ == '__main__':
    main()
    plt.show()


