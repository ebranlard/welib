""" 
Vorticity and tangential velocity for a 2D inviscid vortex patch

See:
 [1] Chapter 33, p.402, Branlard - Wind turbine aerodynamics and vorticity based methods, Springer 2017

"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from welib.tools.colors import fColrs
from welib.tools.clean_exceptions import *
from welib.tools.curves import streamQuiver
from welib.vortilib.elements.InviscidVortexPatch import *
try:
    from numpy import trapezoid
except:
    from numpy import trapz as trapezoid

def main():

    fig,axes = plt.subplots(2, 2, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
    fig.subplots_adjust(left=0.09, right=0.92, top=0.88, bottom=0.09, hspace=0.47, wspace=0.29)

    # --- Plot vorticity distribution for different k values
    r=np.linspace(0,2,100)
    theta=r*0
    ax = axes[0,0]
    for i,k in enumerate([1/2,1,2,4]):
        omega = ivp_omega(r, theta, k=k, polarIn=True)
        ax.plot(r, omega, color=fColrs(i+1), label=r'$k={}$'.format(k))
    ax.set_xlabel('$r$ [m]')
    ax.set_ylabel(r'$\omega_z$ [1/s]')
    ax.autoscale(tight=True)
    ax.legend()
    ax.tick_params(direction='in')
    ax.set_title('Vorticity')


    # --- Plot tangential velocity for different k values
    r=np.linspace(0,2,100)
    theta=r*0
    ax = axes[0,1]
    ax.plot(r, r/2,'k--', label=r'$r/2$ slope')
    for i,k in enumerate([1/2,1,2,4]):
        ur,utheta = ivp_u(r, theta, k=k, polarIn=True, polarOut=True)
        ax.plot(r, utheta, color=fColrs(i+1), label=r'$k={}$'.format(k))
    ax.set_xlabel('$r$ [m]')
    ax.set_ylabel(r'$u_\theta$ [m/s]')
    ax.autoscale(tight=True)
    ax.set_ylim([0,0.35])
    ax.legend()
    ax.tick_params(direction='in')
    ax.set_title('Velocity')


    # --- Plot circulation 
    r     = np.linspace(0,2,100)
    ax = axes[1,0]
    for i,k in enumerate([1/2,1,2,4]):
        Gamma = ivp_Gamma(r, k=k)
        ax.plot(r, Gamma/np.pi, color=fColrs(i+1), label=r'$k={}$'.format(k))

    ax.set_xlabel('$r$ [m]')
    ax.set_ylabel(r'$\Gamma/\pi$ [m$^2$/s]')
    ax.autoscale(tight=True)
    ax.set_ylim([0,0.7])
    ax.legend()
    ax.tick_params(direction='in')
    ax.set_title('Circulation')


    # --- Plot velocity field
    k=1
    Gamma=ivp_Gamma([3],k=k)
    #  Control points
    nX=100
    nY=101
    minSpeed=0
    maxSpeed=1 # Scaled by max
    vX = np.linspace(-4,4,nX)
    vY = np.linspace(-4,4,nY)
    XCP,YCP = np.meshgrid(vX,vY)
    Xcp = XCP.flatten()
    Ycp = YCP.flatten()
    Zcp = Xcp*0
    # Velocity field
    Ux, Uy = ivp_u(Xcp, Ycp, k=k)
    Ux    = Ux.reshape(XCP.shape)
    Uy    = Uy.reshape(XCP.shape)
    Speed = np.sqrt((Ux**2+Uy**2))
    print('min: ',np.min(Speed.ravel()),' - max: ',np.max(Speed.ravel()))
    Speed= Speed/ np.max(Speed) # TODO can easiy be computed analytically

    # Plot
    ax = axes[1,1]
    im = ax.contourf(XCP, YCP, Speed, levels=np.linspace(minSpeed,maxSpeed,250), vmin=minSpeed, vmax=maxSpeed)
    cb=fig.colorbar(im)
    yseed=np.linspace(0.1,3.8,15)
    start=np.array([yseed*0,yseed])
    sp=ax.streamplot(vX,vY,Ux,Uy,color='k',start_points=start.T,linewidth=0.7,density=30,arrowstyle='-')
    qv=streamQuiver(ax,sp,n=7,scale=40,angles='xy')
    ax.set_xlabel('x [m]')
    ax.set_ylabel('y [m]')
    ax.set_aspect('equal','box')
    ax.set_title('Velocity Field')

    fig.suptitle('Vortilib - Inviscid Vorticity Patch')



if __name__ == '__main__':
    main()
    plt.show()
if __name__=="__test__":
    # --- Check that circulation match analytical value
    def circulationSurvey(r, nTheta=100, k=1):
        theta=np.linspace(0,2*np.pi,nTheta+1)
        dTheta=theta[1]-theta[0]
        Xcp=r*np.cos(theta)
        Ycp=r*np.sin(theta)
        Zcp=0*Xcp
        Ux, Uy = ivp_u(Xcp, Ycp, k=k)
        Ut = Uy * np.cos(theta) - Ux * np.sin(theta)
        GammaTheory = ivp_Gamma([r], k=k)[0]
        #GammaCalc   = 2*np.pi * r*Ut[0]
        GammaCalc = r* trapezoid(Ut,theta)
        return  GammaCalc, GammaTheory

    np.testing.assert_almost_equal(*circulationSurvey(0.1))
    print(circulationSurvey(0.1))
    print(circulationSurvey(0.5))
    print(circulationSurvey(0.9))
    print(circulationSurvey(1.0))
    print(circulationSurvey(2.0))
    print(circulationSurvey(3.0))
if __name__=="__export__":
    main()
    from welib.tools.repo import export_figs_callback
    export_figs_callback(__file__)

