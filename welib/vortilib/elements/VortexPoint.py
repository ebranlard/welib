""" 
2D vortex points 
For 3D see VortexParticles

"""
import numpy as np


def vp_u(CP,Pv,Gamma=1,rcore=0,bViscous=False): 
    """ 
    Induced velocity by one 2D vortex point on one Control Point (CP)
    CP: position of control point
    Pv: position of vortex
    """
    DP = CP - Pv
    r2 = DP[0] ** 2 + DP[1] ** 2
    t = np.array([- DP[1],DP[0]])
    if bViscous:
        U = Gamma * t[0] / r2 / (2 * np.pi) * (1 - np.exp(- r2 / rcore ** 2))
        V = Gamma * t[1] / r2 / (2 * np.pi) * (1 - np.exp(- r2 / rcore ** 2))
    else:
        U = Gamma * t[0] / r2 / (2 * np.pi)
        V = Gamma * t[1] / r2 / (2 * np.pi)
    return U,V


def vps_u(CP,XV,Gammas,SmoothModel=0,KernelOrder=2,SmoothParam=None):
    """
    low level interface, for N (2D?) point Vortices

    # Takes matrices as input of size ( : , 3) or (:,2)

    CP: control points where the velocity is computed
# XV: location of point vortex

    [Ui]= fUi_PointVortex2DN([1 0],[0 0],2*np.pi,[],0),[Ui]= fUi_PointVortex2DN([0 1],[0 0],2*np.pi,[],0)
    """

    ncp    = CP.shape[0]
    nv     = XV.shape[0]
    ndim   = XV.shape[1]
    ndimCP = XV.shape[1]
    if ndim != ndimCP:
        raise Exception('Your control points and vortex point do not have the same dimension.')
    if SmoothParam is not None:
        if not hasattr(SmoothParam,'__len__'):
            SmoothParam=np.ones(CP.shape[0])*SmoothParam



    Ui = np.zeros((ncp,ndim))

    # ---Setting up Kernel/smooth param and Exp functions for Smooth model
    if SmoothModel==0:
        fKernel = None
        fE      = None
        frho2   = None
    elif SmoothModel==1:
        frho2 = lambda r2, iv: r2/ SmoothParam[iv]**2
        fE = lambda rho2 : np.exp(- rho2)
        if KernelOrder==2:
            fKernel = lambda rho2 :  1
        elif KernelOrder==4:
            fKernel = lambda rho2 :  1 - rho2
        elif KernelOrder==6:
            fKernel = lambda rho2 :  1 - 2 * rho2 + rho2 ** 2 / 2
        elif KernelOrder==8:
            fKernel = lambda rho2 : 1 - 3 * rho2 + 3 * rho2 ** 2 / 2 - rho2 ** 3 / 6
        else:
            raise Exception('fKernel order not implemented for Majda model')
    elif SmoothModel==2:
        frho2 = lambda r2, iv: r2/ SmoothParam[iv]**2
        fE = lambda rho2 : np.exp(- rho2/2) # NOTE divided by 2
        if KernelOrder==2:
            fKernel = lambda rho2 :  1
            fE = lambda rho2 : np.exp(- rho2) # NOTE not divided by 2
        elif KernelOrder==4:
            fKernel = lambda rho2 :  1 - 1 / 2 * rho2
        elif KernelOrder==6:
            fKernel = lambda rho2 :  1 - rho2 + 1 / 8 * rho2 ** 2
        elif KernelOrder==8:
            fKernel = lambda rho2 : 1 - 3 / 2 * rho2 + 3 / 8 * rho2 ** 2 - 1 / 48 * rho2 ** 3
        elif KernelOrder==10:
            fKernel = lambda rho2 : 1 - 2 * rho2 + 3 / 4 * rho2 ** 2 - 1 / 12 * rho2 ** 3 + 1 / 384 * rho2 ** 4
        else:
            raise Exception('Kernel order not implemented for Gaussian2')
    else:
        raise Exception('Unknown smooth model')


    # loop on control point
    if ndim==2:
        for icp in np.arange(ncp):
            x0 = CP[icp,:]
            # Loop on vortex
            for iv in np.arange(nv):
                xv = XV[iv,:]
                X0XV = x0 - xv
                r2 = X0XV[0] * X0XV[0] + X0XV[1] * X0XV[1]
                if r2 < 1e-15:
                    # We escape No matter the smooth model and order
                    pass
                else:
                    t = np.array([-X0XV[1],X0XV[0]])
                    ## Viscous Model
                    if (SmoothModel == 0 or KernelOrder == 0):
                        Ui[icp,:] = Ui[icp,:] + Gammas[iv] / (2 * np.pi) * t / r2
                    else:
                        Qm = 0
                        if SmoothModel==0:
                            Ui[icp,:] = Ui[icp,:] + Gammas[iv] / (2 * np.pi) * t / r2
                        else:
                            # Smooth model 1 or 2
                            rho2 = frho2(r2,iv)
                            E    = fE(rho2)
                            Qm   = fKernel(rho2)
                            Ui[icp,:] = Ui[icp,:] + Gammas[iv]/(2*np.pi)*t*(1-Qm*E)/r2
    else:
        raise NotImplementedError()
    # fprintf('\n');
    
    return Ui
