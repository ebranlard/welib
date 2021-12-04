import unittest
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# Local 
from welib.CFD.axisym_streamf import *


def testStreamFunctionsForVortexDistribution(distribution='cylinder', plot=False, addU0=True, nr=46, nz=51, verbose=False):
    """ 
    helper function
    """
    from welib.vortilib.elements.VortexAxisymmetric import axisym_predefined_distributions, axisym_u
    from welib.tools.curves import streamQuiver
    #distribution='cylinder'
    # distribution='singular_ring'
    # distribution='regularized_ring'
    #distribution='zero'
    params=dict()
    params['gamma'] = -1
    params['z'] = 0
    params['r'] =1
    params['l'] =100
    params['Gamma'] =-1 # for ring
    params['rc'] =0.3  # core radius for regularized ring
    params['k'] = 2    # patch factor for regularized ring

    zmin = -1.5
    zmax = 1.5
    rmin = 0
    rmax = 3 
    r = np.linspace(rmin,rmax,nr)
    z = np.linspace(zmin,zmax,nz)
    dz=(zmax-zmin)/(nz-1)
    dr=(rmax-rmin)/(nr-1)
    psi = np.zeros((nz,nr))
    om  = np.zeros((nz,nr))
    ur  = np.zeros((nz,nr))
    uz  = np.zeros((nz,nr))
    Rcp,Zcp = np.meshgrid(r,z)

    # --- Baseline
    om, ur1, uz1  = axisym_predefined_distributions(r, z, distribution=distribution, velocity=True, params=params)
    Rcp,Zcp = np.meshgrid(r,z)
    if distribution=='zero' or (distribution =='cylinder' and params['l']>10):
        uz0=uz1
        ur0=ur1
    else:
        ur0, uz0 = axisym_u(Rcp, Zcp, r, z, om)
    if addU0:
        uz0=uz0+1
    if distribution=='zero':
        U0=1
        psi0 = 0.5*U0*Rcp**2
    else:
        psi0=None

    # --- Numerical solution
    BC_Walls={'south':BC_Dir_0, 'west':BC_Dir_uz, 'north':BC_Neu_uz, 'east':BC_Neu_ur}
#     BC_Walls={'south':BC_Dir_0, 'west':BC_Neu_uz, 'north':BC_Neu_uz, 'east':BC_Neu_ur}
#     BC_Walls={'south':BC_Dir_0, 'west':BC_Neu_ur, 'north':BC_Neu_uz, 'east':BC_Neu_ur}

    psi   = streamfunction(om, r, psi, ur0, uz0, dr, dz, BC_Walls=BC_Walls, verbose = verbose)
    ur,uz = velocity_from_streamfunction(psi, r, dr, dz)

    # --- Plot
    if plot:
        vel  = np.sqrt(ur**2 + uz**2)
        vel0 = np.sqrt(ur0**2 + uz0**2)
        minSpeed=np.nanmin(vel0)*0.7
        maxSpeed=np.nanmax(vel0)*1.3
        if minSpeed==maxSpeed:
            minSpeed-=1
            maxSpeed+=1
        levels=np.linspace(minSpeed,maxSpeed,50)
        rStart=np.array([0.25,0.5,0.75,1.25,1.5,1.75])
        start=np.array([rStart*0+0.5,rStart])

        fig,axes = plt.subplots(3, 1, sharex=True, figsize=(6.4,8.8)) # (6.4,4.8)
        fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
        ax=axes[0]
        img = ax.contourf(z, r, vel.T, levels=levels, vmin=minSpeed, vmax=maxSpeed)
        sp=ax.streamplot(z, r, uz.T, ur.T,color='k',start_points=start.T,linewidth=0.7,density=30,arrowstyle='-')
        qv=streamQuiver(ax,sp,n=5,scale=40,angles='xy')
        ax.set_ylabel('r [m]')
        ax.set_title('Finite diff')
        fig.colorbar(img, ax=ax)

        ax=axes[1]
        img = ax.contourf(z, r, vel0.T, levels=levels, vmin=minSpeed, vmax=maxSpeed)
        sp=ax.streamplot(z, r, uz0.T, ur0.T,color='k',start_points=start.T,linewidth=0.7,density=30,arrowstyle='-')
        qv=streamQuiver(ax,sp,n=5,scale=40,angles='xy')
        fig.colorbar(img, ax=ax)
        ax.set_ylabel('r [m]')
        ax.set_xlabel('z [m]')
        ax.set_title('Analytical function')

        ax=axes[2]
        # img = ax.contourf(z, r, psi.T)
        img = ax.pcolormesh(z, r, psi.T, shading='auto')
        fig.colorbar(img, ax=ax)
        ax.set_ylabel('r [m]')
        ax.set_xlabel('z [m]')
        ax.set_title('psi')

        fig,axes = plt.subplots(5, 2, sharey=False, figsize=(6.4,8.8)) # (6.4,4.8)
        fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
        ax=axes[0,0]
        ax.plot(r, uz0[0,:],'k-'    , label='uz Analytical')
        ax.plot(r, uz [0,:],'--'   , label='uz Numerical')
        ax.set_title('uz - Inlet velocity')
        ax=axes[0,1]
        ax.plot(r, ur0[0,:],'k-'    , label='ur Analytical')
        ax.plot(r, ur [0,:],'--'   , label='ur Numerical')
        ax.set_xlabel('r [m]')
        ax.set_ylabel('u [m/s]')
        ax.set_title('ur - Inlet velocity')

        ax=axes[1,0]
        ax.plot(r, uz0[nz-1,:],'k-'    , label='uz Analytical')
        ax.plot(r, uz [nz-1,:],'--'   , label='uz Numerical')
        ax=axes[1,1]
        ax.plot(r, ur0[nz-1,:],'k-'    , label='ur Analytical')
        ax.plot(r, ur [nz-1,:],'--'   , label='ur Numerical')
        ax.set_xlabel('r [m]')
        ax.set_ylabel('u [m/s]')
        ax.set_title('Outlet velocity')

        ax=axes[2,0]
        ax.plot(z, uz0[:,0],'k-'    , label='uz Analytical')
        ax.plot(z, uz [:,0],'--'   , label='uz Numerical')
        ax=axes[2,1]
        ax.plot(z, ur0[:,0],'k-'    , label='ur Analytical')
        ax.plot(z, ur [:,0],'--'   , label='ur Numerical')
        ax.set_xlabel('z [m]')
        ax.set_ylabel('u [m/s]')
        ax.set_title('Axis velocity')

        ax=axes[3,0]
        ax.plot(z, uz0[:,int(nr/3)], 'k-'    , label='uz Analytical')
        ax.plot(z, uz [:,int(nr/3)], '--'   , label='uz Numerical')
        ax=axes[3,1]
        ax.plot(z, ur0[:,int(nr/3)], 'k-'    , label='ur Analytical')
        ax.plot(z, ur [:,int(nr/3)], '--'   , label='ur Numerical')
        ax.set_xlabel('z [m]')
        ax.set_ylabel('u [m/s]')
        ax.set_title('Mid velocity')

        ax=axes[4,0]
        ax.plot(z, uz0[:,nr-1],'k-'    , label='uz Analytical')
        ax.plot(z, uz [:,nr-1],'--'   , label='uz Numerical')
        ax=axes[4,1]
        ax.plot(z, ur0[:,nr-1],'k-'    , label='ur Analytical')
        ax.plot(z, ur [:,nr-1],'--'   , label='ur Numerical')
        ax.set_xlabel('z [m]')
        ax.set_ylabel('u [m/s]')
        ax.set_title('Top velocity')
        ax.legend()

        #fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
        #fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
        #ax.plot(r, 0.5*uz0[0,0]*r**2,'k-'    , label='psi th')
        #ax.plot(r,       psi[0,:]   ,'--'    , label='psi')
        #ax.plot(r[1:], sciint.cumtrapz(r*uz[0,:], r), ':', label='psi2')
        #ax.set_xlabel('r')
        #ax.set_ylabel('psi')
        #ax.legend()

    return ur, uz, ur0, uz0, psi, psi0




# --------------------------------------------------------------------------------}
# --- TEST  
# --------------------------------------------------------------------------------{
class Test(unittest.TestCase):

    def test_velocity_from_psi(self):
        # Uniform flow, psi=1/2 U0 r^2
        U0   = 10
        nr   = 16  # number of nodes
        nz   = 11
        zmin = -1.5
        zmax = 1.5
        rmin = 0
        rmax = 3
        r   = np.linspace(rmin,rmax,nr)
        z   = np.linspace(zmin,zmax,nz)
        dz  = (zmax-zmin)/(nz-1)
        dr  = (rmax-rmin)/(nr-1)
        ur_ref  = np.zeros((nz,nr))
        uz_ref  = np.zeros((nz,nr))
        Rcp,Zcp = np.meshgrid(r,z)
        uz_ref[:,:]=U0

        psi = 0.5*U0*Rcp**2
        ur,uz = velocity_from_streamfunction(psi, r, dr, dz)

        np.testing.assert_almost_equal(ur,ur_ref)
        np.testing.assert_almost_equal(uz,uz_ref)


    def test_streamfunction(self):
        # --- Test on uniform inflow (zero vorticity, constant velocity=1) 
        ur, uz, ur0, uz0, psi, psi0 = testStreamFunctionsForVortexDistribution(distribution='zero', plot=False, addU0=False, nr=7, nz=5, verbose=False)
        np.testing.assert_allclose(psi,psi0)

        # --- Ring
        ur, uz, ur0, uz0, _,_ = testStreamFunctionsForVortexDistribution(distribution='regularized_ring', plot=False, addU0=False, nr=31, nz=15, verbose=False)
        # Inlet velocity
        np.testing.assert_allclose(uz[0,:],uz0[0,:], rtol=0.2)
        # Axis velocity
        np.testing.assert_allclose(uz[:,0],uz0[:,0], rtol=0.1)
        np.testing.assert_allclose(ur[:,0],ur0[:,0], rtol=0.1)
        # Top velocity
        np.testing.assert_allclose(uz[:,-1],uz0[:,-1], rtol=0.1)
        #np.testing.assert_allclose(ur[:,-1],ur0[:,-1], rtol=0.1)

        # --- Cylinder
        ur, uz, ur0, uz0, _,_ = testStreamFunctionsForVortexDistribution(distribution='cylinder', plot=True, addU0=True, nr=31, nz=15, verbose=False)
        # Inlet velocity
        np.testing.assert_allclose(uz[0,:],uz0[0,:], rtol=0.25)
        # Axis velocity
        np.testing.assert_allclose(uz[:,0],uz0[:,0], rtol=0.26)
        np.testing.assert_allclose(ur[:,0],ur0[:,0], rtol=0.1)
        # Top velocity
        #np.testing.assert_allclose(uz[:,-1],uz0[:,-1], rtol=0.21)
        #np.testing.assert_allclose(ur[:,-1],ur0[:,-1], rtol=0.1)
        #plt.show()

if __name__=='__main__':
    unittest.main()


