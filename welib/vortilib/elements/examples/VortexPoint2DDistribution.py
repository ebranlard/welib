""" 
Compute velocity field from a distribution of vorticity in 2D using a set of vortex point

NOTE: see also vorticilib.particles.examples, for methods "vortex particles" method

"""
import numpy as np
import matplotlib.pyplot as plt
# Local 
from vortilib.tools.curves import streamQuiver
# Analytical vorticity distributions
from vortilib.elements.VortexPoint           import * # vps_u
from vortilib.elements.InviscidVortexPatch   import * # ivp_u, ivp_omega
from vortilib.elements.LambOseen             import * # lo_omega
from vortilib.elements.VortexPatch2DGaussian import * # vpg_omega

# --- Parameters
minSpeed=0
nX=30
nY=31
distribution='LambOseen'
# distribution='InviscidVortexPatch'
distribution='VortexPatchGaussian'
distribution='VortexPatchGaussianAsym'
SmoothModel = 0   # Can be inviscid=0 for this example since we evaluate on Vortex points
SmoothParam = 0.1 # as function of diagonal length (should be <1 since no projection of omega)
KernelOrder = 2

# Control points
vX = np.linspace(-8,8,nX)
vY = np.linspace(-8,8,nY)
XCP,YCP = np.meshgrid(vX,vY)
Xcp = XCP.flatten()
Ycp = YCP.flatten()
Area = (vX[1]-vX[0])*(vY[1]-vY[0]) # Area occupied by each particles
dr   = np.sqrt((vX[1]-vX[0])**2 + (vY[1]-vY[0])**2) # Diagonal of cells, used for regularization
# Particles, put at control points...
PartP     = np.column_stack((Xcp, Ycp))
u_th = None
if distribution=='LambOseen':
    omega=lo_omega
    u_th =lo_u
    maxSpeed=0.0508
elif distribution=='InviscidVortexPatch':
    omega=ivp_omega
    u_th =ivp_u
    maxSpeed=0.27
elif distribution=='VortexPatchGaussian':
    omega    = vpg_omega
    u_th     = vpg_u
    maxSpeed = 0.102
elif distribution=='VortexPatchGaussianAsym':
    omega    = vpga_omega
    u_th     = vpga_u
    maxSpeed = 0.057
else:
    raise NotImplementedError('')

PartGammas = omega(Xcp, Ycp)*Area


# --- Computing U
print('Computing...')
U_num = vps_u(PartP,PartP[:,:2],PartGammas,SmoothModel,KernelOrder,SmoothParam=SmoothParam*dr)
print('Done!')
Ux_num = np.reshape(U_num[:,0],(nY,nX))
Uy_num = np.reshape(U_num[:,1],(nY,nX))
Speed_num = np.sqrt(Ux_num**2 + Uy_num**2)
#omega_z_num_calc,__ = matlab_curl_2d(X,Y,U2,V2)
print('min: ',np.min(Speed_num.ravel()),' - max: ',np.max(Speed_num.ravel()))


fig,axes = plt.subplots(1, 2, sharey=False, figsize=(6.2,4.6)) # (6.4,4.8)
fig.subplots_adjust(left=0.12, right=0.96, top=0.96, bottom=0.12, hspace=0.20, wspace=0.30)
ax=axes[0]
im = ax.contourf(XCP, YCP, Speed_num, levels=np.linspace(minSpeed,maxSpeed,250), vmin=minSpeed, vmax=maxSpeed)
yseed=np.linspace(0.1,3.0,10)
start=np.array([yseed*0,yseed])
sp=ax.streamplot(vX,vY,Ux_num,Uy_num,color='k',start_points=start.T,linewidth=0.7,density=30,arrowstyle='-')
qv=streamQuiver(ax,sp,n=7,scale=40,angles='xy')

ax.set_xlabel('x [m]')
ax.set_ylabel('y [m]')
ax.set_aspect('equal','box')
ax.set_title('{} - Numerical'.format(distribution))

# --- Analytical results
if u_th is not None:
    Ux_th, Uy_th = u_th(PartP[:,0],PartP[:,1])
    Ux_th = np.reshape(Ux_th,(nY,nX))
    Uy_th = np.reshape(Uy_th,(nY,nX))
    Speed_th = np.sqrt(Ux_th**2 + Uy_th**2)
    print('min: ',np.min(Speed_th.ravel()),' - max: ',np.max(Speed_th.ravel()))
    ax=axes[1]
    im = ax.contourf(XCP, YCP, Speed_th, levels=np.linspace(minSpeed,maxSpeed,250), vmin=minSpeed, vmax=maxSpeed)
    sp=ax.streamplot(vX,vY,Ux_th,Uy_th,color='k',start_points=start.T,linewidth=0.7,density=30,arrowstyle='-')
    qv=streamQuiver(ax,sp,n=7,scale=40,angles='xy')
    ax.set_xlabel('x [m]')
    ax.set_aspect('equal','box')
    ax.set_title('Theory')
# cb=fig.colorbar(im)



fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.2,4.6)) # (6.4,4.8)
fig.subplots_adjust(left=0.12, right=0.98, top=0.96, bottom=0.12, hspace=0.20, wspace=0.20)
ix=np.argmin(np.abs(vX))
iy=np.argmin(np.abs(vY))
ax.plot(vY,Ux_num[:,ix], label='x=0 num')
ax.plot(vX,Uy_num[iy,:], label='y=0 num')
ax.plot(vY,Ux_th[:,ix] , label='x=0 th')
ax.plot(vX,Uy_th[iy,:] , label='y=0 th')
ax.legend()







plt.show()

