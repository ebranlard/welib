import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from vortilib.elements.InviscidVortexPatch import ivp_omega
from vortilib.elements.VortexPoint import vps_u
from vortilib.mesh import mesh
from vortilib.particles import particles
from vortilib.particles import initialization
from vortilib.particles import projection

from vortilib.maths.vectoranalysis import *

np.set_printoptions(linewidth=500)
# Mesh Params
nx = 30
ny = 33
# nx = 7
# ny = 5
vx = np.linspace(-1.0,1.0,nx)
vy = np.linspace(-1.0,1.0,ny)
XLIM = np.array([-1.5,1.5])
YLIM = XLIM
ZLIM = np.array([0,1])

# --- Setting up Mesh
mesh = mesh.Mesh(vx,vy)
# print(mesh)

# --- Setting up particles
k_patch = 3
fOmega = lambda x,y : ivp_omega(x,y,k_patch,False)
Part = initialization.init_from_mesh_and_fun(fOmega, mesh, location='CellCenter')
Part.removeZeroVorticity()

# --- Selecting Vortex elements for roll up
# CPs_n=Part.P;
# ## Roll-up/convection computation
# [ ConvDx, Utot_n, Grad_n, Uinf_n] = fCompute_dx(CPs_n, Time,  Wind, Algo );
# ## Effective convection
# Part.P=Part.P+ConvDx;

# plt.figure(1)
# ax = plt.axes(projection='3d')
# ax.scatter3D(Part.P[:,0],Part.P[:,1],Part.Intensity / Part.Volume,c=Part.Intensity)#,'+')
# ax.set_xlim(XLIM)
# ax.set_ylim(YLIM)
# ax.set_zlim(ZLIM)
# 
# --- Particle2Mesh
MeshValues,v_p = projection.part2mesh(Part,mesh)
print('MShape',MeshValues.shape)
X,Y = np.meshgrid(vx,vy)
M = np.squeeze(MeshValues[0,:,:])
V = np.squeeze(MeshValues[1,:,:])
M2 = np.transpose(M)
V2 = np.transpose(V)
# omega_proj = M2
# # V2(V2==0)=cell_area;
omega_proj = M2 / V2
# # omega=M2/cell_area;
# 
# --- Computing Ui on mesh
X_flat = X.ravel()
Y_flat = Y.ravel()
print('v_p',v_p.shape)
Gammas = v_p[:,0]
CP = np.column_stack((X_flat,Y_flat))
print('CP',CP.shape)
SmoothModel = 0
SmoothParam = 0
KernelOrder = 0
print('Computing...')
Ui_num_On_Mesh = vps_u(CP,Part.P[:,:2],Gammas,SmoothModel,KernelOrder,SmoothParam)
print('Done!')
U2 = np.reshape(Ui_num_On_Mesh[:,0],(ny,nx))
V2 = np.reshape(Ui_num_On_Mesh[:,1],(ny,nx))
speed2 = np.sqrt(U2 ** 2 + V2 ** 2)
omega_z_num_calc,__ = matlab_curl_2d(X,Y,U2,V2)

plt.figure()
ax = plt.axes(projection='3d')
ax.plot_surface(X,Y,speed2, edgecolor='none',alpha=1,antialiased=False)
plt.title('Velocity norm')


# --- 
plt.figure()
ax = plt.axes(projection='3d')
ax.plot_surface(X,Y,omega_proj, edgecolor='none',alpha=1,antialiased=False)
ax.scatter3D(Part.P[:,0],Part.P[:,1],Part.Intensity / Part.Volume,c=Part.Intensity)
plt.title('Omega projected on grid and particles')

plt.figure()
ax = plt.axes(projection='3d')
ax.plot_surface(X,Y,omega_z_num_calc, edgecolor='none',alpha=1,antialiased=False)
ax.scatter3D(Part.P[:,0],Part.P[:,1],Part.Intensity / Part.Volume,c=Part.Intensity)
plt.title('Omega from Ui computation')

# --- Mesh2Particles
Part2=initialization.init_from_mesh(mesh,location='CellCenter')
Part2.removeZeroVorticity();
# # #  omega_z_th=fOmega_InviscidVortexPatch(Part.P(:,1),Part.P(:,2),k_patch,false);
# 
# 
plt.figure()
ax = plt.axes(projection='3d')
ax.scatter3D(Part.P[: ,0],Part.P[: ,1],Part.Intensity/Part.Volume,marker = 'o',c = 'r',s=3)
ax.scatter3D(Part2.P[:,0],Part2.P[:,1],Part2.Intensity/Part2.Volume,marker = '+',c = 'b')
# 
# 
# 
# # ##
# # figure,
# # contour(X,Y,omega_proj-omega_z_num_calc)
# # colorbar
# # ##
# # # InitCell
# # # # v_p2= m2p_mp4_2d(Part,mesh,xBox,dBox,nx,ny);
# # # #
# # # figure
# # # plot(v_p(1,:)-v_p2(1,:)); title('Intensity')
# # # figure
# # # plot((v_p(2,:)-v_p2(2,:)))/v_p(2,1); title('Volume')
# # # dispatchFigs(1)
# # # ##
# # #
# # # InitCell
# # # PartMesh=zeros(3,nx*ny);
# # # PartMesh(1,:)=X(:);
# # # PartMesh(2,:)=Y(:);
# # # v_p3= m2p_mp4_2d(PartMesh,mesh,xBox,dBox,nx,ny);
# # #
# # # figure, hold on
# # # surf(X,Y,omega_proj)
# # # plot3(PartMesh(1,:),PartMesh(2,:),v_p3(1,:)./v_p3(2,:),'ko')
# # #     figure(2), clf
# # # #     plot(Part.P(:,1),Part.P(:,2),'+')
# # #     plot3(Part.P(:,1),Part.P(:,2),Part.Intensity(:),'+')
# # #     xlim(XLIM);
# # #     ylim(YLIM);
# # #     zlim(ZLIM);
# # #     pause(0.5)
# # figure,
# # surf(squeeze(MeshValues(1,:,:))./squeeze(MeshValues(2,:,:)))
# 
plt.show()
