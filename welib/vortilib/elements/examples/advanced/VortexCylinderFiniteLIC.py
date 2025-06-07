"""
References:
    [1] E. Branlard - Wind Turbine Aerodynamics and Vorticity Based Method, Springer, 2017
"""
# --- General
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable
import os
# --- Local
from welib.vortilib.elements.VortexCylinder import cylinder_tang_u

from welib.tools.colors import darkrainbow as cmap
from welib.tools.colors import darkrainbow, adjust_color_lightness, manual_colorbar
from welib.tools.curves import streamQuiver
from welib.tools.lic    import lic, licImage
import random
random.seed( 30 )
# --- Parameters
R=1
z1=-2*R
z2=2*R
XLIM=[-3*R,3*R] # 
ZLIM=[-5*R,5*R] # 
CLIM=[0.0,1]
gamma_t=-1
nx = 400    # Number of points for velocity evaluation
nz = 201

# LIC params
nLICKernel=31
offset=0.4
spread=1.1
accentuation=1.5


# Movie params
Amplitude=1.0
anim_meth='kernelphaseOffset' # timesfreq, sinepure, kernelphase
movie=False
dpi = 300
freq=5 
TMAX=5/freq
nFrames=220

# --- Flow field and speed

ux0,uy0,uz0 = cylinder_tang_u(0,0,0,gamma_t,R,polar_out=False)
print('uz0',uz0)

zs = np.linspace(ZLIM[0]*1.08,ZLIM[1]*1.11,nz).astype(np.float32)
xs = np.linspace(XLIM[0]*1.08,XLIM[1]*1.11,nx).astype(np.float32)
[Z,X]=np.meshgrid(zs,xs)
Y=X*0
ux,uy,uz = cylinder_tang_u(X,Y,Z,gamma_t,R,z1=z1,z2=z2,polar_out=False)

Speed=np.sqrt((uz**2+ux**2))
print('Speed range',np.min(Speed),np.max(Speed))

# --- Surface and rotor
xr_c = np.array([-R, R , R , -R ,-R])
zr_c = np.array([z1, z1, z2,  z2, z1])
yr_c = 0*xr_c

# fig=plt.figure()
# ax=fig.add_subplot(111)
# im=ax.contourf(Z,X,Speed,30)
# fig.colorbar(im)
# plt.show()
# 
# fig=plt.figure()
# ax=fig.add_subplot(111)
# im=ax.quiver(Z,X,uz,ux)
# plt.show()
# 

# --- Inputs for LIC
n = nLICKernel
texture = np.random.rand(nz,nx).astype(np.float32)
BaseKernel = np.sin(np.arange(n)*np.pi/n)

def lic_step(fig,t,it=0,save=False):
    # --- Plotting
    fig.clf()
#     fig.set_size_inches(12, 12)
    ax=fig.add_subplot(111)

    # --- LIC
    if anim_meth=='timesfreq':
        # NOTE: gives some granularity based on freq and n, not easy to scale
        kernel = BaseKernel*(1+Amplitude*np.sin(2*np.pi*freq*( np.arange(n)/float(n) + t ))) # Time freq
    elif anim_meth=='kernelphase':
        # NOTE: simple 
        kernel = np.sin(np.arange(n)*np.pi/n+ 2*np.pi*freq*t)
    elif anim_meth=='kernelphaseOffset':
        # NOTE: simple 
        kernel = np.sin(np.arange(n)*np.pi/n+ 2*np.pi*freq*t+2*np.pi*t)
    elif anim_meth=='sinpure':
        # NOTE: only changes the intensity, does not give a fow
        kernel = BaseKernel*(1+Amplitude*np.sin(2*np.pi*freq*t))

    MyImage=licImage(zs,xs,uz,ux,texture=texture,kernel=kernel,minSpeed=CLIM[0],maxSpeed=CLIM[1],accentuation=accentuation,offset=offset,spread=spread,axial=False,cmap=cmap)

    # Background
    im=ax.imshow(MyImage,extent=[min(zs),max(zs),max(xs),min(xs)])
    ax.set_xlim(ZLIM[0],ZLIM[1])
    ax.set_ylim(XLIM[0],XLIM[1])
    # Colorbar
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.20)
    manual_colorbar(fig,cmap,cax=cax, norm=mcolors.Normalize(vmin=CLIM[0], vmax=CLIM[1]))
    # Streamlines
    yseed=np.concatenate(( np.linspace(XLIM[0]*0.9,-R*1.3,3),np.linspace(-R*0.5,R*0.5,5), np.linspace(R*1.3,XLIM[1]*0.9,3)))
    start=np.array([yseed*0,yseed])
    sp=ax.streamplot(zs,xs,uz,ux,color='k',start_points=start.T,linewidth=0.7,density=30,arrowstyle='-')
#     # qv=streamQuiver(ax,sp,spacing=0.8,scale=40,angles='xy')
    qv=streamQuiver(ax,sp,n=[3,3,3,5,5,5,5,5,3,3,3],scale=40,angles='xy')
    # Surface
    ax.plot(zr_c,xr_c,'k--')
    ax.set_xlabel('z/R [-]')
    ax.set_ylabel('r/R [-]')

    if save:
        print('i={:04d} t={:.3f} '.format(it,t))
        fig.savefig("_CylinderFiniteLIC\VCF-%04d.png"%it,dpi=dpi)
    return kernel

fig=plt.figure()
# fig.set_size_inches(12, 12)
kernel_t0  =lic_step(fig,t=0,it=0,save=False)
# kernel_tmax=lic_step(fig,t=TMAX,it=0,save=False)
# print(kernel_t0)
# print(kernel_tmax)
# np.testing.assert_almost_equal(kernel_t0,kernel_tmax)
plt.show()

if movie:
    moviename='_'+anim_meth+'_nx{:03d}_nFrames{:04d}_freq{:d}_t{:03d}_A{:04d}_nKernel{:03d}_dpi{:d}.mp4'.format(int(nx),int(nFrames),int(freq*10),int(TMAX*100),int(Amplitude*1000),int(nLICKernel),int(dpi))
    # --- Movie
    for it,t in enumerate(np.linspace(0,TMAX,nFrames)[:-1]):
        lic_step(fig,t,it,save=True)

    os.system('ffmpeg.exe -r 30 -s 1920x1080 -i _CylinderFiniteLIC\VCF-%04d.png   -vcodec libx264 -crf 25 -pix_fmt yuv420p -y '+moviename)


# fig=plt.figure()
# ax=fig.add_subplot(111)
# im=ax.imshow(image,extent=[min(zs),max(zs),max(xs),min(xs)])
# 

plt.show()

# .\ffmpeg.exe -r 30 -s 1920x1080 -i ring-%04d.png   -vcodec libx264 -crf 25 -pix_fmt yuv420p -y movie.mp4

# ffmpeg -r 60 -f image2 -s 1920x1080 -i pic%04d.png -vcodec libx264 -crf 25  -pix_fmt yuv420p test.mp4
# 
#     -r is the framerate (fps)
#     -crf is the quality, lower means better quality, 15-25 is usually good
#     -s is the resolution
#     -pix_fmt yuv420p specifies the pixel format, change this as needed
# 
# lossless codec, e.g. HuffYUV or FFV1:
#     ffmpeg -i frame%04d.png -c:v huffyuv test.avi
#     ffmpeg -i frame%04d.png -c:v ffv1 -qscale:v 0 test.avi
#     ffmpeg -i frame%04d.png -qscale:v 0 test.avi
#     ffmpeg -i frame%04d.png -c:v mjpeg -qscale:v 0 test.avi
