"""
References:
    [1] E. Branlard - Wind Turbine Aerodynamics and Vorticity Based Method, Springer, 2017
"""
# --- General
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable
import os
# --- Local
from welib.vortilib.elements.VortexRing import ring_u
from welib.tools.colors import darkrainbow, adjust_color_lightness, manual_colorbar
from welib.tools.curves import streamQuiver
from welib.tools.lic    import lic

# --- Parameters
XLIM=[-2,2] # 
nx = 600    # Number of points for velocity evaluation
nz = nx 

# LIC params
nLICKernel=31
offset=0.4
spread=1.1
Accentuation=1


# Movie params
Amplitude=1.0
anim_meth='kernelphase' # timesfreq, sinepure, kernelphase
movie=True
dpi = 300
freq=5 
TMAX=5/freq
nFrames=300

# --- Flow field and speed
zs = np.linspace(-2.15,2.15,nx).astype(np.float32)
xs = np.linspace(-2.15,2.15,nx).astype(np.float32)
[Z,X]=np.meshgrid(zs,xs)
Y=X*0
ux,uy,uz = ring_u(X,Y,Z,Gamma=1,R=1,polar_out=False)

Speed=np.sqrt((uz**2+ux**2))
Speed[Speed>1]=1
cmap = darkrainbow
COL=cmap(Speed)
COL=COL[:,:,0:3]

# --- Inputs for LIC
n = nLICKernel
texture = np.random.rand(nz,nz).astype(np.float32)
BaseKernel = np.sin(np.arange(n)*np.pi/n)

# --- DEBUG KERNEL ANIM
# fig=plt.figure()
# ax=fig.add_subplot(111)
# # freq=1
# # Amplitude=0.1
# ax.plot(BaseKernel,'k+',label='BaseKernel')
# for t in np.linspace(0,TMAX,nFrames/10):
# #     kernel = np.sin(np.arange(n)*np.pi/n+ 2*np.pi*freq*t)
#     kernel = BaseKernel*(1+Amplitude*np.sin(2*np.pi*freq*( np.arange(n)/float(n) + t ) ))
# #     kernel = (1+Amplitude*np.sin(2*np.pi*freq*( np.arange(n)/float(n) + t ) ))
#     ax.plot(kernel,label='{}'.format(t))
# plt.legend()
# plt.show()
# raise 
# 
def lic_step(fig,t,it=0,save=False):
    # --- LIC
    if anim_meth=='timesfreq':
        # NOTE: gives some granularity based on freq and n, not easy to scale
        kernel = BaseKernel*(1+Amplitude*np.sin(2*np.pi*freq*( np.arange(n)/float(n) + t ))) # Time freq
    elif anim_meth=='kernelphase':
        # NOTE: simple 
        kernel = np.sin(np.arange(n)*np.pi/n+ 2*np.pi*freq*t)
    elif anim_meth=='sinpure':
        # NOTE: only changes the intensity, does not give a fow
        kernel = BaseKernel*(1+Amplitude*np.sin(2*np.pi*freq*t))
    kernel = kernel*Accentuation
    kernel = kernel.astype(np.float32)

    image=lic(uz,ux,texture=texture,kernel=kernel)
    # image=(image-MIN)/(MAX-MIN)
    image=image-np.mean(image)+nLICKernel/2 # Making sure all images have the same mean
    image=image/nLICKernel # scaling between 0 and 1
    image= offset+spread*image
#     # ---MY IMAGE
    MyImage=adjust_color_lightness(COL, image)
#     MyImage=image

    # --- Plotting
    fig.clf()
#     fig.set_size_inches(12, 12)
    ax=fig.add_subplot(111)
    # Background
    im=ax.imshow(MyImage,extent=[min(zs),max(zs),max(xs),min(xs)])
    ax.set_xlim(-2,2)
    ax.set_ylim(-2,2)
    # Colorbar
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.20)
    manual_colorbar(fig,cmap,cax=cax)
    # Streamlines
    yseed=np.linspace(-0.88,0.88,7)
    start=np.array([yseed*0,yseed])
    sp=ax.streamplot(zs,xs,uz,ux,color='k',start_points=start.T,linewidth=0.7,density=30,arrowstyle='-')
    # qv=streamQuiver(ax,sp,spacing=0.8,scale=40,angles='xy')
    qv=streamQuiver(ax,sp,n=[1,5,5,5,5,5,1],scale=40,angles='xy')
    ax.set_xlabel('z/R [-]')
    ax.set_ylabel('r/R [-]')

    if save:
        print('i={:04d} t={:.3f} '.format(it,t))
        fig.savefig("ring-%04d.png"%it,dpi=dpi)
    return kernel

fig=plt.figure()
# fig.set_size_inches(12, 12)
kernel_t0  =lic_step(fig,t=0,it=0,save=False)
kernel_tmax=lic_step(fig,t=TMAX,it=0,save=False)
print(kernel_t0)
print(kernel_tmax)
np.testing.assert_almost_equal(kernel_t0,kernel_tmax)
plt.show()

if movie:
    moviename='_'+anim_meth+'_nx{:03d}_nFrames{:04d}_freq{:d}_t{:03d}_A{:04d}_nKernel{:03d}_dpi{:d}.mp4'.format(int(nx),int(nFrames),int(freq*10),int(TMAX*100),int(Amplitude*1000),int(nLICKernel),int(dpi))
    # --- Movie
    for it,t in enumerate(np.linspace(0,TMAX,nFrames)[:-1]):
        lic_step(fig,t,it,save=True)

    os.system('ffmpeg.exe -r 30 -s 1920x1080 -i ring-%04d.png   -vcodec libx264 -crf 25 -pix_fmt yuv420p -y '+moviename)


# fig=plt.figure()
# ax=fig.add_subplot(111)
# im=ax.imshow(image,extent=[min(zs),max(zs),max(xs),min(xs)])
# 
# fig=plt.figure()
# ax=fig.add_subplot(111)
# im=ax.contourf(Z,X,Speed,30)

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
