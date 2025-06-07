"""
References:
    [1] E. Branlard - Wind Turbine Aerodynamics and Vorticity Based Method, Springer, 2017
    [2] E. Branlard, M. Gaunaa - Cylindrical vortex wake model: skewed cylinder, application to yawed or tilted rotors - Wind Energy, 2015

Coordinate systems
   c coordinate system used in see [2], rotor in plane z_c=0
   w wind coordinate system where z_w is the wind direction
   theta_yaw : yaw angle, positive around y, pointing upward

   x_c =  x_w cost + z_w sint
   y_c =  y_w
   z_c = -x_w sint + z_w cost

"""
# --- General
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable
import os
# --- Local
from welib.vortilib.elements.VortexCylinderSkewed import svc_tang_u
from welib.tools.colors import darkrainbow as cmap
from welib.tools.colors import darkrainbow, adjust_color_lightness, manual_colorbar
from welib.tools.curves import streamQuiver
from welib.tools.lic    import lic
from welib.tools.tictoc import Timer


# --- Parameters
bWindCoord=True
bAddFlow=True
R         = 1
if bAddFlow:
	CLIM=[0.3,1.1]
else:
    CLIM=[0.0,0.1]
XLIM       = [-3*R,3*R] # 
ZLIM       = [-3*R,3*R] # 
nx = 400    # Number of points for velocity evaluation
nz = nx 
CT        = 0.95
theta_yaw =-30*np.pi/180   # rad
U0        = 1

# LIC params
nLICKernel=31
offset=0.4
spread=1.1
Accentuation=1.5


# Movie params
Amplitude=1.0
anim_meth='kernelphaseOffset' # timesfreq, sinepure, kernelphase
movie=True
dpi = 300
freq=5 
TMAX=5/freq
nFrames=220

# --- Derived params
chi= theta_yaw*(1+0.3*(1-np.sqrt(1-CT))) # rad
if CT==0.4:
    gamma_t = -0.21341 # CT=0.4
elif CT==0.6:
    gamma_t = -0.40 # 
else:
    gamma_t = -0.60414 # CT=0.95
ny = nx
m  = np.tan(chi) 

def Tw2c(x_w,y_w,z_w):
    if bWindCoord:
        x_c =  x_w * np.cos(theta_yaw) + z_w * np.sin(theta_yaw)
        y_c =  y_w
        z_c = -x_w * np.sin(theta_yaw) + z_w * np.cos(theta_yaw)
    else:
        x_c,y_c,z_c = x_w,y_w,z_w
    return x_c,y_c,z_c
def Tc2w(x_c,y_c,z_c):
    if bWindCoord:
        x_w =  x_c * np.cos(theta_yaw) - z_c * np.sin(theta_yaw)
        y_w =  y_c
        z_w =  x_c * np.sin(theta_yaw) + z_c * np.cos(theta_yaw)
    else:
        x_w,y_w,z_w = x_c,y_c,z_c
    return x_w, y_w, z_w



# --- Flow field and speed
zs = np.linspace(ZLIM[0]*1.12,ZLIM[1]*1.08,nx).astype(np.float32)
xs = np.linspace(XLIM[0]*1.08,XLIM[1]*1.08,nx).astype(np.float32)
[Z,X]=np.meshgrid(zs,xs)
Y=X*0
X_c, Y_c, Z_c = Tw2c(X, Y, Z)

with Timer('VelocityField'):
    ux_c,uy_c,uz_c,_,_=svc_tang_u(X_c,Y_c,Z_c,gamma_t,R,m)
    if bAddFlow:
        uz_c=uz_c+U0*np.cos(theta_yaw) # Adding free wind
        ux_c=ux_c+U0*np.sin(theta_yaw)
ux,uy,uz = Tc2w(ux_c,uy_c,uz_c) 

Speed=np.sqrt((uz**2+ux**2))
#Speed=np.sqrt((uz**2))
print('Speed range',np.min(Speed),np.max(Speed))
Speed[Speed>CLIM[1]] = CLIM[1]
Speed[Speed<CLIM[0]] = CLIM[0]

COL=cmap((Speed-CLIM[0])/(CLIM[1]-CLIM[0])) # cmap requires value between 0 and 1
COL=COL[:,:,0:3]


# --- Surface and rotor
xr_c = np.linspace(-R,R,10)
yr_c = 0*xr_c
zr_c = 0*xr_c
xr,yr,zr=Tc2w(xr_c,yr_c,zr_c)

zu_c = np.linspace(0,ZLIM[1]*1.2,10) 
yu_c = 0*zu_c
xu_c = m*zu_c+R
zl_c = np.linspace(0,ZLIM[1]*1.2,10) 
yl_c = 0*zl_c
xl_c = m*zl_c-R
xu,yu,zu=Tc2w(xu_c,yu_c,zu_c)
xl,yl,zl=Tc2w(xl_c,yl_c,zl_c)

 
# fig=plt.figure()
# ax=fig.add_subplot(111)
# im=ax.contourf(Z,X,Speed,30)
# plt.show()
# 
# fig=plt.figure()
# ax=fig.add_subplot(111)
# im=ax.quiver(Z,X,uz,ux)
# plt.show()
# 

# --- Inputs for LIC
n = nLICKernel
texture = np.random.rand(nz,nz).astype(np.float32)
BaseKernel = np.sin(np.arange(n)*np.pi/n)

def lic_step(fig,t,it=0,save=False):
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
    ax.set_xlim(ZLIM[0],ZLIM[1])
    ax.set_ylim(XLIM[0],XLIM[1])
    # Colorbar
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.20)
    manual_colorbar(fig,cmap,cax=cax, norm=mcolors.Normalize(vmin=CLIM[0], vmax=CLIM[1]))
    # Streamlines
    if bAddFlow:
        yseed=np.linspace(XLIM[0]*0.9,XLIM[1]*0.9,17)
    else:
        yseed=np.linspace(-0.88*R,0.88*R,7)
    start=np.array([yseed*0,yseed])
    sp=ax.streamplot(zs,xs,uz,ux,color='k',start_points=start.T,linewidth=0.7,density=30,arrowstyle='-')
#     sp=ax.streamplot(zs,xs,uz,ux,color='k',linewidth=0.7,density=1)
    # qv=streamQuiver(ax,sp,spacing=0.8,scale=40,angles='xy')
    qv=streamQuiver(ax,sp,n=5,scale=40,angles='xy')
    # Surface
    ax.plot(zr,xr,'k--')
    ax.plot(zu,xu,'k--')
    ax.plot(zl,xl,'k--')
    ax.set_xlabel('z/R [-]')
    ax.set_ylabel('r/R [-]')

    if save:
        print('i={:04d} t={:.3f} '.format(it,t))
        fig.savefig("_CylinderSkewedLIC\VCS-%04d.png"%it,dpi=dpi)
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

    os.system('ffmpeg.exe -r 30 -s 1920x1080 -i _CylinderSkewedLIC\VCS-%04d.png   -vcodec libx264 -crf 25 -pix_fmt yuv420p -y '+moviename)


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
