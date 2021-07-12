""" 
Flow field in a plane normal to a vortex segment
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from vortilib.tools.colors import fColrs
from vortilib.tools.clean_exceptions import *
from vortilib.elements.VortexSegment import *
from vortilib.tools.curves import streamQuiver


# --- Parameters
minSpeed=0
maxSpeed=0.6
nX=103
nY=101
L = 2

Epsilon = 0.5*L
EpsilonX = 5*Epsilon
EpsilonY = 1*Epsilon


# --------------------------------------------------------------------------------}
# --- Regular segment 
# --------------------------------------------------------------------------------{
# Vortex segment 
Pa      = np.array([[ 0, 0, -L/2]])
Pb      = np.array([[ 0, 0,  L/2]])
Gamma   = 1
#  Control points
vX = np.linspace(-2*L,2*L,nX)
vY = np.linspace(-2*L,2*L,nY)
XCP,YCP = np.meshgrid(vX,vY)
Xcp = XCP.flatten()
Ycp = YCP.flatten()
Zcp = Xcp*0
# Velocity field
Ux, Uy, Uz = vs_u(Xcp, Ycp, Zcp, Pa, Pb, Gamma=Gamma, RegFunction=3, RegParam=Epsilon)
Ux    = Ux.reshape(XCP.shape)
Uy    = Uy.reshape(XCP.shape)
Speed = np.sqrt((Ux**2+Uy**2))/Gamma *np.pi*L
print('min: ',np.min(Speed.ravel()),' - max: ',np.max(Speed.ravel()))
Ux0=Ux
Uy0=Uy


# --------------------------------------------------------------------------------}
# --- 2D gaussian regularization
# --------------------------------------------------------------------------------{
RegFunction=3
# Velocity field
#nt=[np.sqrt(2)/2,np.sqrt(2)/2,0]
nt=[-np.sqrt(2)/2,np.sqrt(2)/2,0]
nt=[1,0,0]
Ux, Uy, Uz = vs_u(Xcp, Ycp, Zcp, Pa, Pb, Gamma=Gamma, RegFunction=RegFunction, RegParam=EpsilonX, nt=nt, RegParamW=EpsilonY)
Ux    = Ux.reshape(XCP.shape)
Uy    = Uy.reshape(XCP.shape)
Speed = np.sqrt((Ux**2+Uy**2))/Gamma *np.pi*L
print('min: ',np.min(Speed.ravel()),' - max: ',np.max(Speed.ravel()))

def circulationSurvey(r, nTheta):
    theta=np.linspace(0,2*np.pi,nTheta+1)
    dTheta=theta[1]-theta[0]
    Xcp=r*np.cos(theta)
    Ycp=r*np.sin(theta)
    Zcp=0*Xcp
    Ux, Uy, Uz = vs_u(Xcp, Ycp, Zcp, Pa, Pb, Gamma=Gamma, RegFunction=RegFunction, RegParam=EpsilonX, nt=nt, RegParamW=EpsilonY)
    #Ux, Uy, Uz = vs_u(Xcp, Ycp, Zcp, Pa, Pb, Gamma=Gamma, RegFunction=0)
    Ut = Uy * np.cos(theta) - Ux * np.sin(theta)
    U_th = Gamma*(L/2)/(2*np.pi*r * np.sqrt(r**2 + (L/2)**2))
    print(Ut[:3], U_th)
    #GammaCalc = r* np.trapz(Ut,theta)
    GammaCalc = 2*np.pi * r*Ut[0] 
    return  GammaCalc

print(circulationSurvey(4*L, 100))
print(circulationSurvey(2*L, 100))
print(circulationSurvey(L, 100))
print(circulationSurvey(L/2, 100))
print(circulationSurvey(L/5, 100))

print(circulationSurvey(1, 100))
print(circulationSurvey(2, 100))
print(circulationSurvey(4, 100))



def circulationCont(r, nTheta):
    theta=np.linspace(0,2*np.pi,nTheta+1)
    Xcp=r*np.cos(theta)
    Ycp=r*np.sin(theta)
    return Xcp,Ycp


# --- plot
fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.2,4.6)) # (6.4,4.8)
fig.subplots_adjust(left=0.12, right=0.98, top=0.96, bottom=0.12, hspace=0.20, wspace=0.20)

im = ax.contourf(XCP/L, YCP/L, Speed, levels=np.linspace(minSpeed,maxSpeed,250), vmin=minSpeed, vmax=maxSpeed)
cb=fig.colorbar(im)

# yseed=np.linspace(L/5,np.max(vY)*0.85,5)/L
yseed=np.linspace(L/5,L/2,5)/L
start=np.array([yseed*0,yseed])

sp=ax.streamplot(vX/L,vY/L,Ux,Uy,color='k',start_points=start.T,linewidth=0.7,density=30,arrowstyle='-')
qv=streamQuiver(ax,sp,n=7,scale=40,angles='xy')

sp=ax.streamplot(vX/L,vY/L,Ux0,Uy0,color='b',start_points=start.T,linewidth=0.7,density=30,arrowstyle='-')

x0,y0 = circulationCont(L, 100)
ax.plot(x0/L,y0/L,'k--')
ax.set_xlabel('x/L [-]')
ax.set_ylabel('y/L [-]')
ax.set_aspect('equal','box')

# fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.2,4.6)) # (6.4,4.8)
# fig.subplots_adjust(left=0.12, right=0.98, top=0.96, bottom=0.12, hspace=0.20, wspace=0.20)
# ix=np.argmin(np.abs(vX))
# iy=np.argmin(np.abs(vY))
# ax.plot(vY/L,Ux[:,ix]/Gamma *np.pi*L , label='x=0')
# ax.plot(vX/L,Uy[iy,:]/Gamma *np.pi*L , label='y=0')
# ax.legend()
# 


# fig.savefig('VortexFilamentRegularization.pdf')
plt.show()
