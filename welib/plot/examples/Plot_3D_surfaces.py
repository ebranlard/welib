""" 
Examples to plot 3D surfaces
"""
import numpy as np
import matplotlib.pyplot as plt

from welib.plot.surface3d import *



# --- Torus
X,Y,Z= torus(r=0.25, R=1, nTheta=32, nPhi=31, thetaMax=2*np.pi, phiMax=2*np.pi)
fig = plt.figure()
ax = fig.add_subplot(111, projection = '3d')
ax.plot_surface(X, Y, Z, color = 'r', rstride = 1, cstride = 1)
axisEqual3D(ax)


# --- Cylinder
X,Y,Z = cylinder(R1=2, R2=1, z1=0, z2=5, nz=3, nTheta=15, thetaMax=2*np.pi)
fig = plt.figure()
ax = fig.add_subplot(111, projection = '3d')
ax.plot_surface(X, Y, Z, color = 'g', rstride = 1, cstride = 1)
axisEqual3D(ax)


# --- Arbitrary cylinder
theta=np.linspace(0,2*np.pi,30)
x=np.cos(theta)
y=np.sin(2*theta)
X,Y,Z = arbitrary_cylinder(x, y, z1=0, z2=1, nz=3, nTheta=25)
fig = plt.figure()
ax = fig.add_subplot(111, projection = '3d')
ax.plot_surface(X, Y, Z, color = 'b', rstride = 1, cstride = 1)
axisEqual3D(ax)

if __name__=="__main__":
    plt.show()
if __name__=="__test__":
    pass
if __name__=="__export__":
    pass
    #from welib.tools.repo import export_figs_callback
    #export_figs_callback(__file__)
