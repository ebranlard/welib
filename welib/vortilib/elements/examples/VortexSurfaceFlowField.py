import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from InducedVelocities import *

# --------------------------------------------------------------------------------}
# --- Main function
# --------------------------------------------------------------------------------{
def surface_u(x_surf, R_surf, gamma, Xcp, Rcp, nRings =50, includeCylinder=True):
    """ 
    Computes the induced velocity field by a vorticity surface of constant intensity

    INPUTS:
       x_surf: axial points defining the surface (array)
               NOTE: last point will be continued by a semi-inf cylinder
       R_surf: surface radius, for each axial point (array)
       gamma: intensity of the vorticity surface (typically <0) [m/s]
       Xcp : scalar/array/matrix of Control Points axial coordinates
       Rcp : scalar/array/matrix of Control Points radial coordinates (same dimension as Xcp)

       nRings: number of rings used to discretize the surface up until the last value of x_surf
       includeCylinder: if true, a semi-infinite vortex cylinder is used to continue the surface to infinity

    """
    Xcp=np.asarray(Xcp)
    Rcp=np.asarray(Rcp)
    # --- Vortex rings (Quadrature points/discretization of the surface)
    x_rings = np.linspace(np.min(x_surf), np.max(x_surf), nRings) # Equi-spacing the rings 
    R_rings = np.interp(x_rings, x_surf, R_surf) # Interpolating the ring radii based on the input surface distribution 

    dx_rings = x_rings[-1]-x_rings[-2] # Distance between the two last rings
    Gamma_ring = gamma*dx_rings        # Circulation of each vortex ring [m^2/s]
    epsilon_ring  =  dx_rings/2         # Basic epsilon, scaling can be done later (see TODO)

    # --- Vortex cylinder
    x_cyl = x_rings[-1] + dx_rings/2  # Cylinder starts at dx/2 after the last ring
    R_cyl = R_rings[-1]             # Cylinder has same radius as last ring

    # --- Cylinder induced velocity 
    if includeCylinder:
#         ur_cyl, ux_cyl = vc_tang_u_doublet(Rcp, Rcp*0, Xcp-x_cyl, gamma_t=gamma, R=R_cyl, r_bar_Cut=6,polar_out=True)
        ur_cyl, ux_cyl = vc_tang_u(Rcp, Rcp*0, Xcp-x_cyl, gamma_t=gamma, R=R_cyl,polar_out=True)

    # --- Induced velocity from all rings
    ur_rings, ux_rings = np.zeros(Xcp.shape), np.zeros(Xcp.shape)
    for i_ring, (x_ring, R_ring) in enumerate(zip(x_rings, R_rings)):
        if i_ring==0:
            Gamma_ring_scaled=Gamma_ring/2 # first ring represent half a ring..
        else:
            Gamma_ring_scaled   = Gamma_ring   # TODO insert scaling here
        epsilon_ring_scaled = epsilon_ring # TODO insert epsilon hack here
        ur, ux =  ring_u_polar(Rcp, Xcp, Gamma=Gamma_ring_scaled, r0=R_ring, z0=x_ring, epsilon=epsilon_ring_scaled, reg_method='Saffman')
        ur_rings += ur
        ux_rings += ux

    geom = (x_rings, R_rings, x_cyl, R_cyl)

    # --- Total velocity
    if includeCylinder:
        ur = ur_cyl + ur_rings
        ux = ux_cyl + ux_rings
        return ur, ux, geom
    else:
        return ur_rings, ux_rings, geom







# --------------------------------------------------------------------------------}
# --- Main code 
# --------------------------------------------------------------------------------{
# --- Parameters
R               = 1       # Rotor radius
nRings          = 50     # Number of rings used over the vorticity surface region before the semi-inf cylinder
includeCylinder = True    # Include semi-infinite cylinder
U0              = 10      # Freestream
a               = 0.4     # Axial induction
gamma           = -2*a*U0 # intensity of the vorticity surface [m/s]
n_rcp           = 50      # Number of control points in radial direction
n_xcp           = 300     # Number of control points in axial direction
x_max           = 10*R    # Max

# --- Control points used for velocity field
rcp = np.linspace(0,5*R, n_rcp)
xcp = np.linspace(-3*R,x_max*1.5, n_xcp)

# --- Definition of the vorticity surface (NOTE: last point continued by a semi-inf cylinder)

# Example 1: cylinder
# x_surf = np.array([0,x_max ])*R
# R_surf = np.array([1,1] )*R 
# rcp = rcp[ np.abs(rcp-R)>0.01*R ] # hack to remove points close cylinder

# Example 2: expanding cylinder
x_surf = np.linspace(0,x_max, 100)
R_surf = R* np.sqrt( (1-a) / (1-a*(1+x_surf/R/np.sqrt(1+(x_surf/R)**2 )  ))    )

# TODO  Insert your own surface here


# --- Velocity at special points
ur0, ux0, _ = surface_u(x_surf, R_surf, gamma, [0               ], [0], nRings = nRings, includeCylinder=includeCylinder)
urm, uxm, _ = surface_u(x_surf, R_surf, gamma, [np.max(x_surf)/2], [0], nRings = nRings, includeCylinder=includeCylinder)
urw, uxw, _ = surface_u(x_surf, R_surf, gamma, [np.max(x_surf)  ], [0], nRings = nRings, includeCylinder=includeCylinder)
print('Center velocity     :', ur0/U0, ux0/U0)
print('Mid surf velocity   :', urm/U0, uxm/U0)
print('Wake start velocity :', urw/U0, uxw/U0)

# --- Velocity on a grid of points for plotting
Rcp, Xcp = np.meshgrid(rcp,xcp)
ur, ux, geom = surface_u(x_surf, R_surf, gamma, Xcp, Rcp, nRings = nRings, includeCylinder=includeCylinder)
x_rings, R_rings, x_cyl, R_cyl = geom
ux += U0 # Adding free stream velocity


# --- Plotting
# speed = np.sqrt(ur**2 + ux**2)/U0
speed = ux/U0
speedMin=0
speedMax=1.1
levels=np.linspace(speedMin,speedMax,10)

fig,ax = plt.subplots(1,1, figsize=(9,4))
im=ax.contourf(Xcp/R,Rcp/R,speed, levels=levels) # the easiest way to get contourf to comply with cmap is to give levels
fig.colorbar(im)
ax.plot(x_surf/R,R_surf/R,'k-', lw=3, label='Input surface')
ax.plot(x_rings/R, R_rings, 'o'  , label='Rings', ms=2)
ax.plot(np.array([x_cyl,np.max(xcp)])/R, np.array([R_cyl, R_cyl])/R, 'k--'  , label='Vortex cylinder', lw=3)
ax.set_aspect('equal')
ax.legend()


plt.show()
