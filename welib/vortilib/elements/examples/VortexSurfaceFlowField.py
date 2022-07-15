import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from welib.vortilib.elements.VortexAxisymmetric import axisym_surface_u as surface_u

# --------------------------------------------------------------------------------}
# --- Main code 
# --------------------------------------------------------------------------------{
def main():
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
    im=ax.contourf(Xcp/R,Rcp/R, speed, levels=levels) # the easiest way to get contourf to comply with cmap is to give levels
    fig.colorbar(im)
    ax.plot(x_surf/R,R_surf/R,'k-', lw=3, label='Input surface')
    ax.plot(x_rings/R, R_rings, 'o'  , label='Vortex Rings', ms=2)
    ax.plot(np.array([x_cyl,np.max(xcp)])/R, np.array([R_cyl, R_cyl])/R, 'k--'  , label='Vortex cylinder', lw=3)
    ax.set_aspect('equal')
    ax.tick_params(direction='in')
    ax.legend()
    ax.set_title('Vortilib - Flow about an axisymmetric vorticity surface')
    ax.set_xlabel(r'$r/R$ [-]')
    ax.set_ylabel(r'$z/R$ [-]')


if __name__ == '__main__':
    main()
    plt.show()
if __name__=="__test__":
    main()
    pass
if __name__=="__export__":
    main()
    from welib.tools.repo import export_figs_callback
    export_figs_callback(__file__)
