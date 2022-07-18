""" 
Free fall and projectile trajectory

Simulations using either: 
  - one particle under influence of gravity
  - one particle under influence of a force field (corresponding to gravity field)
  - one particle and "Earth" (as a second particle) under influence of gravitational force

"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# Local 
from welib.yams.partdyn.part import*


def freefall(nPart=1, forceField=False):
    # Parameters
    g    = 9.81
    m    = 10
    tmax = 4.52
    z0   = 1/2*g*tmax**2*0

    t  = np.linspace(0,tmax,100)

    # Analytical solution
    zth = -1/2*g*t**2 + z0

    # numerical solution
    p1 = Part(m, r0=(0,0,z0), v0=(0.0,0,0))
    if nPart==1:
        # --- System of one particle
        if not forceField:
            sys = PartSystem([p1], activeForces=['gravity'])
            sys.g = g
        else:
            sys = PartSystem([p1], activeForces=['field'])
            sys.setFieldForce(lambda t,p: np.array([0,0,-g*p.m]))
        res = sys.integrate(t)
    else:
        # --- Using "Earth" as a second particle
        # Hack to get "g"
        mEarth = 5.97219e24
        rEarth = 6.371e6
        G      = g/(mEarth/rEarth**2)
        p2 = Part(mEarth, r0=(0,0,-rEarth), v0=(0.0,0,0))
        sys = PartSystem([p1,p2], activeForces=['gravitational'])
        sys.G=G

    res = sys.integrate(t)
    ##print('Final velocity: ', res.y[5,-1])

    ## --- Plot trajectory
    fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
    fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
    ax.plot(res.t,res.y[2,:] ,'-' , label = 'numerical' )
    ax.plot(res.t,zth        ,'--', label = 'analytical')
    ax.set_xlabel('time [s]')
    ax.set_ylabel('z [m]')
    ax.legend()
    ax.set_title('Free fall')

    #sys.plotEnergy(res)

    return res, zth


def projectile():
    # Parameters
    g    = 9.81 # [m/s^2]
    m    = 10   # [kg]
    v0    = 10  # [m/s]
    angle = 25  # [deg]
    z0    = 10  # [m]

    # Initial velocity
    vx0 = v0*np.cos(angle*np.pi/180)
    vz0 = v0*np.sin(angle*np.pi/180)
    # Time until impact
    tmax = vz0/g+np.sqrt((vz0/g)**2+2*z0/g)
    t  = np.linspace(0,tmax,100)
    # Analytical solution
    xth =  0 + vx0*t
    zth = z0 + vz0*t + -1/2*g*t**2 

    # numerical solution
    p1 = Part(m, r0=(0,0,z0), v0=(vx0,0,vz0))
    sys = PartSystem([p1], activeForces=['gravity'])
    sys.g = g
    res = sys.integrate(t)

    ## --- Plot trajectory
    fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
    fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
    ax.plot(res.y[0,:],res.y[2,:] ,'-' , label = 'numerical' )
    ax.plot(xth,zth        ,'--', label = 'analytical')
    ax.set_xlabel('x [m]')
    ax.set_ylabel('z [m]')
    ax.legend()
    ax.set_title('Projectile')

    #sys.plotEnergy(res)

    return res, (xth,zth)

# --- Freefall 1, or 2 particle, with a force field
ffnum, ffzth   = freefall()
ffnum2, ffzth2 = freefall(2)
ffnum3, ffzth3 = freefall(forceField = True)
# --- Projectile
pjnum, pjth = projectile()

if __name__ == '__main__':
    plt.show()

if __name__ == '__test__':
    # Freefall
    np.testing.assert_almost_equal(ffnum.y[2,:] , ffzth)
    np.testing.assert_almost_equal(ffnum2.y[2,:], ffzth2,3)
    np.testing.assert_almost_equal(ffnum3.y[2,:], ffzth3)
    # Projectile
    np.testing.assert_almost_equal(pjnum.y[0,:], pjth[0])
    np.testing.assert_almost_equal(pjnum.y[2,:], pjth[1])

