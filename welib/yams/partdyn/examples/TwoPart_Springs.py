""" 
Two particles connected by a spring

Second particle has a large mass. 
Gravity can be included or set to 0
Mass 2 can be fixed so that it reduces to a simple oscillator problem

"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# Local 
from welib.yams.partdyn.part import*


def spring(fixP2=True):
    # --- Parameters
    m1 = 10
    m2 = 100
    g  = 9.81   # Gravity
#     g  = 0
    z1 = -2
    z2 = 2
    r_max = 6
    l0 = 2.5
#     l0 = 4

    f = 1
    zeta = 0
    k = (2*np.pi)**2 *m1
    c = 0 # TODO, theory not implemented
    T = 1/f
    t =np.arange(0,5*T,T/100)

    # --- Numerical solution
    p1 = Part(m1, r0=(0,0,z1), v0=(0, 0 ,0))
    p2 = Part(m2, r0=(0,0,z2), v0=(0, 0, 0))
    sys = PartSystem([p1, p2], activeForces=['spring','gravity'])
    sys.connect(p1, p2, 'spring', k=k, l0=l0, c=c) # Connect particles with a spring/damper
    if fixP2:
        sys.fix(p2) # fix second particle so that it doesn't move
    sys.g=g # set gravity
    res = sys.integrate(t, method='Radau') # method 'Radau' works best for this problem

    # --- Analytical solution (NOTE: assumes m2 fixed and no initial velocity)
    zg  = m1*g/k      # equilibrium offset due to gravity
    zl0 = l0 -(z2-z1) # equilibirum offset due to l0
    zeq = z1-zg-zl0   # full equilibrium position
    A   = z1-zeq      # Amplitude is initial position - equilibrium position z~ = A cos(omega t)
    zth = zeq + A * np.cos(2*np.pi*f*res.t)

    # --- Energy
    sys.plotEnergy(res)

    # --- Plot trajectory
    i1,i2=0,2
    fig, ax = sys.plotTrajectories(res, comp=2)
    ax.plot(res.t, zth ,'k:' , label = 'Theory')
    ax.set_title('Springs')
    ax.legend(loc='upper right')

    # --- Plot animation
    if __name__ == '__main__':
        sys.animate2Dtrajectories(res, XLIM=[-r_max, r_max], YLIM=[-r_max, r_max], plane='XZ')

    return res, zth

res, zth = spring()


if __name__ == '__main__':
    plt.show()
if __name__ == '__test__':
    np.testing.assert_almost_equal(res.y[2,:], zth, 3)
