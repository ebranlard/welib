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

def main():
    r_max =80
    m0=10000000
    p1 = Part(m0*1.0, (0 ,0  ,0), (0.0,0,0))
    p2 = Part(m0*1.0, (1 ,3.1,0), (0,0,0))
    p3 = Part(m0*1.0, (5 ,-1.8,0), (0,0,0))

    # --- Case 1 no springs
    #precise=False
    #sys = PartSystem([p1,p2,p3], activeForces=['gravitational','spring'])
    #t   = np.linspace(0,80000,300)

    # --- Case 1 no springs more precise
    #precise=True
    #sys = PartSystem([p1,p2,p3], activeForces=['gravitational','spring'])
    #t   = np.linspace(0,9000,300)

    # --- Case 2 springs precise
    precise=True
    r_max=10
    sys = PartSystem([p1,p2,p3], activeForces=['gravitational','spring'])
    sys.connect(p1, p2, 'spring', k=50)#, l0=r)
    sys.connect(p2, p3, 'spring', k=50)#, l0=r)
    sys.connect(p3, p1, 'spring', k=50)#, l0=r)
    #t   = np.linspace(0,80000,300)
    t   = np.linspace(0,800,300)

    if precise:
        res = sys.integrate(t, method='Radau')
    else:
        res = sys.integrate(t)

    # --- Energy
    sys.plotEnergy(res)

    sys.plot2Dtrajectories(res, plane='XY', XLIM=[-r_max, r_max], YLIM=[-r_max, r_max])

    sys.animate2Dtrajectories(res, plane='XY', XLIM=[-r_max, r_max], YLIM=[-r_max, r_max])



if __name__ == '__main__':
    main()
    plt.show()
