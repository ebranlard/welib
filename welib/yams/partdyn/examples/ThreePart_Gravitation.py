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

def main(test=False, tMax=800):
    m0=10000000
    p1 = Part(m0*1.0, (-2 ,2  ,0), (0.0,0,0))
    p2 = Part(m0*1.0, (-2 ,-5,0), (0,0,0))
    p3 = Part(m0*1.0, (5 ,1,0), (0,0,0))

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
    sys = PartSystem([p1,p2,p3], activeForces=['gravitational','spring'])
    sys.connect(p1, p2, 'spring', k=100)#, l0=r)
    sys.connect(p2, p3, 'spring', k=100)#, l0=r)
    sys.connect(p3, p1, 'spring', k=100)#, l0=r)
    #t   = np.linspace(0,80000,300)
    t   = np.arange(0,tMax,10)

    if precise:
        res = sys.integrate(t, method='Radau')
    else:
        res = sys.integrate(t)

    
    return sys, res



if __name__ == '__main__':
    sys, res = main(test=False, tMax=800)
    r_max=7
    # --- Energy
    sys.plotEnergy(res)
    # --- Trajectories
    fig, ax = sys.plot2Dtrajectories(res, plane='XY', XLIM=[-r_max, r_max], YLIM=[-r_max, r_max])
    fig = sys.animate2Dtrajectories(res, plane='XY', XLIM=[-r_max, r_max], YLIM=[-r_max, r_max], frame=79, fromInit=True, fig=fig)
    #fig = sys.animate2Dtrajectories(res, plane='XY', XLIM=[-r_max, r_max], YLIM=[-r_max, r_max])
    fig.suptitle('PartDyn - Gravitational and spring interactions')
    plt.show()

if __name__=="__export__":
    sys, res = main(test=False, tMax=800)
    r_max=7
    fig, ax = sys.plot2Dtrajectories(res, plane='XY', XLIM=[-r_max, r_max], YLIM=[-r_max, r_max])
    fig = sys.animate2Dtrajectories(res, plane='XY', XLIM=[-r_max, r_max], YLIM=[-r_max, r_max], frame=35, fromInit=True, fig=fig)
    fig.suptitle('PartDyn - Gravitational and spring interactions')

    from welib.tools.repo import export_figs_callback
    export_figs_callback(__file__)

