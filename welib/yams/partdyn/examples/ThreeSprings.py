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
    m1      = 10
    m_      = 100       
    r       = 3
    phi     = 1*np.pi/6
    tension = 1
    r_max = 6

    f = 1
    k = (2*np.pi)**2 *m1
    T = 1/f
    t =np.arange(0, 5*T, T/100)

    # --- Setup system at equilibrium
    p1 = Part(m1, r0=(0,0,0),  v0=(0, 0 ,0))
    #p2 = Point(r=(r*np.cos(2*np.pi/3*0+phi),0,r*np.sin(2*np.pi/3*0+phi)))
    #p3 = Point(r=(r*np.cos(2*np.pi/3*1+phi),0,r*np.sin(2*np.pi/3*1+phi)))
    #p4 = Point(r=(r*np.cos(2*np.pi/3*2+phi),0,r*np.sin(2*np.pi/3*2+phi)))
    p2 = Part(m_, r0=(r*np.cos(2*np.pi/3*0+phi),0,r*np.sin(2*np.pi/3*0+phi)))
    p3 = Part(m_, r0=(r*np.cos(2*np.pi/3*1+phi),0,r*np.sin(2*np.pi/3*1+phi)))
    p4 = Part(m_, r0=(r*np.cos(2*np.pi/3*2+phi),0,r*np.sin(2*np.pi/3*2+phi)))
    sys = PartSystem([p1, p2, p3, p4], activeForces=['spring'])
    sys.connect(p1, p2, 'spring', k=k, l0=r*tension)#, l0=r)
    sys.connect(p1, p3, 'spring', k=k, l0=r*tension)#, l0=r)
    sys.connect(p1, p4, 'spring', k=k, l0=r*tension)#, l0=r)

    # --- Perturb system
    print(sys.barycenter())
    p1.r0 = np.array((0,0,1))
#     p1.v0 = np.array((0.0,0,10))
    print(sys.barycenter())

    res = sys.integrate(t, method='Radau') # method 'Radau' works best for this problem
    R = sys.barycenter(res)

    # --- 
    i1,i2=0,2
    fig,ax = plt.subplots(1,1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
    fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
    ax.plot(res.y[i1  ,:],res.y[i2  ,:] ,'-' , label = 'Position')
    ax.plot(res.y[i1+3,:],res.y[i2+3,:] ,'o' , label = 'Position')
    ax.plot(res.y[i1+6,:],res.y[i2+6,:] ,'o' , label = 'Position')
    ax.plot(res.y[i1+9,:],res.y[i2+9,:] ,'o' , label = 'Position')
    #ax.plot(R[i1,:],R[i2  ,:] ,'k:' , label = 'Barycenter')
    ax.set_xlabel('x [m]')
    ax.set_ylabel('y [m]')
    ax.set_title('')



    # --- Plot animation
    if __name__ == '__main__':
        sys.animate2Dtrajectories(res, XLIM=[-r_max, r_max], YLIM=[-r_max, r_max], plane='XZ')

    return res

res = spring()


if __name__ == '__main__':
    plt.show()
