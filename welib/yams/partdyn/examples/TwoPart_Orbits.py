""" 
Orbit of one particle about another one more massive



Orbital dynamics:

    The two-body problem, with a conservative central force is equivalent to a one body problem.
    with reduced mass mu = m1 m2/(m1+m2), with the force center on the barycenter location

    Quantities conserved:
        L = mu r^2 theta_dot   angular momentum
        E = 1/2 mu rdot^2 +  L^2/(2 mu r^2) - k/r

    k = G m1 m2

    Motion is an ellipse with semi-major axis a, accentricity e, where r is from one of the foci
        E  = -k/(2*a) 
        e = sqrt(1+ 2 E L^2/(mu k^2)) 
        r_min = L^2/(mu k)/(1+e) # min distance from focal point
        r_max = L^2/(mu k)/(1-e) # max distance from focal point

    These equations can be manipulated depending on which parameters are known





"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# Local 
from welib.yams.partdyn.part import*


def orbit(problem = 'moon'):
    # Parameters
    mEarth = 5.97219e24    # 
    rEarth = 6.371e6       # earth radius [m]
    mMoon   = 7.34767309e22 # Moon mass [kg]
    rEMoon  = 384.400e6     # Distance between earth and moon centers [m]
    G      = GRAVITATIONAL_CONSTANT   # Gravitational constant [Nm^2/kg^2[



    if problem=='Moon' :
        # --- Moon about the earth
        # Here we assume a circular orbit.
        # We set the momentum of the earth and moon to be opposite at t=0 so that the 
        # global trajectory is a circle (the earth moves).
        m1 = mMoon
        r1 = rEMoon
        m2 = mEarth
        v1 = np.sqrt(G*m2/r1) # Orbit velocity to obtain a circle
        v2 = -(v1*m1)/m2 # Equating momentum between the two systems to ensure "global" trajectories remain circle
        T = 24*3600*28 # period [s]
        # Ellipse parameters
        a, b, r_min, r_max = r1,r1,r1,r1 # We expect a circular trajectory
        e =0 
        mu = m1*m2/(m1+m2)
        k  = G*m1*m2
        #L = np.sqrt(mu*r1*k)
        #v1 =  L/(mu*r1) # linear velocity at r_min, v1 = r theta_dot
        #v2=0
    elif problem=='Satellite':
        # Here we use the reduce mass for the equations
        # Since the mass of the earth is significantly more than the satellite,
        # it's fair to neglect the motion of the Earth
        m2    = mEarth
        m1    = 295               # Satellite mass [kg]
        r_min = rEarth + 550      # Minimum distance to object 2 [m]
        r_max = rEarth + 20000000 # Maximum distance to object 2 [m]
        # Ellipse parameters
        mu = m1*m2/(m1+m2)
        k  = G*m1*m2
        a  = (r_min+r_max)/2
        E  = -k/(2*a)
        Lrmin  = np.sqrt(2*mu*r_min**2 * (E+k/r_min))  # L = mu r^2 theta_dot, here expressed at r_min using energy equation
        Lrmax  = np.sqrt(2*mu*r_max**2 * (E+k/r_max))  # L = mu r^2 theta_dot, here expressed at r_max using energy equation
        L = (Lrmin+Lrmax)/2
        T = 2*np.pi*np.sqrt(a**3*mu/k)    # Period [s]
        e = np.sqrt(1+2*E*L**2/(mu*k**2)) # eccentricity [-] 
        b = a * np.sqrt(1-e**2)

        # Initial conditions
        r1 = r_min
        v1 =  L/(mu*r_min) # linear velocity at r_min, v1 = r theta_dot
        v2 = 0
        print('Eccentricity',e)

        # --- sanity checks
        #Ermin =  L**2/(2*mu*r_min**2) -k/r_min
        #Ermax =  L**2/(2*mu*r_max**2) -k/r_max
        #r_min2 = L**2/(mu*k)/(1+e)
        #r_max2 = L**2/(mu*k)/(1-e)
    else:
        raise NotImplementedError()

    t = np.arange(0,T,T/100)

    # ---  Numerical solution
    p1 = Part(m1, r0=(r1, 0, 0), v0=(0.0 , v1 ,0))
    p2 = Part(m2, r0=(0 , 0, 0), v0=(0.0 , v2 ,0))
    sys = PartSystem([p1, p2], activeForces=['gravitational'])
    res = sys.integrate(t, method='Radau') # method 'Radau' works best for this problem

    # Barycenter
    R = sys.barycenter(res)


    # --- Energy
    sys.plotEnergy(res)

    ## --- Plot trajectory
    fig,ax = plt.subplots(1,1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
    fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
    ax.plot(res.y[0,:],res.y[1,:] ,'-' , label = '{} Abs. Position'.format(problem) )
    ax.plot(res.y[3,:],res.y[4,:] ,'-' , label = 'Earth Abs. Position' )
    #ax.plot(res.y[0,:]-R[0,:], res.y[1,:]-R[1,:],'-' , label = '{} Rel. Position'.format(problem) )
    #ax.plot(res.y[3,:]-R[0,:], res.y[4,:]-R[1,:],'-' , label = '{} Rel. Position'.format(problem) )
    #ax.plot(R[0,:]          ,  R[1,:],'-' , label = 'Barycenter Position (num.)')
    #ax.plot(res.y[0,:]-res.y[3,:],res.y[1,:]-res.y[4,:],'-' , label = '{} Position wrt. Earth Position'.format(problem) )
    ax.set_xlabel('x [m]')
    ax.set_ylabel('y [m]')
    ax.set_title('{} Orbit'.format(problem))

    # --- Plot main points
    F = np.sqrt(a**2-b**2) # Location of Foci from origin = F=a*e 
    x2,y2 = 0  ,0
    x0,y0 = x2-F,0 # "origin"
    ax.plot(x2          ,y2,'ko') # Focal point 1 ("earth")
    ax.plot(x2+r_min    ,y2,'kd') # periapse
    ax.plot(x2-r_max    ,y2,'kd') # apoapse
    ax.plot(x0          ,y0,'k.') # Origin of ellipse
    ax.plot(x0-F ,y0,'k.')        # Focal point 2
    # Ellipse with coordinate about point 2
    vtheta=np.linspace(0,2*np.pi,100)
    r = a*(1-e**2)/(1+e*np.cos(vtheta))
    xe = r*np.cos(vtheta)+x2
    ye = r*np.sin(vtheta)+y2
    # Ellipse with coordinate about "origin"
    xe2 = a*np.cos(vtheta)+x0
    ye2 = b*np.sin(vtheta)+y0
    ax.plot(xe2,ye2,'k--')
    ax.legend(loc='upper right')
    ax.set_aspect('equal')

    # --- Plot animation
    #sys.animate2Dtrajectories(res, XLIM=[-r_max, r_max], YLIM=[-r_max, r_max])

    return res

res = orbit('Moon')
res = orbit('Satellite')


if __name__ == '__main__':
    plt.show()
