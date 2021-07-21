"""  
 Simulation of a 3D pendulum in python 

 The equation of motions are formulated in different ways:
   - general : suitable for general multbody formulation with constraints
   - minimal : derived using Kane's method, minimal set of coordinates

 The rotational coordinates are implemented using either:
   - Euler angles
   - Bryant angles
   - Euler parameters

"""
import numpy as np
from numpy import cos, sin
import pandas as pd
import matplotlib.pyplot as plt
# Local 
from welib.system.mech_system import MechSystem
from welib.yams.rotations import *

from scipy.integrate import  solve_ivp #odeint


def main(formulation='general', angles='Bryant', tmax=5):
    """ 
    formulation: formulation for equations of motions (general or minimal)
    angles:     choice for angle coordinates (Bryant, Euler, EulerP)
    """
    # --- Parameters
    grav    = 9.82       # acceleration of gravity
    l       = 1.0        # pendulum length [m]
    theta_x0 = 90/180*np.pi # angle wrt x axis at t=0
    omega_z0 = 2*np.pi      # rotational velocity at t=0
    # Inertia
    rho     = 7850       # pendulum density [kg/m^3]
    diam    = 0.02       # pendulum diameter [m]
    s_loc = [0,0,-l/2] # position of COG in body coordinates
    m            = 7850*l*diam**2*np.pi/4  # pendulum mass [kg]
    JxxO         = m*l**2/3                                   # at Origin
    Jxx          = m*l**2/12                                  # at COG
    Jzz          = m*(diam/2)**2/2                            #
    J_loc      = np.array([[Jxx,0,0],[0,Jxx,0], [0,0,Jzz]]) # at COG
    s_loc_skew = np.array([[0,l/2,0],[-l/2,0,0],[0,0,0]])
    J_O  = JxxO
    J_zz = Jzz
    g=grav
    r=l

    # --- Definition of RHS
    x0           = [0,-l/2*np.sin(theta_x0),l/2*np.cos(theta_x0)] # Position of COG at t = 0, in global
    xdot0        = [0,0,0 ]                                       # Linear velocity of COG at t = 0, in global
    omega0_loc = [0,0,omega_z0]                                 # rotational velocity of body at t = 0, in body coordinates
    if angles =='Bryant':
        theta0 = [theta_x0,0,0]
        Itheta     = [3,4,5]
        Ix_dot     = [6,7,8]
        Iomega     = [9,10,11]
        fA = BodyXYZ_A
        fG = BodyXYZ_Gbinv
        stheta =  [r'$\phi_x$', r'$\phi_y$', r'$\phi_z$']
    elif angles =='Euler':
        theta0 = [0,theta_x0,0]
        Itheta     = [3,4,5]
        Ix_dot     = [6,7,8]
        Iomega     = [9,10,11]
        fA = BodyZXZ_A
        fG = BodyZXZ_Gbinv
        stheta =  [r'$\phi$', r'$\theta$', r'$\psi$']
    elif angles =='EulerP':
        theta0 = [np.cos(theta_x0/2), np.sin(theta_x0/2), 0, 0 ];
        Itheta     = [3,4,5,6]
        Ix_dot     = [7,8,9]
        Iomega     = [10,11,12]
        fA = EulerP_A
        fG = EulerP_Gbinv
        stheta =  [r'$e_0$', r'$e_1$', r'$e_2$', r'$e_3$']


#     raise Exception()


    if formulation =='general':
        z0 = np.concatenate((x0,theta0,xdot0,omega0_loc))

        def dzdt(t,z):
            """ 
            z=[x,y,z, theta1..thetan, xdot, ydot, zdot, omegax', omegay', omegaz']
            """
            z  = z.ravel()
            dz = z*0
            x         = z[:3]
            theta     = z[Itheta]
            xdot      = z[Ix_dot]
            omega_loc = z[Iomega]

            # --- Kinematics:
            A     = fA(*theta)
            Gbinv = fG(*theta)

            # --- Newton Euler Equation with constraints
            # Body mass matrix
            M = np.block([[m*np.eye(3)     , np.zeros((3,3))],
                          [np.zeros((3,3)) , J_loc       ]])
            # Jacobian of constraints
            Jac = np.block([np.eye(3) , -A.dot(s_loc_skew)])
            # Full system mass matrix
            MM  = np.block([[M  , Jac.T],
                            [Jac, np.zeros((3,3)) ]])
            # Right hand side
            force  = np.array([0,0, m*grav])
            moment = -np.cross(omega_loc, J_loc.dot(omega_loc))
            gamma  = -A.dot( np.cross(omega_loc, np.cross(omega_loc, s_loc)))
            RHS    = np.concatenate((force, moment, gamma))
            # Solution
            sol = np.linalg.solve(MM, RHS)

            # --- Assign dzdt
            dz[0:3]    = xdot
            dz[Itheta] = Gbinv.dot(omega_loc) # theta_dot
            dz[Ix_dot] = sol[:3]
            dz[Iomega] = sol[3:6]
            return dz

    elif formulation=='minimal':

        z0 = np.concatenate((theta0,omega0_loc)) # TODO theta0dot

        if angles=='Bryant':
            def dzdt(t,z):
                """ 
                z=[theta1..thetan, thetadot1..thetadotn]
                NOTE: J_O at origin!
                """
                z  = z.ravel()
                dz = z*0
                q  = z[0:3] # theta
                qd = z[3:6] # theta_dot

                M = np.zeros((3,3))
                M[0,:] = [J_O*cos (q[1])*cos(q[2]), J_O*sin(q[2]),0]
                M[1,:] = [-J_O*sin(q[2])*cos(q[1]), J_O*cos(q[2]),0]
                M[2,:] = [J_zz*sin(q[1])          ,0  ,        J_zz]

                F = np.zeros(3)
                F[0] = -J_O*qd[0]**2*sin(q[1])*sin(q[2])*cos(q[1])+2*J_O*qd[0]*qd[1]*sin(q[1])*cos(q[2])+J_zz*qd[0]**2*sin(q[1])*sin(q[2])*cos(q[1])-J_zz*qd[0]*qd[1]*sin(q[1])*cos(q[2])+J_zz*qd[0]*qd[2]*sin(q[2])*cos(q[1])-J_zz*qd[1]*qd[2]*cos(q[2])-g*m*r*sin(q[0])*cos(q[2])/2-g*m*r*sin(q[1])*sin(q[2])*cos(q[0])/2
                F[1] = -J_O*qd[0]**2*sin(q[1])*cos(q[1])*cos(q[2])-2*J_O*qd[0]*qd[1]*sin(q[1])*sin(q[2])+J_zz*qd[0]**2*sin(q[1])*cos(q[1])*cos(q[2])+J_zz*qd[0]*qd[1]*sin(q[1])*sin(q[2])+J_zz*qd[0]*qd[2]*cos(q[1])*cos(q[2])+J_zz*qd[1]*qd[2]*sin(q[2])+g*m*r*sin(q[0])*sin(q[2])/2-g*m*r*sin(q[1])*cos(q[0])*cos(q[2])/2
                F[2] = -J_zz*qd[0]*qd[1]*cos(q[1])

                qdd = np.linalg.solve(M, F)
                dz[0:3]=qd
                dz[3:6]=qdd

                return dz
        else:
            raise NotImplementedError()
    else:
        raise NotImplementedError()


    # time integration
    time=np.arange(0,tmax,0.005)

    res = solve_ivp(fun=dzdt, t_span=[time[0], time[-1]], y0=z0, t_eval=time, method='RK45')
    if formulation =='general':
        CG    = res.y[0:3,:]
        theta = res.y[Itheta,:]
        omega = res.y[Iomega,:]

    elif formulation =='minimal':
        theta = res.y[0:3,:]
        q  = res.y[0:3:]
        qd = res.y[3:6,:]
        if angles=='Bryant':
            CG = np.zeros((3,len(time)))
            CG[0,:] =  r*sin(q[1,:])/2
            CG[1,:] = -r*sin(q[0,:])*cos(q[1,:])/2
            CG[2,:] =  r*cos(q[0,:])*cos(q[1,:])/2
            omega = np.zeros((3,len(time)))
            omega[0,:] = qd[0]*cos(q[1])*cos(q[2])+qd[1]*sin(q[2])
            omega[1,:] = -qd[0]*sin(q[2])*cos(q[1])+qd[1]*cos(q[2])
            omega[2,:] = qd[0]*sin(q[1])+qd[2]

    # --- Plot position
    fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
    fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
    ax.plot( time, CG[0,:], label='x')
    ax.plot( time, CG[1,:], label='y')
    ax.plot( time, CG[2,:], label='z')
    ax.set_xlabel('Time [s]')
    ax.set_ylabel('COG position [m]')
    ax.legend()
    ax.tick_params(direction='in')

    ## --- Plot angles
    fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
    fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
    [ax.plot(time, np.mod(theta[i,:],2*np.pi), label=s) for i,s in enumerate(stheta)]
    ax.set_xlabel('Time [s]')
    ax.set_ylabel('Orientation coordinates [rad]')
    ax.legend()
    ax.tick_params(direction='in')

    # --- Plot omega
    fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
    fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
    ax.plot(time, omega[0,:], label=r'$\omega_x$')
    ax.plot(time, omega[1,:], label=r'$\omega_y$')
    ax.plot(time, omega[2,:], label=r'$\omega_z$')
    ax.set_xlabel('Time [s]')
    ax.set_ylabel('Rotational speed [rad/s]')
    ax.legend()
    ax.tick_params(direction='in')

    return time, CG, theta, omega




if __name__=="__main__":
    #time, CG1, theta1, omega1 = main(formulation ='general', angles='Euler')
    #time, CG2, theta2, omega2 = main(formulation ='general', angles='Bryant')
    #time, CG3, theta3, omega3 = main(formulation ='general', angles='EulerP')
    time, CG4, theta4, omega4 = main(formulation ='minimal', angles='Bryant')

    #plt.close(1)
    plt.show()
    #fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
    #fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
    #ax.plot(time, CG1[1] , label='Euler')
    #ax.plot(time, CG2[1] , label='Bryant')
    #ax.plot(time, CG3[1] , label='EulerP')
    #ax.plot(time, CG4[1] , label='Bryant')
    #ax.legend()
    #ax.tick_params(direction='in')
    #plt.show()

if __name__=="__test__":
    time, CG1, theta1, omega1 = main(formulation ='general', angles='Euler' ,tmax=0.1)
    time, CG2, theta2, omega2 = main(formulation ='general', angles='Bryant',tmax=0.1)
    time, CG3, theta3, omega3 = main(formulation ='general', angles='EulerP',tmax=0.1)
    time, CG4, theta4, omega4 = main(formulation ='minimal', angles='Bryant',tmax=0.1)

    np.testing.assert_almost_equal(CG1, CG2, 3)
    np.testing.assert_almost_equal(CG2, CG3, 3)
    np.testing.assert_almost_equal(CG3, CG4, 2)

    np.testing.assert_almost_equal(omega1, omega2, 2)
    np.testing.assert_almost_equal(omega2, omega3, 2)
    np.testing.assert_almost_equal(omega2[0], omega4[0], 1)
    np.testing.assert_almost_equal(omega2[1], omega4[1], 1)
    np.testing.assert_almost_equal(omega2[2], omega4[2], 3)
    pass

if __name__=="__export__":
    time, CG, theta, omega = main(formulation ='minimal', angles='Bryant')
    plt.close('all')
    # --- Plot position
    fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
    fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
    ax.plot( time, CG[0,:], label='x')
    ax.plot( time, CG[1,:], label='y')
    ax.plot( time, CG[2,:], label='z')
    ax.set_xlabel('Time [s]')
    ax.set_ylabel('COG position [m]')
    ax.set_title('3D pendulum - motion')
    ax.legend()
    ax.tick_params(direction='in')

    from welib.tools.repo import export_figs_callback
    export_figs_callback(__file__)

