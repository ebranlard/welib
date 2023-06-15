""" 
Generalized alpha algorithm, for second order dynamical systems.

- ga_step,  ga_integrate, 

   fQ(t, q, qd) = qdd, acceleration 
   fdQdq(t, q, qd) =  Jacobian of acceleration wrt q  ("Instantenous stiffness matrix")
   fdQdqd(t, q, qd) = Jacobian of acceleration wrt qd ("Instantenous damping matrix")

"""
import numpy as np
from scipy.optimize import OptimizeResult as OdeResultsClass 

def ga_step(t, fQ, fdQdq, fdQdqd, q, qd, qdd, a, dt, AlphaM, Beta, Gamma, C, MaxIter=20, SolveTol=1e-5):
    
    # Predict next state (assumes acceleration is constant)
    q += dt*qd + dt**2*(0.5-Beta)*a  # Position (1)
    qd += dt*(1-Gamma)*a    # Velocity (1)
    a = (qdd - AlphaM*a) / (1-AlphaM)
    q += dt**2*Beta*a
    qd += dt*Gamma*a

    # Convergence iterations
    for k in range(MaxIter):

        # Calculate matrices and vector for jacobian calc
        dQdq_n = fdQdq(t, q, qd) 
        dQdqd_n = fdQdqd(t, q, qd)
        Q_n = fQ(t, q, qd)

        # Construct Jacobian
        J = np.eye(len(q)) - dt*Gamma*C*dQdqd_n - dt**2*Beta*C*dQdq_n
        RHS = qdd - Q_n

        # Calculate change in acceleration
        delta_qdd = -np.linalg.solve(J, RHS)

        # If deltas are sufficiently small, exit
        norm = np.linalg.norm(delta_qdd, 2)
        if norm < SolveTol:
            break

        # Update states with deltas
        q += dt**2*Beta*C*delta_qdd     # Displacement
        qd += dt*Gamma*C*delta_qdd      # Velocity
        qdd += delta_qdd                # Acceleration
        a += C*delta_qdd                # Algorithmic Acceleration

    return q, qd, qdd, a, k, norm


def ga_integrate(fQ, fdQdq, fdQdqd, t, x0, qdd0=None, RhoInf=0, MaxIter=20, SolveTol=1e-5):
    nt = len(t)
    nx = len(x0)

    # Time storage
    y  = np.zeros((nx, nt))
    norms = np.zeros(nt)
    iters = np.zeros(nt)

    # Initial conditions
    y[:,0] = x0
    nDOF = int(nx/2)
    q = x0[:nDOF]
    qd = x0[nDOF:]
    if qdd0 is None:
        qdd = fQ(t[0], q, qd)
        #qdd = np.zeros(nDOF) # TODO chose
    else:
        qdd = qdd0

    a = np.zeros(nDOF)  # Initial algorithmic accel, TODO 
    #a = qdd  # Initial algorithmic accel, TODO 

    # Parameters (determined using RhoInf)
    AlphaM = (2.0 * RhoInf - 1.0) / (RhoInf + 1.0)
    AlphaF = RhoInf / (RhoInf + 1.0)
    Gamma = 0.5 - AlphaM + AlphaF
    Beta = (1.0 - AlphaM + AlphaF)**2 / 4.0
    C = (1-AlphaF)/(1-AlphaM)

    # Time stepping
    for k in range(nt-1):
        dt = t[k+1]-t[k]
        q, qd, qdd, a, niter, norm = ga_step(t[k], fQ, fdQdq, fdQdqd, q, qd, qdd, a, dt, AlphaM, Beta, Gamma, C, MaxIter=MaxIter, SolveTol=SolveTol)

        y[:nDOF,k+1] = q
        y[nDOF:,k+1] = qd
        norms[k] = norm
        iters[k] = niter

    res = OdeResultsClass(t=t, y=y) # To mimic result class of solve_ivp
    res['norm'] = norms
    res['iter'] = iters
    return res


