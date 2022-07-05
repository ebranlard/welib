r"""
Set of tools/functions for a second order system expressed in scalar form. 

NOTES: 
    For multiple dimensions, use mech_system
    For a forcing that is a function of states (i.e. f(t,x,xdot)), use mech_system


The characteristic equation considered in this package is:
     
       x_ddot(t) + 2 zeta omega0 x_dot(t) + omega0^2 x(t)  =  b u(t)     (1) 
       
Example for a mass-spring damper: 
     
    omega0=sqrt(k/m), zeta = c/(2*np.sqrt(k*m)), b=1/m if u(t) is the force

"""
import numpy as np
from numpy import sqrt, cos, sin, exp

# --------------------------------------------------------------------------------}
# --- Unforced/homogeneous  responses
# --------------------------------------------------------------------------------{

# --------------------------------------------------------------------------------}
# --- Characteristic transient responses
# --------------------------------------------------------------------------------{
def impulse_response(time, omega0, zeta, b=1, t0=0, A=1, q0=[0,0], both=False):
    """ 
    Impulse response function to a delta dirac (impulse),  u(t)=delta(t-t0) in eq.(1)
      H=x
    """
    t = np.asarray(time)-t0

    H  = np.zeros(time.shape)
    Hp = np.zeros(time.shape)
    bp = t>=0 # boolean for positive times
    t=t[bp]

    if zeta<1:
        omegad = omega0 * sqrt(1-zeta**2)
        H[bp] = b*1/(omegad) * sin(omegad * t) * exp(-zeta * omega0 * t)
        #Hp = 1/(m*omegad)*(-zeta * omega0 * np.sin(omegad * time) + omegad* np.cos(omegad * time)) * np.exp(-zeta * omega0 * time)
        Hp[bp] = -zeta * omega0 * H[bp] + b * cos(omegad * t) * exp(-zeta * omega0 * t)

    elif zeta==1:

        H[bp]  += b * t * exp(-omega0*t)
        Hp[bp]  = b *     exp(-omega0*t)  - b * t * omega0 *  exp(-omega0*t)

    else:
        t1 = -1/(-omega0*zeta + omega0*sqrt(zeta**2-1)) # tau=-1/lambda
        t2 = -1/(-omega0*zeta - omega0*sqrt(zeta**2-1))

        H[bp]  += b * 1/(2*omega0*sqrt(zeta**2-1))*(exp(-t/t1) - exp(-t/t2))
        Hp[bp]  = b * 1/(2*omega0*sqrt(zeta**2-1))*(-1/t1*exp(-t/t1) + 1/t2*exp(-t/t2))

    if both:
        # Also ouputs derivative
        return H*A, Hp*A
    else:
        return H*A

def step_response(time, omega0, zeta, b=1, t0=0, A=1, offset=0, q0=[0,0], both=False):
    r"""
    Response to a unit step u(t)=A*Step(t-t0) in eq.(1)

    NOTE: xdot is 0 at t0, which is one difference with first order systems
    """
    t  = np.asarray(time)-t0
    bp = t>= 0 # boolean for positive times
    t  = t[bp]

    x  = np.zeros(time.shape)+offset
    xd = np.zeros(time.shape)

    if zeta<1: 
        omegad = omega0 * sqrt(1-zeta**2)
        phi    = np.arctan2(zeta, sqrt(1-zeta**2))
        x[bp] += A* b/(omega0**2) * ( 1- exp(-zeta*omega0 *t)/sqrt(1-zeta**2) * cos(omegad*t - phi))
        # Formulation 1 - simple derivation of x
        #xd[bp] = b/(omega0**2) * exp(-zeta*omega0*t)/sqrt(1-zeta**2)*(omegad*sin(omegad*t - phi) +zeta * omega0 * cos(omegad*t - phi))
        # Formulation2 - Simplifying the sum of the cos and sin in formulation 1
        #phi_c = np.arctan2(zeta*omega0*sin(-phi)-omegad*cos(-phi) , zeta*omega0*cos(-phi) + omegad*sin(-phi)  )
        #A_c   = sqrt(omega0**2*zeta**2+omegad**2) = omega0
        #xd[bp] = b/(omega0**2) * exp(-zeta*omega0*t)/sqrt(1-zeta**2)*cos(omegad*t+phi_c)
        # Formulation3  - simplification of formulation 2
        # Phase phi_c is -pi/2, and A_c=omega0
        xd[bp] = A* b/(omega0) * exp(-zeta*omega0*t)/sqrt(1-zeta**2)*sin(omegad*t)

    elif zeta==1: 

        x[bp] += A*b/(omega0**2) * ( 1- exp(-omega0*t)-omega0*t*exp(-omega0*t))
        xd[bp] = A*b * t*exp(-omega0*t)

    else:
        # Time constants tau=-1/lambda
        t1 = -1/(-omega0*zeta + omega0*sqrt(zeta**2-1)) # largest one, drives the response
        t2 = -1/(-omega0*zeta - omega0*sqrt(zeta**2-1)) # smallest, will decay quick

        x[bp] += A*b/(omega0**2) * ( 1- 1/(t2-t1)*(t2*exp(-t/t2)- t1*exp(-t/t1)))
        xd[bp] = A*b/(omega0**2*(t2-t1))*(exp(-t/t2) - exp(-t/t1))
        

    if both:
        return x, xd
    else:
        return x


# --------------------------------------------------------------------------------}
# --- Forced vibrations responses 
# --------------------------------------------------------------------------------{



# --------------------------------------------------------------------------------}
# --- Time integrations 
# --------------------------------------------------------------------------------{
def integrate(t_eval, omega0, zeta, u, b=1, q0=[0,0], method='LSODA', both=False, **options):

    if method in ['LSODA','RK45','RK23','DOP853','Radau','BDF']:
        # NOTE:
        #  - 'DOP853', 'Radau' work good for varying input
        #  - 'LSODA' good for adaptative time step
        from scipy.interpolate import interp1d
        from scipy.integrate import  solve_ivp #odeint

        q0=np.asarray(q0).flatten()

        # Create an interpolant for u
        fU = interp1d(t_eval, u)

        A,B = statematrices(omega0, zeta, b=b)
        odefun = lambda t, q : np.dot(A, q) + np.array([0,b*fU(t)]) # + np.dot(B, [0,fU(t)])

        res = solve_ivp(fun=odefun, t_span=[t_eval[0], t_eval[-1]], y0=q0, t_eval=t_eval, method=method, vectorized=False, **options)   

        if both:
            return res.y[0,:], res.y[1,:]
        else:
            return res.y[0,:]

    elif method.lower() in ['duhamel','convolution']:

        return duhamel(t_eval, omega0, zeta, u, q0=q0, b=b, both=both)

    else:
        raise NotImplementedError('Method {}'.format(method))


def duhamel(time, omega0, zeta, u, b=1, q0=[0,0], both=False):
    r""" 
    Compute time response of a single DOF system using Duhamel's integral

        x(t)    = \int_0^t Phi(t-t') *b* u(t')  dt' + Phi(t-t0)x0
        x(t)    = \int_0^t H(t-t')    *  u(t')  dt' + Phi(t-t0)x0

        x(t)    = \int_0^t H(t-t')   u(t')  dt' 
        xdot(t) = \int_0^t H'(t-t')  u(t')  dt'

    with
     - F: force
     - B: e.g. 1/m for mech system
     - Phi: state transition, Phi(t) 
     - H  : impulse response, H(t)   = Phi * b

     NOTE: impulse response for the time derivative has a discontinuity, 
          making it harder to get the velocity right (more points might be needed).
          One might want to somehow smoothen the impulse response Hp, or ensure that enough points are present.

    INPUTS:
      - time : 1d-array of time stamps
      - tau  : time constant
      - u    : 1d-array of inputs at each time step
    OUTPUTS:
      - x, xd: 1d-array, of position and speed

    TODO: initial conditions
    """
    H, Hp = impulse_response(time, omega0, zeta, b=b, t0=time[0], both=True)
    #dt = t_eval[1]-t_eval[0]
    #x  = np.convolve(F.ravel(), H.ravel()  )[:len(t_eval)]*dt
    #xd = np.convolve(F.ravel(), Hp.ravel() )[:len(t_eval)]*dt
    from welib.tools.signal_analysis import convolution_integral
    x  = convolution_integral(time, u.ravel(), H )
    if both:
        xd = convolution_integral(time, u.ravel(), Hp)
        return x, xd
    else:
        return x



# --------------------------------------------------------------------------------}
# --- System
# --------------------------------------------------------------------------------{
def statematrices(omega0, zeta, b=1):
    """ 
    Return state matrices for first order form of eq.(1)
    """
    A = np.array([[0, 1],
                [-omega0**2 , -2* omega0 *zeta] ])
    
    B = np.array([[0],
                  [b]])
    return A,B




