""" 
Functions for a single DOF, second order system

"""
import numpy as np
import welib.system.secondorder as so

# --------------------------------------------------------------------------------}
# --- Unforced/Homogeneous/Harmonic responses
# --------------------------------------------------------------------------------{
def harmonic_vibration(vt, x0, xdot0, omega0, zeta=0):
    """ 
    harmonic/free vibrations, constants determined based on initial conditions
    """
    omegad = omega0 * np.sqrt(1-zeta**2)
    if zeta<1:
        # --- Underdamped
        if x0==0 and xdot0==0:
            c1=0 # A
            c2=0 # psi
        elif x0==0:
            c2=0 # psi
            c1=xdot0/omegad # A
        else:
            c2=np.arctan((x0 * omegad)/(xdot0+zeta*omega0*x0))  # psi
            c1 = x0/np.sin(c2)
    elif zeta>1:
        pass
    else: # zeta==1:
        # --- Critically damped
        c1 = x0
        c2 = xdot0 + omega0*x0
    x,xdot = _harmonic_vibration_raw(vt, c1, c2, omega0, zeta=zeta)
    return x, xdot, c1, c2


def _harmonic_vibration_raw(vt, c1, c2, omega0, zeta=0):
    """ 
    harmonic/free vibrations
    """
    if zeta<1:
        # --- Underdamped
        x = c1 * np.exp(- zeta * omega0 * vt)*np.sin(omega0 * np.sqrt(1-zeta**2)* vt + c2)
        omegad = omega0 * np.sqrt(1 - zeta ** 2)
        # TODO combine sin and cos
        xdot        = c1 * np.exp(- zeta * omega0 * vt)*( omegad * np.cos(omegad * vt + c2) - zeta * omega0 * np.sin(omegad * vt + c2))
        xdot_approx = c1 * np.exp(- zeta * omega0 * vt) * omegad * np.cos(omegad * vt + c2)

    elif zeta>1:
        # --- Overdamped
        # TODO
        x = (c1*np.cosh(omegac*vt) + c2*np.sinh(omegac*vt))* np.exp(-zeta*omega0*vt)
        xdot=0

    else: # zeta==1:
        # --- Critically damped
        # TODO
        x = [c1 + c2*vt]* np.exp(omega0*vt)
        xdot=0
    return x, xdot



# --------------------------------------------------------------------------------}
# --- Forced damped vibrations
# --------------------------------------------------------------------------------{
def forced_vibration_particular_cst(frat, F0_over_k, zeta):
    """ 
    Constants for the particular solution (H0 and Phi), x=H0 sin(Omega t - Phi)
    INPUTS:
       frat: freqncy ratio Omega/omega0
    """
    H0  = (F0_over_k) / np.sqrt((1 - frat**2)** 2 + (2*zeta*frat)**2)
    phi = np.arctan2(-2 * zeta * frat , (1 - frat ** 2)) # NOTE sign convention for sin(Om*t-phi)
    return H0, phi

def forced_vibration_particular(vt, k, m, F0, Omega, zeta):
    """
    Particualr solution to forced harmonic vibrations, , x=H0 sin(Omega t - Phi)
    """
    omega0       = np.sqrt(k / m)
    H0, phi      = forced_vibration_particular_cst(Omega/omega0, F0/k, zeta)
    x_particular = H0 * np.sin(Omega * vt + phi)
    return x_particular

def forced_vibration_transient(vt, k, m, F0, Omega, zeta, x0, xdot0):
    omega0  = np.sqrt(k / m)
    H0, phi = forced_vibration_particular_cst(Omega/omega0, F0/k, zeta)
    if x0==0 and xdot0==0:
        psi     = np.arctan(np.sqrt(1-zeta**2)/(zeta + Omega / omega0 * np.cos(phi) / np.sin(phi)))
        A       =-H0 * np.sin(phi) / np.sin(psi)
    else:
        # Initial conditions without the particular solution
        x0t    = x0   - H0*np.sin(phi)
        xdot0t = xdot0- H0*Omega*np.cos(phi)
        _,_,A,psi = harmonic_vibration(0, x0t, xdot0t, omega0, zeta)
    x_transient= A * np.exp(- zeta * omega0 * vt) * np.sin(omega0 * np.sqrt(1-zeta**2)*vt + psi)
    return x_transient

def forced_vibration(vt, k, m, F0, Omega, zeta, x0=0, xdot0=0):
    x_particular = forced_vibration_particular(vt, k, m, F0, Omega, zeta)
    x_transient  = forced_vibration_transient(vt, k, m, F0, Omega, zeta, x0, xdot0)
    x = x_transient + x_particular
    return x


# --------------------------------------------------------------------------------}
# --- Characteristic transient responses
# --------------------------------------------------------------------------------{
def impulse_response(time, m, c, k, both=False):
    """ 
    See impulse_response in secondorder.py 
    """
    zeta   = c/(2*np.sqrt(k*m))
    omega0 = np.sqrt(k/m)
    return so.impulse_response(time, omega0, zeta, b=1/m, both=both)

def step_response(time, m, c, k, both=False):
    """ 
    See step_response in secondorder.py 
    """
    zeta   = c/(2*np.sqrt(k*m))
    omega0 = np.sqrt(k/m)
    return so.step_response(time, omega0, zeta, b=1/m, both=both)


# --------------------------------------------------------------------------------}
# --- Time integration 
# --------------------------------------------------------------------------------{

def duhamel(time, F, m, c, k, both=False):
    r""" 
    See duhamel in secondorder.py 

    Compute time response of a single DOF system using Duhamel's integral

        x(t)    = \int_0^t F(t') H(t-t') dt' 
        xdot(t) = \int_0^t F(t') H'(t-t') dt'
    with
     - F: force
     - H: impulse response function
     - H' derivative of impulse function

    NOTE: assumes that initial conditions are (0,0)
    If this is not the case, some transient responses need to be added.

    INPUTS:
      - time : 1d-array of time stamps
      - F    : 1d-array of force at each time step
      - m,c,k: scalar mass, damping, stiffness
    OUTPUTS:
      - x, xd: 1d-array, of position and speed
    """
    zeta   = c/(2*np.sqrt(k*m))
    omega0 = np.sqrt(k/m)
    return so.duhamel(time, omega0, zeta, F, b=1/m, both=both)

if __name__=='__main__':
    pass
