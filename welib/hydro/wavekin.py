import numpy as np
import pandas as pd
from scipy.optimize import fsolve

try:
    from numpy import trapezoid
except:
    from numpy import trapz as trapezoid

def wavenumber(f, h, g=9.81):   
    """ solves the dispersion relation, returns the wave number k
    INPUTS:
      omega: wave cyclic frequency [rad/s], scalar or array-like
      h : water depth [m]
      g: gravity [m/s^2]
    OUTPUTS:
      k: wavenumber
    """
    omega = 2*np.pi*f
    if hasattr(omega, '__len__'): 
        k = np.array([fsolve(lambda k: om**2/g - k*np.tanh(k*h), (om**2)/g)[0] for om in omega])
    else:
        func = lambda k: omega**2/g - k*np.tanh(k*h)
        k_guess = (omega**2)/g 
        k = fsolve(func, k_guess)[0]
    return k

# Functions 
def elevation2d(a, f, k, eps, t, x=0):
    """  wave elevation (eta) 
    INPUTS:
      a : amplitudes,  scalar or array-like of dimension nf
      f : frequencies, scalar or array-like of dimension nf 
      k : wavenumbers, scalar or array-like of dimension nf
      t : time, scalar or array-like of dimension nt
      x : longitudinal position, scalar or array like of dimension (nx)
    OUTPUTS:
      eta: wave elevation
    """ 
    t   = np.atleast_1d(t)
    a   = np.atleast_1d(a)
    f   = np.atleast_1d(f)
    k   = np.atleast_1d(k)
    eps = np.atleast_1d(eps)
    x   = np.atleast_1d(x)
    omega = 2*np.pi * f
    if len(t)==1:
        eta = np.zeros(x.shape)
        for ai,oi,ki,ei in zip(a,omega,k,eps):
            eta += ai * np.cos(oi*t - ki*x + ei)
    elif len(x)==1:
        eta = np.zeros(t.shape)
        for ai,oi,ki,ei in zip(a,omega,k,eps):
            eta += ai * np.cos(oi*t - ki*x + ei)
    else:
        raise NotImplementedError()

    return eta     

def wave_components( eta, time, water_depth=None, g=9.80665, account_for_t0=True, aThreshold=1e-4, cutAboveThreshold=True):
  """Computes wave components, fp, kp, epsp from a wave elevation time series.

  INPUTS:
  - account_for_t0: if True, the wave phases are adjusted to account for the non zero time offset

  OUPUTS:
   - angular frequencies (fp), wave numbers (kp), and phases (epsp)

  """
  eta  = np.asarray(eta)
  time = np.asarray(time)

  n = len(eta)

  dt = (time[-1]-time[0])/(n-1)
  t0 = time[0]

  # Compute one-sided Fast Fourier Transform and associated angular frequencies
  fft_res = np.fft.rfft(eta)
  fp = np.fft.rfftfreq(n, dt)
  omp = 2 * np.pi *fp

  # Calculate wave amplitudes ap matching OpenFAST single-sided scaling
  ap = np.zeros_like(fp)
  if n % 2 == 0:
    ap[0] = np.abs(fft_res[0]) / n
    ap[1:-1] = 2.0 * np.abs(fft_res[1:-1]) / n
    ap[-1] = np.abs(fft_res[-1]) / n
  else:
    ap[0] = np.abs(fft_res[0]) / n
    ap[1:] = 2.0 * np.abs(fft_res[1:]) / n

  # Determine wave phases epsp matching OpenFAST sign conventions
  epsp = np.angle(fft_res) 
  if account_for_t0:
      epsp += - omp * t0
  epsp = np.mod(epsp, 2*np.pi)

  b=np.abs(ap)>aThreshold
  if np.any(b):
      # Find the index of the very last amplitude above the threshold
      last_i = np.where(b)[0][-1]
      ap  [last_i + 1:] = 0.0
      epsp[last_i + 1:] = 0.0
      cutoff_freq = fp[last_i]
      if cutAboveThreshold:
          ap   = ap  [:last_i]
          fp   = fp  [:last_i]
          epsp = epsp[:last_i]

  # Solve the linear dispersion relation omega^2 = g * k * tanh(k * h) for kp
  if water_depth is not None:
      kp = wavenumber(fp, water_depth, g=g) # Wave numbers
  else:
      kp = np.zeros_like(fp) * np.nan
  pSS = {"ap": ap, "fp": fp, "epsp": epsp,  "kp": kp}
  return pSS
# 
def kinematics2d(a, f, k, eps, h, t, z, x=None, Wheeler=False, eta=None): 
    """ 
    2D wave kinematics, longitudinal velocity and acceleration along x 

    z ^
      |
      |--> x   z=0 (sea level)

      -> vel(z,t)

    ~~~~~~      z=-h (sea bed)

    INPUTS:
      a : amplitudes,  scalar or array-like of dimension (nf)
      f : frequencies, scalar or array-like of dimension (nf)
      k : wavenumbers, scalar or array-like of dimension (nf)
      t : time, scalar or array-like of dimension nt
      z : vertical position, scalar or 1d or nd-array-like of dimension(s) (n x ..). NOTE: z=0 sea level, z=-h sea floor
      x : longitudinal position, scalar or 1d or nd-array-like of dimension(s) (n x ..)
    OUTPUTS:
      vel: wave velocity at t,z,x
      acc: wave acceleartion at t,z,x
    """
    t   = np.atleast_1d(t)
    f   = np.atleast_1d(f)
    a   = np.atleast_1d(a)
    eps = np.atleast_1d(eps)
    k   = np.atleast_1d(k)
    z   = np.atleast_1d(z)
    if x is None:
        x=z*0
    else:
        x = np.asarray(x)
    omega = 2 * np.pi * f  # angular frequency

    if Wheeler:
        if eta is None:
            raise Exception('Provide wave elevation (eta), scalar, for Wheeler')

        # User need to provide eta for wheeler stretching
        if len(t)==1:
            z = (z-eta)*h/(h+eta)
        else:
            raise NotImplementedError('Wheeler stretching, need to consider cases where t is function of time')

    z = z+h # 0 at sea bed
        
    if len(t)==1:
        vel = np.zeros(z.shape) 
        acc = np.zeros(z.shape)
        for ai,oi,ki,ei in zip(a,omega,k,eps):
            vel += oi   *ai * np.cosh(ki*z) / np.sinh(ki*h) * np.cos(oi*t-ki*x + ei)
            acc -= oi**2*ai * np.cosh(ki*z) / np.sinh(ki*h) * np.sin(oi*t-ki*x + ei)
    elif len(z)==1:
        vel = np.zeros(t.shape) 
        acc = np.zeros(t.shape)
        for ai,oi,ki,ei in zip(a,omega,k,eps):
            vel += oi   *ai * np.cosh(ki*z) / np.sinh(ki*h) * np.cos(oi*t-ki*x + ei)
            acc -= oi**2*ai * np.cosh(ki*z) / np.sinh(ki*h) * np.sin(oi*t-ki*x + ei)
    else:
        # most likely we have more time than points, so we loop on points
        vel = np.zeros(np.concatenate((z.shape, t.shape)))
        acc = np.zeros(np.concatenate((z.shape, t.shape)))
        for j in np.ndindex(x.shape): # NOTE: j is a multi-dimension index
            for ai,oi,ki,ei in zip(a,omega,k,eps):
                vel[j] += oi   *ai * np.cosh(ki*z[j]) / np.sinh(ki*h) * np.cos(oi*t-ki*x[j] + ei)
                acc[j] -= oi**2*ai * np.cosh(ki*z[j]) / np.sinh(ki*h) * np.sin(oi*t-ki*x[j] + ei)
    return vel, acc



def fcalc(f, h, g, D, ap, CD, CM , rho, t, z, x, k, eps, phi, u_struct): 
    t = np.atleast_1d(t)                                                 
    f = np.asarray([f])
    ap = np.asarray([ap])
    # Cut Phi and Z at water level z = 0 (waves dont act above water height)
    zindex = sum(z<0)
    phi = phi[0:zindex]
    D = D[0:zindex]
    CD = CD[0:zindex]
    CM = CM[0:zindex]
    z = z[0:zindex]
    u_struct = u_struct[0:zindex]
    
    # Wave kinematics
    u,du = kinematics2d(ap, f, k, eps, h, t, z, x=x, Wheeler=False)
    u=u.reshape(len(z),1) - u_struct.reshape(len(z),1)
    du=du.reshape(len(z),1)

    # Morison
    P = (.5 * rho * CD * D * u * np.abs(u)) + (rho * CM * np.pi*(D**2)/4 * du) #N/m inline Force from wave
    GF = trapezoid((P*phi.reshape(len(phi), 1)), z , axis = 0 ) #work - N       #Generalized force
    M = P * (z.reshape(len(z), 1) + h) #Nm/m [N]

    F_total = trapezoid(P, z, axis=0)                                           
    M_total = trapezoid(M, z, axis=0) #[Nm]

    return u, du, GF, P, M, F_total, M_total


if __name__ == '__main__':
    pass
