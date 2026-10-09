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





def wave_components(eta, time, water_depth=None, g=9.80665, account_for_t0=True, 
                    cutAboveThreshold=True, aThreshold=1e-4, 
                    minimize_nComp=True, target_R2=0.001, target_eps=0.02,
                    eta_out = False,
                    test=False,
                    plot=False,
                    verbose=False,
                    ):
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
  
    # --- Cut above a threshold
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

    # Reduce number of components 
    if minimize_nComp:
        if water_depth is None:
            raise Exception('Provide water_depth for component minimization testing')
        if verbose:
            print('Reducing number of wave components')
        ap, fp, epsp, kp = reduce_wave_components(ap, fp, epsp, kp, time, eta, target_R2=target_R2, target_eps=target_eps, verbose=verbose)

    pSS = {"ap": ap, "fp": fp, "epsp": epsp,  "kp": kp}


    # --- Testing reconstruction
    if test or plot or eta_out:
        if water_depth is None:
            raise Exception('Provide water_depth')
        eta_sim = elevation2d(pSS['ap'], pSS['fp'], pSS['kp'], pSS['epsp'], time, x=0)
        pSS['time_sim'] = time
        pSS['eta_sim'] = eta_sim
    if plot:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8))
        fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
        ax.plot(time, eta, label='Ref')
        ax.plot(time, eta_sim, label='Sim')
        ax.set_xlabel('Time [s]')
        ax.set_ylabel('Wave elevation [m]')
        ax.legend()
        pSS['ax'] = ax
    if test:
        from welib.tools.stats import rsquare, mean_rel_err
        R2,_ = rsquare(y_ref=eta, y_sim=eta_sim)
        eps  = mean_rel_err(y1=eta, y2=eta_sim, method='meanabs')

        #np.testing.assert_almost_equal(eta_sim, eta, test_decimal)
        np.testing.assert_array_less(1-R2, target_R2)
        np.testing.assert_array_less(eps , target_eps)

    return pSS

def reduce_wave_components(ap, fp, epsp, kp, time, eta, target_R2, target_eps, verbose=False):
    """Reduces wave components using binary search to minimize count while satisfying error tolerances."""
    from welib.tools.stats import rsquare, mean_rel_err
    if target_R2 is None and target_eps is None:
        return ap, fp, epsp, kp

    # Sort components by amplitude descending to prioritize high energy terms
    idx_sort = np.argsort(ap)[::-1]
    ap_s     = ap[idx_sort]
    fp_s     = fp[idx_sort]
    epsp_s   = epsp[idx_sort]
    kp_s     = kp[idx_sort]


    # Binary search to find the minimum number of components meeting tolerance
    low = 1
    high = len(ap_s)
    best_M = high

    while low <= high:
        mid = (low + high) // 2
        eta_test = elevation2d(ap_s[:mid], fp_s[:mid], kp_s[:mid], epsp_s[:mid], time, x=0)

        r2_val, _ = rsquare(y_ref=eta, y_sim=eta_test)
        err_val = mean_rel_err(y1=eta, y2=eta_test, method='meanabs')
        if verbose:
            print(f'ncomp: {mid:6d}  r2= {r2_val:7.5f}  eps={err_val:4.3f}%')

        met_r2 = target_R2 is None or r2_val >= target_R2
        met_err = target_eps is None or err_val <= target_eps

        if met_r2 and met_err:
            best_M = mid
            high = mid - 1
        else:
            low = mid + 1

    # Retain the optimal subset and resort chronologically by frequency
    ap_out = ap_s[:best_M]
    fp_out = fp_s[:best_M]
    epsp_out = epsp_s[:best_M]
    kp_out = kp_s[:best_M]

    idx_freq = np.argsort(fp_out)
    return ap_out[idx_freq], fp_out[idx_freq], epsp_out[idx_freq], kp_out[idx_freq]




def kinematics2d(a, f, k, eps, h, t, z, x=None, Wheeler=False, eta=None): 
    r""" 
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

    cosh(k z)/ \sinh(k h) = e^{-k(h - z)} + e^{-k(h + z)} + e^{-k(3h - z)} + ...

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

    if f[0]==0:
        f   = f[1:]
        a   = a[1:]
        eps = eps[1:]
        k   = k[1:]
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
#         np.seterr(all='raise')
        vel = np.zeros(z.shape) 
        acc = np.zeros(z.shape)
        for ai,oi,ki,ei in zip(a,omega,k,eps):
            #exponent = -ki* (h*z)
            hyp_ratio = np.cosh(ki*z) / np.sinh(ki*h)
            hyp_ratio[np.isnan(hyp_ratio)] = 0
            #hyp_ratio = np.exp(-ki*(h-z))
            #hyp_ratio = np.where(exponent < -500.0, 0.0, np.exp(exponent)) # BUGGY
#             except:
#                 import pdb; pdb.set_trace()
            vel += oi   *ai * hyp_ratio * np.cos(oi*t-ki*x + ei)
            acc -= oi**2*ai * hyp_ratio * np.sin(oi*t-ki*x + ei)
    elif len(z)==1:
        vel = np.zeros(t.shape) 
        acc = np.zeros(t.shape)
        for ai,oi,ki,ei in zip(a,omega,k,eps):
            hyp_ratio = np.cosh(ki*z) / np.sinh(ki*h)
            vel += oi   *ai * hyp_ratio * np.cos(oi*t-ki*x + ei)
            acc -= oi**2*ai * hyp_ratio * np.sin(oi*t-ki*x + ei)
    else:
        # most likely we have more time than points, so we loop on points
        vel = np.zeros(np.concatenate((z.shape, t.shape)))
        acc = np.zeros(np.concatenate((z.shape, t.shape)))
        for j in np.ndindex(x.shape): # NOTE: j is a multi-dimension index
            for ai,oi,ki,ei in zip(a,omega,k,eps):
                hyp_ratio = np.cosh(ki*z[j]) / np.sinh(ki*h)
                vel[j] += oi   *ai * hyp_ratio * np.cos(oi*t-ki*x[j] + ei)
                acc[j] -= oi**2*ai * hyp_ratio * np.sin(oi*t-ki*x[j] + ei)
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
