""" 
Tools for stochastic process, in particular, the methods:

- sample_from_autospectrum: to generate a time series based on a spectrum
- autocovariance_* : compute the autocoviariance
- autospectrum_* : compute the autospectrum



# TODO for autocovariance look at
# http://statsmodels.sourceforge.net/devel/generated/statsmodels.tsa.stattools.acovf.html http://statsmodels.sourceforge.net/devel/generated/statsmodels.tsa.stattools.ccovf.html
"""



import numpy as np
from numpy.random import uniform


# --------------------------------------------------------------------------------
# --- Helper functions for Time series generation from auto spectrum
# --------------------------------------------------------------------------------
def _getPhasesAmplitudesFactors(N1, phases=None, amplitudesFactors=None, randomAmplitudes=True):
    if randomAmplitudes:
        if amplitudesFactors is not None:
            if len(amplitudesFactors)!=N1:
                raise Exception('Length of amplitudes should be {}'.format(N1))
            if amplitudesFactors[0]!=1:
                raise Exception('First amplitudeFactor (for zero frequency) must be 1 by convention.')
            U1 = amplitudesFactors
        else:
            U1 = uniform(0, 1, N1)  # amplitude
            U1[0] = 1
    else:
        U1 = np.ones(N1)*np.exp(-1)

    if phases is not None:
        if len(phases)!=N1:
            raise Exception('Length of phases should be half the length of time vector')
        if phases[0]!=0:
            raise Exception('First phase (for zero frequency) must be zero by convention.')
        phi1 = phases
    else:
        phi1 = uniform(0, 2*np.pi, N1)  # random phases between 0 and 2pi
        phi1[0] = 0
    return phi1, U1

def _getX1_sumcos(N, Sf_or_So, df_or_dom, phases=None, amplitudesFactors=None, randomAmplitudes=True):
    N1 = len(Sf_or_So)
    phi1, U1 = _getPhasesAmplitudesFactors(N1, phases=phases, amplitudesFactors=amplitudesFactors, randomAmplitudes=randomAmplitudes)
    A1 = np.sqrt(2*Sf_or_So * df_or_dom) 
    A2 = np.sqrt(-2*np.log(U1))/np.sqrt(2)
    A2[0]=1
    a1 = A1 * A2 * np.exp(  1j * phi1)
    X1 = N/2*a1
    X1[0] = np.sqrt(2)*X1[0]
    return X1, A1, A2, a1, phi1

def _getX1_boxmuller(N, Sf_or_So, df_or_dom, phases=None, amplitudesFactors=None, randomAmplitudes=True):
    # TODO This is the same as sumcos
    # TODO I'm not sure whether we should use N or N1 different random numbers here
    N1 = len(Sf_or_So)
    phi1, U1 = _getPhasesAmplitudesFactors(N1, phases=phases, amplitudesFactors=amplitudesFactors, randomAmplitudes=randomAmplitudes)
    # Decomposition 1
    A1 = np.sqrt(2 * df_or_dom * Sf_or_So) 
    A2 = np.sqrt(-2*np.log(U1))/np.sqrt(2);
    A2[0]=1
    a1 = A1 * A2 * np.exp(1j * phi1)
    X1 = N/2*a1
    # Decomposition 2
    A1_ = np.sqrt(N * df_or_dom * Sf_or_So/2)
    W1_ = np.sqrt(N/2) * np.sqrt(-2*np.log(U1))*np.exp(1j * phi1)
    X1_ = A1_*W1_
    if np.max(np.abs(X1-X1_))>1e-10:
        raise Exception('Problem in formulation')
    X1[0] = np.sqrt(2)*X1[0]
    return  X1, A1, A2, a1, phi1

def _getX1_basic(N, Sf_or_So, df_or_dom):
    N1 = len(Sf_or_So)
    A1 = np.sqrt(2*Sf_or_So * df_or_dom) 
    A2 = 1
    phi1=0
    a1 = A1 * A2 
    X1 = N/2*a1
    X1[0] = np.sqrt(2)*X1[0]
    return X1, A1, A2, a1, phi1


# --------------------------------------------------------------------------------
# --- Main function for time series generation from auto spectrum 
# --------------------------------------------------------------------------------
def sample_from_autospectrum(tMax, dt, f_S, angularFrequency=True, 
        method='sumcos-irfft', 
        seed=None, phases=None, amplitudesFactors=None, randomAmplitudes=True):
    """ 
    Generate a time series x, based on an single-sided autospectrum

    INFO:
        single-sided =  2 double-sided
                S_X  = 2 S_XX

    Fundamental relation of FFT/IFFT:
        N = 1/ (dt * df)

        exp ( i 2 pi k n / N  )  = exp ( i 2 pi  f_k t_n )   for k in [0..N-1], n in [0..N-1]


    INPUTS:
     - tMax: maximum time of time series
     - dt  : time increament 
     - f_S: single-sided autospectrum function with interface:
         - S_X_om = f_S(omega)   if angularFrequency is True
         or
         - S_X_f  = f_S(f)       if angularFrequency is False
     - method: method used to generate the time series:
         - 'sumcos-manual': 
         - 'sumcos-irfft': 
         - 'sumcos-idft-manual': 
         - 'sumcos-idft-vectorized': 
         - 'sumcos-idft-ifft': 
         - 'boxmuller-*': 
         - 'ifft': inverse FFT (fastest)
     - randomAmplitudes:  if True, randomized the amplitudes (not only the phases) 
     - phases: vector of phases (between 0 and 2pi) to use for all positive frequencies. 
     - amplitudesFactors:  vectors of amplitudes (0,1) to use for all positive frequencies

    OUTPUTS:
     - t: time vector
     - x: time series
     - f_or_om: frequency in [Hz] if not angularFrequency, or in [rad/s] if angularFrequency
     - Sf_or_So: Spectrum corresponding to f_or_om

    """
    methods = method.lower().strip().split('-')
    # -- Define time vector
    N  = int(np.round(tMax/dt)) + 1
    t     = np.arange(0,tMax+dt/2,dt)
    # -- Define frequencies and spectrum for single sided approaches
    f1     = np.fft.rfftfreq(N, dt) # NOTE: when N is even, we have an extra frequency
    fMax   = f1[-1]                 # approx (1/dt)/2  Highest frequency
    df1    = 1/(dt*N)               # approx  1./tMax  Frequency resolution
    N1     = len(f1)
    omega1 = 2*np.pi*f1
    # -- Define frequencies and spectrum
    if angularFrequency:
        Sf_or_So  = np.array([f_S(omi)     for omi in omega1])
        df_or_dom = 2*np.pi*df1
        f_or_om   = omega1
    else:
        Sf_or_So  = np.array([f_S(fi)      for fi in f1])
        df_or_dom = df1
        f_or_om   = f1
    if seed is not None:
        # Initialized random number generator
        np.random.seed(seed)
    
    if methods[0] in ['sumcos', 'boxmuller']:
        # --------------------------------------------------------------------------------
        # --- SUM OF COSINES / Box-Muller
        # --------------------------------------------------------------------------------
        X1, A1, A2, a1, phi1 = _getX1_sumcos(N, Sf_or_So, df_or_dom, phases=phases, amplitudesFactors=amplitudesFactors, randomAmplitudes=randomAmplitudes)

        if methods[1] == 'manual':
            A1 = A1*A2
            x = np.zeros(t.shape)
            for ai,oi,ei in zip(A1,omega1,phi1):
                x += ai * np.cos(oi*t + ei)

        elif methods[1] == 'irfft' :
            x = np.fft.irfft(X1, N)

        elif methods[1] == 'idft':
            from welib.tools.spectral import double_sided_DFT_real, check_DFT_real, IDFT

            if np.mod(N,2)==0:
                print('[WARN] IDFT introduces small error for now when N is even')
                # Keep me: what we use to do:
                #   f1 = np.fft.rfftfreq(N, dt)[:-1] # NOTE the last frequency should not be used
            X = double_sided_DFT_real(X1, N=N)
            check_DFT_real(X)
            x = IDFT(X, method=methods[2])
            x = np.real(x)

        else:
            NotImplementedError('Method {}'.format(method))

    elif methods[0]=='randn':
        # --------------------------------------------------------------------------------
        # --- RANDOM NORMAL complex numbers 
        # --------------------------------------------------------------------------------
        X1_bas, A1, A2, a1, phi1 = _getX1_basic(N, Sf_or_So, df_or_dom)

        a = np.random.randn(N1)
        b = np.random.randn(N1)
        # First part of the complex vector
        X1_rand = (a + b * 1j )/np.sqrt(2)
        X1_rand[0] = 1
        X1 = X1_bas * X1_rand

        x = np.fft.irfft(X1, N)

        #X = double_sided_DFT_real(X1, N=N)
        #x=np.fft.ifft(X)
        #x = np.real(x)

    else:
        raise NotImplementedError('Method {}'.format(method))

    return t, x, f_or_om, Sf_or_So










def autocovariance_num(x, **kwargs):
    from welib.tools.signal_analysis import correlation
    sigma2 = np.var(x)
    rho_XX, tau = correlation(x, **kwargs)
    kappa_XX = rho_XX*sigma2
#     autoCov = 0
#     for i in np.arange(0, N-k):
#         autoCov += ((Xi[i+k])-Xs)*(Xi[i]-Xs)
#     return (1/(N-1))*autoCov
# def lagged_auto_cov(Xi,t):
#     """
#     for series of values x_i, length N, compute empirical auto-cov with lag t
#     defined: 1/(N-1) * \sum_{i=0}^{N-t} ( x_i - x_s ) * ( x_{i+t} - x_s )
#     """
#     N = len(Xi)
# 
#     # use sample mean estimate from whole series
#     Xs = np.mean(Xi)
# 
#     # construct copies of series shifted relative to each other, 
#     # with mean subtracted from values
#     end_padded_series = np.zeros(N+t)
#     end_padded_series[:N] = Xi - Xs
#     start_padded_series = np.zeros(N+t)
#     start_padded_series[t:] = Xi - Xs
# 
#     auto_cov = 1./(N-1) * np.sum( start_padded_series*end_padded_series )
#     return auto_cov

    return kappa_XX, tau


def autocovariance_fft(omega, S_X, onesidedIn=True):
    """ 
     Return Kappa_XX autocovariance from oneside (S_X) or double sided (S_XX) spectrum
       - Using FFT
       - Based on discrete values of the S_X

      S_X(\omega) = 2 S_XX(\omega)

    """
    if onesidedIn:
        k = np.fft.irfft(S_X)
        om_max = omega[-1]
        n      = len(omega)
        domega = (omega[-1]-omega[0])/(n-1)
        Fs     = 1/domega       
        tau    = np.arange(n)*Fs/(n)*np.pi
        K_X    = k[:n]*om_max #*2
    else:
        raise NotImplementedError()

    return tau, K_X

def autocovariance_int(tau, f_S_X, omega_max):
    """
    Return autocovariance function
    """
    from scipy.integrate import nquad
    K_XX = np.zeros(tau.shape) 
    for i, taui in enumerate(tau):
        integrand = lambda om : np.cos(om * taui) * f_S_X(om)
        K_XX[i], _, _ = nquad(integrand, [[0, omega_max]] , full_output=1)
    return K_XX


def autocorrcoeff_num(x, **kwargs):
    from welib.tools.signal_analysis import correlation
    rho_XX, tau = correlation(x, **kwargs)
    return rho_XX, tau

def autospectrum_fft(tau, kappa_XX, onesided=True, method='fft_wrap'):
    """ 
     Return one-sided (S_X) or double-sided (S_XX) autospectrum:
       - Using FFT
       - Based on discrete values of the auto-covariance kappa_XX.

      S_X(\omega) = 2 S_XX(\omega)

    """
    if method=='fft_wrap':
        from welib.tools.spectral import fft_wrap
        f, S_X, Info = fft_wrap(tau, kappa_XX, output_type='amplitude', averaging='None', detrend=False)
        omega =2*np.pi * f
    else:
        print('[WARN] autospectrum_FFT might need work')
        k = np.fft.rfft(kappa_XX)
        tau_max = tau[-1]
        n      = int(len(tau)/2)
        dtau   = (tau[-1]-tau[0])/(n-1)
        Fs     = 1/dtau       
        omega  = np.arange(n)*Fs/(n) * 2 *np.pi
        S_X    = np.real(k[:n])/n * tau_max / (np.pi)
#     print('>> k', np.real(k[:5]))
#     print('>> k', np.abs(k[:5]))
#     print('>> S', S_X[:5])
#     print('>>> tau+max',tau_max)
#     print('>>> len(k)',len(k))
#     print('>>> len(S_X)',len(S_X))
#     print('>>> len(om  )',len(om))

    return omega, S_X

def autospectrum_int(omega, f_kappa, tau_max):
    """
    Return one-sided (S_X) or double-sided (S_XX) autospectrum:
      - Using numerical quadrature integration
      - Based on a "continuous" function kappa_XX(tau)
    """
    from scipy.integrate import nquad
    S_X = np.zeros(omega.shape) 
    # for i, omi in enumerate(om_th):
    #     f_S_X_integrand = lambda tau : 1/np.pi * np.exp(-1j * omi * tau) * f_kappa(tau)
    #     S_X_i[i], _ = nquad(f_S_X_integrand, [[0, tau_max]] )
    for i, omi in enumerate(omega):
        f_S_X_integrand = lambda tau : 2/np.pi * np.cos(omi * tau) * f_kappa(tau)
        S_X[i], _, _  = nquad(f_S_X_integrand, [[0, tau_max]] , full_output=1)
    return S_X

