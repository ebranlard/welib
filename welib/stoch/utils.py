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
def _getPhasesAmplitudesFactors(Nf, phases=None, amplitudesFactors=None, randomAmplitudes=True):
    """ return uniformly distributed phases and ampitudes """
    if randomAmplitudes:
        if amplitudesFactors is not None:
            if len(amplitudesFactors)!=Nf:
                raise Exception('Length of amplitudes should be {}'.format(Nf))
            if amplitudesFactors[0]!=np.exp(-1) :
                raise Exception('First amplitudeFactor (for zero frequency) must be exp-1 by convention.')
            AU = amplitudesFactors
            AU[0] = np.exp(-1) # safety
        else:
            AU = uniform(0, 1, Nf)  # amplitude
            AU[0] = np.exp(-1) # used to be 1
    else:
        AU = np.ones(Nf)*np.exp(-1)

    if phases is not None:
        if len(phases)!=Nf:
            raise Exception('Length of phases should be half the length of time vector')
        if phases[0]!=0:
            raise Exception('First phase (for zero frequency) must be zero by convention.')
        phiU = phases
    else:
        phiU = uniform(0, 2*np.pi, Nf)  # random phases between 0 and 2pi
        phiU[0] = 0
    return phiU, AU

def _getX_sumcos(N, Sf_or_So, df_or_dom, phases=None, amplitudesFactors=None, randomAmplitudes=True):
    """ Note same a Box Muller method"""
    Nf = len(Sf_or_So)
    phi, U = _getPhasesAmplitudesFactors(Nf, phases=phases, amplitudesFactors=amplitudesFactors, randomAmplitudes=randomAmplitudes)
    A1 = np.sqrt(2 * Sf_or_So * df_or_dom) 
    A2 = np.sqrt(-np.log(U)) #A2[0]=1
    A = A1 * A2  
    X = N/2 * A * np.exp(1j * phi)
    X[0] = np.sqrt(2)*X[0]
    return X, A, phi

def _getX_basic(N, Sf_or_So, df_or_dom):
    Nf = len(Sf_or_So)
    A = np.sqrt(2 * Sf_or_So * df_or_dom) 
    phi=0
    X = N/2 * A * np.exp(1j * phi)
    X[0] = np.sqrt(2)*X[0]
    return X, A, phi


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
         - 'sumcos-irfft':  fastest
         - 'sumcos-idft-manual': 
         - 'sumcos-idft-vectorized': 
         - 'sumcos-idft-ifft': 
         - 'boxmuller-*': 
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
    #print('N', N, 'len(t)',len(t))
    # -- Define frequencies and spectrum for single sided approaches
    fp     = np.fft.rfftfreq(N, dt) # NOTE: when N is even, we have an extra frequency
    fMax   = fp[-1]                 # approx (1/dt)/2  Highest frequency
    df1    = 1/(dt*N)               # approx  1./tMax  Frequency resolution
    Nf     = len(fp)
    omega1 = 2*np.pi*fp
    # -- Define frequencies and spectrum
    if angularFrequency:
        Sf_or_So  = np.array([f_S(omi)     for omi in omega1])
        df_or_dom = 2*np.pi*df1
        f_or_om   = omega1
    else:
        Sf_or_So  = np.array([f_S(fi)      for fi in fp])
        df_or_dom = df1
        f_or_om   = fp
    if seed is not None:
        # Initialized random number generator
        np.random.seed(seed)
    
    if methods[0] in ['sumcos', 'boxmuller']:
        # --------------------------------------------------------------------------------
        # --- SUM OF COSINES / Box-Muller
        # --------------------------------------------------------------------------------
        X, A1, phi1 = _getX_sumcos(N, Sf_or_So, df_or_dom, phases=phases, amplitudesFactors=amplitudesFactors, randomAmplitudes=randomAmplitudes)

        if methods[1] == 'manual':
            A1 = A1
            x = np.zeros(t.shape)
            for ai,oi,ei in zip(A1,omega1,phi1):
                x += ai * np.cos(oi*t + ei)

        elif methods[1] == 'irfft' :
            x = np.fft.irfft(X, N)

        elif methods[1] == 'idft':
            from welib.tools.spectral import double_sided_DFT_real, check_DFT_real, IDFT

            if np.mod(N,2)==0:
                print('[WARN] IDFT introduces small error for now when N is even')
                # Keep me: what we use to do:
                #   fp = np.fft.rfftfreq(N, dt)[:-1] # NOTE the last frequency should not be used
            X = double_sided_DFT_real(X, N=N)
            check_DFT_real(X)
            x = IDFT(X, method=methods[2])
            x = np.real(x)

        else:
            NotImplementedError('Method {}'.format(method))

    elif methods[0]=='randn':
        # --------------------------------------------------------------------------------
        # --- RANDOM NORMAL complex numbers 
        # --------------------------------------------------------------------------------
        X_bas, A1, phi1 = _getX_basic(N, Sf_or_So, df_or_dom)

        a = np.random.normal(size=Nf)
        b = np.random.normal(size=Nf)
        #mag = np.sqrt(a**2+b**2)
        #phi = np.atan2(b,a)
        #X_rand = mag*np.exp(1j*phi) # TODO try that and compare to BoxMuller
        # First part of the complex vector
        X_rand = (a + b * 1j )/np.sqrt(2)
        X_rand[0] = 1
        X = X_bas * X_rand

        x = np.fft.irfft(X, N)

        #X = double_sided_DFT_real(X, N=N)
        #x=np.fft.ifft(X)
        #x = np.real(x)

    else:
        raise NotImplementedError('Method {}'.format(method))

    return t, x, f_or_om, Sf_or_So










def autocovariance_num(x, **kwargs):
    from welib.tools.signal_analysis import autoCorrCoeff
    sigma2 = np.var(x)
    rho_XX, tau = autoCorrCoeff(x, **kwargs)
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
    r""" 
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
    from welib.tools.signal_analysis import autoCorrCoeff
    rho_XX, tau = autoCorrCoeff(x, **kwargs)
    return rho_XX, tau

def autospectrum_fft(tau, kappa_XX, onesided=True, method='fft_wrap', verbose=False):
    r""" 
     Return one-sided (S_X) or double-sided (S_XX) autospectrum:
       - Using FFT
       - Based on discrete values of the auto-covariance kappa_XX.

      S_X(\omega) = 2 S_XX(\omega)

    """
    if onesided is False:
        raise NotImplementedError()
    if method=='fft_wrap':
        from welib.tools.spectral import fft_wrap
        f, S_X, Info = fft_wrap(tau, kappa_XX, output_type='amplitude', averaging='None', detrend=False)
        omega =2*np.pi * f
    else:
        if verbose:
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



def plot_pdf(y, method='histogram', n=100, 
        ax=None, label=None, sty='-',# plot options
        **kwargs):            # pdf options
    from welib.tools.stats import pdf
    import matplotlib.pyplot as plt
    xi,yi = pdf(y, method=method, n = n, **kwargs)
    if ax is None:
        fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
        fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
    ax.plot(xi, yi, sty, label=label)
#     ax.set_xlabel('')
    ax.set_ylabel('Probability density function')
#     ax.legend()   
    return ax
