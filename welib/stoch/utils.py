""" 



# TODO Look at
# http://statsmodels.sourceforge.net/devel/generated/statsmodels.tsa.stattools.acovf.html http://statsmodels.sourceforge.net/devel/generated/statsmodels.tsa.stattools.ccovf.html
"""



import numpy as np
from scipy.integrate import nquad

from welib.tools.signal_analysis import correlation
from welib.tools.spectral import fft_wrap

def autocovariance_num(x, **kwargs):
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


def autocorrcoeff_num(x, **kwargs):
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
    S_X = np.zeros(omega.shape) 
    # for i, omi in enumerate(om_th):
    #     f_S_X_integrand = lambda tau : 1/np.pi * np.exp(-1j * omi * tau) * f_kappa(tau)
    #     S_X_i[i], _ = nquad(f_S_X_integrand, [[0, tau_max]] )
    for i, omi in enumerate(omega):
        f_S_X_integrand = lambda tau : 2/np.pi * np.cos(omi * tau) * f_kappa(tau)
        S_X[i], _, _  = nquad(f_S_X_integrand, [[0, tau_max]] , full_output=1)
    return S_X

def autocovariance_int(tau, f_S_X, omega_max):
    """
    Return autocovariance function
    """
    K_XX = np.zeros(tau.shape) 
    for i, taui in enumerate(tau):
        integrand = lambda om : np.cos(om * taui) * f_S_X(om)
        K_XX[i], _, _ = nquad(integrand, [[0, omega_max]] , full_output=1)
    return K_XX
