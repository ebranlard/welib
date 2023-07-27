""" 



# TODO Look at
# http://statsmodels.sourceforge.net/devel/generated/statsmodels.tsa.stattools.acovf.html http://statsmodels.sourceforge.net/devel/generated/statsmodels.tsa.stattools.ccovf.html
"""



import numpy as np
from scipy.integrate import nquad

from welib.tools.signal_analysis import correlation
from welib.tools.spectral import fft_wrap




def sample_from_autospectrum(tMax, dt, f_S, angularFrequency=True, 
        method='ifft', fCutInput=None):
    """ 
    Generate a time series x, based on an single-sided autospectrum

    INPUTS:
     - tMax: maximum time of time series
     - dt  : time increament 
     - f_S: single-sided autospectrum function with interface:
         - S_X_om = f_S(omega)   if angularFrequency is True
         or
         - S_X_f  = f_S(f)       if angularFrequency is False
     - method: method used to generate the time series:
         - 'ifft': inverse FFT (fastest)
         - 'sum':  sum of cosines
     - fCutInput: maximum frequency. If None, set based on dt/2

    OUTPUTS:
     - x: time series

    """
    from numpy.random import uniform
    # -- Define frequencies and spectrum
    fCut = (1./dt)/2      # Highest frequency    [Hz]
    if fCutInput is not None:
        fMax = fCutInput
    else:
        fMax = fCut
 
    df = 1./tMax                  # Frequency resolution [Hz]
    f = np.arange(0,fMax+df/2,df) # array of frequencies [Hz]
    omega = 2*np.pi*f             # array of angular frequencies [rad/s]
    t = np.arange(0,tMax+dt/2,dt)

    # NOTE: this is evaluated here, because the ifft evaluates it elsewhere..
    if angularFrequency:
        S = np.array([f_S(omi) for omi in omega])
    else:
        S = np.array([f_S(fi) for fi in f])

    if method=='sum':
        # --- Method 1
        if angularFrequency:
            ap = np.sqrt(2*S*df*2*np.pi)      # amplitude  # TODO see if needed to adjust if angularFrequency
        else:
            ap = np.sqrt(2*S*df)      # amplitude  # TODO see if needed to adjust if angularFrequency
        eps = uniform(0,2*np.pi,len(f))  # random phases between 0 and 2pi
        u = np.zeros(t.shape)
        for ai,oi,ei in zip(ap,omega,eps):
            u += ai * np.cos(oi*t + ei)
    elif method=='ifft':
        # NOTE: The time series from the FFT will repeat itself at the middle. 
        # So we need twice the number of frequencies! 
        # Also, we use a double side spectrum
        if angularFrequency:
            # TODO Somehow the variance of first and last points are too high
            # So we increas the target time vector
            Ioffset=20
            #tMax2 = tMax                   # TODO sort it out
            #fMax=fMax*4
            tMax2 = tMax+dt*(2*Ioffset)                 # TODO sort it out
            fMax=fMax*4
        else:
            tMax2 = 4*tMax                     # TODO sort it out
        df2 = 1/(tMax2)                 
        f2 = np.arange(0,fMax+df2/2,df2)/2  # TODO sort it out
        # --- Method 2
        #  Value sets of two random normal distributed variables
        N = int(tMax2* fMax + 1)
        Nh= int(N/2) 
        #print('N',N)
        a = np.random.randn(Nh+1)
        b = np.random.randn(Nh+1)
        # First part of the complex vector
        xlhalf1 = a + b * 1j  # TODO TODO the b does nto seem to have an impact
        xlhalf1[0] = 0 # zero mean

        omega = 2.*np.pi*f2[:Nh+1] 

        #print('[WARN] Need to sort out ifft method for sample generation')
        if angularFrequency:
            Sw = np.array([f_S(omi) for omi in omega]) # TODO WHAT SHOULD BE DONE
            #Sw = np.array([f_S(omi*2*np.pi) for omi in omega])
        else:
            #Sw = np.array([f_S(omi/(2*np.pi)) for omi in omega]) # TODO WHAT SHOULD BE DONE!
            Sw = np.array([f_S(omi) for omi in omega])
        if angularFrequency:
            #sig = np.sqrt(tMax2/(2*np.pi)* Sw) # see a= np.sqrt(2*S*df)
            sig = np.sqrt(tMax2*Sw) # see a= np.sqrt(2*S*df)
            #sig = Sw #np.sqrt(2*df2*Sw) # see a= np.sqrt(2*S*df)
            # TODO somehow the variance of the first and last points are too high
            #sig[0] =sig[0]/20
        else:
            sig = np.sqrt(tMax2/(2*np.pi)* Sw) # see a= np.sqrt(2*S*df)
#         sig = np.sqrt(tMax2/(1)* Sw) # see a= np.sqrt(2*S*df)
#         sig = np.sqrt(2*Sw*df2) 
        xlhalf1 *= sig

        xlhalf2 = np.flipud(xlhalf1[1:])

        # Total complex vector
        xl=np.concatenate((xlhalf1,xlhalf2))



        # Fourier inverse
        u=np.fft.ifft(xl)

#         import matplotlib.pyplot as plt
#         fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
#         fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
#         ax.plot(np.real(xl )  , label='')
# #         ax.plot(np.real(sig)   , label='sig')
# #         ax.plot(np.real(a  ) , label='a')
# #         ax.plot(np.real(Sw ) , label='Sw')
#         ax.set_xlabel('')
#         ax.set_ylabel('')
#         ax.legend()
#         fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
#         fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
#         ax.plot(np.real(u )  , label='re')
#         ax.plot(np.imag(u )  , label='im')
#         ax.plot(np.abs(u )  , label='ab')
#         ax.plot(np.conjugate(u )*u,'--'  , label='z z*')
# #         ax.plot(np.real(sig)   , label='sig')
# #         ax.plot(np.real(a  ) , label='a')
# #         ax.plot(np.real(Sw ) , label='Sw')
#         ax.set_xlabel('')
#         ax.set_ylabel('')
#         ax.legend()
#         plt.show()

        # Remove zero imaginairy part and normalize
        if angularFrequency:
            Nh= len(t)-1 # TODO check consistency with Nh above
            # NOTE: imaginary or real works equivalently
            u=np.real(u)*fMax*np.sqrt(np.pi/2)
            #u=np.imag(u)*fMax*np.sqrt(np.pi/2)
            #u=u[20:Nh+21] 
            #u=u[20:Nh+21] # TODO somehow the variance of the first and last points are too high

#             import matplotlib.pyplot as plt
#             fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
#             fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
#             I=np.arange(len(u))
#             I2=I[Ioffset:Nh+Ioffset+1]
#             ax.plot(I,u      , label='')
#             ax.plot(I2,u[I2] , label='')
#             ax.set_xlabel('')
#             ax.set_ylabel('')
#             plt.show()

            # TODO somehow the variance of the first and last points are too high
            #u=u[0:Nh+1]
            #u=u[3:Nh+4]
            u=u[Ioffset:Nh+Ioffset+1]
        else:
            u=np.sqrt(2)*np.pi*fMax*np.real(u)
            #u=fMax*np.real(u)/2
            u=u[0:Nh+1]

    u=u-np.mean(u)
    return t, u, f, S















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
