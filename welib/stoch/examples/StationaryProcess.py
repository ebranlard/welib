import numpy as np
import matplotlib.pyplot as plt

from welib.stoch.distribution import*
from welib.stoch.stationary_process import *


# --------------------------------------------------------------------------------}
# --- Rational auto-spectrum 
# --------------------------------------------------------------------------------{
def rational(nDiscr=6000, verbose=True, nSamples=200):
    # Distribution parameters
    m = 1
    S0 = 1/np.pi
    omega0 = 2*np.pi
    zeta = 0.1 # NOTE: as zeta increase the spectral peak spreads, S_X(0) increases and the accuracy decreases somehow
    # Derived parameters
    sigma2= np.pi*S0/(2*zeta*omega0**3*m**2)
    omegad = np.sqrt(1-zeta**2)*omega0
    # Process integration parameters
    omega_max= 10*omega0
    tau_max  = 10*2*np.pi/(omega0) * 1/zeta # the smaller zeta, the longer we need
    time_max = 4*tau_max
    # --- Initialize Process
    proc = SampledStochasticProcess(name='rational', omega_max=omega_max, tau_max=tau_max, nDiscr=nDiscr, verbose=verbose, time_max=time_max)
    # Theory
    proc._f_S_X_th = lambda omega: 2*S0/m**2  / (omega**4 + 4*omega**2*omega0**2*zeta**2 - 2*omega**2*omega0**2+omega0**4)
    #proc._f_S_X_th = lambda omega: np.real(2*S0/m**2  /( ( -omega**2 + 2*zeta*omega0*omega*1j +omega0**2) * ( -omega**2 - 2*zeta*omega0*omega*1j +omega0**2) ))
    proc._f_kappa_XX_th = lambda tau:  sigma2 * np.exp(-zeta*omega0*np.abs(tau))*(np.cos(omegad*tau)+ zeta/(1-zeta**2)*np.sin(omegad*tau))
    proc._var_th = sigma2
    proc._mu_th = 0


    # --- Samples
    proc.generate_samples_from_autospectrum(nSamples=nSamples, method=gen_method)
    #proc.generate_samples_from_autospectrum(nSamples=100, method='sum')
    proc.samples_stats(stats=stats) #nTau = None)
    proc.plot_samples()
    proc.plot_var()

    # --- Computation from analytical variance
    S_X_i = proc.autospectrum_int  (method='fft', nDiscr=50*nDiscr) #, tau_max=tau_max) # still need work
    K_X_i = proc.autocovariance_int(method='fft') #, nDiscr=5000, omega_max=omega_max)
    S_X_i = proc.autospectrum_int  (method='quad')
    K_X_i = proc.autocovariance_int(method='quad')

    # --- Autocovariance
    axK = proc.plot_autocovariance()
    # --- Autospectrum
    axS = proc.plot_autospectrum()
    #axS.plot([omega0,omega0],[0, 2*S0/m**2 / (4*omega0**4*zeta**2)], 'k--')
    axS.axvline(omega0, ls=':', c='k') # Remember, max is not exactly at omega0
    #axS.set_xlim([0,2*omega0])

# --------------------------------------------------------------------------------}
# --- Harmonic process 
# --------------------------------------------------------------------------------{
def harmo(nDiscr=500, verbose=True, nSamples=500):

    # Distribution parameters
    s      = 4
    omega0 = 2*np.pi

    # Derived parameters
    T = 1/(omega0/(2*np.pi))
    time=np.linspace(0,40*T, 1000)
    dt = (time[-1]-time[0])/(len(time)-1)
    tau_max   = 10*T     # TODO
    omega_max = 3*omega0

    # Process
    proc = HarmonicProcess(s=s, omega0=omega0, verbose=verbose, nDiscr=nDiscr, tau_max=tau_max, omega_max=omega_max)

    # --- Samples
    # NOTE: no autospectrum
    xi = proc.generate_samples(nSamples=nSamples, time=time)
    proc.samples_stats(stats=stats) #nTau = None)
    proc.plot_samples()
    proc.plot_var()
    # proc.plot_mean()

    # --- Computations
    proc.autospectrum_int(method='quad')
    proc.autospectrum_int(method='fft', nDiscr=20*nDiscr, tau_max=tau_max*3)
    #proc.autocovariance_int(method='fft')

    proc.plot_autocovariance()
    # proc.plot_autocorrcoeff()

    #
    ax = proc.plot_autospectrum()
    ax.plot([omega0,omega0], [0,proc._var_th] ,'k--'   , label='Theory', lw=2)
    ax.set_xlabel(r'$\omega$ [rad/s]')
    ax.set_ylabel(r'$S_X(\omega)$')
    ax.set_xlim([0, 2*omega0])
    ax.legend()

    
def expo(nDiscr=2000, verbose=True, nSamples=200):
    # Distribution parameters
    lambd = 10
    sigma = 3
    # Process integration parameters
    omega_max = - 1/lambd* np.log( 1e-3/ (lambd*sigma**2))
    tau_max   = np.sqrt( lambd**2 * (sigma**2/1e-3-1))
    time_max  = 2*tau_max

    # --- Initialize Process
    proc = ExponentialProcess(lambd, sigma, omega_max=omega_max, tau_max=tau_max, nDiscr=nDiscr, verbose=verbose, time_max=time_max)

    # --- Samples
    proc.generate_samples_from_autospectrum(nSamples=nSamples, method=gen_method)
    #proc.generate_samples_from_autospectrum(nSamples=100, method='sum')
    proc.samples_stats(stats=stats) #nTau = None)
    proc.plot_samples()
    proc.plot_var()
    proc.plot_mean()

    # --- Computations
    S_X_i = proc.autospectrum_int(method='fft', nDiscr=30*nDiscr, tau_max=tau_max) # still need work
    K_X_i = proc.autocovariance_int(method='fft', nDiscr=5*nDiscr, omega_max=omega_max)
    S_X_i = proc.autospectrum_int(method='quad')
    K_X_i = proc.autocovariance_int(method='quad')

    # --- Autocovariance
    axK = proc.plot_autocovariance()
    # --- Autospectrum
    axS = proc.plot_autospectrum()
#     axS.set_yscale('log')
#     axS.set_xlim([0,5])
#     axS.set_ylim([1e-4,125])


# --------------------------------------------------------------------------------}
# --- BAND-LIMITED WHITE NOISE 
# --------------------------------------------------------------------------------{
def bandedwhite(nDiscr=2000, verbose=True, nSamples=200):
    # Distribution parameters
    omega0 = 2*np.pi
    zeta = 0.5
    sigma = 1 
    # Derived parameters
    B    = 2*zeta*omega0
    omega_1 = omega0-B/2
    omega_2 = omega0+B/2
    T = 1/(omega0/(2*np.pi))
    tau_max = 20*T
    time_max = 2*tau_max

    # --- Initialize Process
    proc = BandedWhiteNoiseProcess(omega0, zeta, sigma, omega_max=2*omega_2, tau_max=tau_max, nDiscr=nDiscr, verbose=verbose, time_max=time_max)
    print(proc)

    # --- Samples
    proc.generate_samples_from_autospectrum(nSamples=nSamples, method=gen_method)
    #proc.generate_samples_from_autospectrum(nSamples=200, method='sum')
    proc.samples_stats(stats=stats) #nTau = None)
    proc.plot_samples()
    proc.plot_var()
    # # proc.plot_mean(

    # --- Moments
    moments = proc.autospectrum_moments(method='num')
    moments2= proc.autospectrum_moments(method='quad')
    if verbose:
        print('Moments: ', [sigma**2 /( B*(j+1))*( omega_2**(j+1) - omega_1**(j+1)) for j in range(4) ])
        print('Moments: ', moments)
        print('Moments: ', moments2)

    # --- Computations from provided functions
    S_X_i = proc.autospectrum_int(method='fft', nDiscr=10*nDiscr, tau_max=30)
    K_X_i = proc.autocovariance_int(method='fft', omega_max=20*omega0, nDiscr=5*nDiscr) # needs high res
    S_X_i = proc.autospectrum_int(method='quad')
    K_X_i = proc.autocovariance_int(method='quad')
    # 
#     proc.plot_autocorrcoeff()
    # --- Autocovariance
    axK = proc.plot_autocovariance()
    # --- Autospectrum
    axS = proc.plot_autospectrum()

#     tau = np.linspace(0, 2*tau_max, 1000000)
#     kappa_XX = np.array([proc._f_kappa_XX_th(t) for t in tau])
#     k = np.fft.rfft(kappa_XX)
#     tau_max = tau[-1]
#     n      = int(len(tau)/2)
#     dtau   = (tau[-1]-tau[0])/(n-1)
#     Fs     = 1/dtau       
#     om    = np.arange(n)*Fs/(n) * 2 *np.pi
#     S_X    = k[:n]/n * tau_max / (np.pi)
# #     print('>> k', np.real(k[:5]))
# #     print('>> k', np.abs(k[:5]))
# #     print('>> S', S_X[:5])
# #     print('>>> tau+max',tau_max)
# #     print('>>> len(k)',len(k))
# #     print('>>> len(S_X)',len(S_X))
# #     print('>>> len(om  )',len(om))
#     axS.plot(om, np.real(S_X))
#     axS.set_xlim([0,10])

def bandedwhiteSym(nDiscr=2000, verbose=True, nSamples=200):
    # Distribution parameters
    sigma = 1 
    # Derived parameters
    B = 2*np.pi*10  # [rad/s]
    T = 1/(B/(2*np.pi))
    tau_max = 30*T        # TODO
    time_max = 2*tau_max  # TODO
    omega_max = B

    # --- Initialize Process
    proc = CutOffWhiteNoiseProcess(B, sigma, omega_max=omega_max, tau_max=tau_max, nDiscr=nDiscr, verbose=verbose, time_max=time_max)
    print(proc)

    # --- Samples
    proc.generate_samples_from_autospectrum(nSamples=nSamples, method=gen_method)
    #proc.generate_samples_from_autospectrum(nSamples=200, method='sum')
    proc.samples_stats(stats=stats) #nTau = None)
    proc.plot_samples()
    proc.plot_var()
    # # proc.plot_mean(

    # --- Moments
    moments = proc.autospectrum_moments(method='num')
    moments2= proc.autospectrum_moments(method='quad')
    if verbose:
        print('Moments: ', [sigma**2 /( B*(j+1))*( (B/2)**(j+1) - (-B/2)**(j+1)) for j in range(4) ])
        print('Moments: ', moments)
        print('Moments: ', moments2)

    # --- Computations from provided functions
    S_X_i = proc.autospectrum_int(method='fft', nDiscr=10*nDiscr, tau_max=30)
    K_X_i = proc.autocovariance_int(method='fft', omega_max=10*omega_max, nDiscr=5*nDiscr) # needs high res
    S_X_i = proc.autospectrum_int(method='quad')
    K_X_i = proc.autocovariance_int(method='quad')
    # 
#     proc.plot_autocorrcoeff()
    # --- Autocovariance
    axK = proc.plot_autocovariance()
    # --- Autospectrum
    axS = proc.plot_autospectrum()


if __name__ == '__main__':
    from numpy.random import seed
    seed(12)
    stats=[]
    stats+=['avg_spectra']
#     stats+=['correlation']
    gen_method ='sumcos-irfft'

    bandedwhiteSym()
#     bandedwhite()
#     harmo()
#     expo()
#     rational()
    plt.show()

if __name__ == '__test__':
    stats=[]
    stats+=['avg_spectra']
    stats+=['correlation']
    gen_method ='sumcos-irfft'

    bandedwhite(nDiscr = 100, verbose = False, nSamples=5)
    harmo      (nDiscr = 100, verbose = False, nSamples=5)
    expo       (nDiscr = 100, verbose = False, nSamples=5)
    rational   (nDiscr = 500, verbose = False, nSamples=5)
