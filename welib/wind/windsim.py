import numpy as np
from numpy.random import uniform
from welib.stoch.utils import sample_from_autospectrum

def pointTS(f, S, **kwargs):
    """ 
    Generate a wind time series at a single point based on a spectrum
    INPUTS:
     - f: array of frequencies [Hz]
     - S(f): spectrum at f 
    OUTPUTS:
     - u: wind speed
    """
    if hasattr(S,'__len__'):
        f_S = lambda f_: np.interp(f_, f, S)
    else:
        f_S = S

    t, u, f, S = sample_from_autospectrum(tMax=tMax, dt=dt, f_S=f_S, angularFrequency=False, **kwargs)

    return t, u, f, S


def pointTSKaimalWrap(time, U0=8, sigma=1, L=300, **kwargs):
    tMax = time[-1]
    dt = (time[-1]-time[0])/(len(time)-1)
    t, u, f, Sf = pointTSKaimal(tMax, dt, U0=8, sigma=1, L=300, angularFrequency=False, **kwargs)
    return u

def pointTSKaimal(tMax, dt, U0=8, sigma=1, L=300, angularFrequency=False, **kwargs):
    """
    Generate a wind time series at a single point based on the Kaimal spectrum.

    INPUTS:
    ( for the generation of a sample, see welib.stoch.utils.sample_from_autospectrum:)
     - tMax   : maximum time for time vector 
               NOTE: different conventions exists for tMax. here we use:
               tMax = time[-1] = (N-1) dt  

     - dt     : time step
     - method : 
    (for the Kaimal spectrum, see welib.wind.spectra)
     - U0    : mean wind speed
     - sigma : standard deviation
     - L     : length scale
     OUTPUTS:
      - t: time vector
      - u: velocity time series
      - f : frequencies in [Hz]
      - Sf : spectrum as function of frequency
    """
    from welib.wind.spectra import kaimal

    f_S = lambda f_or_om: kaimal(f_or_om, U0, sigma, L, angularFrequency=angularFrequency)
    t, u, f_or_omega, Sf_or_So = sample_from_autospectrum(tMax=tMax, dt=dt, f_S=f_S, angularFrequency=angularFrequency, **kwargs)

    if angularFrequency:
        f   = f_or_om/(2*np.pi)
        Sf  = Sf_or_So*2*np.pi
    else:
        f   = f_or_omega
        Sf  = Sf_or_So

    # Add the influence of the mean
    u += U0
    tSuperMax = tMax + dt
    Sf[0] = U0**2 * tSuperMax
    return t, u, f, Sf

if __name__ == '__main__':
    import matplotlib.pyplot as plt
    from welib.tools.spectral import fft_wrap, DFT, IDFT
    from welib.tools.tictoc import Timer
    from welib.wind.spectra import kaimal

    U0    = 8      # Wind speed [m/s], for Kaimal spectrum
    I     = 0.14   # Turbulence intensity [-], for Kaimal spectrum
    L     = 340.2  # Length scale [m], for Kaimal spectrum
    tMax  = 800    # Maximum time for time series [s]
    dt    = 0.01   # Time step [s]
    tMax  = 200+0*dt    # Maximum time for time series [s]
    angularFrequency = False
    angularFrequency = True
    randomAmplitudes = False # Only for Box Muller method for now
    #randomAmplitudes = True # Only for Box Muller method for now

    us=[]
    ts=[]
    Methods = []
    Methods += ['sumcos-manual']
    Methods += ['sumcos-idft-ifft']
    Methods += ['sumcos-irfft']
#     Methods += ['boxmuller-irfft']

    fig,axes = plt.subplots(2, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
    fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.30, wspace=0.20)

    vSty=['-','--',':','-.']

    for iMethod, method in enumerate(Methods): #,'sum']:
        # NOTE: we use the same seed number for all methods for reproductibility
        with Timer(method):
            t, u, freq, Sf = pointTSKaimal(tMax, dt, U0, U0*I, L, method=method, randomAmplitudes=randomAmplitudes, seed=11)

        # --- Compute FFT of wind speed
        f_fft, S_fft, Info   = fft_wrap(t, u, output_type='PSD', averaging='none')

        # --- Plots
        ax=axes[0]
        ax.plot(t, u, vSty[iMethod], label='Generated '+method)
        ax.tick_params(direction='in')
        ax.autoscale(enable=True, axis='both', tight=True)
        ax.set_xlabel('Time [s]')
        ax.set_ylabel(r'Wind speed [m/s]')
        ax.set_title('Wind - wind generation at point')

        ax=axes[1]
        ax.plot(f_fft , S_fft, vSty[iMethod], label='Generated '+method)
        if iMethod==len(Methods)-1:
            ax.plot(freq  , Sf    , 'k--', label='Kaimal')
        ax.set_yscale('log')
        ax.set_xscale('log')
        ax.legend()
        ax.set_xlabel('Frequency [Hz]')
        ax.set_ylabel(r'Spectral density [m$^2$ s]')
        ax.tick_params(direction='in')
        ax.autoscale(enable=True, axis='both', tight=True)

        print('Integration Kaimal:', np.trapezoid(Sf, freq))
        print('Integration Gener1:', np.trapezoid(S_fft, f_fft))
        Sf[0]    = 0
        S_fft[0] = 0
        print('Integration Kaimal:', np.trapezoid(Sf, freq))
        print('Integration Gener1:', np.trapezoid(S_fft, f_fft))
        print('Sigma2            :', np.var(u))
    
    plt.show()
