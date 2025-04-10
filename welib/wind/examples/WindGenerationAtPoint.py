""" 
Generate random wind time series at a given point based on the Kaimal spectrum
Compute the FFT to verify that the spectrum is retrieved.

"""
from welib.wind.windsim import *

from welib.tools.spectral import fft_wrap
import matplotlib.pyplot as plt
try:
    from numpy import trapezoid
except:
    from numpy import trapz as trapezoid

def generateTSPlot(method='sumcos-irfft', seed=12, randomAmplitudes=True):
    U0    = 8      # Wind speed [m/s], for Kaimal spectrum
    I     = 0.14   # Turbulence intensity [-], for Kaimal spectrum
    L     = 340.2  # Length scale [m], for Kaimal spectrum
    tMax  = 600    # Maximum time for time series [s]
    dt    = 0.01   # Time step [s]

    # --- Generate time series based on Kaimal spectrum
    t, u, freq, S =pointTSKaimal(tMax, dt, U0, U0*I, L, method=method, seed=seed, randomAmplitudes=randomAmplitudes)

    # --- Compute FFT of wind speed
    f_fft, S_fft, Info = fft_wrap(t, u, output_type='PSD', averaging='none')

    # --- Plots
    fig,axes = plt.subplots(2, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
    fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.30, wspace=0.20)

    ax=axes[0]
    ax.plot(t, u)
    ax.tick_params(direction='in')
    ax.autoscale(enable=True, axis='both', tight=True)
    ax.set_xlabel('Time [s]')
    ax.set_ylabel(r'Wind speed [m/s]')
    ax.set_title('Wind - wind generation at point')

    ax=axes[1]
    ax.plot(f_fft, S_fft, '-'  , label='Generated')
    ax.plot(freq, S     , 'k--', label='Kaimal')
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.legend()
    ax.set_xlabel('Frequency [Hz]')
    ax.set_ylabel(r'Spectral density [m$^2$ s]')
    ax.tick_params(direction='in')
    ax.autoscale(enable=True, axis='both', tight=True)


    return t, u, freq, S, f_fft, S_fft



if __name__=='__main__':
    from welib.tools.tictoc import Timer
    methods = []
    methods += ['sumcos-irfft']
    methods += ['sumcos-manual']
    for method in methods:
        with Timer(method):
            t, u, freq, S, f_fft, S_fft = generateTSPlot(method=method, randomAmplitudes=False)
        S[0]     = 0
        S_fft[0] = 0
        print('Integration Kaimal:', trapezoid(S, freq))
        print('Integration       :', trapezoid(S_fft, f_fft))
        print('Sigma^2           :', np.var(u))
    plt.show()

if __name__=='__main__':
    from welib.tools.tictoc import Timer
    methods = []
    methods += ['sumcos-irfft']
    methods += ['sumcos-manual']
    for method in methods:
        t, u, freq, S, f_fft, S_fft = generateTSPlot(method=method, randomAmplitudes=False)

if __name__=='__export__':
    generateTSPlot(method='sumcos-irfft', randomAmplitudes=False)

    from welib.tools.repo import export_figs_callback
    export_figs_callback(__file__)

