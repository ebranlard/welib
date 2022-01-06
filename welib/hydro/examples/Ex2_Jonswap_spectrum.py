"""
Plot the JONSWAP spectrum for a given sea state
"""
import numpy as np
from welib.hydro.spectra import jonswap
import matplotlib.pyplot as plt

# --- Parameters
t  = np.arange(0,3600.1,1)    # time vector  [s]
dt = t[1]-t[0]                # timestep [s] 
Hs = 8.1                      # Significant wave height [m]
Tp = 12.7                     # Peak period [s]

# --- Derived parameters
df   = 1./np.max(t) # Step size for frequency
fMax = (1./dt)/2    # Highest frequency 
freq = np.arange(df, fMax+df/2, df)

# --- Spectrum and amplitude
S = jonswap(freq, Hs, Tp) # Spectral density [m^2.s]
ap = np.sqrt(2*S*df)      # Wave amplitude [m]

# Find location of maximum energy
iMax = np.argmax(S)

# --- Plots
fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
ax.plot(freq, S)
ax.plot(freq[iMax], S[iMax], 'ko')
ax.set_xlabel('Frequency [Hz]')
ax.set_ylabel(r'Spectral density [m^2 s]')
ax.set_title('Hydro - Jonswap spectrum')
ax.tick_params(direction='in')


if __name__=="__main__":
    plt.show()
if __name__=="__test__":
    np.testing.assert_almost_equal(S[iMax], 113.8770176)
if __name__=="__export__":
    from welib.tools.repo import export_figs_callback
    export_figs_callback(__file__)


