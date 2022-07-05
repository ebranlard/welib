import numpy as np
import matplotlib.pyplot as plt
from welib.tools.signal_analysis import *
from numpy.random import seed; seed(0)

# --- Parameters
dt    = 1
n     = 10000
coeff = 0.95 # 1:full corr, 0: no-corr
nMax  = 180
# --- Create a correlated time series
tvec = np.arange(0,n)*dt
ts   = correlated_signal(coeff, n)
# --- Compute correlation coefficient
R, tau = correlation(ts, nMax=nMax, dt=dt)

# --- Plot
fig,axes = plt.subplots(2, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
ax=axes[0]
# Plot time series
ax.plot(tvec,ts)
ax.set_xlabel('t [s]')
ax.set_ylabel('u [m/s]')
ax.tick_params(direction='in')
# Plot correlation
ax=axes[1]
ax.plot(tau,  R              ,'-o', label='Computed')
ax.plot(tau, coeff**(tau/dt) ,'--' ,label=r'Theoretical -  c$^{\tau/dt}$') # analytical coeff^n trend
ax.set_xlabel(r'$\tau$ [s]')
ax.set_ylabel(r'$R(\tau)$ [-]')
ax.set_title('Signal - Correlation coefficient')
ax.legend()


if __name__=='__main__':
    plt.show()
if __name__ == '__test__':
    pass
if __name__=="__export__":
    from welib.tools.repo import export_figs_callback
    export_figs_callback(__file__)
