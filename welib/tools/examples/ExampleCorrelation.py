import numpy as np
import matplotlib.pyplot as plt
from pybra.signal import *

if __name__=='__main__':
    # --- Parameters
    dt    = 1
    n     = 10000
    coeff = 0.95 # 1:full corr, 0: no-corr
    nMax  = 180
    # --- Create a correlated time series
    tvec = np.arange(0,n)*dt
    ts   = correlated_signal(coeff, n)
    # --- Compute correlation coefficient
    R, tau = correlation(x, nMax=nMax, dt=dt)

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
    ax.plot(tau,  R              ,'b-o', label='computed')
    ax.plot(tau, coeff**(tau/dt) , 'r--' ,label='coeff^{tau/dt}') # analytical coeff^n trend
    ax.set_xlabel(r'$\tau$ [s]')
    ax.set_ylabel(r'$R(\tau)$ [-]')
    ax.legend()
    plt.show()






