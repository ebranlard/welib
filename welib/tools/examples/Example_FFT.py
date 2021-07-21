import numpy as np
import matplotlib.pyplot as plt
from welib.tools.spectral import fft_wrap

Tmax =5
dt  = 0.01
t  = np.arange(0,Tmax,dt)
f1 = 1; 
f2 = 5; 
A1 = 5; 
A2 = 10; 
x  = A1*np.sin(2*np.pi*f1*t) + A2*np.sin(2*np.pi*f2*t) + np.random.randn(len(t))*2.0



# --- Amplitude
f1, S1, Info  = fft_wrap(t,x,output_type='amplitude',averaging='Welch', nExp=7) # <<< TODO Play with nExp
f2, S2, Info  = fft_wrap(t,x,output_type='amplitude',averaging='none')

fig,ax = plt.subplots(2, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.30, wspace=0.20)
ax[0].plot(t,x)
ax[1].plot(f1,S1,'r',label ='Welch')
ax[1].plot(f2,S2,'k--',label ='None')
ax[0].set_xlabel('Time [s]')
ax[1].set_xlabel('Frequency [Hz]')
ax[0].set_ylabel('Signal')
ax[1].set_ylabel('FFT Amplitude')
# ax[1].set_yscale("log", nonposy='clip')
ax[1].legend()
ax[0].set_title('Signal - FFT')


# --- PSD, log y
f1, S1, Info  = fft_wrap(t,x,output_type='PSD',averaging='Welch', nExp=9) # <<< TODO Play with nExp
f2, S2, Info  = fft_wrap(t,x,output_type='PSD',averaging='none')

fig2,ax = plt.subplots(2, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
fig2.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.30, wspace=0.20)
ax[0].plot(t,x)
ax[1].plot(f1,S1,'r',label ='Welch')
ax[1].plot(f2,S2,'k--',label ='None')
#ax[1].set_yscale("log", nonposy='clip')
ax[0].set_xlabel('Time [s]')
ax[1].set_xlabel('Frequency [Hz]')
ax[0].set_ylabel('Signal')
ax[1].set_ylabel('PSD Amplitude')
ax[1].legend()
ax[0].set_title('Signal - PSD')


if __name__=="__main__":
    plt.show()
if __name__=="__test__":
    pass
if __name__=="__export__":
    plt.close(fig2)
    from welib.tools.repo import export_figs_callback
    export_figs_callback(__file__)


