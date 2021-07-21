import numpy as np
from numpy.random import uniform

def pointTS(f, S):
    """ 
    Generate a wind time series at a single point based on a spectrum
    INPUTS:
     - f: array of frequencies [Hz]
     - S: spectrum at f
    OUTPUTS:
     - u: wind speed
    """
    pass


def pointTSKaimal(tMax, dt, U0, sigma, L, method='ifft', fCutInput=None):
    """
    Generate a wind time series at a single point based on a s
    
    """
    from welib.wind.spectra import kaimal

    # -- Define frequencies and spectrum
    fCut = (1./dt)/2      # Highest frequency 
    df = 1./tMax          # Frequency resolution
    f = np.arange(0,fCut+df/2,df)
    t = np.arange(0,tMax+dt/2,dt)
    S = kaimal(f, U0, sigma, L)

    if method=='sum':
        # --- Method 1
        ap = np.sqrt(2*S*df)      # amplitude 
        omega = 2*np.pi*f
        eps = uniform(0,2*np.pi,len(f))  # random phases between 0 and 2pi
        u = np.zeros(t.shape)
        for ai,oi,ei in zip(ap,omega,eps):
            u += ai * np.cos(oi*t + ei)
    elif method=='ifft':
        tMax2 = 4*tMax                     # TODO sort it out
        if fCutInput is not None:
            fMax = fCutInput
        else:
            fMax = fCut
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
        xlhalf1 = a + b * 1j
        xlhalf1[0] = 0 # zero mean

        omega = 2.*np.pi*f2[:Nh+1] 
        Sw = kaimal(omega, U0, sigma, L) # S(omega)
        sig = np.sqrt(tMax2/(2*np.pi)* Sw) # see a= np.sqrt(2*S*df)
        xlhalf1 *= sig
        #print(N)
        #print(sig[0:5], sig[-1])

        xlhalf2 = np.flipud(xlhalf1[1:])

        # Total complex vector
        xl=np.concatenate((xlhalf1,xlhalf2))

        # Fourier inverse
        u=np.fft.ifft(xl)

        # Remove zero imaginairy part and normalize
        u=np.sqrt(2)*np.pi*fMax*np.real(u)
        u=u[0:Nh+1]
# ###########################################################
#  # First approach
#  ###########################################################
#  generateWSFromSpectrum=function(S,N, fs ,U){
#  T=N/fs ;
#  X=numeric(N)
#  for ( l in 1:(N/2) ){
#  X[ l ]=rnorm(1 ,mean = 0 , sd = sqrt(T/(2∗pi )∗S[ l ]) ) ;
#  }
#  x= Re(( fft (c(T∗U/(2∗pi ) ,X) ,inverse=T)/N ) )∗2∗pi∗fs
#  print(c(mean(x) ,sd(x) )
#  return(x)
#  }
# 
#  ###########################################################
#  # Second approach
#  ###########################################################
#  generateWSFromSpectrum2=function(S,N, fs ,U){
#  T=N/fs ;
#  I=complex( real=0,imaginary=1)
#  a=numeric(N/2)
#  b=numeric(N/2)
#  for ( l in 1:(N/2) ){
#  a [ l ]=rnorm(1 ,mean = 0 , sd = sqrt(T/(2∗pi )∗S[ l ])/sqrt (2) ) ;
#  b[ l ]=rnorm(1 ,mean = 0 , sd = sqrt(T/(2∗pi )∗S[ l ])/sqrt (2) ) ;
#  }
#  X=numeric(N)
#  X[2:(N/2)]=(a+I∗b) [1:(N/2−1) ]/2
#  X[((N/2+2) :N)]=(a−I∗b) [(N/2−1) :1]/2
#  X[1]=T∗U/(2∗pi ) ;
#  X[N/2+1]=a [N/2]
#  x= Re(( fft (X, inverse=T)/N ) )∗2∗pi∗fs
#  print(mean(x) )
#  print(sd(x) )
#  return(x)
#  }

    # Setting mean to U0
    u=u-np.mean(u) + U0
    #print('u', u.shape)
    #print('t', t.shape)
    #print('f', f.shape)
    #print('S', S.shape)

    return t, u, f, S

if __name__ == '__main__':

    from welib.tools.spectral import fft_wrap
    from welib.tools.tictoc import Timer

    seed(11)

    U0    = 8      # Wind speed [m/s], for Kaimal spectrum
    I     = 0.14   # Turbulence intensity [-], for Kaimal spectrum
    L     = 340.2  # Length scale [m], for Kaimal spectrum
    tMax  = 800    # Maximum time for time series [s]
    dt    = 0.01   # Time step [s]

    for method in ['sum','ifft']:
        with Timer(method):
            t, u, freq, S =pointTSKaimal(tMax, dt, U0, U0*I, L, method=method)

        # --- Compute FFT of wind speed
        f_fft, S_fft, Info = fft_wrap(t, u, output_type='PSD', averaging='none')

        # --- Plots
        import matplotlib.pyplot as plt
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
        ax.plot(f_fft, S_fft, '-', label='Generated')
        ax.plot(freq, S     , 'k', label='Kaimal')
        ax.set_yscale('log')
        ax.set_xscale('log')
        ax.legend()
        ax.set_xlabel('Frequency [Hz]')
        ax.set_ylabel(r'Spectral density [m$^2$ s]')
        ax.tick_params(direction='in')
        ax.autoscale(enable=True, axis='both', tight=True)
    
    plt.show()
