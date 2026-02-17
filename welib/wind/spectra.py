import numpy as np
import pandas as pd

def kaimal_omega(omega, V0, sigma, L):
    r""" 
    Single-sided Kaimal spectrum S_X(\omega)

    INPUTS:
     - omega: frequency, in [rad/s]
     - V0: reference/hub velocity velocity
     - sigma: wind standard deviation, turbulence intensity I= sigma/V0. 
     - L: turbulence length scale. 
    OUTPUTS:
      - S_X: Single-sided Kaimal spectrum as function of omega
    """
    omega = np.asarray(omega)
    S = 1/(2*np.pi)*(4*sigma**2 * L / V0) / (1 + 6 * omega/(2*np.pi) * L/V0) ** (5./3) 
    S[omega==0] = 0 # V0**2 * (tMax+dt), but we will add the mean at the end
    return S

def kaimal_f(f, V0, sigma, L):
    r""" 
    Single-sided Kaimal spectrum S_X(f)

    INPUTS:
     - f: frequency, in [Hz]
     - V0: reference/hub velocity velocity
     - sigma: wind standard deviation, turbulence intensity I= sigma/V0. 
     - L: turbulence length scale. 
    OUTPUTS:
      - S_X: Single-sided Kaimal spectrum as function of f
    """
    f = np.asarray(f)
    S = 1*(4*sigma**2 * L / V0) / (1 + 6 * f * L/V0) ** (5./3) 
    S[f==0] = 0 # V0**2 * (tMax+dt), but we will add the mean at the end
    return S

def kaimal(f, V0, sigma, L, angularFrequency=False):
    """ 
    Single-sided Kaimal spectrum S_X

    See also welib.standards.IEC for definition according to classes, ETM/NTM

    INPUTS:
     - f: frequency, in [Hz] if angularFrequency=False, in [rad/s] if angularFrequency is True
             Sf(f) df = So(om)dom   =>   So(om) = Sf(f) / (2pi)

     - V0: reference/hub velocity velocity
     - sigma: wind standard deviation, turbulence intensity I= sigma/V0. 
     - L: turbulence length scale. 
              IEC: L = 0.7*min(60, zhub)* [8.10, 2.70, 0.66]
    OUTPUTS:
      - S: kaimal spectrum
    """
    if angularFrequency:
        # Input "f" is actually omega
        omega = f
        f = omega/(2*np.pi)
        scale=1/(2*np.pi)
    else:
        scale=1
    S = scale * (4 * sigma**2 * L / V0) / (1 + 6 * f * L / V0) ** (5./3) 
    f = np.asarray(f)
    S = np.asarray(S)
    b=f==0
    S[b] = 0
    return S




if __name__ == '__main__':
    import matplotlib.pyplot as plt
    import numpy as np
    omega = np.linspace(0,10, 100)
    S= kaimal_omega(omega, V0=10, sigma=1.5, L=300)
    S2= kaimal(omega, V0=10, sigma=1.5, L=300, angularFrequency=True)
    fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
    fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
    ax.plot(omega, S , label='')
    ax.plot(omega, S2,'--', label='')
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.legend()
    plt.show()

