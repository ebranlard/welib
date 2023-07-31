import numpy as np
import pandas as pd


def kaimal(f, V0, sigma, L, angularFrequency=False):
    """ 
    Kaimal spectrum

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
    pass
