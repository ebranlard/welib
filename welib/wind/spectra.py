import numpy as np
import pandas as pd


def kaimal(f, V0, sigma, L):
    """ 
    Kaimal spectrum

    See also welib.standards.IEC for definition according to classes, ETM/NTM

    INPUTS:
     - f: frequency (NOTE: can be in Hz or rad/s)
     - V0: reference/hub velocity velocity
     - sigma: wind standard deviation, turbulence intensity I= sigma/V0. 
     - L: turbulence length scale. 
              IEC: L = 0.7*min(60, zhub)* [8.10, 2.70, 0.66]
    OUTPUTS:
      - S: kaimal spectrum
    """
    return (4 * sigma**2 * L / V0) / (1 + 6 * f * L / V0) ** (5./3)



if __name__ == '__main__':
    pass
