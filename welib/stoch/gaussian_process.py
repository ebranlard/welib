""" 
Gaussian process

"""

import numpy as np
# local
from welib.stoch.stationary_process import StationaryProcess



# --------------------------------------------------------------------------------}
# --- White noise 
# --------------------------------------------------------------------------------{
""" 
White noise, is a homogeneous (stationary) Gaussian process with mean value function mu_X=0 and an autospectra density function for all angular freqeuncies
S_XX = S0

"""
class WhiteNoise(StationaryProcess):
    pass




if __name__ == '__main__':
    from welib.essentials import *
    from welib.stoch.utils import plot_pdf
    from welib.stoch.variable import *
    from welib.stoch.distribution import *
    from welib.stoch.generator import *
    import matplotlib.pyplot as plt

    
    plt.show()
