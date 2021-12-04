import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# Local 
from welib.airfoils.Polar import Polar
import welib.weio as weio

MyDir=os.path.dirname(__file__)


def main_correction3D(test=False):
    polarFile_in = os.path.join(MyDir,'../data/DU21_A17.csv')

    r_over_R     = 0.2
    chord_over_r = 3./5.
    tsr          = 10

    polar = Polar(polarFile_in, compute_params=True)
    #ADpol = polar.toAeroDyn(polarFile_AD)
    polar3D= polar.correction3D(r_over_R, chord_over_r, tsr)

    return polar, polar3D


if __name__ == '__main__':
    polar,polar3D = main_correction3D()

    import matplotlib.pyplot as plt
    plt.plot(polar.alpha, polar.cl       , label= 'cl')
    plt.plot(polar3D.alpha, polar3D.cl   , label= 'cl 3d')
    plt.plot(polar.alpha, polar.cl_inv   , label= 'cl inv')
    plt.ylim([-1.5,2])
    plt.legend()
    plt.show()

if __name__ == '__test__':
    polar,polar3D = main_correction3D()
