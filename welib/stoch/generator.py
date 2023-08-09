import numpy as np
from numpy.random import uniform
from numpy.random import normal


def generateGaussianNoiseBoxMuller(mu=0, sigma=1, size=1):

    u1 = uniform(0, 1      , size)  # random phases between 0 and 2pi
    u2 = uniform(0, 2*np.pi, size)  # random phases between 0 and 2pi

    #W = sigma  * np.sqrt(-2*np.log(u1))*np.exp(1j * u2)

    mag = sigma * np.sqrt(-2.0 * np.log(u1))
    z0  = mag * np.cos(u2) + mu
    #z1  = mag * np.sin(u2) + mu
    return z0
