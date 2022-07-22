import numpy as np



def dqdt_lorenz(t, q, sigma, beta, rho):
    """ 
    Lorenz system
    """
    x, y, z = q
     
    dx = sigma * (y - x)
    dy = x * (rho - z) - y
    dz = x * y - beta * z
     
    return [dx, dy, dz]
