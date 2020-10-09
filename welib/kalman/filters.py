import numpy as np

def moving_average(a, n=3) :
    """ 
    perform moving average, return a vector of same length as input
    """
    a=a.ravel()
    a = np.concatenate(([a[0]]*(n-1),a)) # repeating first values
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    ret=ret[n - 1:] / n
    return ret

