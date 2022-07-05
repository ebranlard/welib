import numpy as np


def Fspring3D(r1,r2,k,l0):
    """ return linear spring force between two points """
    r1 = np.asarray(r1).flatten()
    r2 = np.asarray(r2).flatten()
    l  = np.sqrt((r2-r1).dot(r2-r1))
    dl = (l-l0)
    return -k * dl *(r2-r1)/l

def Fdamp3D(r1,r2,r1d,r2d,c):
    """ return linear damping force between two points """
    r1  = np.asarray(r1).flatten()
    r2  = np.asarray(r2).flatten()
    r1d = np.asarray(r1d).flatten()
    r2d = np.asarray(r2d).flatten()
    l  = np.sqrt(r1.dot(r2))
    e  = (r2-r1)/l
    dv = (r2d-r1d).dot(e)
    return -c * dv * e
