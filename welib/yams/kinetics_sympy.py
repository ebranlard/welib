import sympy as sp


def Fspring3D(r1,r2,k,l0):
    """ return linear spring force between two points """
    l  = sp.sqrt((r2-r1).dot(r2-r1))
    return -k * (l-l0)*(r2-r1)/l

def Fdamp3D(r1,r2,r1d,r2d,c):
    """ return linear damping force between two points """
    l  = sp.sqrt(r1.dot(r2))
    e  = (r2-r1)/l
    dv = (r2d-r1d).dot(e)
    return -c * dv * e



def rigidTransformationTwoPoints_Loads(Ps, Pd):
    """ 
    Relate loads at source node to destination node:

    See also welib.FEM.utils
      fd = T.dot(fs)
       T =[ I3           0  ] =  [ I3         0  ]
          [ skew(Ps-Pd)  I3 ]    [ skew(r0)   I3 ]
    """
    from welib.yams.yams_sympy import skew
    T = sp.eye(6) # 1 on the diagonal
    T[3:6,0:3] = skew(Ps-Pd)
    return T

