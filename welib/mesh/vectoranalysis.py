import numpy as np


# --------------------------------------------------------------------------------}
# --- Functions to reproduce Matlab's gradient, div, and curl: axis 0 and 1 inverted!
# --------------------------------------------------------------------------------{
def matlab_gradient_1d(F,*args,**kwargs):
    """ Computes gradient of a 1d vector field in python like matlab gradient function. 
    call syntax:
       gradient(F) or gradient(F,dx) or gradient(F,x)
    """
    shape_in=F.shape
    G=np.gradient(F.ravel(),*args,**kwargs)
    return G.reshape(shape_in)

def matlab_gradient_2d(F,*args,**kwargs):
    """ Computes gradient of a 1d vector field in python like matlab gradient function. 

    call syntax:
       gradient(F) or gradient(F,dx,dy) or gradient(F,x,y)
    """
    G = np.gradient(F.T,*args,**kwargs)
    G = (G[0].T, G[1].T)
    return G

def matlab_div_2d(*args):
    """ Computes divergence of a 2D vector field in python like matlab div function

    call syntax:
       div(U,V) or div(X,Y,U,V)
    """
    if len(args)==2:
        sz = args[0].shape
        dx = np.arange(1,sz[1]+1) # NOTE: matlab, x is along the 2nd dimension..
        dy = np.arange(1,sz[0]+1)
        U=args[0]
        V=args[1]
    elif len(args)==4:
        dx = args[0][0,:];
        dy = args[1][:,0];
        U=args[2]
        V=args[3]
    else:
        raise Exception('Input 2 or 4 arguments')
    #retval  = matlab_gradient_2d (U, dx, dy)[0]
    #retval += matlab_gradient_2d(V.T, dy, dx)[0].T;
    return np.gradient(U, dx, axis=1) +  np.gradient(V,  dy, axis=0)

def matlab_curl_2d(*args):
    """ Computes curl of a D vector field in python like matlab curl function
    
    call syntax:
       curl(U,V) or curl(X,Y,U,V)
    """
    if len(args)==2:
        sz = args[0].shape
        dx = np.arange(1,sz[1]+1) # NOTE: matlab, x is along the 2nd dimension..
        dy = np.arange(1,sz[0]+1)
        U=args[0]
        V=args[1]
    elif len(args)==4:
        dx = args[0][0,:];
        dy = args[1][:,0];
        U=args[2]
        V=args[3]
    else:
        raise Exception('Input 2 or 4 arguments')
    dFx_dy = matlab_gradient_2d(U.T, dy, dx)[0].T
    dFy_dx = matlab_gradient_2d(V, dx, dy)[0]
    rot_z  = dFy_dx - dFx_dy                    
    av     = rot_z / 2
    return rot_z, av


# --------------------------------------------------------------------------------}
# --- Pythonic versions 
# --------------------------------------------------------------------------------{
def div(f):
    """
    TODO
    Computes the divergence of the vector field f, corresponding to dFx/dx + dFy/dy + ...
    :param f: List of ndarrays, where every item of the list is one dimension of the vector field
    :return: Single ndarray of the same shape as each of the items in f, which corresponds to a scalar field
    """
    num_dims = len(f)
    return np.ufunc.reduce(np.add, [np.gradient(f[i], axis=i) for i in range(num_dims)])







