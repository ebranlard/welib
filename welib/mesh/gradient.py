import numpy as np
import math


def differ_stencil ( x0, o, p, x ):
    """ Determine coefficients C to approximate the derivative at X0
         of order O and precision P, using finite differences, so that 
          d^o f(x)/dx^o (x0) = sum ( 0 <= i <= o+p-1 ) c(i) f(x(i)) + O(h^(p))
        where H is the maximum spacing between X0 and any X(I).
    Based on a Fortran code from John Burkardt licensed under the GNU LGPL license.
    """
    n = o + p
    assert(len(x)==n)
    dx = x - x0
    b = np.zeros(n)
    c = np.zeros(n)
    b[o] = 1
    c = vm_solve(n, dx, b)
    fact = math.factorial(o)
    c *= fact
    return c

def vm_solve (n, a, b):
    """ 
    Solves a linear  system Ax = b  (Vandermonde) for x
    Based on a Fortran code from John Burkardt licensed under the GNU LGPL license.

    The R8VM storage format is used for an M by N Vandermonde matrix.
    An M by N Vandermonde matrix is defined by the values in its second
    row, which will be written here as X(1:N).  The matrix has a first 
    row of 1's, a second row equal to X(1:N), a third row whose entries
    are the squares of the X values, up to the M-th row whose entries
    are the (M-1)th powers of the X values.  The matrix can be stored
    compactly by listing just the values X(1:N).
   
    Vandermonde systems are very close to singularity.  The singularity
    gets worse as N increases, and as any pair of values defining
    the matrix get close.  Even a system as small as N = 10 will
    involve the 9th power of the defining values.
    INPUTS:
     - n: number of rows and columns of A
     - A(n): input matrix
     - B(n): right hand side
    OUTPUTS:
     - x(n): solution of the linear system
    """
    #  Check for explicit singularity.
    for j in range(n-1):
        for i in range(j+1, n):
            if a[i] == a[j]:
              raise Exception('System is singular')
    # Solve
    x = b.copy()
    for j in range(n-1):
        for i in  range(n-1, j, -1):
            x[i] = x[i] - a[j] * x[i-1]
    for j in  range(n-1, -1, -1):
        for i in  range(j+1, n):
            x[i] = x[i] / ( a[i] - a[i-j-1] )
        for i in  range(j, n-1):
            x[i] = x[i] - x[i+1]
    return x

def gradient_regular(f,dx=1,order=4):
    # Compute gradient of a function on a regular grid with different order
    # INPUTS
    #   f    : 1d field
    #   order: differentiation order
    # Optional argument
    #   dx   : cell spacing in each dimensions, by default d = 1 

    n    = len(f)  # Grid Dimensions

    gradf=np.zeros(len(f))
    d2=2*dx
    if order==2:
        # --------------------------------------------------------------------------------
        # --- ORDER 2 
        # --------------------------------------------------------------------------------
        # left boundary: forward difference   (a=0, b=3)&
        # right boundary: backward difference (a=3, b=0)&
        # elsewhere: centered difference (a=1,b=1)&
        for i in range(n):
            if (i==0):
                df_dx     =  ( - 3 * f[i  ] + 4 * f[i+1] - 1 * f[i+2])/d2;
            elif (i==n-1):
                df_dx     =  (   1 * f[i-2] - 4 * f[i-1] + 3 * f[i  ])/d2;
            else:
                df_dx     =  ( - 1 * f[i-1] + 1 * f[i+1])/d2;
            gradf[i]=df_dx;
    elif (order==4) :
        # --------------------------------------------------------------------------------
        # --- ORDER 4 
        # --------------------------------------------------------------------------------
        # left boundary: forward difference   (a=0, b=5) and (a=1, b=3)&
        # right boundary: backward difference (a=3, b=1) and (a=0, b=5)&
        # elsewhere: centered difference (a=2,b=2)&
        for i in range(n):
            if (i==0) :
                df_dx = ( - 25/6 * f[ i]+ 8      * f[ i+1]- 6    * f[ i+2] +8/3 * f[ i+3] -1/2   * f[ i+4])/d2;
            elif (i==1):
                df_dx= ( - 1/2   * f[ i-1] - 5/3 * f[ i] +3      * f[ i+1]-1    * f[ i+2]+ 1/6   * f[ i+3])/d2;
            elif (i==n-2):
                df_dx= ( - 1/6   * f[ i-3] +1    * f[ i-2] -3    * f[ i-1]+ 5/3 * f[ i]+ 1/2     * f[ i+1])/d2;
            elif (i==n-1):
                df_dx = (1/2     * f[ i-4]-8/3   * f[ i-3]+ 6    * f[ i-2] - 8  * f[ i-1] + 25/6 * f[ i])/d2;
            else:
                df_dx= ( 1/6     * f[ i-2] - 4/3 * f[ i-1] + 4/3 * f[ i+1]- 1/6 * f[ i+2])/d2;
            gradf[i]=df_dx;
    else:
       raise Exception('Order not implemented: ',order);


    return gradf

if __name__ == '__main__':
    dx=0.1
    x     = np.arange(0,1,dx)
    phi   = x**2
    dphi  = gradient_regular(phi,dx)
    ddphi = gradient_regular(dphi,dx)
    print(phi, len(phi))
    print(dphi, len(dphi))
    print(ddphi, len(ddphi))
