import numpy as np

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
