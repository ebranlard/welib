""" 
Rankine half body / Rankine nose

Rankine oval
"""
import numpy as np

# --------------------------------------------------------------------------------}
# --- Rankine nose  / Half body
# --------------------------------------------------------------------------------{
""" 
Rankine nose:

- The stagnation point is at (if the source is at (0,0)):
    x0 = -Sigma/(2*np.pi*U0)   (upstream)

- Maximum extent far downstream (large x) is:
    y_max = +/- Sigma/(2*U0)

- Extent right above the source is:
   y = +/- Sigma/(4*U0)

- Separating streamline: psi = Sigma/2

- Cp on surface: 2 sin(theta) * cos(theta) /(pi-theta) + (sin(theta)/(pi-theta))**2

- Minimmum pressure 
      theta_m = atan ( pi-theta_m / (pi-theta_m -1))  -> theta_m approx 92.9569deg
      Cp_min = -0.5866

"""
def rn_stag(Sigma=1, U0=0):
    """ Stagnation point (relevant if freestream)"""
    if abs(U0)>0:
        ys=0
        xs=Sigma/(2*np.pi*U0)
    else:
        ys=0
        xs=0
    return xs, ys

def rn_coord(Sigma=1, U0=0, x_max=None, ne=100):
    x0 = Sigma/(2*np.pi*U0)
    if x_max is None:
        x_max = 8*x0
    x = np.linspace(-x0, x_max, ne)
    y = rn_coords_yx(x, Sigma=Sigma, U0=U0)
    return x, y

def rn_coord_theta(Sigma=1, U0=0, x_max=4, ne=100, bodyAt0=False, both=True):
    """ 
    The separating streamline delimiting the inner and outer flow is given by the value Sigma/2
    """
    if U0==0:
        return np.nan
    y_max = Sigma/(2*U0)
    x0 = Sigma/(2*np.pi*U0)
    #theta_max = np.arctan2(y_max,x_max)
    # Obtained by using y(x_max)  = x0 (pi-theta_m) ->  tan(theta_m) = ym/xm ~ theta_m
    # Alternative, nonlinear solve for theta_max
    theta_max = np.pi * x0/(x0+x_max) 
    if both:
        eps=0.001
        theta1 = np.linspace(theta_max, np.pi-eps, ne)          # upper surface
        theta2 = np.linspace(np.pi+eps, -theta_max+2*np.pi, ne) # lower surface
        theta = np.concatenate( (theta1, theta2) )
    else:
        theta = np.linspace(theta_max, np.pi, ne)          # upper surface
    # Equation for r obtained from the streamline equation = Sigma/2
    r = Sigma/(2*np.pi*U0)*(np.pi-theta)/np.sin(theta)
    x = r*np.cos(theta)
    y = r*np.sin(theta)
    if bodyAt0:
        x += x0
    return x, y

def rn_coords_yx(x, Sigma=2*np.pi, U0=1):
    """ Return y(x) (using intepolation for now...)
    TODO make it dimensionless
    """
    ne = len(x)*10
    x_max = np.max(x)+1
    y_max = Sigma/(2*U0)
    eps=0.001
    x0 = Sigma/(2*np.pi*U0)
    #theta_max = np.arctan2(y_max,x_max)
    #print(theta_max)
    # Obtained by using y(x_max)  = x0 (pi-theta_m) ->  tan(theta_m) = ym/xm ~ theta_m
    # Alternative, nonlinear solve for theta_max
    theta_max = np.pi * x0/(x0+x_max) 
    theta = np.linspace(np.pi, theta_max, ne)          # upper surface
    r = Sigma/(2*np.pi*U0)*(np.pi-theta)/np.sin(theta)
    x1 = r*np.cos(theta)
    y1 = r*np.sin(theta)
    x1[0] = -x0
    y1[0] = 0
    y = np.interp(x, x1, y1)
    return y

def rn_rtheta(theta, Sigma=2*np.pi, U0=1):
    """ r(theta) for rankine nose surface"""
    r = Sigma/(2*np.pi*U0)*(np.pi-theta)/np.sin(theta)
    return r

def rn_u(x, y, Ps=[0,0], Sigma=1,  U0=1, **kwargs):
    from welib.vortilib.elements.SourcePoint import sp2d_u
    return sp2d_u(x, y, Ps, Sigma=Sigma, U0=U0, **kwargs)

def rn_CpTheta(theta):
    """
     Cp on surface: 2 sin(theta) * cos(theta) /(pi-theta) + (sin(theta)/(pi-theta))**2
    """
    bOK = np.abs(theta-np.pi)>1e-8
    Cp =np.zeros_like(theta)
    th = theta[bOK]
    Cp[bOK]=- (2 *np.sin(th) * np.cos(th) /(np.pi-th) + (np.sin(th)/(np.pi-th))**2)
    Cp[~bOK]=1
    return Cp

def rn_Cpx(x, Sigma=2*np.pi, U0=1):
    """
        Cp(x) 
    TODO make it dimensionless
    """
    y = rn_coords_yx(x, Sigma=Sigma, U0=U0)
    x0 = Sigma/(2*np.pi*U0)
    theta = np.arctan2(y, x)
    Cp = rn_CpTheta(theta)
    return Cp



# --------------------------------------------------------------------------------}
# ---  Rankine Oval 
# --------------------------------------------------------------------------------{
""" 
Rankine Oval:


"""

def ro_params(R, t, U0=1, method='least_squares', verbose=True):
    """ 
    find Sigma and b such that the oval has length "2R" and thickness "2t"
     Need to solve for:
       tan (2pi U_0 t_0  / Sigma ) = 2 b t_0 / (t_0^2 - b^2)
       Sigma = pi U0 (R^2 - b^2) / b 
     
    """
    from scipy.optimize import fsolve
    from scipy.optimize import least_squares

    def ro_height_1(Sigma, b, t):
        h = (t**2 - b**2)/(2*b) * np.tan(2*np.pi*U0*t/Sigma)
        return h

    def equations(p):
        Sigma, b = p
        #eq1 = np.tan(2*np.pi*U0*t  / Sigma ) - 2*b*t / (t**2 - b**2)
        eq1 = U0* t  - Sigma/(2*np.pi) * np.arctan2(2*b*t , (t**2-b**2))
        eq2 = Sigma - np.pi* U0*(R**2 - b**2) / b 
        return (eq1, eq2)

    if method=='least_squares':
        smin=0
        smax=25
        bmin=0
        bmax=1
#         if t/R<0.1:
#             smin = 0
#             smax = 0.1
#             bmin = 0.95
#             bmax = 1.0001
#         elif t/R<0.6:
#             smin = 0
#             smax = 0.7
#             bmin = 0.7
#             bmax = 1.0001
#         elif t/R<0.98:
#             smin = .6
#             smax = 6
#             bmin = 0.1
#             bmax = 0.8
#         elif t/R>=0.98:
#             smin = 5.5
#             smax = 25
#             bmin = 0
#             bmax = 0.2
#         else:
#             raise NotImplementedError()
        # --- Try 1
        res = least_squares(equations, (np.pi*U0*R*(smin+smax)/2, R*(bmin+bmax)/2), bounds = ((0, bmin*R), (smax*np.pi*U0*R, bmax*R)))
        Sigma, b = res.x
        h = ro_height_1(Sigma, b, t)
        if np.abs(h-t)/R*100<1:
            return Sigma, b
        if verbose:
            print('[FAILED] Found height:',h, 'target:', t)

    else:
        Sigma, b =  fsolve(equations, (np.pi*U0*R, R*0.8))
    return Sigma, b

def ro_curve(x,y, b, Sigma, U0):
    """" Implicit curve defining the oval"""
    residual = U0*y - Sigma/(2*np.pi) * np.arctan2( 2* b * y, (x**2 + y**2 -a**2))
    return residual

def ro_coord(b, Sigma, U0, ne=100):
    x0, y0= ro_stag(b=b, Sigma=Sigma, U0=U0)
    xs = np.linspace(x0[0], x0[1], ne)
    ys = ro_coord_yx(xs, b=b, Sigma= Sigma, U0=U0)
    return xs, ys

def ro_coord_yx(x, b, Sigma, U0=1, method='least_squares', verbose=True):
    """ 
    """
    from scipy.optimize import least_squares
    from scipy.optimize import fsolve

    vy = np .zeros_like(x)
    y0 = 0.2*b
    xstag = np.sqrt(b**2 + Sigma*b / (np.pi*U0))
    for ix, x in enumerate(x):
        if np.abs(np.abs(x)-xstag)<0.00001*xstag:
            vy[ix]=0
        else:
            def equation(p):
                y = np.abs(p)
                return U0*y - Sigma/(2*np.pi) * np.arctan2( 2* b * y, (x**2 + y**2 -b**2))
                #return U0*y + Sigma/(2*np.pi) *( np.arctan2( y, x+b) - np.arctan2(y, x-b))
            #res = least_squares(equation, (y0), bounds = (0, ))
            vy[ix] = np.asarray(fsolve(equation, y0))[0]
            y0 = max(vy[ix],0.2*b)
        #print(y0)
    return vy


def ro_stag(b, Sigma, U0=1):
    xstag = np.sqrt(b**2 + Sigma*b / (np.pi*U0))
    xs=[-xstag,xstag]
    ys=[0,0]
    return xs, ys

def ro_psi(x, y, Sigma=1, b=1, U0=1, ):
    """ 
    The stagnation streamline is given for psi =0
    """
    r2= y**2 + y**2
    theta1 = np.arctan2(y,x+b)
    theta2 = np.arctan2(y,x-b)
    #T = np.arctan2( 2*b*y , (r2-b**2)) 
    #psi = U0*y + Sigma/(2*np.pi) * T
    psi = U0*y + Sigma/(2*np.pi) * (theta1-theta2)
    return psi

def ro_u(x, y, Sigma=1, b=1, U0=1):
    from welib.vortilib.elements.SourcePoint import sp2d_u
    U1, V1= sp2d_u(x, y, [-b,0], Sigma=Sigma)
    U2, V2= sp2d_u(x, y, [b,0], Sigma=-Sigma)
    U = U1 + U2 + U0
    V = V1 + V2 
    return U, V


def ro_Cp(x, y, U0=1, Sigma=1, b=0):
    U, V =ro_u(X, Y, U0=U0, Sigma=Sigma, b=b)
    U2 = U**2 + V**2 
    Cp = 1-U2/U0**2
    return Cp

def ro_Cpx(x, U0=1, Sigma=1, b=0):
    y = ro_coord_yx(x, U0=U0, Sigma=Sigma, b=b)
    U, V =ro_u(x, y, U0=U0, Sigma=Sigma, b=b)
    U2 = U**2 + V**2 
    Cp = 1-U2/U0**2
    return Cp









if __name__ == '__main__':
    import matplotlib.pyplot as plt

    # --------------------------------------------------------------------------------}
    # --- Rankine nose Cp  
    # --------------------------------------------------------------------------------{
    #Sigma =1
    #U0 =1 
    #y_max = Sigma/(2*U0)
    #eps=0.001
    #theta_max = np.arctan2(y_max, 4)
    #theta = np.linspace(theta_max, np.pi-eps, 100)          # upper surface
    #r = rn_rtheta(theta, Sigma=Sigma, U0=U0 )
    #x = r*np.cos(theta)
    #Cp = rn_CpTheta(theta)
    #x0 = -Sigma/(2*np.pi*U0)
    #x2 = np.linspace(x0, 4, 30)
    #Cpx = rn_Cpx(x2, Sigma=Sigma, U0=U0)
    #fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
    #fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
    #ax.plot(x, Cp, label='')
    #ax.plot(x2, Cpx, 'k.', label='')
    #ax.set_xlabel('')
    #ax.set_ylabel('')
    ##ax.legend()

    #x0 = -Sigma/(2*np.pi*U0)
    #x = np.linspace(x0, 4, 150)
    #y = rn_coords_yx(x, Sigma=Sigma, U0=U0)
    #fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
    #fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
    #ax.plot(x, y, label='')
    #ax.plot(x,-y, label='')
    #ax.set_xlabel('')
    #ax.set_ylabel('')



    # --------------------------------------------------------------------------------}
    # --- Rankine Oval parameters variation
    # --------------------------------------------------------------------------------{
    U0=1
    vt = np.linspace(0.0001,0.9999, 100)
    sig  = vt*0
    bb   = vt*0
    sig2 = vt*0
    bb2  = vt*0
    for it, t in enumerate(vt):
        Sigma, b = ro_params(R=1, t=t*1, U0=U0) #, method='ajkj')
        sig[it] = Sigma
        bb[it] = b
        Sigma, b = ro_params(R=2, t=t*2, U0=U0) #, method='ajkj')
        sig2[it] = Sigma
        bb2[it] = b
    fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
    fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
    ax.plot(vt, sig /(U0*np.pi*1)    , label='Sigma')
    ax.plot(vt, bb                 , label='b')
    ax.plot(vt, sig2/(U0*np.pi*2),'--'   , label='Sigma')
    ax.plot(vt, bb2/2            ,'--'   , label='b')
    ax.set_xlabel('t/R [-]')
    ax.set_ylabel('')
    ax.legend()


# 
    plt.show()
