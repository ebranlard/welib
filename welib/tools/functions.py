import numpy as np

import matplotlib.pyplot as plt

def delta(x, epsilon=None):
    """ Dirac function"""
    if epsilon is not None and epsilon>=1e-5:
        rng=(-epsilon,epsilon)
        d = smooth_delta(x, e=1, rng=rng, method='gaussian')
        return d

    else:
        d = np.zeros(x.shape)
        vx=np.unique(x)
        if np.min(vx)>0:
            return d
        elif np.max(vx)<0:
            return d
        i = np.argmin(np.abs(vx))
        x0_new = vx[i]
        dx=vx[1]-vx[0]
        d[x==x0_new] = 1/dx
        return d

def Pi(x, epsilon=None):
    """ Step function"""
    if epsilon is not None and epsilon>=1e-5:
        rng=(-epsilon,epsilon)
        p = smooth_heaviside(x, k=1, rng=rng, method='exp')
        return p
    else:
        p = np.zeros(x.shape)
        p[x>0] =1
        return p

def smooth_heaviside(x, k=1, rng=(-np.inf, np.inf), method='exp'):
    """ 
    Smooth approximation of Heaviside function where the step occurs between rng[0] and rng[1]:
       if rng[0]<rng[1]: then  f(<rng[0])=0, f(>=rng[1])=1
       if rng[0]>rng[1]: then  f(<rng[1])=1, f(>=rng[0])=0

    exp: 
       rng=(-inf,inf):  H(x)=[1 + exp(-2kx)            ]^-1
       rng=(-1,1):      H(x)=[1 + exp(4kx/(x^2-1)      ]^-1
       rng=(0,1):       H(x)=[1 + exp(k(2x-1)/(x(x-1)) ]^-1

    INPUTS: 
        x  : scalar or vector of real x values \in ]-infty; infty[ 
        k  : float >=1, the higher k the "steeper" the heaviside function
        rng: tuple of min and max value such that f(<=min)=0  and f(>=max)=1. 
             Reversing the range makes the Heaviside function from 1 to 0 instead of 0 to 1
        method: smooth approximation used (e.g. exp or tan)

    """
    if k<1:
        raise Exception('k needs to be >=1')
    eps = 1e-2
    mn,mx = rng
    x     = np.asarray(x)
    H     = np.zeros(x.shape)
    if mn<mx:
        H[x<=mn] = 0
        H[x>=mx] = 1
        b        = np.logical_and(x>mn, x<mx)
    else:
        H[x<=mx] = 1
        H[x>=mn] = 0
        b        = np.logical_and(x<mn, x>mx)
    x=x[b]
    if method=='exp':
        if np.abs(mn)==np.inf and np.abs(mx)==np.inf:
            # Infinite support
            x[k*x>100 ]  = 100./k
            x[k*x<-100] = -100./k
            if mn<mx:
                H[b] = 1 / ( 1+np.exp(-  k * x))
            else:
                H[b] = 1 / ( 1+np.exp(   k * x))
        elif np.abs(mn)!=np.inf and np.abs(mx)!=np.inf:
            n=4.
            # Compact support
            s= 2./(mx -mn) * (x-(mn+mx)/2.) # transform compact support into ]-1,1[ 
            x = -n*s/(s**2-1.)             # then transform   ]-1,1[  into ]-inf,inf[
            x[k*x>100 ]  = 100./k
            x[k*x<-100] = -100./k
            H[b] = 1./(1+np.exp(-k*x))
        else:
            raise NotImplementedError('Heaviside with only one bound infinite')
    else:
        # TODO tan approx
        raise NotImplementedError()
    return H


def smooth_delta_inf(xx, e, method='gaussian'):
    """ 
    Smooth approximation of delta function between -inf and inf

    INPUTS: 
      - x  : scalar or vector of real x values \in ]-infty; infty[ 
      - e  : defining the "epsilon". NOTE: this value has different meaning depending on the method.
      - method: method used to approximate the delta function. 
           - 'sin':        1/(pi x) sin(x/e)
           - 'gaussian':   1./(2*np.sqrt(np.pi*e)) *np.exp( -xx**2/(4.*e))
           - 'frac':       1/pi e/(x^2-e^2)
           - 'exp-heaviside': derivative of the heaviside function with exp smoothening
    OUTPUTS:
      - value(s) at x
    """
    if method=='frac':
        return 1./np.pi * (e/(xx**2+e**2))
    elif method=='gaussian':
        return 1./(2*np.sqrt(np.pi*e)) *np.exp( -xx**2/(4.*e))
    elif method=='sin':
        return 1./(np.pi*xx)*np.sin(xx/e)
    elif method=='exp-heaviside':
        # NOTE: this is buggy..
        k=e
        xx[k*xx >100]  = 100./k
        xx[k*xx<-100] = -100./k
        E        = np.exp(-1*k*xx)
        return  1.*k*E / (1.+E)**2
    else:
        raise NotImplementedError('Method {}'.format(method))


def smooth_delta(x, e=1, rng=(-np.inf, np.inf), method='gaussian'):
    """ 
    Smooth approximation of delta function between rng[0] and rng[1]
    The integral over the range should be 1.

    INPUTS: 
      - x  : scalar or vector of real x values \in ]-infty; infty[ 
      - e  : defining the "epsilon". NOTE: this value has different meaning depending on the method.
      - rng: tuple of min and max value such that f(<=min)=0  and f(>=max)=0
      - method: method used to approximate the delta function. 
           - 'sin':        1/(pi x) sin(x/e)
           - 'gaussian':   1./(2*np.sqrt(np.pi*e)) *np.exp( -xx**2/(4.*e))
           - 'frac':       1/pi e/(x^2-e^2)
           - 'exp-heaviside': derivative of the heaviside function with exp smoothening
    OUTPUTS:
      - value(s) at x
    """
    eps = 1e-2
    mn,mx = rng
    if mn>=mx:
        raise Exception('Range needs to be increasing')
    x     = np.asarray(x)
    delta = np.zeros(x.shape)
    delta[x<=mn] = 0
    delta[x>=mx] = 0
    b  = np.logical_and(x>mn, x<mx)
    x=x[b]

    if np.abs(mn)!=np.inf and np.abs(mx)!=np.inf:
        n=4. 
        # Compact support
        s     = 2./(mx -mn) * (x-(mn+mx)/2) # first transform ]-mn,mx[ into ]-1,1[
        x_new = -n*s/(s**2-1)              # then transform   ]-1,1[  into ]-inf,inf[
        # NOTE: the first transformation is linear, so the integral scales accordingly
        #       the second is not, for now we compute the scaling numerically...
        scale=2./(mx-mn) # scaling due to linear transformation
        s0 = np.linspace(-1-1e-5,1-1e-5,100)
        x0 = -n*s0/(s0**2-1)
        if method=='frac': # further scaling adjustements
            scale/=np.trapz(1/np.pi * (e/(x0**2+e**2)), s0)
        elif method=='gaussian':
            scale/=np.trapz(1./(2*np.sqrt(np.pi*e)) *np.exp( -x0**2/(4.*e)), s0)
        elif method=='sin':
            scale*=n
        elif method=='exp-heaviside':
            pass
        else:
            raise NotImplementedError('Method {}'.format(method))

        if method=='exp-heaviside':
            k=e
            E =  np.exp( k*4.*s/(s**2-1.-eps)) # NOTE: adding an epsilon to limit overflow 
            delta[b] = 4.*k*(s**2+1) * E/ ( (s**2-1)**2 * (1.+E)**2)*scale
        else:
            delta[b] = smooth_delta_inf(x_new, e, method=method) * scale

    elif np.abs(mn)==np.inf and np.abs(mx)==np.inf:
        delta[b] = smooth_delta_inf(x, e, method=method)
#         k=e
#         if np.abs(mn)!=np.inf and np.abs(mx)!=np.inf:
#         else:
#             # Avoid exp overflow
#             x[x>10]  = 10
#             x[x<-10] = -10
#             E        = np.exp(-1*k*x)
#             delta[b] =  1*k*E / (1+E)**2
# 
#         elif method=='exp-heaviside':
#             x0[x0>10]  = 10
#             x0[x0<-10] = -10
#             k=e
#             E = np.exp(-1*k*x0)
#             scale/=np.trapz(1*k*E / (1+E)**2, s0)

#             s= 2/(mx -mn) * (x-(mn+mx)/2) # transform compact support into ]-1,1[ 
    else:
        # TODO tan approx
        raise NotImplementedError()
    return delta



if __name__=='__main__':

    # --- Smooth heaviside
#     delta= 1
#     eps = 1e-15
#     eps2 = 1e-2
#     x = np.linspace(-10,10,100000)
#     R2=(-7,7)
#     k2=2
#     H1=smooth_heaviside(x, k=1, rng=(-np.inf, np.inf))
#     H2=smooth_heaviside(x, k=k2, rng=R2)
#     H3=smooth_heaviside(x, k=1, rng=(0, 3))
#     H4=smooth_heaviside(x, k=1, rng=(-0.5, 0.5))
#     H5=smooth_heaviside(x, k=1, rng=(np.inf,-np.inf))
#     D2_num=np.diff(H2,prepend=0)/(x[1]-x[0])
#     fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
#     fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
# #     ax.plot(x, H1, '-', label='H1')
# #     ax.plot(x, H2, '--', label='H2')
# #     ax.plot(x, H3, '-', label='H3')
# #     ax.plot(x, H4, ':', label='H4')
# #     ax.plot(x, H5, ':', label='H5')
#     ax.set_xlabel('x ')
#     ax.set_ylabel('')
#     ax.legend()
#     ax.tick_params(direction='in')
#     plt.show()

    # --- Smooth Dirac between -inf and inf, check different methods
    epsilon=10.1
    x = np.linspace(-1e2,1e2,100001)
    fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
    fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
    for method in ['frac','gaussian']:
        D=smooth_delta_inf(x, e=epsilon, method=method)
        print('Method {:15s}:  integral:{}'.format(method, np.trapz(D,x)))
        ax.plot(x, D/np.max(D), '-', label=method)
    ax.set_xlabel('x ')
    ax.set_ylabel('')
    #ax.set_yscale('log')
    ax.legend()
    ax.tick_params(direction='in')
    ax.set_title('Smooth Dirac - infinite support - different methods')


    # --- Smooth Dirac - Check different methods
    epsilon=1.1
    rng=(-1,1)
    x = np.linspace(-2,2,100001)
    fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
    fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
    for method in ['frac','exp-heaviside','gaussian']:#, 'sin']:
        D=smooth_delta(x, e=epsilon, rng=rng, method=method )
        print('Method {:15s}:  integral:{}'.format(method, np.trapz(D,x)))
        ax.plot(x, D/np.max(D), '-', label=method)
    ax.set_xlabel('x ')
    ax.set_ylabel('')
    ax.legend()
    ax.tick_params(direction='in')
    ax.set_title('Smooth Dirac - compact support - different methods')

    # --- Smooth Dirac - Check different epsilon
    method='frac'; epsilons = [0.3,0.7,1.1,1.7]
    #m='exp-heaviside'; epsilons=[1,1.1,1.4,1.6]
    rng=(-1,1)
    x = np.linspace(-1.1,1.1,100000)
    fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
    fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
    for epsilon in epsilons:
        D=smooth_delta(x, e=epsilon, rng=rng, method=method )
        print('Epsilon {:8.5f}  integral:{}'.format(epsilon, np.trapz(D,x)))
        ax.plot(x, D/np.max(D), '-', label=f'e={epsilon}')
    ax.set_xlabel('x ')
    ax.set_ylabel('')
    #ax.set_yscale('log')
    ax.legend()
    ax.tick_params(direction='in')
    ax.set_title(f'Smooth Dirac - method {method} - different epsilon')

    # --- Smooth Dirac Sample
    method='frac'
    epsilon=1.6
    rng=(-1,1)
    x = np.linspace(-1,1,100001)
    nSamples=[3,5,7]
    fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
    fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
    D=smooth_delta(x, e=epsilon, rng=rng, method=method )
    vmax=np.max(D)
    ax.plot(x, D/vmax, 'k-')
    for n in nSamples:
        x2 = np.linspace(-1, 1, n)
        D2 = smooth_delta(x2, e=epsilon, rng=rng, method=method )
        print('n {}: values: {}'.format(n, D2/vmax))
        ax.plot(x2, D2/vmax, '-o', label= f'n={n}')
    ax.set_xlabel('x ')
    ax.set_ylabel('')
    ax.legend()
    ax.tick_params(direction='in')
    ax.set_title(f'Smooth Dirac - method {method} - different number of points')


    #from welib.tools.figure import *
    #export2png()

    plt.show()

    

