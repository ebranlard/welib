""" 
pdf: probability density function

Reference:
   [1] 

"""


import numpy as np

# --------------------------------------------------------------------------------
# --- Uniform
# --------------------------------------------------------------------------------
def uniform_pdf(x, a, b):
    f = np.zeros_like(x)
    bn = np.logical_and(x>=a,x<=b)
    f[bn]=1/(b-a)
    return f

# --------------------------------------------------------------------------------
# --- Normal/Gaussian 
# --------------------------------------------------------------------------------
def gaussian_pdf(x, mu, sig):
    return  1/( sig * np.sqrt(2*np.pi) ) * np.exp(-0.5 * (x-mu)**2/sig**2)


def gaussian2d_covar(sigmas, rho):
    C = np.array(
       [[ sigmas[0]**2             , sigmas[0]*sigmas[1]*rho ],
        [ sigmas[0]*sigmas[1]*rho  , sigmas[1]**2           ]] )
    return C


def gaussian2d_pdf(x1, x2, mus, sigmas, rho=0):
    """ """
    scale = 1/(2*np.pi*sigmas[0]*sigmas[1]*np.sqrt(1-rho**2))
    xi1 = (x1 - mus[0])/sigmas[0]
    xi2 = (x2 - mus[1])/sigmas[1]
    f = scale * np.exp(-0.5* (xi1**2-2*rho*xi1*xi2+xi2**2)/(1-rho**2) )
    return f


def gaussian2d_pdf_meshgrid(x1, x2, mus, sigmas, rho=0, method='analytical'):
    """ """
    C = gaussian2d_covar(sigmas, rho)

    X1, X2 = np.meshgrid(x1,x2)
    if method=='analytical':
        F = gaussian2d_pdf(X1, X2, mus, sigmas, rho=rho)
    else:
        F = np.zeros(X1.shape)
        for i1,xx1 in enumerate(x1):
            for i2,xx2 in enumerate(x2):
                x = np.array([xx1,xx2])
                ff=gaussianNd_pdf(x, mus, C)
                F[i2,i1] = ff
    return X1, X2, F

def gaussianNd_pdf(x, mu, C):
    """ 
    Probability density function of a normally distributed stochastic vector of n-dimension.

    x: array of shape (ndim x nsamp) where the pdf is sampled
    mu: means, array of length ndim
    C: Covariance matrix, ndim x ndim
    """
    x = np.asarray(x)
    if not hasattr(mu, '__len__'):
        # Standard 1d Gaussian
        f = gaussian_pdf(x, mu, np.sqrt(C))
        return f
    else:
        ndim = len(mu)
        assert(C.shape == (ndim, ndim)) 
        assert(x.shape[0] == (ndim)) 
        Cdet = np.linalg.det(C)
        Cinv = np.linalg.inv(C)
        if len(x.shape)==1:
            # nsamp =1
            xp = (x.T-mu).T
            f = 1/( (2*np.pi)**(ndim/2) * np.sqrt(Cdet) ) * np.exp(-0.5 * (xp.T).dot(Cinv.dot(xp)))
        else:
            f = np.zeros(x.shape)
            Xp = (x.T-mu)
            for xp in Xp:
                f = 1/( (2*np.pi)**(ndim/2) * np.sqrt(Cdet) ) * np.exp(-0.5 * (xp.T).dot(Cinv.dot(xp)))
        return f


# --------------------------------------------------------------------------------
# --- Rayleigh 
# --------------------------------------------------------------------------------
def rayleigh_pdf(x, sig):
    x = np.asarray(x)
    f = np.zeros_like(x)
    bn = x>=0
    x = x[bn]
    f[bn]= x/sig**2 * np.exp(-x**2/(2*sig**2))
    return f


# --------------------------------------------------------------------------------
# --- Gamma 
# --------------------------------------------------------------------------------
def gamma_pdf(x, alpha, beta):
    import math
    x = np.asarray(x)
    f = np.zeros_like(x)
    bn = x>=0
    x = x[bn]
    f[bn]= beta**(alpha+1)/math.gamma(alpha+1) * x**alpha * np.exp(-beta*x)
    return f


# --------------------------------------------------------------------------------
# --- Exponential 
# --------------------------------------------------------------------------------
def exponential_pdf(x, beta):
    x = np.asarray(x)
    f = np.zeros_like(x)
    bn = x>=0
    x = x[bn]
    f[bn]= beta * np.exp(-beta * x)
    return f


# --------------------------------------------------------------------------------
# --- Weibull
# --------------------------------------------------------------------------------
def weibull_pdf(x, A, k, x0):
    x = np.asarray(x)
    f = np.zeros_like(x)
    bn = x>=x0
    x=x[bn]-x0
    f[bn]= k/A * (x/A)**(k-1) * np.exp(-(x/A)**k)
    return f

# --------------------------------------------------------------------------------
# --- Lognormal
# --------------------------------------------------------------------------------
def lognormal_pdf(x, mu, sig):
    x = np.asarray(x)
    f = np.zeros_like(x)
    bn = x>0
    x=x[bn]
    f[bn]= 1/(np.sqrt(2*np.pi) * sig) * 1/x * np.exp( - ( np.log(x) - mu )**2 / (2*sig**2))
    return f




if __name__ == '__main__':
    pass

