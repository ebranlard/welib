""" 

Reference:
   [1] 

"""


import numpy as np



def gaussian_pdf(x, mu, C):
    """ 
    Probability density function of a normally distributed stochastic vector of n-dimension.

    x: array of shape (ndim x nsamp) where the pdf is sampled
    mu: means, array of length ndim
    C: Covariance matrix, ndim x ndim
    """
    x = np.asarray(x)
    if not hasattr(mu, '__len__'):
        # Standard 1d Gaussian
        ndim = 1
        var = C
        sig = np.sqrt(C)
        f = 1/( sig * np.sqrt(2*np.pi) ) * np.exp(-0.5 * (x-mu)**2/var)
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

def gaussian2d_pdf(x1, x2, mus, sigmas, rho=0):
    """ """
    scale = 1/(2*np.pi*sigmas[0]*sigmas[1]*np.sqrt(1-rho**2))
    xi1 = (x1 - mus[0])/sigmas[0]
    xi2 = (x2 - mus[1])/sigmas[1]
    f = scale * np.exp(-0.5* (xi1**2-2*rho*xi1*xi2+xi2**2)/(1-rho**2) )
    return f


def gaussian2d_pdf_meshgrid(x1, x2, mus, sigmas, rho=0, method='analytical'):
    """ """
    C = np.array(
       [[ sigmas[0]**2             , sigmas[0]*sigmas[1]*rho ],
        [ sigmas[0]*sigmas[1]*rho  , sigmas[1]**2           ]] )

    X1, X2 = np.meshgrid(x1,x2)
    if method=='analytical':
        F = gaussian2d_pdf(X1, X2, mus, sigmas, rho=rho)
    else:
        F = np.zeros(X1.shape)
        for i1,xx1 in enumerate(x1):
            for i2,xx2 in enumerate(x2):
                x = np.array([xx1,xx2])
                ff=gaussian_pdf(x, mus, C)
                F[i2,i1] = ff
    return X1, X2, F




if __name__ == '__main__':
    sigmas = [5,6]
    mus    = [0,3]
    rho=0.7

    # --- joint probability density function
    x1 = np.linspace(-20,20,100)
    x2 = np.linspace(-20,20,80)
    X1, X2, F     = gaussian2d_pdf_meshgrid(x1, x2, mus, sigmas, rho)
    _, _, F_gen = gaussian2d_pdf_meshgrid(x1, x2, mus, sigmas, rho, method='general')
    print('Error norm:', np.linalg.norm(F-F_gen))

    # --- Marginal distributions
    F1 = np.trapz(F.T, x2)
    F2 = np.trapz(F, x1)
    F1th = gaussian_pdf(x1, mus[0], sigmas[0]**2)
    F2th = gaussian_pdf(x2, mus[1], sigmas[1]**2)

    # Joint disitribution if variables were independent
    _, _, F_indpt = gaussian2d_pdf_meshgrid(x1, x2, mus, sigmas, rho=0)
#     F_indpt = (np.atleast_2d(F1th).T).dot(np.atleast_2d(F2th)).T
    F1_indpt = np.trapz(F_indpt.T, x2)
    F2_indpt = np.trapz(F_indpt, x1)


    # --- Contour
    import matplotlib.pyplot as plt
    fig,axes = plt.subplots(2, 2, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
    ax = axes[0,0]
    ax.contourf(X1,X2,F)
    ax.plot(mus[0], mus[1], 'ko')
    #ax.plot(x1, mus[1]+(x1-mus[0])*rho, 'k--') # need influence of sigma
    ax.axis('equal')
    ax.set_title('Correlated')

    ax = axes[0,1]
    ax.contourf(X1,X2,F_indpt)
    ax.plot(mus[0], mus[1], 'ko')
    ax.axis('equal')
    ax.set_title('Independent')


    ax = axes[1,0]
    ax.plot(x1, F1th, 'k-')
    ax.plot(x2, F2th, 'k-')
    ax.plot(x1, F1, '--', label='F1')
    ax.plot(x2, F2, '--', label='F2')
    ax.set_xlabel('x')
    ax.set_ylabel('')
    ax.legend()
    

    ax = axes[1,1]
    ax.plot(x1, F1th, 'k-')
    ax.plot(x2, F2th, 'k-')
    ax.plot(x1, F1_indpt, '--', label='F1')
    ax.plot(x2, F2_indpt, '--', label='F2')
    ax.set_xlabel('x')
    ax.set_ylabel('')
    ax.legend()
    

    plt.show()



