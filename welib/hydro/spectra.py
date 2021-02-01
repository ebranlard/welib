import numpy as np

def jonswap(freq, Hs, Tp=None, Tz=None, g=9.81):
    """
    Jonswap spectrum for a set of frequencies and spectral parameters.

    INPUTS:
      - freq : frequencies [Hz], scalar or array, typically freq=np.arange(df, fMax, df)
      - Hs : significant wave height [m], scalar
      either:
      - Tp  : peak period [s], scalar
      or:
      - Tz  : average zero crossing period [s], scalar
      - g: acceleration of gravity

    OUTPUTS:
      - S : spectral amplitude
    
    """
    def jonswap_gamma(Hs,Tp):
        """ return peak enhancement factor gamma"""
        p=Tp/np.sqrt(Hs);
        if p<=3.6:
            return 5
        elif  3.6<p and p<=5:
            return np.exp(5.75-1.15*p)
        else:
            return 1

    if Tp is not None:
        gamma = jonswap_gamma(Hs,Tp)
    else:
        # Iterate to find Tp and gamma based on Tz
        gamma0  = 1     # initial guess
        tol_rel = 0.005
        tol_abs = 0.1
        it      = 0
        gamma = None
        while it<1000:
            it += 1
            Tp      = Tz*np.sqrt((11+gamma1)/(5+gamma0))
            gamma1  = gamma_jonswap(Hs,Tp)
            err_abs = abs(gamma1-gamma0)
            err_rel = err_abs/gamma0 
            if err_abs<tol_abs and err_rel<tol_rel:
                gamma=gamma1
                break
            gamma0=gamma0+(gamma1-gamma0)/2
        if gamma is None:
            Exception('Maximum number of iteration reached while searching Tp and Gamma')

    freq = np.asarray(freq)
    fp = 1/Tp
    sigma = np.ones(freq.shape)*0.09
    sigma[freq<=fp]=0.7
    beta  = np.exp(-0.5*(((freq/fp)-1.)*(1./sigma))**2)
    alpha = 5*(Hs**2*fp**4/g**2)*(1-0.287*np.log(gamma))*np.pi**4
    S     = alpha*g**2/(2*np.pi)**4*freq**(-5)*np.exp(-1.25*(freq/fp)**(-4))*gamma**beta
    # Alternative expression for g=9.81)
    # S = 0.3125*Hs**2*Tp*((f/fp)**(-5))* np.exp(-1.25*(f/fp)**(-4))* (1-0.287*np.log(gamma))* gamma**beta
    return S
        
