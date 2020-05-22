""" 3rd edition of IEC standard 61400-1 """
def fVref(WT_class='I'):
    if WT_class== 'I':
          Vref = 50.0
    elif WT_class== 'II':
          Vref = 42.5
    elif WT_class=='III':
        Vref = 37.5
    else:
        raise Exception('Unknown wind turbine class '+IEC_class)
    return Vref


def fIref(IEC_class='A'):
    if IEC_class=='A':
        Iref=0.16
    elif IEC_class=='B':
        Iref=0.14
    elif IEC_class=='C':
        Iref=0.12
    else:
        raise Exception('Unknown class '+IEC_class)
    return Iref

def sigma1_ETM(WS, IEC_class='A'): 
    # "sigma1" from extreme turbulence model
    Iref = fIref(IEC_class)
    c    = 2
    Vave = 10
    return  c * Iref * (0.072 * (Vave / c + 3) * (WS / c - 4) + 10)

def sigma1_NTM(WS, IEC_class='A'):
    # Normal turbulence model
    # Function sigma, offshore
    #   sigma1 = fE_sigma(U_10_hub) + 1.28*fD_sigma(U_10_hub)
    Iref = fIref(IEC_class)
    return Iref * (0.75 * WS + 5.6) 
    
def ETM(WS, IEC_class='A'): 
    # Turbulence intensity from extreme turbulence model
    return sigma1_ETM(WS, IEC_class)* WS

def NTM(WS, IEC_class='A'):
    # Turbulence intensity from normal turbulence model
    return sigma1_ETM(WS, IEC_class)* WS

def Lambda(zhub, IECedition=3):
    """
    IEC length scale.  Lambda = 0.7*min(Zhat,zhub)
    Where: Zhat = 30,60 for IECedition = 2,3, respectively.
    """
    if IECedition <= 2:
        return 0.7 * min(30, zhub)
    return 0.7 * min(60, zhub)


def KaimalSpectrum(f, zhub, uhub, Model='ETM', IEC_class='A' ):
    """IEC Kaimal spectrum model.
    See IEC 61400-3 Appendix B.2
        S_k(f) = \frac{4 \sigma_k^2 L_k/V }{(1+6 f L_k/V)^{5/3}} 
        k=0,1,2 = longi, lateral, upward
    """
    if Model == 'ETM':
        sigma1 = sigma1_ETM(uhub, IEC_class)
    elif Model == 'NTM':
        sigma1 = sigma1_NTM(uhub, IEC_class)

    S_k = np.zeros(3,len(f))
    sigma_k = sigma1       * np.array([1.0,  0.8,  0.5 ])
    L_k     = Lambda(zhub) * np.array([8.10, 2.70, 0.66])
    for k in [0,1,2]: # loop on components
        S_k[k] = (4 * sigma[k]**2 * L_k[k]/uhub)/( 1 + 6 * f * L_k[k]/uhub)**(5/3)
    return S_k

def VonKarmanSpectrum(f, zhub, uhub, Model='ETM', IEC_class='A'):
    r"""IEC Von-Karman spectral model
    The form of this model is,
    .. math::
          S_u(f) = \frac{4 \sigma^2/\hat{f}}{(1+71(f/\hat{f})^2)^{5/6}}
          S_v(f) = S_w(f) = (1+189(f/\hat{f})^2\frac{2\sigma^2/\hat{f}}{(1+71 (f/\hat{f})^2)^{11/6}}
    Where,
      :math:`\hat{f} = \bar{u}_{hub}/\Lambda`
    """
    if Model == 'ETM':
        sigma1 = sigma1_ETM(uhub, IEC_class)
    elif Model == 'NTM':
        sigma1 = sigma1_NTM(uhub, IEC_class)
    
    spec    = np.zeros(3,len(f))
    sig2 = 4 * sigma1 ** 2 
    L_u  = 3.5 * Lambda(zhub)/uhub
    dnm = 1 + 71 * (f * L_u) ** 2 # TODO 71 or 70.8...?
    S_k[0] =     sig2 * L_u/(dnm)**(5./6)
    S_k[1] = sig2 / 2 * L_u/(dnm)**(11./6) * (1 + 189 * (f * L_u) ** 2)
    S_k[2] = S_k[1]
    return spec


if __name__=='__main__':
    import numpy as np
    import matplotlib.pyplot as plt
    ws=np.linspace(0,30,100)
    fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
    fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
    ax.plot(ws,ETM(ws,'A')*100, '-', label='ETM A', c='k')
    ax.plot(ws,ETM(ws,'B')*100, '--', label='ETM B', c='k')
    ax.plot(ws,ETM(ws,'C')*100, '.-', label='ETM C', c='k')
    ax.plot(ws,NTM(ws,'A')*100, '-', label='NTM A',c=[0.5,0.5,0.5])
    ax.plot(ws,NTM(ws,'B')*100, '--', label='NTM B',c=[0.5,0.5,0.5])
    ax.plot(ws,NTM(ws,'C')*100, '-.', label='NTM C',c=[0.5,0.5,0.5])
    ax.set_ylim([0,88])
    ax.set_xlim([0,30])
    ax.set_xlabel('Wind speed [m/s]')
    ax.set_ylabel('Turbulence intensity [%]')
    ax.legend()
    ax.tick_params(direction='in')
    plt.show()
