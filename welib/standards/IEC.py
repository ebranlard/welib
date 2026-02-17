""" 3rd edition of IEC standard 61400-1 """
import numpy as np

# --------------------------------------------------------------------------------}
# --- Reference values 
# --------------------------------------------------------------------------------{
def fVref(WT_class='I'):
    """ Reference velocity """
    if WT_class== 'I':
        Vref = 50.0
    elif WT_class== 'II':
        Vref = 42.5
    elif WT_class=='III':
        Vref = 37.5
    elif WT_class == 'IV':
        V_ref = 30.
    else:
        raise Exception('Unknown wind turbine class '+turbulence_class)
    return Vref

def fIref(turbulence_class='A'):
    """ Reference turbulence intensity """
    if turbulence_class=='A+':
        Iref=0.18
    elif turbulence_class=='A':
        Iref=0.16
    elif turbulence_class=='B':
        Iref=0.14
    elif turbulence_class=='C':
        Iref=0.12
    else:
        raise Exception('Unknown class '+turbulence_class)
    return Iref


def Lambda(zhub, IECedition=3):
    """
    IEC length scale.  Lambda = 0.7*min(Zhat,zhub)
    Where: Zhat = 30,60 for IECedition = 2,3, respectively.
    """
    if IECedition <= 2:
        return 0.7 * min(30, zhub)
    return 0.7 * min(60, zhub)

# --------------------------------------------------------------------------------}
# --- Turbulence models 
# --------------------------------------------------------------------------------{
def sigma1_ETM(WS, turbulence_class='A'): 
    """ "sigma1" from extreme turbulence model """
    Iref = fIref(turbulence_class)
    c    = 2.
    Vave = 10.
    return  c * Iref * (0.072 * (Vave / c + 3.) * (WS / c - 4.) + 10.)

def sigma1_NTM(WS, turbulence_class='A'):
    # Normal turbulence model
    # Function sigma, offshore
    #   sigma1 = fE_sigma(U_10_hub) + 1.28*fD_sigma(U_10_hub)
    Iref = fIref(turbulence_class)
    return Iref * (0.75 * WS + 5.6) 

def Sigma1_zHub(z_hub): 
    if z_hub > 60:
       return 42
    else:
        return 0.7*z_hub
    
def ETM(WS, turbulence_class='A'): 
    # Turbulence intensity from extreme turbulence model
    return sigma1_ETM(WS, turbulence_class)/ WS

def NTM(WS, turbulence_class='A'):
    # Turbulence intensity from normal turbulence model
    return sigma1_NTM(WS, turbulence_class)/ WS

def NTM_TI(WS, turb_class='A'):
    """ Normal Turbulence intensity model according to the IEC standards
    INPUT: 
       -  WS: array of wind speed [m/s]
       -  turb_class: turbine class, "A", "B" or "C"
    OUTPUT:
       - TI: turbulence intensity at WS
    """
    Iref={'A+':0.18, 'A':0.16, 'B':0.14, 'C':0.12}[turb_class]
    c    = 2.
    Vave = 10.
    sigma =  c * Iref * (0.072 * (Vave / c + 3.) * (WS / c - 4.) + 10.)
    TI = sigma/WS
    return TI


def EWM(WS, WT_class='I'):
    # Extreme wind speed model: 6.3.2.1
    V_ref = fVref(WT_class=WT_class)
    # Steady
    V_e50 = 1.4*V_ref
    V_e1  = 0.8*V_e50
    # Turbulent
    V_50    = V_ref
    V_1     = 0.8*V_50
    sigma_1 = 0.11*WS
    return sigma_1, V_e50, V_e1, V_50, V_1


# --------------------------------------------------------------------------------}
# --- Gusts 
# --------------------------------------------------------------------------------{
def EOG(WS, D, z_hub, WT_class='I',turbulence_class='A', vert_slope=0, tStart=0, filename=None, shearExp=None):
    """
    D    : rotor diameter [m]
    z_hub: hub height [m]
    WS   : hub velocity [m/s]
    vert_slope: vertical slope [deg]
    """
    # Extreme operating gust: 6.3.2.2
    T  = 10.5
    dt = T/100
    t  = np.linspace(tStart, tStart+T, int((T/dt)+1))

    # constants from standard
    if shearExp is None:
        alpha = 0.2
    else:
        alpha = shearExp

    # Flow angle adjustments
    V_hub      = WS*np.cos(vert_slope*np.pi/180)
    V_vert_mag = WS*np.sin(vert_slope*np.pi/180)

    sigma_1 = NTM(V_hub, turbulence_class=turbulence_class)
    __, __, V_e1, __, __ = EWM(V_hub, WT_class=WT_class)

    # Contant variables
    V              = np.zeros_like(t)+V_hub
    V_dir          = np.zeros_like(t)
    V_vert         = np.zeros_like(t)+V_vert_mag
    shear_horz     = np.zeros_like(t)
    shear_vert     = np.zeros_like(t)+alpha
    shear_vert_lin = np.zeros_like(t)
    V_gust         = np.zeros_like(t)

    V_gust = min([ 1.35*(V_e1 - V_hub), 3.3*(sigma_1/(1+0.1*(D/Sigma1_zHub(z_hub)))) ])

    V_gust_t = np.zeros_like(t)
    V_gust_t = 0. - 0.37*V_gust*np.sin(3*np.pi*(t-tStart)/T)*(1-np.cos(2*np.pi*(t-tStart)/T))
    V_gust_t[t<tStart] = 0
    V_gust_t[t>tStart+T] = 0

    data = np.column_stack((t, V, V_dir, V_vert, shear_horz, shear_vert, shear_vert_lin, V_gust_t))

    if filename is not None:
        from welib.weio.fast_wind_file import FASTWndFile
        import pandas as pd
        f = FASTWndFile()
        f.header[0] = '! EOG wind file for OpenFAST. V_hub0 = {:.2f}m/s - sigma_1: {:.2f}m/s - vert_slope={:.2f}deg'.format(V_hub, sigma_1, vert_slope)
        f.data = pd.DataFrame(data = data, columns=f.colNames)
        f.write(filename)

    return t, V_gust_t+WS, V_vert_mag, shear_vert


    

# --------------------------------------------------------------------------------}
# --- Turbulence spectra 
# --------------------------------------------------------------------------------{


def KaimalSpectrum(f, zhub, uhub, Model='ETM', turbulence_class='A' ):
    """IEC Kaimal spectrum model.
    See IEC 61400-3 Appendix B.2
        S_k(f) = \frac{4 \sigma_k^2 L_k/V }{(1+6 f L_k/V)^{5/3}} 
        k=0,1,2 = longi, lateral, upward
    """
    if Model == 'ETM':
        sigma1 = sigma1_ETM(uhub, turbulence_class)
    elif Model == 'NTM':
        sigma1 = sigma1_NTM(uhub, turbulence_class)

    S_k = np.zeros(3,len(f))
    sigma_k = sigma1       * np.array([1.0,  0.8,  0.5 ])
    L_k     = Lambda(zhub) * np.array([8.10, 2.70, 0.66])
    for k in [0,1,2]: # loop on components
        S_k[k] = (4 * sigma[k]**2 * L_k[k]/uhub)/( 1 + 6 * f * L_k[k]/uhub)**(5/3)
    return S_k

def VonKarmanSpectrum(f, zhub, uhub, Model='ETM', turbulence_class='A'):
    r"""IEC Von-Karman spectral model
    The form of this model is,
    .. math::
          S_u(f) = \frac{4 \sigma^2/\hat{f}}{(1+71(f/\hat{f})^2)^{5/6}}
          S_v(f) = S_w(f) = (1+189(f/\hat{f})^2\frac{2\sigma^2/\hat{f}}{(1+71 (f/\hat{f})^2)^{11/6}}
    Where,
      :math:`\hat{f} = \bar{u}_{hub}/\Lambda`
    """
    if Model == 'ETM':
        sigma1 = sigma1_ETM(uhub, turbulence_class)
    elif Model == 'NTM':
        sigma1 = sigma1_NTM(uhub, turbulence_class)
    
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
