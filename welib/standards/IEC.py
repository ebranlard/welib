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
    
def ETM(WS, IEC_class='A'): 
    # Extreme turbulence model

    Iref = fIref(IEC_class)
    c    = 2
    Vave = 10
    TI = (c * Iref * (0.072 * (Vave / c + 3) * (WS / c - 4) + 10)) / WS
    
    return TI

def NTM(WS, IEC_class='A'):
    # Normal turbulence model
    #Function sigma, offshore
    # sigma1 = fE_sigma(U_10_hub) + 1.28*fD_sigma(U_10_hub)

    Iref = fIref(IEC_class)
    TI = Iref * (0.75 * WS + 5.6) / WS
    return TI
