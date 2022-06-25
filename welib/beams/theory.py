import numpy as np
import scipy.optimize as sciopt

def UniformBeamBendingModes(Type,EI,rho,A,L,w=None,x=None,Mtop=0,norm='tip',nModes=4):
    """
    returns Mode shapes and frequencies for a uniform beam in bending

    References:
      [1] Inman : Engineering vibration
      [2] Nielsen
      [3] Thomsen: Vibration and stability
    
    Author: E. Branlard"""
    if x is None or len(x)==0:
        x = np.linspace(0,L,101)
    if np.amax(x) != L:
        raise Exception('Max of x should be equal to L')

    # Useful Intermediate functions [3] p374
    F = lambda u: np.sinh(u) + np.sin(u)
    G = lambda u: np.cosh(u) + np.cos(u)
    H = lambda u: np.sinh(u) - np.sin(u)
    J = lambda u: np.cosh(u) - np.cos(u)
    Jp  = F
    Jpp = G
    Hp  = J
    Hpp = F
    Gp  = H
    Gpp = J
    Fp  = G
    Fpp = H


    # Dimensionless spanwise position
    x0 = x / L
    Type = Type.lower().replace('fixed','clamped')
    s = Type.split('-')
    if s[0] =='unloaded':

        modesU = lambda x0, l: x0*0
        modesV = lambda x0, l: x0*0
        modesK = lambda x0, l: x0*0

        # See [1] p.335  [2] p129 [3] p375
        if   'unloaded-clamped-clamped' == Type:
            freq_function = lambda x: -1 + np.cosh(x) * np.cos(x)
            freq_guess    = lambda i: (2*(i+1)+1)*np.pi/2
            modesU        = lambda x0, l: H(l) * J  (l*x0) - J(l) * H  (l*x0)
            modesV        = lambda x0, l: H(l) * Jp (l*x0) - J(l) * Hp (l*x0)
            modesK        = lambda x0, l: H(l) * Jpp(l*x0) - J(l) * Hpp(l*x0)
        elif 'unloaded-clamped-hinged' == Type:
            freq_function = lambda x: np.tan(x) - np.tanh(x)
            freq_guess    = lambda i: (4*(i+1)+1)*np.pi/4
            modesU        = lambda x0, l: H(l) * J  (l*x0) - J(l) * H  (l*x0)
            modesV        = lambda x0, l: H(l) * Jp (l*x0) - J(l) * Hp (l*x0)
            modesK        = lambda x0, l: H(l) * Jpp(l*x0) - J(l) * Hpp(l*x0)
        elif 'unloaded-clamped-free' == Type:
            # NOTE: cosh(beta_n)cos(beta_n) =-1
            #  sigma_n = [ np.sinh(beta_n) - sin(beta_n) ]/[cosh(beta_n) + cos(beta_n)]
            #  for j>5, a good approx is B(j) = (2*j-1)np.pi/2  and S(j)=1;
            #B  = [1.87510407, 4.69409113, 7.85475744,10.99554073,14.13716839, (2*6-1)*np.pi/2];
            #S  = [0.734095514 1.018467319 0.999224497 1.000033553 0.999998550 1];
            freq_function = lambda x: 1 + np.cosh(x) * np.cos(x)
            freq_guess    = lambda i: (2*(i+1)-1)*np.pi/2
            modesU        = lambda x0, l: F(l) * J  (l*x0) - G(l) * H  (l*x0)
            modesV        = lambda x0, l: F(l) * Jp (l*x0) - G(l) * Hp (l*x0)
            modesK        = lambda x0, l: F(l) * Jpp(l*x0) - G(l) * Hpp(l*x0)
        elif 'unloaded-clamped-guided' == Type:
            freq_function = lambda x: np.tan(x) + np.tanh(x)
            freq_guess    = lambda i: (4*(i+1)-1)*np.pi/4
            modesU        = lambda x0, l: J(l) * J  (l*x0) - F(l) * H  (l*x0)
            modesV        = lambda x0, l: J(l) * Jp (l*x0) - F(l) * Hp (l*x0)
            modesK        = lambda x0, l: J(l) * Jpp(l*x0) - F(l) * Hpp(l*x0)
        elif 'unloaded-topmass-clamped-free' == Type:
            # The geometrical stiffning is not accounted for here
            if Mtop is None:
                raise Exception('Please specify value for Mtop for %s',Type)
            M = rho * A * L
            freq_function = lambda x: 1+np.cosh(x)*np.cos(x)-x*Mtop/M*(np.sin(x)*np.cosh(x)-np.cos(x)*np.sinh(x))
            freq_guess    = lambda i: (2*(i+1)-1)*np.pi/2
            modesU        = lambda x0, l: F(l) * J  (l*x0) - G(l) * H  (l*x0)
            modesV        = lambda x0, l: F(l) * Jp (l*x0) - G(l) * Hp (l*x0)
            modesK        = lambda x0, l: F(l) * Jpp(l*x0) - G(l) * Hpp(l*x0)
        elif 'unloaded-hinged-hinged' == Type: # simply-supported
            freq_function = lambda x: np.sin(x) 
            freq_guess    = lambda i: (i+1)*np.pi # NOTE: exact..
            modesU        = lambda x0, l:  np.sin(l*x0)
            modesV        = lambda x0, l:  np.cos(l*x0)
            modesK        = lambda x0, l: -np.sin(l*x0)
        elif 'unloaded-hinged-guided' == Type:
            freq_function = lambda x: np.cos(x) 
            freq_guess    = lambda i: (2*(i+1)-1)*np.pi/2 # NOTE: exact..
            modesU        = lambda x0, l:  np.sin(l*x0)
            modesV        = lambda x0, l:  np.cos(l*x0)
            modesK        = lambda x0, l: -np.sin(l*x0)
        elif 'unloaded-guided-guided' == Type:
            freq_function = lambda x: np.sin(x) 
            freq_guess    = lambda i: (i+1)*np.pi  # NOTE: exact..
            modesU        = lambda x0, l:  np.cos(l*x0)
            modesV        = lambda x0, l: -np.sin(l*x0)
            modesK        = lambda x0, l: -np.cos(l*x0)
        elif 'unloaded-free-free' == Type:
            freq_function = lambda x: -1 + np.cosh(x) * np.cos(x)
            freq_guess    = lambda i: (2*(i+1)+1)*np.pi/2
            modesU        = lambda x0, l: H(l) * G  (l*x0) - J(l) * F  (l*x0)
            modesV        = lambda x0, l: H(l) * Gp (l*x0) - J(l) * Fp (l*x0)
            modesK        = lambda x0, l: H(l) * Gpp(l*x0) - J(l) * Fpp(l*x0)
        elif 'unloaded-free-hinged' == Type:
            freq_function = lambda x: np.tan(x) - np.tanh(x)
            freq_guess    = lambda i: (4*(i+1)+1)*np.pi/4
            modesU        = lambda x0, l: F(l) * G  (l*x0) - G(l) * F  (l*x0)
            modesV        = lambda x0, l: F(l) * Gp (l*x0) - G(l) * Fp (l*x0)
            modesK        = lambda x0, l: F(l) * Gpp(l*x0) - G(l) * Fpp(l*x0)
        elif 'unloaded-free-guided' == Type:
            freq_function = lambda x: np.tan(x) + np.tanh(x)
            freq_guess    = lambda i: (4*(i+1)-1)*np.pi/4
            modesU        = lambda x0, l: F(l) * G  (l*x0) - H(l) * F  (l*x0)
            modesV        = lambda x0, l: F(l) * Gp (l*x0) - H(l) * Fp (l*x0)
            modesK        = lambda x0, l: F(l) * Gpp(l*x0) - H(l) * Fpp(l*x0)
        else:
            raise Exception('unknown type %s',Type)

        # Solve for "lambda"
        B = np.zeros(nModes)
        for i in np.arange(nModes):
            B[i] = sciopt.fsolve(freq_function, freq_guess(i))
        # Frequency
        freq = (B / L) ** 2 / (2 * np.pi) * np.sqrt(EI / (rho * A))
        # --- Mode shapes
        ModesU = np.zeros((len(B),len(x0)))
        ModesV = np.zeros((len(B),len(x0)))
        ModesK = np.zeros((len(B),len(x0)))
        for i in np.arange(nModes):
            l = B[i]
            ModesU[i,:] = modesU(x0,l)
            ModesV[i,:] = modesV(x0,l) * l
            ModesK[i,:] = modesK(x0,l) * l**2

    elif s[0] == 'loaded':
        if 'loaded-clamped-free' == Type:
            if w is None:
                w = A * rho
            if L==0:
                raise Exception('Please specify value for L for %s',Type)
            B = np.array([1.875,4.694])
            freq = (B / L) ** 2 / (2 * np.pi) * np.sqrt(EI / w)
        else:
            raise Exception('unknown type %s',Type)
    else:
        raise Exception('Unknown %s'^Type)
    ## Going back to physical dimension
    x = x0 * L
    ModesV = ModesV/L
    ModesK = ModesK/L**2
    ## Normalization of modes
    if norm=='tip':
        for i in np.arange(nModes):
            tipVal = ModesU[i,-1]
            fact = 1 / tipVal
            ModesU[i,:] = ModesU[i,:] * fact
            ModesV[i,:] = ModesV[i,:] * fact
            ModesK[i,:] = ModesK[i,:] * fact
    elif norm=='max':
        for i in np.arange(nModes):
            maxVal = np.max(np.abs(ModesU[i,:]))
            iMaxVal = np.argmax(np.abs(ModesU[i,:]))
            fact = 1 / ModesU[i, iMaxVal]
            ModesU[i,:] = ModesU[i,:] * fact
            ModesV[i,:] = ModesV[i,:] * fact
            ModesK[i,:] = ModesK[i,:] * fact
    else:
        raise Exception('Norm not implemented or incorrect: `%s`'%norm)

    return freq,x,ModesU,ModesV,ModesK


    
# --------------------------------------------------------------------------------}
# --- Longitudinal modes 
# --------------------------------------------------------------------------------{
def UniformBeamLongiModes(Type,E,rho,A,L,x=None,nModes=4,norm='tip_norm'):
    """
    Returns longitudinals modes for a uniform beam
    """
    if x is None:
        x = np.linspace(0,L,101)
    if np.amax(x) != L:
        raise Exception('Max of x should be equal to L')
    # Dimensionless spanwise position
    x0 = x / L
    if Type.lower()=='unloaded-clamped-free':
        c = np.sqrt(E / rho)
        freq = np.zeros(nModes)
        ModesU = np.full([nModes,len(x0)],np.nan)
        for j in np.arange(nModes):
            omega_j     = c/L*(np.pi/2 + j*np.pi)
            freq[j]     = omega_j/(2*np.pi)
            ModesU[j,:] = np.sin(omega_j/c*x)
    else:
        raise Exception('Unknown %s'%Type)
    
    ## Computation of derivatives if no analytical functions
    # V=fgradient_regular(U(i,:),4,dx);
    # K=fgradient_regular(V(i,:),4,dx);
    ## Going back to physical dimension
    # x=x0*L;
    # ModesV = ModesV/L;
    # ModesK = ModesK/L^2;
    
    ## Normalization of modes
    if norm=='tip_norm':
        for i in np.arange(nModes):
            fact = 1 / ModesU[i,-1]
            ModesU[i,:] = ModesU[i,:] * fact
    
    return freq,x,ModesU

def UniformBeamTorsionModes(Type,G,Kt,Ip,rho,A,L,x=None,nModes=4,norm='tip_norm'):
    """
    Returns torsional modes for a uniform beam
    """
    if x is None:
        x = np.linspace(0,L,101)
    if np.amax(x) != L:
        raise Exception('Max of x should be equal to L')
    # Dimensionless spanwise position
    x0 = x / L
    if Type.lower()=='unloaded-clamped-free':
        c = np.sqrt(G*Kt/(rho*Ip))
        freq = np.zeros(nModes)
        ModesV = np.full([nModes,len(x0)],np.nan)
        # NOTE: equations are the same as longi
        for j in np.arange(nModes):
            omega_j     = c/L*(np.pi/2 + j*np.pi)
            freq[j]     = omega_j/(2*np.pi)
            ModesV[j,:] = np.sin(omega_j/c*x)
    else:
        raise Exception('Unknown %s'%Type)
    
    ## Normalization of modes
    if norm=='tip_norm':
        for i in np.arange(nModes):
            fact = 1 / ModesV[i,-1]
            ModesV[i,:] = ModesV[i,:] * fact
    
    ## Computation of derivatives if no analytical functions
    #if exist('fgradient_regular'):
    #    dx = x(2) - x(1)
    #    ModesK = np.zeros((ModesV.shape,ModesV.shape))
    #    for i in np.arange(1,ModesV.shape[1-1]+1).reshape(-1):
    #        ModesK[i,:-1] = fgradient_regular(ModesV(i,:),4,dx)
    ModesK = []
    
    return freq,x,ModesV,ModesK


if __name__=='__main__':
    pass
