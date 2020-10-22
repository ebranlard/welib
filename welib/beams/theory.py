import numpy as np
import scipy.optimize as sciopt

def UniformBeamBendingModes(Type,EI,rho,A,L,w=None,x=None,Mtop=0,norm='tip_norm',nModes=4):
    """
    returns Mode shapes and frequencies for a uniform beam in bending

    References:
      Inman : Engineering variation
    
    Author: E. Branlard"""
    if x is None or len(x)==0:
        x = np.linspace(0,L,101)
    if np.amax(x) != L:
        raise Exception('Max of x should be equal to L')

    # Dimensionless spanwise position
    x0 = x / L
    s = Type.split('-')
    if s[0].lower()=='unloaded':
        # --- "Theory" (clamped-free, vertical, no weight)
        # See Inman, p.335 or Nielsen1 p129
        if 'unloaded-clamped-free' == (Type.lower()):
            # NOTE: cosh(beta_n)cos(beta_n) =-1
            #  sigma_n = [ np.sinh(beta_n) - sin(beta_n) ]/[cosh(beta_n) + cos(beta_n)]
            #  for j>5, a good approx is B(j) = (2*j-1)np.pi/2  and S(j)=1;
            #B  = [1.87510407, 4.69409113, 7.85475744,10.99554073,14.13716839, (2*6-1)*np.pi/2];
            #S  = [0.734095514 1.018467319 0.999224497 1.000033553 0.999998550 1];
            B = np.zeros(nModes)
            for i in np.arange(nModes):
                B[i] = sciopt.fsolve(lambda x: 1 + np.cosh(x) * np.cos(x), (2*(i+1)-1)*np.pi/2)
        elif 'unloaded-topmass-clamped-free' == (Type.lower()):
            # The geometrical stiffning is not accounted for here
            if Mtop is None:
                raise Exception('Please specify value for Mtop for %s',Type)
            M = rho * A * L
            B = np.zeros(nModes)
            for i in np.arange(nModes):
                B[i] = sciopt.fsolve(lambda x: 1+np.cosh(x)*np.cos(x)-x*Mtop/M*(np.sin(x)*np.cosh(x)-np.cos(x)*np.sinh(x)),(2*(i+1)-1)*np.pi/2)
        else:
            raise Exception('unknown type %s',Type)
        #S  = ( sinh(B)-sin(B) ) ./ ( cosh(B) + cos(B));  # Sigma
        #C  = ( cosh(B)+cos(B) ) ./ ( sinh(B) + sin(B));  # Sigma
        SS = np.sinh(B) + np.sin(B)
        CC = np.cosh(B) + np.cos(B)
        # Frequency
        freq = (B / L) ** 2 / (2 * np.pi) * np.sqrt(EI / (rho * A))
        # --- Mode shapes
        ModesU = np.zeros((len(B),len(x0)))
        ModesV = np.zeros((len(B),len(x0)))
        ModesK = np.zeros((len(B),len(x0)))
        for i in np.arange(nModes):
            ModesU[i,:] =            SS[i] * (np.cosh(B[i]*x0) - np.cos(B[i] * x0)) - CC[i] * (np.sinh(B[i] * x0) - np.sin(B[i] * x0))
            ModesV[i,:] = B[i]    * (SS[i] * (np.sinh(B[i]*x0) + np.sin(B[i] * x0)) - CC[i] * (np.cosh(B[i] * x0) - np.cos(B[i] * x0)))
            ModesK[i,:] = B[i]**2 * (SS[i] * (np.cosh(B[i]*x0) + np.cos(B[i] * x0)) - CC[i] * (np.sinh(B[i] * x0) + np.sin(B[i] * x0)))
            #  ModesU(i,:)  =        cosh(B[i]*x0)-cos(B[i]*x0) - S[i]*(sinh(B[i]*x0)-sin(B[i]*x0)) ;
            #  ModesV(i,:) = B[i]  *(sinh(B[i]*x0)+sin(B[i]*x0) - S[i]*(cosh(B[i]*x0)-cos(B[i]*x0)));
            #  ModesK(i,:) = B[i]^2*(cosh(B[i]*x0)+cos(B[i]*x0) - S[i]*(sinh(B[i]*x0)+sin(B[i]*x0)));
            #  ModesU(i,:)  =        cosh(B[i]*x0)-cos(B[i]*x0) - C[i]*(sinh(B[i]*x0)-sin(B[i]*x0)) ;
            #  ModesV(i,:) = B[i]  *(sinh(B[i]*x0)+sin(B[i]*x0) - C[i]*(cosh(B[i]*x0)-cos(B[i]*x0)));
            #  ModesK(i,:) = B[i]^2*(cosh(B[i]*x0)+cos(B[i]*x0) - C[i]*(sinh(B[i]*x0)+sin(B[i]*x0)));
    elif s[0].lower()=='loaded':
        if 'loaded-clamped-free' == (Type.lower()):
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
    if norm=='tip_norm':
        for i in np.arange(nModes):
            fact = 1 / ModesU[i,-1]
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
