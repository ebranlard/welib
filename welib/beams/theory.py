import numpy as np
import unittest
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


# --------------------------------------------------------------------------------}
# --- TEST  
# --------------------------------------------------------------------------------{
class Test(unittest.TestCase):
    def test_beam_theory_bending(self):
        # --- Uniform beam, testing shape functions and frequencies
        freq_ref = [0.81404393,5.10152622 ,14.28442115,27.99176434,46.27239213,69.12294146]
        Umid_ref = [0.33952311,-0.71366583,0.01968759 ,0.70711864 ,0.00085144 ,-0.70710680]
        Vmid_ref = [0.01163054,0.00453142 ,-0.05551999,0.00045037 ,0.09996480 ,0.00003058]
        Kmid_ref = [0.00011938,0.00157253 ,0.00012147 ,-0.00854920,0.00001702 ,0.02111106]
        L  = 100                  ;
        EI = 1.868211939147334e+12;
        m  = 8.828201296825122e+03;
        freq,x,U,V,K = UniformBeamBendingModes('unloaded-clamped-free',EI,m,A=1,L=L,nModes=6)
        #import matplotlib.pyplot as plt
        #plt.figure
        #plt.plot(x,U[0,:])
        #plt.plot(x,U[1,:])
        #plt.show()
        np.testing.assert_almost_equal(freq,freq_ref)
        np.testing.assert_almost_equal(U[:,50],Umid_ref)
        np.testing.assert_almost_equal(V[:,50],Vmid_ref)
        np.testing.assert_almost_equal(K[:,50],Kmid_ref)

    def test_beam_theory_longi(self):
        # --- Longitudinal beam
        freq_ref=[ 12.93048538, 38.79145615, 64.65242691, 90.51339768]
        Umid_ref=[ 0.70710678, -0.70710678,  -0.70710678,  0.70710678]
        L   = 100
        D   = 8
        t   = 0.045
        A   = np.pi*((D/2)**2-(D/2-t)**2)
        E   = 210e9                       # Young modulus [Pa] [N/m^2]
        rho = 7850
        freq,x,U =  UniformBeamLongiModes('unloaded-clamped-free',E,rho,A,L,nModes=4,norm='tip_norm')

        np.testing.assert_almost_equal(freq,freq_ref)
        np.testing.assert_almost_equal(U[:,50],Umid_ref)
        #import matplotlib.pyplot as plt
        #plt.figure()
        #for i in np.arange(4):
        #    plt.plot(x,U[i,:])
        #print(freq)
        #plt.show()
    
    def test_beam_theory_torsion(self):
        # --- Torsion of a uniform beam
        freq_ref=([5.61858268, 16.85574804, 28.09291340, 39.33007876])
        Vmid_ref=[ 0.70710678,-0.70710678, -0.70710678,   0.70710678]

        L   = 100
        G   = 79.3e9                   #% Shear modulus. Steel: 79.3  [Pa] [N/m^2]
        D   = 8
        t   = 0.045
        A   = np.pi * ((D / 2) ** 2 - (D/ 2 - t) ** 2)
        rho = 7850
        Ip  = np.pi / 32 * (D ** 4 - (D - 2 * t) ** 4)
        Kt  = np.pi / 64 * (D ** 4 - (D - 2 * t) ** 4)
        freq,x,V,_ = UniformBeamTorsionModes('unloaded-clamped-free',G,Kt,Ip,rho,A,L)
        #import matplotlib.pyplot as plt
        #plt.figure()
        #for i in np.arange(4):
        #    plt.plot(x,V[i,:])
        #print(freq)
        #plt.show()
        np.testing.assert_almost_equal(freq,freq_ref)
        np.testing.assert_almost_equal(V[:,50],Vmid_ref)

if __name__=='__main__':
    unittest.main()
