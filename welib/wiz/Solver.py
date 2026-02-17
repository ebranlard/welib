"""
References:
    [1] E. Branlard and M. Gaunaa, Superposition of vortex cylinders for steady and unsteady simulation of rotors of finite tip-speed ratio - Wind Energy, 2014

    [X] E. Branlard, M. Gaunaa - Cylindrical vortex wake model: right cylinder - Wind Energy, 2014
    [X] E. Branlard, M. Gaunaa - Cylindrical vortex wake model: skewed cylinder, application to yawed or tilted rotors - Wind Energy, 2015
    [X] E. Branlard - Wind Turbine Aerodynamics and Vorticity Based Method, Springer, 2017
"""
# --- General
import unittest
import numpy as np
import scipy.optimize as sciopt
try:
    from scipy.integrate import cumulative_trapezoid as cumtrapz
except:
    from scipy.integrate import cumtrapz


def Ct_const_cutoff(CT0,r_bar_cut,vr_bar,r_bar_tip=None):
    """ Returns an almost constant Ct, 
       linearly dropping to zero below a cut-off radius
       linearly dropping to zero after a tip radius
    """
    Ct=np.ones(vr_bar.shape)*CT0
    I=vr_bar<r_bar_cut
    Ct[I]=CT0*vr_bar[I]/r_bar_cut
    if r_bar_tip is not None:
        I=vr_bar>r_bar_tip
        Ct[I]=CT0*(1-(vr_bar[I]-r_bar_tip)/(1-r_bar_tip))
    return Ct 


def InductionsFromPrescribedCtCq_ST(vr_bar,Ct,Cq,Lambda,bSwirl):
    """ 
    Returns the stream tube theory inductions based on a given Ct and Cq.
    Based on script fGetInductions_Prescribed_CT_CQ_ST
    """
    lambda_r=Lambda*vr_bar
    # --- Stream Tube theory
    a_ST       = 1/2*(1-np.sqrt(1-Ct))
    if bSwirl:
        a_prime_ST = Cq/(4*(1-a_ST)*lambda_r)
    #     a_prime_ST = 0.5*(sqrt(1+(4*a_ST.*(1-a_ST)./lambda_r.^2))-1);
    else:
        a_prime_ST =0
    return a_ST,a_prime_ST

def CirculationFromPrescribedCt(vr_bar,vCt,Lambda,bSwirl): 
    """ Finds k to match a given Ct 
    Solves equation (32) of [1]
    Based on script fGetInductions_Prescribed_CT_CQ2.py
    """
    vk = np.zeros(vr_bar.shape)
    vlambda_r = Lambda * vr_bar
    if bSwirl:
        for ir,(r_bar,lambda_r,Ct) in enumerate(zip(vr_bar,vlambda_r,vCt)):
            if Ct>0 and r_bar>0:
                res=sciopt.minimize_scalar(lambda k:np.abs(k*(1+k/(4*lambda_r**2))-Ct), bounds=[0,1.8], method='bounded')
                vk[ir] = res.x
    else:
        vk = vCt
    return vk


def InductionsFromCirculation_VC_Cont(vr_bar,lambda_r,k,bSwirl,method='analytical',bHighThrustCorr=True,Cq=None):
    """ 
    Based on "Continuous formulation". Based on script fInductionCoefficientsCylinders.py
    If Cq is provided, an attempt of using it is made.
         Not recommended, results in lower Tangential inductions
         The circulation obtained from Cq is not smooth

    For now the continuous results seems a bit off on the axial induction compared to the discrete ones
    The tangential induction seem to be the same for the continuous and discrete case
    """
    # Critical values for high-thrust Spera correction
    ac   = 0.34
    Ct_c = 4*ac*(1-ac)
    misc={}
    # Safety
    if lambda_r[0]==0:
        lambda_r[0]=lambda_r[1]*0.01;
    # --- Finding Cylinders intensities
    if method.lower()=='analytical':
        if bSwirl:
            a_prime = k/(4*lambda_r**2) # NOTE: entirely determined from Gamma! Eq.(38)
            Ir=np.arange(len(k)-1,-1,-1)
            Ct_rot = -8*cumtrapz(lambda_r[Ir]**2*a_prime[Ir]**2/vr_bar[Ir], vr_bar[Ir])# Eq.(6)
            Ct_rot = np.concatenate(([0],Ct_rot))
            Ct_rot=Ct_rot[Ir];
        else:
            a_prime = k*0
            Ct_rot = k*0
        Ct_KJ  = k*(1+a_prime)
        Ct_eff = Ct_KJ-Ct_rot 
        a=np.zeros(vr_bar.shape)
        if bHighThrustCorr:
            Icorr=Ct_eff> Ct_c
            Inorm=Ct_eff<=Ct_c
            # Spera correction
            a[Icorr]=(Ct_eff[Icorr]-4*ac**2)/(4*(1-2*ac));
            a[Inorm]=1/2*(1-np.sqrt(1-Ct_eff[Inorm]))
        else:
            a=1/2*(1-np.sqrt(1-Ct_eff));

        # Using Cq
        if Cq is not None:
            a_2     = 1 - lambda_r*Cq/k
            k2      = Cq*lambda_r/(1-a)
            a_prime = k2/(4*lambda_r**2)
            misc['k_Cq'] = k2
    else:
        raise NotImplementedError('')
    misc['Ct_KJ']  = Ct_KJ
    misc['Ct_eff'] = Ct_eff
    misc['Ct_rot'] = Ct_rot
    return a,a_prime,misc

def WakeVorticityFromCirculation_Cont(r_cp,Gamma_cp,R,U0,Omega,bSwirl,method='analytical',bHighThrustCorr=True):
    """
    Uses continuous formulation to find wake vorticity that matches a given circulation distribution
    (need to solve for pitch angle, and hence induction)
    """
    r_cp     = np.asarray(r_cp).ravel()
    Gamma_cp = np.asarray(Gamma_cp).ravel()
    if r_cp[0]==0:
        r_cp[0]=r_cp[1]*0.5;
    # Non dimensional parameters
    k        = Omega*Gamma_cp/(np.pi*U0**2)
    vr_bar   = r_cp/R              
    lambda_r = Omega*r_cp/U0
    # Finding inductions
    a,a_prime,misc= InductionsFromCirculation_VC_Cont(vr_bar,lambda_r,k,bSwirl,method=method,bHighThrustCorr=bHighThrustCorr)
    # Computing convection
    misc['Vz'] = U0*(1-2*a)
    misc['h']  = misc['Vz']*2*np.pi/(Omega*(1+2*a_prime))
    misc['a'] = a
    misc['a_prime'] = a_prime
    misc['Gamma_cp'] = Gamma_cp
    misc['r_cp'] = r_cp
    # Vortex intensities
    Gamma_tilde = Gamma_cp - np.concatenate((Gamma_cp[1:],[0])) #Gamma_tilde = Gamma_i-Gamma_{i+1}
    gamma_t     = - Gamma_tilde/misc['h']
    if bSwirl:
        gamma_l     =   Gamma_tilde/(2*np.pi*r_cp)
        Gamma_r     = - Gamma_cp[0]
    else:
        gamma_l     = 0
        Gamma_r     = 0
    return gamma_t,gamma_l,Gamma_r,misc




def WakeVorticityFromCirculation_Discr(r_cp,Gamma_cp,R,U0,Omega,nB,bSwirl,method='analytical',bTanInd_k=True,bHighThrustCorr=True,r_cyl=None):
    """ 
    Implements discrete algorithm from [1] to find gamma_t, and the inductions
    """
    r_cp     = np.asarray(r_cp).ravel()
    Gamma_cp = np.asarray(Gamma_cp).ravel()
    if r_cyl is None:
        r_cyl    = np.zeros(r_cp.shape)
        if len(r_cyl)==1:
            r_cyl=np.array([R])
        else:
            r_cyl[:-1] = (r_cp[1:] + r_cp[:-1])/2
            r_cyl[-1]  = r_cp[-1] + (r_cp[-1]-r_cp[-2])/2
    else:
        r_cyl    = np.asarray(r_cyl).ravel()
    misc={}
    # --- Checks
    if np.amax(r_cyl) < np.amax(r_cp):
        warnings.warn('Cylinders are assumed to enclose the control points')
        import pdb; pdb.set_trace()
    if r_cyl[0] == 0:
        raise Exception('First cylinder at 0, impossible')

    ##  Interpolating Circulation on the control points
    Gamma_cp   = Gamma_cp * nB # IMPORTANT
    lambda_rcp = Omega*r_cp/U0
    k_r        = Gamma_cp*Omega/(np.pi*U0**2)      # Equation (27) of [Superp]
    Ct_KJ_cp   = k_r * (1+k_r / (4*lambda_rcp**2)) # Equation (32) of [Superp]
    ## Some init

    n = len(r_cyl)
    # Critical values for high-thrust Spera correction
    ac   = 0.34
    Ct_c = 4*ac*(1-ac)
    # ---  Finding Cylinders intensities
    if method.lower()=='analytical':
        # ---  Case 1: Using k / ct to determine gamma
        # Using a_prime to correct Ct
        # Tangential convection velocity
        if not bTanInd_k:
            ap_cyl_conv = 0 * Gamma_cp
        else:
            # Tangential convection velocity as the mean of velocity on both side of the cylinder)
            # See Equation (26) of [1]:  a'_c,i = (\Gamma_i+\Gamma_i+1)/(4 pi Omega R_i^2)
            ap_cyl_conv      = np.zeros(Gamma_cp.shape)
            ap_cyl_conv[:-1] = (Gamma_cp[:-1] + Gamma_cp[1:])/(4*np.pi*Omega*r_cyl[:-1])**2
            ap_cyl_conv[-1]  = Gamma_cp[-1]/(4*np.pi*Omega*r_cyl[-1])**2
        # --- Allocation
        b           = np.zeros(n) # Cummulative sum of gamma_bar
        S           = np.zeros(n) # Term under the square root
        Sbis        = np.zeros(n) # Term under the square root, account
        Ctrot       = np.zeros(n)
        Ctrot_loc   = np.zeros(n)
        Cti         = np.zeros(n)
        Cteff       = np.zeros(n)
        gamma_t_bar = np.zeros(n)
        # --- Init
        C        = Gamma_cp*Omega/(np.pi*U0**2)*(1 + ap_cyl_conv)
        lambda_R = Omega * r_cyl / U0
        # Looping through radii in reversed order
        for i in np.arange(n-1,-1,-1):
            if i==n-1:
                b[i]         = 1
                S[i]         = 1 - C[i]
                Sbis[i]      = 1 - C[i]
                Ctrot[i]     = 0
                Ctrot_loc[i] = 0
            else:
                b[i]         = 1 + np.sum(gamma_t_bar[i+1:])                         # Eq.(25)
                S[i]         = b[i]**2 - (k_r[i] - k_r[i+1])*(1 + ap_cyl_conv[i])    # Eq.(28)
                Ctrot_loc[i] = (k_r[i+1]/2)**2*(1/lambda_R[i]**2-1/lambda_R[i+1]**2) # Eq.(33)
            if not bSwirl:
                Ctrot[i] = 0
                Cti[i]   = k_r[i]
            else:
                Ctrot[i] = np.sum(Ctrot_loc[i+1:])
                Cti[i]   = k_r[i]*(1+k_r[i]/(4*lambda_R[i]**2)) # Eq.(32)
            Cteff[i] = Cti[i] - Ctrot[i] #Eq.(38)
            Sbis[i] = 1 - Cteff[i]
            if bSwirl :
                gamma_t_bar[i] = - b[i] + np.sqrt(Sbis[i]) #Eq.(28)
            else:
                gamma_t_bar[i] = - b[i] + np.sqrt(S[i]) 
#             print('b',b,'S',S)
            if bHighThrustCorr:
                if (Cteff[i] > Ct_c):
                    aeff = (Cteff[i] - 4*ac**2)/(4*(1-2*ac)) # Eq.(40)
                    gamma_t_bar[i] = - b[i] + (1-2*aeff)     # Eq.(41)
            #if np.abs(imag(gamma_t_bar(i))) > 0:
            #    gamma_t_bar[i] = - b[i] + 0
        Vc      = U0 * (b + gamma_t_bar / 2)
        h       = 2*np.pi*Vc/(Omega*(1+ap_cyl_conv))
        #print('Vc',Vc,'b',b,'gamma_t_bar',gamma_t_bar,'ap',ap_cyl_conv,'Omega',Omega,'h',h)


        # Analytical 1
        a_1=-np.cumsum(gamma_t_bar[-1::-1]/2);
        a_1=a_1[-1::-1]
        #if isequal(Algo.VCYL.PitchMethod,'Analytical')
        a=a_1;
        #elif isequal(Algo.VCYL.PitchMethod,'WrongAnalytical')
        #    # Analytical 2
        #    a_2=-(gamma_t_hat/2)/U0;
        #    a=a_2;
        #else
        #    a=a_1;
        #end
        if bSwirl:
            a_prime=Gamma_cp/(4*np.pi*Omega*r_cp**2)
        else:
            a_prime=a*0

        misc['a'] = a
        misc['a_prime'] = a_prime
        misc['Gamma_cp'] = Gamma_cp
        misc['r_cp'] = r_cp
        misc['h']     = h

        # --- Vortex intensities
        Gamma_tilde = Gamma_cp - np.concatenate((Gamma_cp[1:],[0])) #Gamma_tilde = Gamma_i-Gamma_{i+1}
        #gamma_t     = - Gamma_tilde/h
        gamma_t = gamma_t_bar * U0
        if bSwirl:
            gamma_l     =   Gamma_tilde/(2*np.pi*r_cp)
            Gamma_r     = - Gamma_cp[0]
        else:
            gamma_l     = None
            Gamma_r     = None

        return gamma_t,gamma_l,Gamma_r, misc

# --------------------------------------------------------------------------------}
# --- Main interfaces 
# --------------------------------------------------------------------------------{
def WakeVorticityFromCt(r,Ct,R,U0,Omega):
    """ Returns the wake vorticity intensity for a given Ct"""
    if Omega==np.inf:
        Omega=800 # TODO, implement infinite lambd formulae directly
    Lambda=Omega*R/U0
    bSwirl=Lambda<20
    vk=CirculationFromPrescribedCt(r/R,Ct,Lambda,bSwirl)
    Gamma = vk*np.pi*U0**2/Omega
    return WakeVorticityFromCirculation_Discr(r,Gamma,R,U0,Omega,nB=1,bSwirl=bSwirl)

def WakeVorticityFromGamma(r,Gamma,R,U0,Omega):
    """ Returns the wake vorticity intensity for a given circulation """
    if Omega==np.inf:
        Omega=800 # TODO, implement infinite lambd formulae directly
    Lambda=Omega*R/U0
    bSwirl=Lambda<20
    return WakeVorticityFromCirculation_Discr(r,Gamma,R,U0,Omega,nB=1,bSwirl=bSwirl)


# --------------------------------------------------------------------------------}
# --- TESTS
# --------------------------------------------------------------------------------{
class TestSolver(unittest.TestCase):

    def test_Solver_gamma(self):
        # NOTE: TODO, difference in discrete and continuous axial induction
        bSwirl    = True
        U0        = 7
        R         = 50
        r_bar_cut = 0.11
        Lambda    = 6
        CT0       = 0.95
        nCyl      = 50
        nB=1
        Omega = Lambda * U0 / R
        if nCyl==1:
            vr_bar = np.array([0.995])
        else:
            vr_bar = np.linspace(0.005,0.995,nCyl)
        # --- Cutting CT close to the root
        Ct_AD = Ct_const_cutoff(CT0,r_bar_cut,vr_bar)
        ## Finding k to match a given Ct
        vk=CirculationFromPrescribedCt(vr_bar,Ct_AD,Lambda,bSwirl)
        r_cp     = vr_bar * R
        Gamma_cp = vk*np.pi*U0**2/Omega
        gamma_t_d,gamma_l_d,Gamma_r_d,misc_d= WakeVorticityFromCirculation_Discr(r_cp,Gamma_cp,R,U0,Omega,nB,bSwirl)
        gamma_t_c,gamma_l_c,Gamma_r_c,misc_c= WakeVorticityFromCirculation_Cont (r_cp,Gamma_cp,R,U0,Omega,bSwirl)

        lambda_r = Lambda * vr_bar
        a_VC_c, a_prime_VC_c, misc_c = InductionsFromCirculation_VC_Cont(vr_bar,lambda_r,vk,bSwirl)

        a_VC_d, a_prime_VC_d = misc_d['a'], misc_d['a_prime']
        # NOTE:
        # gamma_t more or less match along the blade except at extremities
        #print(Gamma_r_d,Gamma_r_c)
        #import matplotlib.pyplot as plt
        #if bSwirl:
        #    plt.figure()
        #    plt.plot(vr_bar,gamma_l_d,'-d', label='discrete formulation')
        #    plt.plot(vr_bar,gamma_l_c,'-o', label='continuous formulation')
        #    plt.legend()
        #    plt.ylabel('$\gamma_l$ [-]')
        #    plt.xlabel('$r/R$ [-]')
        #    plt.title('SuperpCylindersvsADk')
        #plt.figure()
        #plt.plot(vr_bar,2*lambda_r*a_prime_VC_d,'-d',label='KJ Discr')
        #plt.plot(vr_bar,2*lambda_r*a_prime_VC_c,'-o', label='KJ Cont')
        #plt.legend()
        #plt.ylabel('v_t/U_0 = 2a \lambda_r [-]')
        #plt.xlabel('r/R [-]')
        #plt.title('SuperpCylindersvsADTangentialVelocity')

        #plt.figure()
        #plt.plot(vr_bar,1-a_VC_d,'-d',label = 'KJ Disc')
        #plt.plot(vr_bar,1-a_VC_c,'-o' ,label = 'KJ Cont')
        #plt.legend()
        #plt.ylabel('v_a/U_0 = 1-a [-]')
        #plt.xlabel('r/R [-]')
        #plt.title('SuperpCylindersvsADAxialVelocity')

        #plt.figure()
        #plt.plot(vr_bar,gamma_t_d,'-d', label='discrete formulation')
        #plt.plot(vr_bar,gamma_t_c,'-o', label='continuous formulation')
        #plt.legend()
        #plt.ylabel('$\gamma_t$ [-]')
        #plt.xlabel('$r/R$ [-]')
        #plt.title('SuperpCylindersvsADk')
        #plt.show()

        ###
        np.testing.assert_almost_equal(a_VC_d[10],0.3035,decimal=4)
        np.testing.assert_almost_equal(a_VC_c[10],0.2925,decimal=4)
        np.testing.assert_almost_equal(a_prime_VC_d[10],0.1355,decimal=4)
        np.testing.assert_almost_equal(a_prime_VC_c[10],0.1355,decimal=4)
        # Continuous and discrete are identical for a_prime!
        np.testing.assert_almost_equal(a_prime_VC_c,a_prime_VC_d,decimal=4)

#     def test_Solver_main(self):
#         TODO TODO This was started from Task7_Contours_Lambda
#         import matplotlib.pyplot as plt
#         U0     = 1
#         R      = 1
#         Lambda = 2
#         CT0    = 0.95
#         Omega=Lambda*U0/R
#         nB=1
# 
#         nCyl=30;
# 
#         # Prescribed CT distribution
#         r_bar_cut=0.11
#         r_min=0.0 # ???
# 
#         if Lambda>20:
#             bSwirl=False 
#         else:
#             bSwirl=True 
# 
#         vr_bar_cp  = np.linspace(r_min,1,nCyl) ;
#         # Cylinders higher than CP - TODO, messy
#         if nCyl>9:
#             dx=vr_bar_cp[1]-vr_bar_cp[0]
#             vr_bar_cyl = vr_bar_cp+dx/2
#         else:
#             if nCyl==1:
#                 vr_bar_cyl = vr_bar_cp
#             else:
#                 vr_bar_cyl = vr_bar_cp+0.03;
#         # Prescribing a Ct
#         Ct_Prescr = Ct_const_cutoff(CT0,r_bar_cut,vr_bar_cp)
#         Cq_Prescr = Ct_Prescr*0
# 
#         # Getting inductions and cylinder strength
#         lambda_r = Lambda * vr_bar_cp
#         k = CirculationFromPrescribedCt(vr_bar_cp,Ct_Prescr,Lambda,bSwirl)
#         Gamma_cp = k*np.pi*U0**2/Omega
#         Ct_KJ = k*(1+k/(4*lambda_r**2))
#         r_cp = vr_bar_cp*R
# 
# #         gamma_t,gamma_l,Gamma_root, h = WakeVorticityFromCirculation_Discr(r_cp,Gamma_cp,vr_bar_cyl*R,R,U0,Omega,nB,bSwirl)
#         gamma_t, gamma_l, Gamma_root, h, Vz, a, a_prime = WakeVorticityFromCirculation(r_cp,Gamma_cp,R,U0,Omega,nB,bSwirl)
# 
# 
#         #r_in=vr_bar_cp*R;
#         #r_cp =r_in;
#         #r_cyl=vr_bar_cyl*R;
#         #a_in=a_VC;
#         #a_prime_in=a_prime_VC;
#         #Algo.bSwirl=bSwirl;
#         #Algo.VCYL.PitchMethod='Analytical'; % Analytical or Pitch
#         #Algo.VCYL.bHighThrustCorr=true;
# #         #Algo.VCYL.bTI_k = bSwirl; % tangential induction in k
#         plt.figure()
#         plt.plot(vr_bar_cp,Ct_Prescr,label='Ct prescr')
#         plt.plot(vr_bar_cp,k, label='k')
#         plt.plot(vr_bar_cp,Ct_KJ,'--',label='Ct_KJ')
#         plt.plot(vr_bar_cp,a,'-',label='a')
#         plt.plot(vr_bar_cp,a_prime,'-',label='a_prime')
#         plt.legend()
#         plt.ylabel('$C_t$ [-]')
#         plt.xlabel('$r/R$ [-]')
#         plt.title('SuperCylindersCtDistribution')
# 
#         plt.figure()
#         plt.plot(vr_bar_cp,gamma_t,'-',label='gamma_t')
#         plt.plot(vr_bar_cp,gamma_l,'-',label='gamma_l')
#         plt.legend()
# 
#         plt.show()

if __name__ == "__main__":
    unittest.main()
