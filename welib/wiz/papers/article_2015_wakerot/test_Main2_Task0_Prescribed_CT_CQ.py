"""

Reproduces the figures from 
    [1] Branlard, Troldborg, Gaunaa (2015) A vortex based BEM-like algorithm accounting for wake rotation

"""

# --- General
import matplotlib.pyplot as plt
import unittest
import numpy as np
# --- Local
from welib.wiz.Solver import *


def main(test=False):
    try:
        import welib.weio as weio
    except:
        print('Skipping test')
        return
    bSwirl    = True
    U0        = 7
    R         = 50
    r_bar_cut = 0.11
    Lambda   = 6
    CT0       = 0.95
    nCyl      = 50
    nB=1
    ## Setting up the Loading
    Omega = Lambda * U0 / R
    vr_bar = np.linspace(0.005,0.995,nCyl)
    # -- Setting up Cq
    df=weio.read('MadsenAD.csv').toDataFrame()
    Cq_AD = df['Cq']
    va_AD = df['Va']
    vt_AD = df['Vt']
    # --- Cutting CT close to the root
    Ct_AD = Ct_const_cutoff(CT0,r_bar_cut,vr_bar)
    ## Finding k to match a given Ct
    vk=CirculationFromPrescribedCt(vr_bar,Ct_AD,Lambda,bSwirl)
    lambda_r=Lambda*vr_bar

    ## Determination of inductions
    a_ST ,a_prime_ST          = InductionsFromPrescribedCtCq_ST(vr_bar,Ct_AD,Cq_AD,Lambda,bSwirl)
    a_VC , a_prime_VC , misc  = InductionsFromCirculation_VC_Cont(vr_bar,lambda_r,vk,bSwirl,Cq=Cq_AD)
    a_VCc, a_prime_VCc, miscc = InductionsFromCirculation_VC_Cont(vr_bar,lambda_r,vk,bSwirl)
    vk2    = misc['k_Cq']
    Ct_eff = miscc['Ct_eff']
    Ct_rot = miscc['Ct_rot']
    ##  Using vortex cylinders function (based on k/Gamma)
    r_cp     = vr_bar * R
    Gamma_cp = vk*np.pi*U0**2/Omega
    _,_,_,misc= WakeVorticityFromCirculation_Discr(r_cp,Gamma_cp,R,U0,Omega,nB,bSwirl)
    a_VCd, a_prime_VCd=misc['a'], misc['a_prime']

    ##
    plt.figure()
    plt.plot(vr_bar,vk, label='From C_t')
    plt.plot(vr_bar,vk2, label='From C_q')
    plt.legend()
    plt.ylabel('k [-]')
    plt.xlabel('r/R [-]')
    plt.title('SuperpCylindersvsADk')
    ##
    plt.figure(110)
    plt.plot(vr_bar,Ct_AD,'k-',label='C_t AD')
    plt.plot(vr_bar,Ct_rot,'--',label='C_t_rot')
    plt.plot(vr_bar,Ct_eff,'-.',label='C_t eff')
    plt.plot(vr_bar,Cq_AD,'-', label='C_q AD')
    plt.legend()
    plt.ylabel('Load coefficients [-]')
    plt.xlabel('r/R [-]')
    plt.title('Ct')
    plt.ylim(np.array([0,1]))
    ##
    plt.figure(111)
    plt.plot(vr_bar,va_AD  ,'k-',label = 'AD'  )
    #plt.plot(vr_bar,1-a_VC ,'-' ,label = 'KJ Cont, with Cq')
    plt.plot(vr_bar,1-a_VCd,'-d',label = 'KJ 2')
    plt.plot(vr_bar,1-a_VCc,'o',label = 'KJ Cont')
    plt.plot(vr_bar,1-a_ST ,'--',label = 'ST'  )
    plt.legend()
    plt.ylabel('v_a/U_0 = 1-a [-]')
    plt.xlabel('r/R [-]')
    plt.title('SuperpCylindersvsADAxialVelocity')

    plt.figure(112)
    plt.plot(vr_bar,vt_AD,'k-',label='AD')
    plt.plot(vr_bar,2*lambda_r*a_prime_VC, label='KJ Cont')
    plt.plot(vr_bar,np.multiply(2 * lambda_r,a_prime_VCd),'-d',label='KJ Discr')
    plt.plot(vr_bar,2*lambda_r*a_prime_VCc,'o', label='KJ Cont')
    plt.plot(vr_bar,np.multiply(2 * lambda_r,a_prime_ST),'--', label='ST')
    plt.legend()
    plt.ylabel(r'v_t/U_0 = 2a \lambda_r [-]')
    plt.xlabel(r'r/R [-]')
    plt.title('SuperpCylindersvsADTangentialVelocity')

    if not test:
        plt.show()
    else:
        plt.close('all')

class Test(unittest.TestCase):
    def test_Article_WakeRot_Prescribed_CT_CQ(self):
#         import sys
#         if sys.version_info >= (3, 0):
        main(test=False)

if __name__ == "__main__":
    main()
