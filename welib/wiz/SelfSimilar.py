"""
References:
    [1] N. Troldborg, A.R. Meyer Fortsing, Wind Energy, 2016
"""
import numpy as np
import unittest
    
def ss_u(Xcp, Ycp, Zcp, gamma_t, R=1, lambda_= 0.587, eta=1.32): 
    """
    Self similar approximation of the induced velocity field upstream of a rotor
    """
    # --- Constants
    beta  = np.sqrt(2)
    alpha = 8./9.

    # 
    ui0=gamma_t/2

    r = np.asarray( np.sqrt(Xcp**2 + Ycp**2) )
    z = np.asarray( Zcp                      )

    r12 = np.sqrt(lambda_ * (eta + (z / R) ** 2)) * R        # Eq. (13) from [1]

    epsilon = r / r12
    FEPS = (1/np.cosh(beta * epsilon)) ** alpha                  # Eq. (6) from [1]

    ui0_x = ui0*(1 + z/np.sqrt(z**2 + R**2))          # Eq. (7) from [1]
    uz    = ui0_x * FEPS
    bWake=Zcp>0
    uz[bWake]=0
    return uz

class TestSimilarity(unittest.TestCase):
    def test_Similarity_paper(self):
        # --- Paper check
        R       = 100
        U0      = 10
        CT0     = 0.8
        fact    = 1.1                     
        a       = 0.5*(1-np.sqrt(1-fact*CT0))# <<< NOTE gamma factor here!!!
        gamma_t = -2*U0*a                 
        Z0 = [-R, -2*R, -3*R]
        uz_list=[]
        for z0 in Z0:
            rcp = np.linspace(0,5,50)*R
            zcp = rcp*0+ z0
            uz_list.append(ss_u (rcp, rcp*0, zcp, gamma_t, R))

        #--- Verifying that mid induction at r=0 is consistent with axis scaling
        ui0=gamma_t/2
        for z0,uz in zip(Z0,uz_list):
            np.testing.assert_almost_equal(uz[0], ui0*(1 + z0/np.sqrt(z0**2 + R**2)), 2)

        # ---- Plot
        #import matplotlib.pyplot as plt
        #fig,ax = plt.subplots(1,1)
        #for z0,uz in zip(Z0,uz_list):
        #    ax.plot(rcp/R, (U0+uz)/U0    , label='x={:.0f}R'.format(z0/R))
        #ax.set_xlim([0  , 5])
        #ax.set_ylim([0.9, 1])
        #ax.set_xlabel('r/R [-]')
        #ax.set_ylabel('U/U_0 [-]')
        #ax.legend()
        #plt.show()

if __name__ == "__main__":
#     TestSimilarity().test_Similarity_example()
    unittest.main()
