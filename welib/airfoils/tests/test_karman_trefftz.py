import unittest
import numpy as np    
import os
from welib.airfoils.karman_trefftz import * 

scriptDir=os.path.dirname(__file__)

class TestKT(unittest.TestCase):
    def test_one(self):
        import matplotlib.pyplot as plt

        ## Parameters for Karman-Trefftz airfoil
        U0 = 1
        alpha = 5*np.pi/180;
        xc = - 0.2
        yc = 0.1
        tau_deg = 10
        # ---  Display parameters
        n = 100
        # --- 
        l = 2 - tau_deg / 180 # l = 2 - tau/np.pi
        rc, beta, Gamma =  cyl_params(xc, yc, a=1, alpha=alpha, U0=U0)
        nx = 100
        ny = 100
        XLIM = np.array([- 3.5,3.5])
        YLIM = np.array([- 3,3])
        Xg,Yg = np.meshgrid(np.linspace(XLIM[0],XLIM[1],nx), np.linspace(YLIM[0], YLIM[1], ny))

        # --- Main functions
        Xp, Yp = KT_shape(xc, yc, l=l, n=n) # ,Xg,Yg)
        XCp, theta, Cp, U_circ, V_circ = KT_wall(xc, yc, l=l, n=n, U0=U0, alpha=alpha)
        Ug, Vg, CP = KT_flow(Xg, Yg, xc, yc, l=l, U0=U0, alpha=alpha)
        # --- Renaming..
        U, V, X, Y = Ug, Vg, Xg, Yg
        xa, ya = Xp, Yp
        chord     = np.max(Xp)-np.min(Xp)
        thickness = np.max(Yp)-np.min(Yp)

        # --- 
        np.testing.assert_almost_equal(chord, 4.0119868, 5)
        np.testing.assert_almost_equal(thickness, 1.006537, 5)
        np.testing.assert_almost_equal(np.max(Cp), 0.99947, 5)
        np.testing.assert_almost_equal(np.min(Cp),-2.0695, 5)



#         # --- Plot
#         ## Plotting pressure distribution about the airfoil
#         fig,axes = plt.subplots(1, 3, sharey=False, figsize=(12.8,4.8)) # (6.4,4.8)
#         fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
#         ax = axes[0]
#         ax.plot(Xp, Yp    ,'k.-')
#         ax.plot(Xp[0], Yp[0] , 'bo')
#         ax.plot(Xp[1], Yp[1] , 'rd')
#         ax.set_xlabel('')
#         ax.set_ylabel('')
#         ax.legend()
#         ax.set_aspect('equal','box')
# 
#         ax = axes[1]
#         ax.plot((XCp - np.amin(XCp)) / (np.amax(XCp) - np.amin(XCp)), -Cp    , label='')
#         ax.set_xlabel('')
#         ax.set_ylabel('-Cp')
#     #     plt.title(sprintf('Karman-Trefftz C_p \alpha = %.1f deg.',alpha_deg))
#     #     plt.xlim(np.array([0,1]))
#     #     plt.axis('ij')
# 
#         ax=axes[2]
#         im =ax.contourf(Xg, Yg, CP, 15)
#         ax.fill(Xp, Yp  ) #  ,'k.-')
#         ax.set_aspect('equal','box')
#         fig.colorbar(im)
            
#         plt.show()

if __name__ == '__main__':
    unittest.main()
