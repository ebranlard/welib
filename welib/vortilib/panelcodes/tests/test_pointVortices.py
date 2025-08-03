import unittest
import numpy as np    
import os as os
import matplotlib.pyplot as plt
from welib.CFD.flows2D import *
from welib.vortilib.panelcodes.pointVortices import *

scriptDir = os.path.dirname(__file__)

class Test(unittest.TestCase):

    def test_cylinder(self, plot=False):
        # Lifting cylinder
        # --- Parameters
        alpha=0 # NOTE: alpha is 0 but we have lift
        ge=None
        Cp_theory = None
        Uth_theory=None
        theta_mid=None
        curv_method='zero'
        curv_method='Menger'
        #curv_method='Lewis'
        backDiagCorr=False
        backDiagCorr=True
        coords_filename =None

        # Cylinder
        m=151
        U0    = 1
        R=2
        Gamma =-4*np.pi*U0*R*0.5
        theta   =-np.linspace(0,2*np.pi, m+1)
        theta_TE = np.arcsin(Gamma/(4*np.pi*U0*R))
        theta   += theta_TE
        iTE=-1        
        iTE=0
        XP = R*np.cos(theta)
        YP = R*np.sin(theta)
        dtheta=theta[1]-theta[0]
        theta_mid   = theta[:-1]+dtheta/2
        Uth_theory = -2.0*U0*np.sin(theta_mid) + Gamma/(2*np.pi*R)# Utheta at CP
        ge         = -Uth_theory
        Cp_theory  = 1-(Uth_theory)**2/U0**2
        XCp_th  = (XP[0:-1]+XP[1:])/2
        hasLift   = abs(Gamma)>0

        # --- Derived parameters
        Vinf_x = U0*np.cos(alpha)
        Vinf_y = U0*np.sin(alpha)
        fU =lambda X, Y : (X*0+Vinf_x, X*0+Vinf_y) # External velocity function

        # --- Panel method
        gammas, out = VPts_panel_solve(XP, YP, fU=fU, hasLift=hasLift, iTE=iTE, curv_method=curv_method, backDiagCorr=backDiagCorr)
        CP = out['CP']

        ds_mean = np.mean(out['ds'])
        VP = out['VP']

        if plot:
            plot_airfoil(XP, YP, Uwall=out['Vwall'], ntScale=0.1, UScale=0.8)
            #plot_airfoil(XP, YP, nt=True, ntScale=0.1, UScale=0.8)

            fig,axes = plt.subplots(1, 3, sharey=False, figsize=(12.8,4.8)) # (6.4,4.8)
            fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
            ax=axes[0]
            if ge is not None:
                ax.plot(VP[:,0], ge    , 'k-', label='Theory')
            ax.plot(out['VP'][:,0], gammas, '--', label='gammas')
            ax=axes[1]
            if theta_mid is not None:
                ax.plot(theta_mid, Uth_theory, 'k-', label='Utheta')
                ax.plot(theta_mid, out['Un'], '--', label='Un')
                ax.plot(theta_mid, out['Ur'], ':', label='Ur')
                ax.plot(theta_mid, out['Ut'],  '.', label='Ut')
                ax.plot(theta_mid, out['Uth'], '--', label='Utheta')

                ax.legend()

            ax=axes[2]
            if Cp_theory is not None:
                ax.plot(XCp_th, Cp_theory, 'k-', label='Theory', lw=2)
            ax.plot(out['VP'][:,0], out['Cp'], '--', label='Cp')
            ax.legend()
            plt.show()

        # --- Test
        #print(XCp_th, len(XCp_th))
        #print(CP[:,0], len(CP[:,0]))
        np.testing.assert_allclose(CP[:,0], XCp_th)
        np.testing.assert_array_almost_equal(out['Cp'][:-1], Cp_theory[1:], 2) # TODO weird offset

    def test_ellipse(self, plot=False):
        # --- Parameters
        alpha=0
        ge=None
        Cp_theory = None
        Uth_theory=None
        theta_mid=None
        curv_method='zero'
        curv_method='Menger'
    #     curv_method='Lewis'
        backDiagCorr=False
        backDiagCorr=True
        coords_filename =None

        # Ellipse
        m=31
        U0    = 1
        major=1
        ratio=0.05
        alpha = 30*np.pi/180
        theta   =-np.linspace(0,2*np.pi, m+1)
        abyr= np.sqrt((1.-ratio)/(1.+ratio))# ! see Lewis p 50
        XP = 0.5*major       *(np.cos(theta))
        YP = 0.5*major*ratio*  np.sin(theta)
        dtheta=theta[1]-theta[0]
        theta_mid   = theta[:-1]+dtheta/2
        f = np.sqrt(1.0+ (abyr**4) - 2.*(abyr**2)*np.cos(2.*theta_mid) )
        Uth_theory =-(-2.0*U0*np.sin(theta_mid-alpha)-2*np.sin(alpha))/f          # See Lewis p50
        ge         =-(-2.0*U0*np.sin(theta_mid-alpha)-2*np.sin(alpha))/f          # See Lewis p50
        iTE=0
        Cp_theory  = 1-(Uth_theory)**2/U0**2
        XCp_th  = (XP[0:-1]+XP[1:])/2
        hasLift =abs(alpha)>0

        # --- Derived parameters
        Vinf_x = U0*np.cos(alpha)
        Vinf_y = U0*np.sin(alpha)
        fU =lambda X, Y : (X*0+Vinf_x, X*0+Vinf_y) # External velocity function

        # --- Panel method
        gammas, out = VPts_panel_solve(XP, YP, fU=fU, hasLift=hasLift, iTE=iTE, curv_method=curv_method, backDiagCorr=backDiagCorr)
        CP = out['CP']

        ds_mean = np.mean(out['ds'])
        VP = out['VP']

        if plot:
            plot_airfoil(XP, YP, Uwall=out['Vwall'], ntScale=0.1, UScale=0.8)
            #plot_airfoil(XP, YP, nt=True, ntScale=0.1, UScale=0.8)

            fig,axes = plt.subplots(1, 3, sharey=False, figsize=(12.8,4.8)) # (6.4,4.8)
            fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
            ax=axes[0]
            if ge is not None:
                ax.plot(VP[:,0], ge    , 'k-', label='Theory')
            ax.plot(out['VP'][:,0], gammas, '--', label='gammas')
            ax=axes[1]
            if theta_mid is not None:
                ax.plot(theta_mid, Uth_theory, 'k-', label='Utheta')
                ax.plot(theta_mid, out['Un'], '--', label='Un')
                ax.plot(theta_mid, out['Ur'], ':', label='Ur')
                ax.plot(theta_mid, out['Ut'],  '.', label='Ut')
                ax.plot(theta_mid, out['Uth'], '--', label='Utheta')

                ax.legend()

            ax=axes[2]
            if Cp_theory is not None:
                ax.plot(XCp_th, Cp_theory, 'k-', label='Theory', lw=2)
            ax.plot(out['VP'][:,0], out['Cp'], ':', label='Cp')
            ax.legend()
            plt.show()

        # --- Test
        #np.testing.assert_allclose(CP[:,0], XCp_th)
        #np.testing.assert_array_almost_equal(out['Cp'][:-1], Cp_theory[1:], 1)
        #np.testing.assert_array_almost_equal(out['Cp'], Cp_theory, 1)


if __name__ == '__main__':
    unittest.main()
    #Test().test_cylinder(plot=True)
    #Test().test_ellipse(plot=True)
    #plt.show()
