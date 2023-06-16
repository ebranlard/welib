import unittest
import numpy as np
import os
MyDir=os.path.dirname(__file__)

from scipy.integrate import solve_ivp

from welib.airfoils.Polar import Polar
from welib.airfoils.DynamicStall import * 

# --------------------------------------------------------------------------------}
# ---  
# --------------------------------------------------------------------------------{
class TestDynamicStall(unittest.TestCase):
    def assertNaN(self,x):
        self.assertTrue(np.isnan(x))

    def test_oye(self):
        #FFA-W3-241 airfoil Dyna Stall
        P=Polar(os.path.join(MyDir,'../data/FFA-W3-241-Re12M.dat'), compute_params=True)

        omega       = 12.57
        T           = 2*np.pi/omega
        tau         = 0.08
        alpham      = 20
        dt          = 0.01                   # time step
        # 
        fs_prev = P.fs_interp(alpham) # init with steady value
        Cl0 = P.cl_interp(alpham) # init with steady value
        Cl_new,fs_prev_new = P.dynaStallOye_DiscreteStep(alpham,tau,fs_prev,dt)

        # Testing that value at t=0 is equal to the steady state cl
        np.testing.assert_almost_equal(Cl_new, Cl0, decimal=4)
        self.assertEqual(fs_prev_new,fs_prev)

        # An increase of alpha from the steady value should have dCl/dt>0
        Cl_new,fs_prev_new = P.dynaStallOye_DiscreteStep(alpham+1,tau,fs_prev,dt)
        self.assertEqual( (Cl_new-Cl0)>0 ,True)
        self.assertEqual( (fs_prev_new-fs_prev)<0 ,True)

        # A decrease of alpha from the steady value should have dCl/dt<0
        Cl_new,fs_prev_new = P.dynaStallOye_DiscreteStep(alpham-1,tau,fs_prev,dt)
        self.assertEqual( (Cl_new-Cl0)<0 ,True)
        self.assertEqual( (fs_prev_new-fs_prev)>0 ,True)


    def test_convergence(self):
        # Starting from a wrong set point, the Cl value should converge to the steady Cl value
        # Script params, reading polar
        radians=True
        P=Polar(os.path.join(MyDir,'../data/FFA-W3-241-Re12M.dat'), compute_params=True, radians=radians)
        U0, chord = 10, 0.1591
        alpha_st  = 3 * P._alpha0 
        tau_t     = np.linspace(0,40,30)
        vt        = chord * tau_t / (2*U0)
           
        # Oye's Parameters
        p_oye = dynstall_oye_param_from_polar(P, tau_chord=chord/U0)
        p_mhh = dynstall_mhh_param_from_polar(P, chord, constants='OpenFAST')
        # Inputs
        u=dict()
        u['U']         = lambda t: U0
        u['U_dot']     = lambda t: 0 
        u['alpha']     = lambda t: alpha_st
        u['omega']     = lambda t: 0
        u['alpha_34']  = u['alpha']

        # Init values, off
        y0_oye = [0]
        y0_mhh = [0,0,0,0]

        Cl_mhh  = np.zeros(len(vt))
        Cl_oye  = np.zeros(len(vt))
        ## Integration using solve_ivp
        np.seterr(under='ignore')
        sol_mhh = solve_ivp(lambda t,x: dynstall_mhh_dxdt(t,x,u,p_mhh), t_span=[0, max(vt)], y0=y0_mhh, t_eval=vt)
        for it,t in enumerate(vt):
            Cl_mhh[it],_,_ = dynstall_mhh_outputs(t,sol_mhh.y[:,it],u,p_mhh)

        ## Integration using solve_ivp
        sol_oye = solve_ivp(lambda t,x: dynstall_oye_dxdt(t,x,u,p_oye), t_span=[0, max(vt)], y0=y0_oye, t_eval=vt)
        for it,t in enumerate(vt):
            Cl_oye[it] = dynstall_oye_output(vt[it],sol_oye.y[0,it],u,p_oye)

        ## Steady values
        Cl_st  = P.cl_interp(alpha_st)
        fs_st  = P.fs_interp(alpha_st) 

        ## --- Test that the last value is the steady state one
        np.testing.assert_almost_equal(Cl_mhh[-1], Cl_st, decimal=3)
        np.testing.assert_almost_equal(Cl_oye[-1], Cl_st, decimal=3)
        np.testing.assert_almost_equal(sol_oye.y[0,-1], fs_st, decimal=3)

        # --- Plot, keep me
        #import matplotlib.pyplot as plt
        #fig=plt.figure()
        #ax = fig.add_subplot(111)
        #ax.plot(tau_t,Cl_mhh[:]/Cl_st ,'--',label = 'Cl dynamic (MHH)')
        #ax.plot(tau_t,Cl_oye[:]/Cl_st ,'-' ,label = 'Cl dynamic (Oye)')
        #ax.set_xlabel('Dimensionless time [-]')
        #ax.set_ylabel('Cl [-]')
        #plt.legend()
        #plt.show()

        # 
        #y0_mhh = dynstall_mhh_steady(0,u,p_mhh)


    def test_mhh_wagner_step(self):
        # Step from alpha0 to alpha0+2, testing the circulatory response (history), 
        # The Cl result is compared to Wagner's function
        radians=True # <<<
        P=Polar(os.path.join(MyDir,'../data/FFA-W3-241-Re12M.dat'), compute_params=True, radians=radians)

        U0, chord = 10, 0.1591
        alpha1  = P._alpha0 
        alpha2  = alpha1+2*np.pi/180
        tau_t   = np.linspace(0,30,100)
        vt      = chord * tau_t / (2*U0)

        ## MHH Parameters and Inputs
        np.seterr(under='ignore')
        p = dynstall_mhh_param_from_polar(P, chord, constants='Jones')
        u=dict()
        u['U']         = lambda t: U0
        u['U_dot']     = lambda t: 0 
        u['alpha']     = lambda t: alpha1 if t<=0 else alpha2 
        u['omega']     = lambda t: 0
        u['alpha_34']  = u['alpha']
        ## Steady values
        Cl_st2 = P.cl_interp(alpha2)
        y0_mhh = dynstall_mhh_steady(0,u,p)

        Cl_mhh  = np.zeros(len(vt))
        # Integration using solve_ivp
        sol_mhh = solve_ivp(lambda t,x: dynstall_mhh_dxdt(t,x,u,p), t_span=[0, max(vt)], y0=y0_mhh, t_eval=vt)
        for it,t in enumerate(vt):
            Cl_mhh[it],_,_ = dynstall_mhh_outputs(t,sol_mhh.y[:,it],u,p)

        Cl_wag_Jones=wagner(tau_t, constants='Jones')

        np.testing.assert_almost_equal(Cl_mhh[1:]/Cl_st2,Cl_wag_Jones[1:],decimal=4)

        # --- Plot, keep me
        #import matplotlib.pyplot as plt
        #fig=plt.figure()
        #ax = fig.add_subplot(111)
        #ax.plot(tau_t    ,Cl_wag_Jones,'k'  ,label='Wagner function (Jones approx.)')
        #ax.plot(tau_t[1:],Cl_mhh[1:]/Cl_st2 ,'--',label = 'Cl dynamic (MHH)')
        #ax.set_xlabel('Dimensionless time 2 U_0 t/c [-]')
        #ax.set_ylabel('Cl/Cl_ref [-]')
        #plt.ylim([0.3,1.1])
        #plt.title('Response to an angle of attack change')
        #plt.legend()
        #plt.show()



if __name__ == '__main__':
    unittest.main()
