import unittest
import numpy as np    
from welib.BEM.unsteadyBEM import *
from numpy import cos, sin, arctan2, pi, arccos, exp, abs, min, sqrt
import os
from welib.BEM.unsteadyBEM import _fInductionCoefficients

MyDir=os.path.dirname(__file__)

class Test(unittest.TestCase):
    def test_relaxation_issue(self):
        # NOTE: 
        # Without relaxation, high thrust correction will lead to large hysteris 
        # in axial induction
        bSwirl=True
        CTcorrection='AeroDyn'
        swirlMethod='HAWC2'
        relaxation=0.5
        a_last=np.array([0])

        R_a2g = np.array([[ 0.99999829,  0.00185005,  0.        ],
                         [-0.00185005,  0.99999829,  0.        ],
                         [ 0.        ,  0.        ,  1.        ]])
        R_p2g = np.array([[ 1.,  0.,  0.],
                          [ 0., -1.,  0.],
                          [ 0.,  0.,  1.]])

        lambda_r = np.array([6.454224083216526])
        sigma    = np.array([0.010992800851847919])
        r        = np.array([61.6333])
        R=63
        nB=3
        V0=np.array([10])

        # Velocity in global
        Vwnd_g = np.array([10,0,0])
        Vstr_g = np.array([0, -64.54224083,0])
        Vind_g = np.array([0,0,0]) # init

        a_store=np.zeros(100)
        for it in np.arange(len(a_store)):
            Vrel_g = Vwnd_g+Vind_g-Vstr_g

            # Airfoil coordinates
            Vrel_a = (R_a2g.T).dot(Vrel_g)
            # Polar coordinates
            Vstr_p = (R_p2g.T).dot(Vstr_g) # Structural velocity in polar coordinates
            Vrel_p = (R_p2g.T).dot(Vrel_g)
            Vwnd_p = (R_p2g.T).dot(Vwnd_g) # Wind Velocity in polar coordinates

            Vrel_norm = np.array([sqrt(Vrel_a[0]**2 + Vrel_a[1]**2)])
            phi_p     = np.array([arctan2(Vrel_p[0],-Vrel_p[1])])  # NOTE: using polar grid for phi
            F = 2./pi*arccos(exp(-(nB *(R-r))/(2*r*sin(phi_p))))

            alpha = np.arctan2(Vrel_a[0],Vrel_a[1]) # angle of attack [rad]
            Cl=2*np.pi*(alpha+4.432*np.pi/180)
            Cd=0
            C_xa       ,C_ya        = Cl*cos(alpha)+ Cd*sin(alpha  )   ,  -Cl*sin(alpha)+ Cd*cos(alpha)
            C_g = R_a2g.dot(np.array([C_xa, C_ya, 0]))
            C_p = (R_p2g.T).dot(C_g)
            cnForAI = np.array([C_p[0]])
            ctForTI = np.array([C_p[1]])
            a,aprime,CT = _fInductionCoefficients(Vrel_norm, V0, F, cnForAI, ctForTI, lambda_r, sigma, phi_p, 
                    bSwirl=bSwirl, CTcorrection=CTcorrection, swirlMethod=swirlMethod,
                    relaxation=relaxation, a_last=a_last)
            a_last = a
            a_store[it]=a

            Vind_qs_p = np.array([-a[0]*V0[0], -aprime[0]*Vstr_p[1], 0])
            Vind_g    = R_p2g.dot(Vind_qs_p) # global


        np.testing.assert_almost_equal(a_store[-1], a_store[-2], 5)
        np.testing.assert_almost_equal(a_store[-2], a_store[-3], 5)
        np.testing.assert_almost_equal(a_store[-3], a_store[-4], 5)

    def test_equilibirumDynWake(self):
        # Test that wake reaches the same equilibrium with or without the dynamic wake
        dt       = 0.1
        # Read a FAST model to get structural parameters for blade motion
        motion = PrescribedRotorMotion()
        motion.init_from_FAST(os.path.join(MyDir,'../../../data/NREL5MW/Main_Onshore_OF2.fst'), tilt=0)
        motion.setType('constantRPM', RPM=10.0)
        # Read a FAST model to get Aerodynamic parameters to initialze unstady BEM code
        BEM = AeroBEM()
        BEM.init_from_FAST(os.path.join(MyDir,'../../../data/NREL5MW/Main_Onshore_OF2.fst'))


        # --- Simulation 1, no dynamic inflow, starting at equilibrium 
        tmax=5*dt
        BEM.bDynaWake = False # dynamic inflow model
        BEM.timeStepInit(0,tmax,dt) 
        xdBEM = BEM.getInitStates()
        for it,t in enumerate(BEM.time):
            motion.update(t)
            xdBEM = BEM.timeStep(t, dt, xdBEM, motion.psi,
                    motion.origin_pos_gl, motion.omega_gl, motion.R_b2g, 
                    motion.R_ntr2g, motion.R_bld2b,
                    motion.pos_gl, motion.vel_gl, motion.R_s2g, motion.R_a2g,
                    firstCallEquilibrium=it==0
                    )
        a     = BEM.AxInd.copy()
        aprime= BEM.TnInd.copy()

        # --- Simulation 2, dynamic inflow, starting at equilibrium 
        tmax=5*dt
        BEM.bDynaWake = True # dynamic inflow model
        BEM.timeStepInit(0,tmax,dt) 
        xdBEM = BEM.getInitStates()
        for it,t in enumerate(BEM.time):
            motion.update(t)
            xdBEM = BEM.timeStep(t, dt, xdBEM, motion.psi,
                    motion.origin_pos_gl, motion.omega_gl, motion.R_b2g, 
                    motion.R_ntr2g, motion.R_bld2b,
                    motion.pos_gl, motion.vel_gl, motion.R_s2g, motion.R_a2g,
                    firstCallEquilibrium=it==0
                    )
        a2     = BEM.AxInd.copy()
        aprime2= BEM.TnInd.copy()



        # --- Simulation 3, no dynamic inflow
        tmax     = 6
        BEM.bDynaWake = False # dynamic inflow model
        BEM.timeStepInit(0,tmax,dt) 
        xdBEM = BEM.getInitStates()
        for it,t in enumerate(BEM.time):
            motion.update(t)
            xdBEM = BEM.timeStep(t, dt, xdBEM, motion.psi,
                    motion.origin_pos_gl, motion.omega_gl, motion.R_b2g, 
                    motion.R_ntr2g, motion.R_bld2b,
                    motion.pos_gl, motion.vel_gl, motion.R_s2g, motion.R_a2g
                    )
        a3     = BEM.AxInd.copy()
        aprime3= BEM.TnInd.copy()


        # --- Simulation 4 Dynamic Wake
        dt  = 0.5
        tmax=40
        BEM.bDynaWake = True # dynamic inflow model
        BEM.timeStepInit(0,tmax,dt) 
        xdBEM = BEM.getInitStates()
        for it,t in enumerate(BEM.time):
            motion.update(t)
            xdBEM = BEM.timeStep(t, dt, xdBEM, motion.psi,
                    motion.origin_pos_gl, motion.omega_gl, motion.R_b2g, 
                    motion.R_ntr2g, motion.R_bld2b,
                    motion.pos_gl, motion.vel_gl, motion.R_s2g, motion.R_a2g
                    )
        a4     = BEM.AxInd.copy()
        aprime4= BEM.TnInd.copy()

        # --- Tests
        np.testing.assert_almost_equal(a2[-1], a[-1], 5)
        np.testing.assert_almost_equal(aprime2[-1], aprime[-1], 5)
        np.testing.assert_almost_equal(a3[-1], a[-1], 4)
        np.testing.assert_almost_equal(aprime3[-1], aprime[-1], 4)
        np.testing.assert_almost_equal(a4[-1], a[-1], 3)
        np.testing.assert_almost_equal(aprime4[-1], aprime[-1], 3)


        #import matplotlib.pyplot as plt
        #fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
        #fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
        #ax.plot(a [:,0,10],'k-', label='Eq. no dyn')
        #ax.plot(a2[:,0,10],'--', label='Eq. w/ dyn')
        #ax.plot(a3[:,0,10], label='no dyn')
        #ax.plot(a4[:,0,10], label='dyn')
        #ax.set_xlabel('')
        #ax.set_ylabel('')
        #ax.legend()
        #ax.tick_params(direction='in')
        #plt.show()



if __name__ == '__main__':
    unittest.main()
