import unittest
import numpy as np
from welib.system.singledof import *

# --------------------------------------------------------------------------------}
# --- TESTS
# --------------------------------------------------------------------------------{
class Test(unittest.TestCase):
    def test_harmonic_vibration(self):
        # --- underdamped 
        zeta   = 0.1/(2*np.pi)
        psi    = - np.pi / 2
        omega0 = 5

        x, xdot, A, psi =  harmonic_vibration(0, x0=0, xdot0=0, omega0=omega0, zeta=zeta)
        np.testing.assert_almost_equal((x,xdot),(0,0))

        x, xdot, A, psi =  harmonic_vibration(0, x0=0, xdot0=10, omega0=omega0, zeta=zeta)
        np.testing.assert_almost_equal((x,xdot),(0,10))

        x, xdot, A, psi =  harmonic_vibration(0, x0=10, xdot0=0, omega0=omega0, zeta=zeta)
        np.testing.assert_almost_equal((x,xdot),(10,0))
        np.testing.assert_almost_equal((A,psi),(10,np.pi/2),1)

        x, xdot, A, psi =  harmonic_vibration(0, x0=3, xdot0=10, omega0=omega0, zeta=zeta)
        np.testing.assert_almost_equal((x,xdot),(3,10))

        # Keep me
        #import matplotlib.pyplot as plt
        #T=2*np.pi/omega0
        #vt = np.linspace(0, 10*2*np.pi/omega0,300)
        #x, xdot, A, psi =  harmonic_vibration(vt, x0=-1, xdot0=0, omega0=omega0, zeta=zeta)
        #print('A',A,'psi',psi)
        #fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
        #fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
        #ax.plot(vt/T ,x/A,'k')
        #ax.plot(vt/T, np.exp(-zeta*omega0*vt),'k--');
        #ax.plot(vt/T,-np.exp(-zeta*omega0*vt),'k--');
        #ax.set_xlabel(r'$t/T_0$ [-]')
        #ax.set_ylabel(r'$x/A$ [-]')
        #ax.set_title('SingleDOFDampedVibrations')
        #ax.tick_params(direction='in')
        #fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
        #fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
        #ax.plot(vt/T ,xdot,'k')
        #ax.set_xlabel(r'$t/T_0$ [-]')
        #ax.set_ylabel(r'$\dot{x}$ [-]')
        #ax.set_title('SingleDOFDampedVibrationsVelocity')
        #ax.tick_params(direction='in')
        #plt.show()

    def test_forced_vibrations(self):
        # TODO, could test:
        #    frequency of max amplitide
        #    max amplitude
        zeta = 0.5 / (2 * np.pi)
        k = 50
        m = 2
        omega0 = np.sqrt(k / m)
        Omega = 1.0 * omega0
        F0 = 1
        T = 2 * np.pi / Omega
        vt = np.linspace(0,10 * T,1000)
        H0, phi = forced_vibration_particular_cst(Omega/omega0, F0/k, zeta)
        x= forced_vibration(vt, k, m, F0, Omega, zeta, x0=0, xdot0=0)
        # ---
        # vt_peak   = 0:(pi/Omega):vt(end)+phi/Omega                                          ;
        # x1_peak   = A*exp(-zeta*omega0*vt_peak) .* sin(omega0 * sqrt(1-zeta^2) * vt_peak + psi);
        # x2_peak   = H_0 * sin(Omega*vt_peak - phi)                                         ;
        # x_th_peak = x1_peak+x2_peak                                                            ;
#         import matplotlib.pyplot as plt
#         fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
#         fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
#         ax.plot(vt/T, x/H0,'k')
# #         ax.plot(vt/T, x,'k')
#         ax.set_xlim(np.array([0,10]))
#         ax.set_ylim(np.array([-1,1]))
#         ax.set_xlabel(r'$t/T_0$ [-]')
#         ax.set_ylabel(r'$x/H_0$ [-]')
#         ax.set_title('SingleDOFForcedVibrations')
        # --- Showing amplitude
#         vzeta = np.array([0,0.13,0.18,0.3,0.5,1])
#         vfrat = np.linspace(0,3,200) # frequenc ratio
#         F0_over_k = 1
#         fig,axes = plt.subplots(1, 2, sharey=False, figsize=(12.8,4.8)) # (6.4,4.8)
#         fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
#         for iz,zeta in enumerate(vzeta):
#             H0,phi = forced_vibration_particular_cst(vfrat, F0_over_k, zeta) 
#             axes[0].plot(vfrat,H0 , label=r'$\zeta ={:.2f}$'.format(zeta))
#             axes[1].plot(vfrat,phi*180/np.pi, label=r'$\zeta ={:.2f}$'.format(zeta))
#         ax=axes[0]
#         ax.plot(np.array([0,3]),np.array([1,1]),'k:',lw=0.5)
#         ax.plot(np.array([1,1]),np.array([0,4]),'k:',lw=0.5)
#         axes[0].set_xlabel(r'$\Omega/\omega_0$ [-]')
#         axes[0].set_ylabel(r'$H_0/(F_0/k)$ [m]')
#         axes[0].set_ylim(np.array([0,4]))
#         axes[0].set_xlim(np.array([0,3]))
#         axes[1].set_xlabel(r'$\Omega/\omega_0$ [-]')
#         axes[1].set_ylabel(r'$\phi$ [deg]')
#         axes[1].set_ylim(np.array([-180,0]))
#         axes[1].set_xlim(np.array([0,3]))
#         axes[0].legend()
#         axes[1].legend()
#         ax.tick_params(direction='in')
#         ax.set_title('SingleDOFResonanceDamping')
        # ---
#         plt.show()

if __name__=='__main__':
    unittest.main()
