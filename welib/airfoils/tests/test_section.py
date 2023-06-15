import matplotlib.pyplot as plt
import unittest
import numpy as np    
import os as os
from welib.airfoils.section import Section

scriptDir = os.path.dirname(__file__)

class TestSection(unittest.TestCase):
    def test_section_x(self):
        # Structure
        m = 1
        c = 0
        k = 61.7
        x0   = -0.02 # initial position
        xd0  = 0     # initial velocity

        # Inputs
        Uy=2    # Incoming wind speed "Omega r" [m/s]
        Ux=0    # Incoming wind speed "U0" [m/s]
        theta_p = -15


        # TODO Interface likely to change
        sec = Section()
        sec._chord = 0.2
        sec._rho = 1.225
        sec._tau = 0.08
        sec._y_AQ = 0
        sec._y_AT = 0
        M = sec.setMassMatrix(m=m)
        C = sec.setDampMatrix(cxx=c, cyy=c)
        K = sec.setStifMatrix(kxx=k, kyy=k)
        sec.polarFromCSV(polarFilename=os.path.join(scriptDir,'../data/tjaere11_ds.csv'))
        sec.setDynStall(ds_model='oye')
        sec.setParameters(sx_sim='x')
        sec.setConstantInputs(Ux, Uy, theta_p=theta_p*np.pi/180)

        # --- Time integration
        fs_i= 0.88792530768895 # TODO
        t = np.linspace(0, 2, 7*60)
        q0 = sec.setInitialConditions(x=x0, fs=fs_i)
        res, df = sec.integrate(t, y0=q0)
        np.testing.assert_almost_equal(df['x1'].values[-1], 0.035058, 5)
        np.testing.assert_almost_equal(df['x2'].values[-1], 0.069730, 5)
        np.testing.assert_almost_equal(df['x3'].values[-1], 0.779248, 5)
        #sec.plot_states()
        #import matplotlib.pyplot as plt
        #plt.show()

        # --- Frequencies
        u0 = [Ux, Uy, theta_p*np.pi/180]
        xop = sec.equilibrium(x0=q0, u0=u0)
        freq_d, zeta, Q, freq_0 = sec.eigA(x0=xop, u0=u0)
        #print('Frequencies',freq_0)
        #print('Damping    ',zeta*100)
        # print('Q\n',np.abs(Q))
        # TODO TODO TODO FIGURE OUT which one, adapt _dx
        #np.testing.assert_almost_equal(freq_0[0]  , 1.207931, 5)
        #np.testing.assert_almost_equal(zeta[0]*100, 1.080974, 5)

        np.testing.assert_almost_equal(freq_0[0]  , 1.24125, 5)
        np.testing.assert_almost_equal(zeta[0]*100, 5.64522, 5)


    def test_section_dyninflow(self):
        from welib.dyninflow.DynamicInflow import dyninflow_oye_sim, tau1_oye, tau2_oye
        # Parameters
        U0     = 10 # mean wind speed [m/s]
        R      = 65  # rotor radius [m]
        r_bar  = 0.5 # spanwise station investigated [-]
        a1    = 0.1 # axial induction used before the step
        a2    = 0.2 # axial induction used after the step
        tstep = 10  # time at which the step of axial induction occurs
        tmax  = 20  # simulation time
        time = np.linspace(0,tmax,1000)

        # --- DYNINFLOW STANDALONE 
        p=dict()
        p['tau1'] = tau1_oye(a1, R, U0)
        p['tau2'] = tau2_oye(r_bar ,p['tau1'])
        p['k']    = 0.6

        # Define functions of time (NOTE: these functions would be set within an unsteady BEM algorithm)
        # Look at the function dyninflow_oye_sim for lower level interfaces
        u=dict()
        u['Vqs']  = lambda t: a1*U0 if t<tstep else a2*U0  
        df_c = dyninflow_oye_sim(time, u, p, x0=None, prefix='', method='continuous')

        # --- USING SECTION
        # Inputs
        uop = [U0, 0, 0]
        # TODO interface might change
        sec = Section()
        sec._chord = 0.2 
        sec._rho = 1.225
        sec._r_bar = r_bar
        sec._tau = 0.08 # DS Oye TODO
        sec._R = R
        sec._y_AQ = 0
        sec._y_AT = 0
        sec.polarFromCSV(polarFilename=os.path.join(scriptDir,'../data/tjaere11_ds.csv'))
        sec.setDynInflow(di_model='oye', a0=a1, ap0=0)
        sec.setConstantInputs(uop[0], uop[1], uop[2])
        sec.u_sim['Ux'] = lambda t : U0 if t<tstep else U0*a2/a1
        sec.setParameters(sx_sim='')
        # Integration
        q0 = sec.setInitialConditions(wxr=0, wx=0, equilibrium=False, di_eq=True, ds_eq=True, uop=uop, verbose=False)
        res, df = sec.integrate(time, y0=q0, method='LSODA')

        ## Plot
        #fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
        #fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
        #ax.plot(time, df_c['Vdyn_[m/s]']/U0, 'k-' , label=r'$a_{dyn}$ (continuous formulation)')
        #ax.plot(time, df_c['Vqs_[m/s]']/U0 , 'k-', label=r'$a_{qs}$')
        #ax.plot(time,-df  ['Wqs_x_AC']/U0 , '--', label=r'$a_{qs}$ (section)')
        #ax.plot(time,-df  ['Wdyn_x_AC']/U0 , ':', label=r'$a_{dyn}$ (section)')
        #ax.legend()
        #ax.set_xlabel('Time [s]')
        #ax.set_ylabel('Axial induction [-]')
        #ax.set_ylim([0.09,0.21])
        #ax.set_title('Dynamic Inflow (Oye) - induction step')
        #plt.show()

        np.testing.assert_array_almost_equal(-df['Wdyn_x_AC'], df_c['Vdyn_[m/s]'], 5)
        np.testing.assert_array_almost_equal(-df['Wqs_x_AC'],  df_c['Vqs_[m/s]'], 5)




if __name__ == '__main__':
    unittest.main()
