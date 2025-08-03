import unittest
import numpy as np    
import os as os
import matplotlib.pyplot as plt

from welib.CFD.flows2D import flowrate2D, flowfield2D, flowfield2D_plot
from welib.vortilib.panelcodes.panel_tools import plot_line
from welib.vortilib.panelcodes.panel_examples import getCase
from welib.vortilib.panelcodes.pointSources import PS_solve, PS_velocity
from welib.vortilib.panelcodes.pointSources_discontinuous import PS_disc_solve

scriptDir = os.path.dirname(__file__)

class Test(unittest.TestCase):

    def test_cylinder_CpUt(self, plot=False):
    
        # --- Case/ Geometry
        for case in ['cylinder', 'ellipse']:
            for method in ['CSP', 'DSP']:

                # --- Parameters
                offset = 2.0
                cas = getCase(case, U0=1.2)
                R = cas['R']

                # --- Panel method
                if method=='CSP':
                   out = PS_solve(cas['XP'], cas['YP'], Uxy=cas['Uxy'], offset=offset)
                else:
                    SP1 = np.column_stack([cas['XP'][0:-1], cas['YP'][0:-1]])
                    SP2 = np.column_stack([cas['XP'][1:]  , cas['YP'][1:]])
                    out = PS_disc_solve(SP1, SP2, Uxy=cas['Uxy'], offset=offset)
           
                # --- Tests
                #print('MaxError ',np.max(np.abs(Cp_theory-out['Cp'])))
                np.testing.assert_almost_equal(out['Cp'], cas['Cp'], 5)
                np.testing.assert_almost_equal(out['Ut'], cas['Ut'], 5)
            
                # --- Plots   
                if plot:
                    fig, axes = plt.subplots(1, 4, figsize=(10, 3.5))
                    fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.40)

                    ax = plot_line(cas['XP'], cas['YP'], Uwall=out['Vwall'], ntScale=0.5, UScale=0.8, ax=axes[0])
                    ax.plot(out['SP'][:,0], out['SP'][:,1], '+')

                    ax=axes[1]
                    ax.plot(out['theta_CP'], out['Un'], label='Un')
                    ax.plot(out['theta_CP'], out['Ut'], label='Ut')
                    ax.plot(out['theta_CP'], cas['Ut'], 'k.', label='Theory')

                    ax.legend()
                    ax=axes[2]
                    ax.plot(out['theta_CP'], out['Cp'], label='Cp')
                    ax.plot(out['theta_CP'], cas['Cp'], 'k.', label='Theory')
                    ax.legend()
                    
                    # --- Flow field
                    vel = lambda X, Y : PS_velocity(X, Y, out['SP'], out['Sigmas'])
                    X, Y, U, V =  flowfield2D(vel, xmax=3.5, nx=50, U0x=cas['U0'], L=R, rel=True)
                    if case=='cylinder':
                        RR =np.sqrt(X**2+Y**2)
                        U[RR<1]=0
                        V[RR<1]=0
                    ax =  flowfield2D_plot(X, Y, U, V, ax=axes[3], minVal=0, maxVal=cas['maxVal'], bounded=False, rel=True)
                    ax.plot(cas['XP']/R, cas['YP']/R, 'k-',lw=3)
                    #ax.set_title('Source panel, cylinder')

        if plot:
            plt.show()




if __name__ == '__main__':
    #Test().test_cylinder_CpUt(plot=True)
    unittest.main()
