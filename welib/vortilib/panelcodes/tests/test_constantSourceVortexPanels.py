import unittest
import numpy as np    
import os as os
from welib.CFD.flows2D import *
from welib.vortilib.panelcodes.constantSourceVortexPanels import *
from welib.vortilib.panelcodes.panel_examples import getCase
from welib.vortilib.panelcodes.panel_tools import *

scriptDir = os.path.dirname(__file__)

class Test(unittest.TestCase):

    def test_naca2412(self, plot=False):

        # --- Case/ Geometry
        cas = getCase('NACA2412_lift')

        # --- Panel method
        out = CSPN_CVP1_solve(cas['XP'], cas['YP'], Uxy=cas['Uxy'], alpha=cas['alpha'])

        np.testing.assert_allclose(cas['M']     ,out['M'])
        np.testing.assert_allclose(cas['rhs']   ,out['rhs'])
        np.testing.assert_allclose(cas['sigmas'],out['sigmas'])
        np.testing.assert_allclose(cas['gamma'] ,out['gamma'])
        np.testing.assert_allclose(cas['Cp']    ,out['Cp'])
        np.testing.assert_allclose(cas['Cl']    ,out['Cl'])
        np.testing.assert_allclose(cas['Cm']    ,out['Cm'])

        if plot:
            plot_pressure_force_bars(out['CP'], out['Cp'], out['n'], scale=0.1, ax=None)
            plot_Cp(out['CP'], out['Cp'], Cp_ref=cas['Cp'])

            # --- Flow field
            #PP  = np.column_stack((cas['XP'], cas['YP']))
            #vel = lambda X, Y : CSPN_CVP1_u(X, Y, PP, out['sigmas'], -out['gamma']) # TODO gamma sign issue
            #X, Y, U, V =  flowfield2D(vel, xmax=1.5, xmin=-0.5, ymin=-0.3, ymax=0.3, nx=140, ny=120, Uxy=cas['Uxy'], rel=False)
            ##bIn = points_inside(XP, YP, X, Y)
            ##U[bIn] = 0
            ##V[bIn] = 0

            #ax =  flowfield2D_plot(X, Y, U, V, ax=None, minVal=0, maxVal=2, bounded=False, rel=False)
            #ax.fill(cas['XP'], cas['YP'], 'k')
            #ax.set_aspect('equal')

            ## --- Cp field
            #Speed = np.sqrt(U**2 + V**2) 
            #CpXY = 1 - Speed**2/cas['U0']**2 
            #ax =  flowfield2D_plot(X, Y, U, V, Speed=CpXY, ax=None, nLevels=31, minVal=-2, maxVal=1, bounded=False, rel=False, cmap='jet')
            #ax.fill(cas['XP'], cas['YP'], 'k')
            #ax.set_aspect('equal')

            plt.show()
        

if __name__ == '__main__':
    Test().test_naca2412(plot=True)
    unittest.main()
