import unittest
import numpy as np    
import os
from welib.CFD.flows2D import *
from welib.vortilib.panelcodes.linearVortexPanels import *
from welib.vortilib.panelcodes.panel_tools import *
from welib.vortilib.panelcodes.panel_examples import *


scriptDir = os.path.dirname(__file__)

class Test(unittest.TestCase):

    def test_vondevooren(self, plot=False):
        # --- Parameters
        cas = getCase('VonDeVooren_lift', solver='LVP')

        # --- Panel method
        out = LVP_solve(cas['XP'], cas['YP'], Uxy=cas['Uxy'], closed=True, alpha=cas['alpha'])

        if plot:
            plot_pressure_force_bars(out['CP'], out['Cp'], out['n'], scale=0.1, ax=None)
            plot_Cp(out['CP'][:,0], out['Cp'], Cp_ref=cas['Cp'])
            plt.show()

        # --- Test
        np.testing.assert_allclose(out['Cp'], cas['Cp'], 3)
        

if __name__ == '__main__':
    #Test().test_vondevooren(plot=True)
    unittest.main()
