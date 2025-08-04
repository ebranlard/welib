""" 
Test the different panel codes on lifting geometries
"""

import numpy as np
import matplotlib.pyplot as plt

from welib.vortilib.panelcodes.panel_tools import *
from welib.vortilib.panelcodes.panel_examples import getCase

from welib.vortilib.panelcodes.pointVortices import PV_solve
from welib.vortilib.panelcodes.constantSourceVortexPanels import CSPN_CVP1_solve
from welib.vortilib.panelcodes.linearVortexPanels import LVP_solve

def compare(cases = None):
    if cases is None:
        cases=['VonDeVooren_lift', 'NACA2412_lift']

    for case in cases:

        cas = getCase(case, mid_out=True)
        XP, YP = cas['XP'], cas['YP']
        Uxy = cas['Uxy']

        outPV = PV_solve(XP, YP, Uxy=Uxy)
        outCSPN= CSPN_CVP1_solve(XP, YP, Uxy=Uxy)
        outLVP = LVP_solve(XP, YP, Uxy=Uxy)

        ax = None
        ax = plot_Cp(outCSPN['CP'][:,0], outCSPN['Cp'], label='Constant Source and Vortex Panels', ax=ax, sty='d', Cp_ref=cas['Cp'] )
        ax = plot_Cp(outLVP['CP'][:,0], outLVP['Cp'], label='Linear Vortex Panel', ax=ax, sty='^' )
        ax = plot_Cp(outPV['CP'][:,0], outPV['Cp'], label='Point Vortex', ax=ax, sty='.')


if __name__ == '__main__':
    #compare(cases=['VonDeVooren'])
    compare(cases=None)
    plt.show()


if __name__ == '__test__':
    compare(cases=None)
