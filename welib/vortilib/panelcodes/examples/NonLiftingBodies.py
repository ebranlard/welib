""" 
Test the different panel codes on non-lifting geometries
"""

import numpy as np
import matplotlib.pyplot as plt

from welib.vortilib.panelcodes.panel_tools import *
from welib.vortilib.panelcodes.panel_examples import getCase

from welib.vortilib.panelcodes.pointSources import PS_solve
from welib.vortilib.panelcodes.pointVortices import PV_solve
from welib.vortilib.panelcodes.constantSourcePanels import CSP_solve
from welib.vortilib.panelcodes.constantVortexPanels import CVP_solve
from welib.vortilib.panelcodes.constantSourceVortexPanels import CSPN_CVP1_solve
from welib.vortilib.panelcodes.linearVortexPanels import LVP_solve

def compare(cases = None):
    if cases is None:
        cases=['cylinder', 'ellipse']

    for case in cases:

        cas = getCase(case)
        XP, YP = cas['XP'], cas['YP']
        Uxy = cas['Uxy']

        outPS = PS_solve(XP, YP, Uxy=Uxy, offset=2)
        outPV = PV_solve(XP, YP, Uxy=Uxy, hasLift=False)
        outCSP = CSP_solve(XP, YP, Uxy=Uxy)
        outCVP = CVP_solve(XP, YP, Uxy=Uxy)
        outCSPN= CSPN_CVP1_solve(XP, YP, Uxy=Uxy)
        outLVP = LVP_solve(XP, YP, Uxy=Uxy)

        ax = None
        ax = plot_Cp(outCSP['CP'][:,0], outCSP['Cp'], label='Constant Source Panels', ax = ax, sty='o', Cp_ref=cas['Cp'] )
        ax = plot_Cp(outPS['CP'][:,0], outPS['Cp'], label='Point Sources', ax=ax, sty='>')
        ax = plot_Cp(outCSPN['CP'][:,0], outCSPN['Cp'], label='Constant Source and Vortex Panels', ax=ax, sty='d' )
        ax = plot_Cp(outCVP['CP'][:,0], outCVP['Cp'], label='Constant Vortex', ax=ax, sty='d')
        ax = plot_Cp(outLVP['CP'][:,0], outLVP['Cp'], label='Linear Vortex Panel', ax=ax, sty='^' )
        ax = plot_Cp(outPV['CP'][:,0], outPV['Cp'], label='Point Vortex', ax=ax, sty='.')
        ax.set_title(case)


if __name__ == '__main__':
    #compare(cases=['ellipse'])
    compare(cases=None)
    plt.show()


if __name__ == '__test__':
    compare(cases=None)
