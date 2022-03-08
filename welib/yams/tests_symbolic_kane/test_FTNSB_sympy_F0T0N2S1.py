""" 
Analytical equation of motions for a wind turbine using 3 degrees of freedom:
 - nacelle yaw 
 - nacelle tilt
 - shaft rotation psi

Equations should match a gyroscope/top


NOTE: test not ready

""" 
import numpy as np
import unittest
from sympy import Symbol
from sympy.parsing.sympy_parser import parse_expr
from welib.yams.models.FTNSB_sympy import get_model
from welib.yams.models.FTNSB_sympy_symbols import *

def main(unittest=False):

    model = get_model('F0T0N2S1', mergeFndTwr=True)
    N=model.nac
    R=model.rot
    N.noInertia()
    N.noMass()
    #R.origin.set_pos(N.origin, 0*N.frame.x)

    model.kaneEquations()

    #print(model.kane.mass_matrix[0,:])
    #print(model.kane.mass_matrix[1,:])
    #print(model.kane.mass_matrix[2,:])
    #model.smallAngleApprox(l.smallAngles, extraSubs)
    return model



class TestF0T0N2S1(unittest.TestCase):
    def test_F0T0N2S1(self):
        #Test expression of mass and forcing
        model=main(unittest=True)
        print('Gyroscope test not ready')

if __name__=='__main__':
    np.set_printoptions(linewidth=500)
    #main(unittest=False)
    unittest.main()
    
