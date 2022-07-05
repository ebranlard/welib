""" 
Analytical equation of motions for a wind turbine using 1 degree of freedom:
 - flexible fore-aft tower

""" 
import numpy as np
import unittest

from sympy import Symbol
from sympy.physics.mechanics import dynamicsymbols
from sympy.parsing.sympy_parser import parse_expr
from welib.yams.models.FTNSB_sympy import get_model
from welib.yams.models.FTNSB_sympy_symbols import *

def main(unittest=False):

#     'rot_elastic_type':'Body',     #<<< Very important, SmallRot, or Body, will affect the rotation matrix
#     'rot_elastic_subs':True,       #<<< Very important, will substitute alpha_y with nuy q. Recommended True
#     'rot_elastic_smallAngle':False,#<<< Very important, will perform small angle approx: sin(nu q) = nu q and nu^2=0 !!! Will remove all nu^2 and nu^3 terms!! might not be recommended
#     'orderMM':2, #< order of taylor expansion for Mass Matrix
#     'orderH':2,  #< order of taylor expansion for H term


    model = get_model('F0T1RNA', mergeFndTwr=True, yaw='zero', tilt='fixed', tiltShaft=False,
            rot_elastic_type='SmallRot',
            rot_elastic_smallAngle=True,
            aero_torques=True
            )
    model.kaneEquations()
    extraSubs=model.shapeNormSubs

    model.smallAngleApprox(model.smallAngles, extraSubs)

    return model



class TestF0T1RNA(unittest.TestCase):
    def test_F0T1RNA(self):
        #Test expression of mass and forcing

        model=main(unittest=True)
        twr=model.twr


        import sys
        if sys.version_info.major < 3 :
            print('>>>> Skipping test for python 2, due to sign error')
            return

        # --- Mass matrix
        MM = model._sa_mass_matrix
        MM_noTwr= parse_expr('Matrix([[-2*M_RNA*v_yT1c*(x_RNAG*sin(theta_tilt) - z_RNAG*cos(theta_tilt)) + M_RNA]])')# M_e^0_T_11 + M_e^1_1_T_11*q_T1(t)
        MM_twr  = model.twr.MM 
        #print('MM',MM)
        #print('MM_noT',MM_noTwr)
        #print('MM_twr',MM_twr)
        #MM_twr  = twr.Me.eval(twr.q)
        DMM = MM-MM_twr-MM_noTwr
        DMM.simplify()
        #print('DMM')
        self.assertEqual(DMM[0,0],0)

        # forcing
        FF = model._sa_forcing[0,0].subs([(dynamicsymbols('g'),Symbol('g')),(dynamicsymbols('T_a'),Symbol('T_a')), ])
        #forcing1=parse_expr('-D_e^0_T_11*Derivative(q_T1(t), t) - K_e^0_T_11*q_T1(t)')
        forcing1 = twr.bodyElasticForce(twr.q, twr.qdot)[6,0]
        forcing2 = parse_expr('T_a*cos(theta_tilt) + v_yT1c*(M_RNA*g*x_RNAG*cos(theta_tilt) + M_RNA*g*z_RNAG*sin(theta_tilt) + T_a*z_NR - T_a*q_T1(t)*sin(theta_tilt) + M_y_a(t))')
        #print('ff       ',FF)
        #print('forcing 1',forcing1)
        #print('forcing 2',forcing2)

        DFF = FF-forcing2+forcing1
        DFF = DFF.simplify()
        #print('DFF',DFF)

        self.assertEqual(DFF,0)



if __name__=='__main__':
    np.set_printoptions(linewidth=500)
    #main(unittest=False)
    unittest.main()
    
