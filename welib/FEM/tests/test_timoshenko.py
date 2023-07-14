import unittest
import os
import numpy as np
from welib.FEM.timoshenko import *
from welib.FEM.frame3d import *


MyDir=os.path.dirname(__file__)

class Test(unittest.TestCase):

    def test_compare_frame(self):
        # Test that frame3d and timoshenko are consistent (shear=False, kappa=0, and crossterms 0 in mass matrix)
        # By default frame3d has main_axis x
        # By default timosh  has main_axis z
        kappa=0

        L    = 100                       # Beam Length [m]
        G    = 79.3e9                    # Shear modulus. Steel: 79.3  [Pa] [N/m^2]
        E    = 210e9                     # Young modulus [Pa] [N/m^2]
        D    = 8
        t    = 0.045
        A   = np.pi*( (D/2)**2 - (D/2-t)**2)
        rho  = 7850
        m0   = rho*A
        EA  = E*A
        # Data when main axis is x
        Ixx_x = np.pi/32*(D**4-(D-2*t)**4)# Polar second moment of area [m^4]
        Iyy_x = np.pi/64*(D**4-(D-2*t)**4)# Planar second moment of area [m^4]
        Izz_x = Iyy_x
        Kt   = Ixx_x # Kt = Polar for circular sections
        # Data when main axis is z
        Izz_z = np.pi/32*(D**4-(D-2*t)**4)# Polar second moment of area [m^4]
        Ixx_z = np.pi/64*(D**4-(D-2*t)**4)# Planar second moment of area [m^4]
        Iyy_z = Ixx_z


        np.set_printoptions(linewidth=300, precision=3)
        # Frame 3d
        ke_f_x, me_f_x, _ = frame3d_KeMe(E, G, Kt, EA, E*Ixx_x, E*Iyy_x, E*Izz_x, L, A, rho*A*L, main_axis='x') 
        ke_f_z, me_f_z, _ = frame3d_KeMe(E, G, Kt, EA, E*Ixx_z, E*Iyy_z, E*Izz_z, L, A, rho*A*L, main_axis='z') 
        # Timoshenko
        ke_t_x, me_t_x    = timoshenko_KeMe(L, A, Ixx_x, Iyy_x, Izz_x, kappa, E=E, G=G, rho=rho, main_axis='x', shear=False, crossterms=False)
        ke_t_z, me_t_z    = timoshenko_KeMe(L, A, Ixx_z, Iyy_z, Izz_z, kappa, E=E, G=G, rho=rho, main_axis='z', shear=False, crossterms=False)

        # Compare stiffnesses in both coordinate systems
        np.testing.assert_almost_equal(ke_t_x/1e7, ke_f_x/1e7, 7)
        np.testing.assert_almost_equal(ke_t_z/1e7, ke_f_z/1e7, 7)

        np.testing.assert_almost_equal(me_t_x/1e5, me_f_x/1e5, 7)
        np.testing.assert_almost_equal(me_t_z/1e5, me_f_z/1e5, 7)

#         print('me_f_z\n',me_f_z)
#         print('me_t_z\n',me_t_z)
#         print('me_f_z-me_t_z\n',me_f_z-me_t_z)
#         print('me_f_z\n',me_f_z)
#         print('me_t_z\n',me_t_z)
#         print('me_f_z-me_t_z\n',me_f_z-me_t_z)


if __name__=='__main__':
    unittest.main()
