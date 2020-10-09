import unittest
import sys
#sys.path.insert(0,'..')
from welib.system.system import *
from welib.system.statespacelinear import lti_state_space_function, lti_output_function

import random


# --------------------------------------------------------------------------------}
# --- TESTS
# --------------------------------------------------------------------------------{
class Test(unittest.TestCase):
    def get_LTI(self, nx=2, nu=3):
        random.seed(0)
        p=dict()
        p['A'] = np.random.randint(5, size=(nx, nx)).astype(float)
        p['B'] = np.random.randint(5, size=(nx, nu)).astype(float)
        p['C'] = np.random.randint(5, size=(nx, nx)).astype(float)
        p['D'] = np.random.randint(5, size=(nx, nu)).astype(float)
        
        x0 = np.random.randint(5, size=(nx, 1)).astype(float)
        u0 = np.random.randint(5, size=(nu, 1)).astype(float)

        return p, x0, u0



    def test_system_implicit(self):
        # --- test that implicit function returns 0 at a solved operating point from explicit model
        # Create some simple state space 
        nx = 5
        p, x0, u0 = self.get_LTI(nx)
        # Initial system
        sys = System(Fx=lti_state_space_function, Y=lti_output_function, interface='xup', param = p )

        # Compute RHS
        xdot0=lti_state_space_function(0,x0,u0,p) 

        # Generate implicit function
        F=sys.implicit_function

        # Check that residual is zeros 
        residual = F(0,xdot0,x0,u0,p)
        np.testing.assert_almost_equal(residual, np.zeros((nx,1)))

    def test_system_linearization_lti_statespace(self):
        # --- Test the general "system" linearization using a LTI state space system
        
        # Create some simple state space 
        p, x0, u0 = self.get_LTI()
        sys = System(Fx=lti_state_space_function, Y=lti_output_function, interface='xup', param = p )

        delta_x  = np.zeros(x0.shape)+0.001
        delta_xd = np.zeros(x0.shape)+0.001
        delta_u  = np.zeros(u0.shape)+0.001

        # --- Test implicit linearization feature of an explicit system
        op=(0,x0,u0)
        A,B,C,D=sys.linearize(op, dx=delta_x, dxp=delta_xd, du=delta_u, use_implicit=True)
        np.testing.assert_almost_equal(A, p['A'])
        np.testing.assert_almost_equal(B, p['B'])
        np.testing.assert_almost_equal(C, p['C'])
        np.testing.assert_almost_equal(D, p['D'])

        # --- Test explicit linearization feature of an explicit system
        A,B,C,D=sys.linearize(op, dx=delta_x, du=delta_u)
        np.testing.assert_almost_equal(A, p['A'])
        np.testing.assert_almost_equal(B, p['B'])
        np.testing.assert_almost_equal(C, p['C'])
        np.testing.assert_almost_equal(D, p['D'])
        #print('A')
        #print(A)
        #print('B')
        #print(B)
        #print('C')
        #print(C)
        #print('D')
        #print(D)


if __name__=='__main__':
    unittest.main()
