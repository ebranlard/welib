import unittest
import sys
#sys.path.insert(0,'..')
from welib.system.linearization import *


    
# --------------------------------------------------------------------------------}
# --- TESTS
# --------------------------------------------------------------------------------{
class Test(unittest.TestCase):
    def test_numerical_jacobian(self):
        # --- Setup of same functions and their analytical jacobians
        def f_test(x, u):
            xdot=np.zeros(len(x))
            xdot[0] = x[0]**2 + x[2] + 3*u[0] 
            xdot[1] = x[1]**3            
            xdot[2] =             u[2]**3
            return xdot
        def f_testp(x, u, p):
            xdot=np.zeros(len(x))
            xdot[0] = x[0]**2 + x[2] + 3*u[0] 
            xdot[1] = x[1]**3                 + p[0]
            xdot[2] =             u[2]**3   + p[1]**2
            return xdot
        def jacx_test(x, u):
            jacx=np.zeros((x.shape[0],x.shape[0]))
            jacx[0,0] = 2*x[0];  jacx[0,2] = 1 
            jacx[1,1] = 3*x[1]**2; 
            return jacx
        def jacu_test(x, u):
            jacu=np.zeros((x.shape[0],u.shape[0]))
            jacu[0,0] = 3
            jacu[2,2] = 3*u[2]**2
            return jacu
        # Define an operating point
        p=np.zeros(2)
        x0=np.zeros((3,1))
        u0=np.zeros((3,1))
        x0[:,0]=[1,2,3]
        u0[:,0]=[10,20,30]
        dx=[0.01]*3
        du=[0.01]*3
        p[0]=100
        p[1]=2

        # --- Test of x jacobian
        jacx_ref = (jacx_test(x0,u0))
        jacx_num = numerical_jacobian(f_test , (x0,u0), 0, dx   )
        jacx_nump= numerical_jacobian(f_testp, (x0,u0), 0, dx, p)
        np.testing.assert_almost_equal(jacx_ref,jacx_num , 4)
        np.testing.assert_almost_equal(jacx_ref,jacx_nump, 4)
        #print(jacx_ref)
        #print(jacx_num)
        #print(jacx_nump)

        # --- Test of u jacobian
        jacu_ref = jacu_test(x0,u0)
        jacu_num = numerical_jacobian(f_test , (x0,u0), 1, du )
        jacu_nump= numerical_jacobian(f_testp, (x0,u0), 1, du, p)
        np.testing.assert_almost_equal(jacu_ref,jacu_num , 4)
        np.testing.assert_almost_equal(jacu_ref,jacu_nump, 4)
        #print(jacu_num)
        #print(jacu_nump)

    def test_linearizeFunction(self):
        # Test wrapper funciton which linearize with respect to any arguments
        A=np.array([[0,1],[2,3]])
        B=np.array([[5], [6]])

        def F(x,u,p=None):
            return A.dot(x) + B.dot(u)

        x0=[0,0]
        u0=[0]
        dx=[0.01]*2
        du=[0.01]*1
        A0,B0  = linearize_function(F, (x0, u0), [0,1], (dx, du))
        np.testing.assert_almost_equal(A,A0, 4)
        np.testing.assert_almost_equal(B,B0, 4)
 
        A0,B0 = linearize_Fxu(F, x0, u0, dx, du)
        np.testing.assert_almost_equal(A,A0, 4)
        np.testing.assert_almost_equal(B,B0, 4)
 

        def F(x,p=None):
            return A.dot(x)
        A0 = linearize_Fx(F, x0, dx)
        np.testing.assert_almost_equal(A,A0, 4)




if __name__=='__main__':
    unittest.main()
