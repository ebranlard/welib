import unittest

import numpy as np
from welib.yams.yams_sympy_tools import *
from sympy import symbols, diff, cos, sin, exp
from sympy.physics.vector import dynamicsymbols

# --------------------------------------------------------------------------------}
# --- TESTS
# --------------------------------------------------------------------------------{
class TestYAMSSPTools(unittest.TestCase):
    def test_subs_no_diff(self):
        # --- Test substitution without replacing time dderivatives
        x,y = dynamicsymbols('x, y')
        z   = symbols('z')
        a,b = symbols('a,b')
        time = dynamicsymbols._t

        # Function of x
        expr      = a*cos(x) + b*cos(x)*diff(x,time)+ cos(x)**2 * diff(diff(x,time),time)**2
        expr0_ref =    a +     b*       diff(x,time)+      1    * diff(diff(x,time),time)**2
        expr0     = subs_no_diff(expr, [(x,0),(y,0)] )
        self.assertEqual(expr0, expr0_ref)


        # One parameter is not a function
        expr      = a*cos(x) + b*cos(x)*diff(x,time)+ cos(x)**2 * diff(diff(x,time),time)**2
        expr0_ref =    a +                          +      1    * diff(diff(x,time),time)**2
        expr0     = subs_no_diff(expr, [(x,0),(b,0)] )
        self.assertEqual(expr0, expr0_ref)

        # two parameters as function of time
        expr      = a*diff(y,time)*cos(x) + y*cos(x)*diff(x,time)+ cos(x)**2 * diff(diff(x,time),time)**2
        expr0_ref = a*diff(y,time)        + b       *diff(x,time)+      1    * diff(diff(x,time),time)**2
        expr0     = subs_no_diff(expr, [(x,0),(y,b)] )
        self.assertEqual(expr0, expr0_ref)

    def test_linearization(self):
        # --- Test linearization without replacing time derivatives
        x,y = dynamicsymbols('x, y')
        z   = symbols('z')
        a,b,c,d = symbols('a,b,c,d')

        # Linearization 0th order of sin
        expr      = a*sin(b*x)
        expr0_ref =    0
        expr0     = linearize(expr, [(x,0)], order=0)
        self.assertEqual(expr0, expr0_ref)

        # Linearization 1st order of sin
        expr      = a*sin(b*x)
        expr0_ref =    a*b*x 
        expr0     = linearize(expr, [(x,0)], order=1)
        self.assertEqual(expr0, expr0_ref)

        # Linearization 2nd order of sin
        expr      = a*sin(b*x)
        expr0_ref =    a*b*x  - a*b**3*x**3/6
        expr0     = linearize(expr, [(x,0)], order=3)
        self.assertEqual(expr0, expr0_ref)

        # Linearization 0th order of cos
        expr      = a*cos(b*x)
        expr0_ref =    a
        expr0     = linearize(expr, [(x,0)], order=0)
        self.assertEqual(expr0, expr0_ref)

        # Linearization 1st order of cos
        expr      = a*cos(b*x)
        expr0_ref =    a
        expr0     = linearize(expr, [(x,0)], order=1)
        self.assertEqual(expr0, expr0_ref)

        # Linearization 2nd order of cos
        expr      = a*cos(b*x)
        expr0_ref =    a       - a*b**2*x**2/2
        expr0     = linearize(expr, [(x,0)], order=2)
        self.assertEqual(expr0, expr0_ref)


        # Linearization 1st order, mix of sin and cos
        expr      = a*cos(x) + b*sin(x) + c*sin(x)**2 + d*cos(y)*sin(y)
        expr0_ref =    a     + b*x                    + d*y
        expr0     = linearize(expr, [(x,0),(y,0)] )
        self.assertEqual(expr0, expr0_ref)

        # Linearization polynomial vars is polynomial
        expr      = a * (x +y)**2
        expr0_ref = a * x**2 + 2* a*x*y + a* y**2
        expr0     = linearize(expr, [(x,0),(y,0)], order=3 )
        self.assertEqual(expr0, expr0_ref)

        # Linearization polynomial vars is polynomial
        expr      = a * exp(x +y)
        expr0_ref = a*x + a*y + a
        expr0     = linearize(expr, [(x,0),(y,0)], order=1 )
        self.assertEqual(expr0, expr0_ref)

    def test_linearization_with_diff(self):

        x,y      = dynamicsymbols('x, y')
        a,b,c,d  = symbols('a,b,c,d')
        theta    = dynamicsymbols('theta')
        T0  = symbols('T0')
        TD0  = symbols('TD0')
        time = dynamicsymbols._t
        xdot = diff(x,time)
        thetadot = diff(theta,time)

        # Linearization 1st order with time derivatives
        expr      = a*cos(x) + b*cos(x)*diff(x,time)+ cos(x)**2 * diff(diff(x,time),time)**2
        expr0_ref =    a +     b*       diff(x,time)+      1    * diff(diff(x,time),time)**2
        expr0     = linearize(expr, [(x,0)] )
        self.assertEqual(expr0, expr0_ref)

        expr      = a*sin(y) + b*sin(x)*diff(x,time)+ sin(x)**2 * diff(diff(y,time),time)**2
        expr0_ref = a* y     + b*x     *diff(x,time)+   0
        expr0     = linearize(expr, [(x,0),(y,0)], order=1 )
        self.assertEqual(expr0, expr0_ref)

        # Linearization polynomial vars is polynomial
        expr      = 3*thetadot*theta
        expr0_ref = 3
        expr0_ref = 3*T0*TD0 + 3*T0*( thetadot-TD0) + 3*TD0*(theta-T0)
        expr0     = linearize(expr, [(theta, T0), (thetadot, TD0)], order=1 )
        self.assertEqual(expr0, expr0_ref)

        # ISSUE when you do a substitution within a derivative, the substitution of derivative does not work
        expr      = thetadot**3
        expr0_ref = TD0**3 + 3*TD0**2* (thetadot-TD0)
        expr02_ref =-TD0**3 - 3*TD0**2* (thetadot-TD0)
        expr2     = (thetadot**3).subs([(theta, -theta)])
        expr0  =  linearize(expr  ,  [(thetadot,TD0)], order=1 )
        expr02 =  linearize(expr2 ,  [(thetadot,TD0)], order=1 )
        self.assertEqual(expr0, expr0_ref)
        self.assertEqual(expr02, expr02_ref)



    def test_jacobian(self):

        from welib.yams.yams_sympy_tools import myjacobian
        from sympy import Symbol, Matrix
        x,y   = symbols('x, y')
        x0,y0 = symbols('x0, y0')
        a,b   = symbols('a, b')

        f1= a*x + b*y**2
        f2= a*x + b*x**2
        f = Matrix([[f1],[f2]])
        J = myjacobian(f, [x,y], [x0,y0])
        self.assertEqual(J[0,0], a)
        self.assertEqual(J[1,0], a+2*b*x0)

        J = myjacobian(f1, [x,y], [x0,y0])
        self.assertEqual(J[0,0], a)


if __name__=='__main__':
    #TestYAMSSPTools().test_jacobian()
    #TestYAMSSPTools().test_linearization_with_diff()
    unittest.main()
