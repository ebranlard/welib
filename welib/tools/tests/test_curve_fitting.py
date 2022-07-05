import unittest
import numpy as np
from welib.tools.curve_fitting import *
from welib.tools.curve_fitting import sinusoid

class TestFitting(unittest.TestCase):

    # --------------------------------------------------------------------------------}
    # --- Testing of high level functions  
    # --------------------------------------------------------------------------------{
    def test_fit_gaussian(self):
        mu,sigma,y0=0.5,1.2,10
        x=np.linspace(0,1,10)
        y=gaussian_w_offset(x,(mu,sigma,y0))
        y_fit, pfit, fitter = fit_gaussian(x, y, offset=True)
        #fitter.plot(); import matplotlib.pyplot as plt; plt.show()
        np.testing.assert_array_almost_equal(y,y_fit)
        np.testing.assert_almost_equal(mu   ,fitter.model['coeffs']['mu'])
        np.testing.assert_almost_equal(sigma,fitter.model['coeffs']['sigma'])
        np.testing.assert_almost_equal(y0   ,fitter.model['coeffs']['y0'])

    def test_fit_gaussian(self):
        mu,sigma,y0=0.5,1.2,10
        x=np.linspace(0,1,10)
        y=gaussian_w_offset(x,(mu,sigma,y0))
        y_fit, pfit, fitter = fit_gaussian(x, y, offset=True)
        #fitter.plot(); import matplotlib.pyplot as plt; plt.show()
        np.testing.assert_array_almost_equal(y,y_fit)
        np.testing.assert_almost_equal(mu   ,fitter.model['coeffs']['mu'])
        np.testing.assert_almost_equal(sigma,fitter.model['coeffs']['sigma'])
        np.testing.assert_almost_equal(y0   ,fitter.model['coeffs']['y0'])

    def test_fit_sinusoid(self):
        A     = 101
        B     = -200.5
        omega = 10
        phi   = np.pi/3
        x=np.linspace(0,10,300)
        y=sinusoid(x,(A,omega,phi,B))  #+ np.random.normal(0, 0.1, len(x))
        y_fit, pfit, fitter = fit_sinusoid(x, y)
        #fitter.plot(); import matplotlib.pyplot as plt; plt.show()
        np.testing.assert_array_almost_equal(y,y_fit)
        np.testing.assert_almost_equal(A    ,fitter.model['coeffs']['A'],1)
        np.testing.assert_almost_equal(B    ,fitter.model['coeffs']['B'])
        np.testing.assert_almost_equal(omega,fitter.model['coeffs']['omega'])
        np.testing.assert_almost_equal(phi  ,fitter.model['coeffs']['phi'])

    def test_fit_polynomial(self):
        x=np.linspace(0,1,10)
        y=x**2
        exponents=[0,1,2]
        y_fit, pfit, fitter = fit_polynomial(x, y, exponents=exponents)
        np.testing.assert_array_almost_equal(y,y_fit)
        np.testing.assert_almost_equal(1 , fitter.model['coeffs']['c'])
        np.testing.assert_almost_equal(0 , fitter.model['coeffs']['a'])

        y_fit, pfit, fitter = fit_polynomial(x, y, order=3)
        np.testing.assert_array_almost_equal(y,y_fit)
        np.testing.assert_almost_equal(1 , fitter.model['coeffs']['b'])
        np.testing.assert_almost_equal(0 , fitter.model['coeffs']['a'])

    # --------------------------------------------------------------------------------}
    # --- Testing of predefined models
    # --------------------------------------------------------------------------------{
    def test_gaussian(self):
        mu,sigma=0.5,1.2
        x=np.linspace(0,1,10)
        y=gaussian(x,(mu,sigma))
        y_fit, pfit, fitter = model_fit('predef: gaussian', x, y)
        #fitter.plot(); import matplotlib.pyplot as plt; plt.show()
        np.testing.assert_array_almost_equal(y,y_fit)
        np.testing.assert_almost_equal(mu   ,fitter.model['coeffs']['mu'])
        np.testing.assert_almost_equal(sigma,fitter.model['coeffs']['sigma'])

    def test_gaussian_w_offset(self):
        mu,sigma,y0=0.5,1.2,10
        x=np.linspace(-0,1,10)
        y=gaussian_w_offset(x,(mu,sigma,y0))
        y_fit, pfit, fitter = model_fit('predef: gaussian-yoff', x, y)
        #fitter.plot(); import matplotlib.pyplot as plt; plt.show()
        np.testing.assert_array_almost_equal(y,y_fit)
        np.testing.assert_almost_equal(mu   ,fitter.model['coeffs']['mu'])
        np.testing.assert_almost_equal(sigma,fitter.model['coeffs']['sigma'])
        np.testing.assert_almost_equal(y0   ,fitter.model['coeffs']['y0'])

    def test_powerlaw_alpha(self):
        u_ref,z_ref,alpha=20,12,0.12
        x = np.linspace(0,1,10)
        y=powerlaw_all(x,(alpha,u_ref,z_ref))

        fun_kwargs = {'u_ref':u_ref,'z_ref':z_ref}
        y_fit, pfit, fitter = model_fit('predef: powerlaw_alpha', x, y, p0=(0.1), **fun_kwargs)
        #fitter.plot(); import matplotlib.pyplot as plt; plt.show()
        np.testing.assert_array_almost_equal(y,y_fit)
        np.testing.assert_almost_equal(alpha ,fitter.model['coeffs']['alpha'])

    def test_powerlaw_u_alpha(self):
        u_ref,z_ref,alpha=10,12,0.12
        x = np.linspace(0,1,10)
        y=powerlaw_all(x,(alpha,u_ref,z_ref,alpha))

        fun_kwargs = {'z_ref':z_ref}
        y_fit, pfit, fitter = model_fit('predef: powerlaw_u_alpha', x, y, **fun_kwargs)
        np.testing.assert_array_almost_equal(y,y_fit)
        np.testing.assert_almost_equal(alpha ,fitter.model['coeffs']['alpha'])
        np.testing.assert_almost_equal(u_ref ,fitter.model['coeffs']['u_ref'])

#     def test_powerlaw_all(self):
#         u_ref,z_ref,alpha=10,12,0.12
#         x = np.linspace(0,1,10)
#         y=powerlaw_all(x,(alpha,u_ref,z_ref,alpha))
# 
#         y_fit, pfit, fitter = model_fit('predef: powerlaw_all', x, y)
#         np.testing.assert_array_almost_equal(y,y_fit)
#         np.testing.assert_almost_equal(alpha ,fitter.model['coeffs']['alpha'])
# # NOTE: cannot test for u_ref or z

    def test_expdecay(self):
        A,k,B=0.5,1.2,10
        x=np.linspace(0,1,10)
        y=expdecay(x,(A,k,B))
        y_fit, pfit, fitter = model_fit('predef: expdecay', x, y)
        #fitter.plot(); import matplotlib.pyplot as plt; plt.show()
        np.testing.assert_array_almost_equal(y,y_fit)
        np.testing.assert_almost_equal(A,fitter.model['coeffs']['A'])
        np.testing.assert_almost_equal(k,fitter.model['coeffs']['k'])
        np.testing.assert_almost_equal(B,fitter.model['coeffs']['B'])

    def test_weibull(self):
        A, k = 10, 2.3,
        x=np.linspace(0.01,1,10)
        y=weibull_pdf(x,(A,k))
        y_fit, pfit, fitter = model_fit('predef: weibull_pdf', x, y)
        np.testing.assert_array_almost_equal(y,y_fit)
        np.testing.assert_almost_equal(A,fitter.model['coeffs']['A'],5)
        np.testing.assert_almost_equal(k,fitter.model['coeffs']['k'])

    # --------------------------------------------------------------------------------}
    # --- Testing of Predefined fitters 
    # --------------------------------------------------------------------------------{
    def test_secondorder_impulse(self):
        A, omega0, zeta, B, t0 = 2, 2, 0.01, 100, 10
        x=np.linspace(0,50,1000)
        y=secondorder_impulse(x,(A,omega0, zeta, B, t0)) #0* + np.random.normal(0, 0.1, len(x))
        #y_fit, pfit, fitter = model_fit('predef: secondorder_impulse', x, y)
        y_fit, pfit, fitter = model_fit('fitter: secondorder_impulse', x, y)
        #fitter.plot(); import matplotlib.pyplot as plt; plt.show()
        np.testing.assert_array_almost_equal(y,y_fit)
        np.testing.assert_almost_equal(omega0   ,fitter.model['coeffs']['omega'])
        np.testing.assert_almost_equal(zeta     ,fitter.model['coeffs']['zeta'])

    def test_secondorder_step(self):
        A, omega0, zeta, B, t0 = 2, 2, 0.8, 100, 10
        x=np.linspace(0,50,1000)
        y=secondorder_step(x,(A,omega0, zeta, B, t0)) #0* + np.random.normal(0, 0.1, len(x))
        #y_fit, pfit, fitter = model_fit('predef: secondorder_step', x, y)
        y_fit, pfit, fitter = model_fit('fitter: secondorder_step', x, y)
        #fitter.plot(); import matplotlib.pyplot as plt; plt.show()
        np.testing.assert_array_almost_equal(y,y_fit)
        np.testing.assert_almost_equal(omega0   ,fitter.model['coeffs']['omega'])
        np.testing.assert_almost_equal(zeta     ,fitter.model['coeffs']['zeta'])

    def test_polycont(self):
        k = 2.0
        x = np.linspace(0,1,10)
        y = k * x**3
        y_fit, pfit, fitter = model_fit('fitter: polynomial_continuous', x, y, order=3)
        #fitter.plot(); import matplotlib.pyplot as plt; plt.show()
        np.testing.assert_array_almost_equal(y,y_fit)
        np.testing.assert_almost_equal(k ,fitter.model['coeffs']['a'])
        np.testing.assert_almost_equal(0 ,fitter.model['coeffs']['b'])
        np.testing.assert_almost_equal(0 ,fitter.model['coeffs']['c'])
        np.testing.assert_almost_equal(0 ,fitter.model['coeffs']['d'])

    def test_polydisc(self):
        exponents=[0,3,5]
        a,b,c = 2.0, 3.0, 4.0
        x = np.linspace(0,1,10)
        y = a + b*x**3 + c*x**5
        y_fit, pfit, fitter = model_fit('fitter: polynomial_discrete', x, y, exponents=exponents)
        #fitter.plot(); import matplotlib.pyplot as plt; plt.show()
        np.testing.assert_array_almost_equal(y,y_fit)
        np.testing.assert_almost_equal(a ,fitter.model['coeffs']['a'])
        np.testing.assert_almost_equal(b ,fitter.model['coeffs']['b'])
        np.testing.assert_almost_equal(c ,fitter.model['coeffs']['c'])

    def test_sinusoid(self):
        A     = 101
        B     = -200.5
        omega = 10
        phi   = np.pi/3
        x=np.linspace(0,10,300)
        y=sinusoid(x,(A,omega,phi,B)) # + np.random.normal(0, 0.1, len(x))
        y_fit, pfit, fitter = model_fit('fitter: sinusoid', x, y, physical=False)
        #fitter.plot(); import matplotlib.pyplot as plt; plt.show()
        np.testing.assert_array_almost_equal(y,y_fit)
        np.testing.assert_almost_equal(A    ,fitter.model['coeffs']['A'],1)
        np.testing.assert_almost_equal(B    ,fitter.model['coeffs']['B'])
        np.testing.assert_almost_equal(omega,fitter.model['coeffs']['omega'])
        np.testing.assert_almost_equal(phi  ,fitter.model['coeffs']['phi'])

    def test_gentorque(self):
        pass # TODO
#         GBRatio= 27.5647     #; % Gearbox ratio (-)
#         SpdGenOn  = 14*GBRatio#
#         RtGnSp = 1207.61    # % Rated generator speed for simple variable-speed generator control (HSS side) (rpm) 
#         RtTq   = 1790.49    # % Rated generator torque/constant generator torque in Region 3 for simple variable-speed generator control (HSS side) (N-m) 
#         Rgn2K  = 0.0004128  # % Generator torque constant in Region 2 for simple variable-speed generator control (HSS side) (N-m/rpm^2) 
#         SlPc   = 6          # % Rated generator slip percentage in Region 2 1/2 for simple variable-speed generator control (%) 
# #         x=np.linspace(300,1500,100)
#         x=np.linspace(300,1000,100)
#         y=gentorque(x, (RtGnSp, RtTq  , Rgn2K , SlPc , SpdGenOn))
# 
#         bounds='RtGnSp=(1200,1300) , RtTq=(1500,1800), Rgn2K=(0.0,0.01) ,SlPc=(0,20) , SpdGenOn=(10,500)'
#         p0 = [1250, 1700,0.001, 10, 50]
#         y_fit, pfit, fitter = model_fit('fitter: gentorque', x, y)
# 
#         y_fit, pfit, fitter = model_fit('predef: gentorque', x, y, bounds=bounds, p0=p0)
# #         np.testing.assert_array_almost_equal(y,y_fit)
#         print(fitter)
#         import matplotlib.pyplot as plt
# 
#         fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
#         fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
#         ax.plot(x, y     ,'o', label='')
#         ax.plot(x, y_fit ,'-', label='')
#         ax.plot(x, fitter.model['fitted_function'](x) ,'.', label='')
#         ax.set_xlabel('')
#         ax.set_ylabel('')
#         ax.legend()
#         ax.tick_params(direction='in')
#         plt.show()


    # --------------------------------------------------------------------------------}
    # --- Testing of user models
    # --------------------------------------------------------------------------------{
    def test_evalpoly(self):
        exponents=[0,3,5]
        a,b,c = 2.0, 3.0, 4.0
        x = np.linspace(0,1,10)
        y = a + b*x**3 + c*x**5
        y_fit, pfit, fitter = model_fit('eval: {a} + {b}*x**3 + {c}*x**5', x, y)
        np.testing.assert_array_almost_equal(y,y_fit)
        np.testing.assert_almost_equal(a ,fitter.model['coeffs']['a'])
        np.testing.assert_almost_equal(b ,fitter.model['coeffs']['b'])
        np.testing.assert_almost_equal(c ,fitter.model['coeffs']['c'])

    def test_evalpowerlaw(self):
        u_ref,z_ref,alpha=10,12,0.12
        x = np.linspace(0,1,10)
        y=powerlaw_all(x,(alpha,u_ref,z_ref))
        y_fit, pfit, fitter = model_fit('eval: {u_ref}*(x/{z_ref})**{alpha}', x, y, p0=(8,9,0.1), bounds=(0.001,100))
        np.testing.assert_array_almost_equal(y,y_fit)

    def test_lowlevelpoly(self):
        x=np.linspace(0,1,10)
        y=x**2
        exponents=[0,1,2]
        y_fit, pfit, model = fit_polynomial_discrete(x, y, exponents)
        np.testing.assert_array_almost_equal(y,y_fit)
        np.testing.assert_almost_equal(1 , model['coeffs']['c'])
        np.testing.assert_almost_equal(0 , model['coeffs']['a'])

        y_fit, pfit, model = fit_polynomial_continuous(x, y, 3)
        np.testing.assert_array_almost_equal(y,y_fit)
        np.testing.assert_almost_equal(1 , model['coeffs']['b'])
        np.testing.assert_almost_equal(0 , model['coeffs']['a'])

    def test_lowlevelpowerlaw(self):
        u_ref,z_ref,alpha=10,12,0.12
        x = np.linspace(0,1,10)
        y=powerlaw_all(x,(alpha,u_ref,z_ref))

        y_fit, pfit, model = fit_powerlaw_u_alpha(x, y, z_ref=z_ref, p0=(9,0.1))
        np.testing.assert_array_almost_equal(y,y_fit)
        np.testing.assert_almost_equal(alpha , model['coeffs']['alpha'])
        np.testing.assert_almost_equal(u_ref , model['coeffs']['u_ref'])

#     def test_debug(self):
#         # --- Try Gaussian
#         x=np.linspace(0,1,10)
#         y=gaussian(x,(0.5,1.2))
#         y_fit, pfit, fitter = model_fit('predef: gaussian', x, y) #, p0=(0,1))
# #         fitter = ModelFitter('eval: {a}*(1.0/{b}+2/0)**{c}', x, y, p0=(8,9,0.1))
# #         fitter = ModelFitter('eval: {a}/x', x, y, p0=(8,9,0.1))
# 
#         # --- Plot 
#         y_fit=fitter.data['y_fit']
#         print(fitter)
# 
#         import matplotlib.pyplot as plt
#         fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
#         fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
#         ax.plot(x, y     ,'o', label='')
#         ax.plot(x, y_fit ,'-', label='')
#         ax.plot(x, fitter.model['fitted_function'](x) ,'.', label='')
#         ax.set_xlabel('')
#         ax.set_ylabel('')
#         ax.legend()
#         ax.tick_params(direction='in')
#         plt.show()

    def test_extract_var(self):
        var, _ = extract_variables('{a}*x + {b}')
        self.assertEqual(var,['a','b'])

        var, _ = extract_variables('{BB}*x + {a}*{BB}')
        self.assertEqual(var,['BB','a'])

        var, _ = extract_variables('{a}*x + {{b}}') #< TODO Won't work
        #self.assertEqual(var,['a','b'])

    def test_key_tuples(self):
        self.assertEqual(extract_key_tuples('a=(1,2)'),{'a':(1,2)})

        self.assertEqual(extract_key_tuples('a=(1, 2),b =(inf,0),c= ( -inf , 0.3e+10)'),{'a':(1,2),'b':(inf,0),'c':(-inf,0.3e+10)})

    def test_key_num(self):
        self.assertEqual(extract_key_num('a=2'),OrderedDict({'a':2}))
        self.assertEqual(extract_key_num('all=0.1,b =inf, c= -0.3e+10'),OrderedDict({'all':0.1,'b':inf,'c':-0.3e+10}))

    def test_key_misc(self):
        self.assertEqual(extract_key_miscnum('a=2'),{'a':2})

        #np.testing.assert_almost_equal(d['a'],(2,3))
        d=extract_key_miscnum('a=(2,3)')
        self.assertEqual(d['a'],(2,3))
        d=extract_key_miscnum('a=[2,3]')
        np.testing.assert_almost_equal(d['a'],[2,3])

        d=extract_key_miscnum('a=[2,3],b=3,c=(0,)')
        np.testing.assert_almost_equal(d['a'],[2,3])
        self.assertEqual(d['b'],3)
        self.assertEqual(d['c'],(0,))


if __name__ == '__main__':
    #TestFitting().test_sinusoid()
    #TestFitting().test_secondorder_step()
    unittest.main()




    ## --- Generate file for testing
    #import welib.system.firstorder as fo
    #import welib.system.secondorder as so
    #m      = 250.0              # system mass
    #k      = 40.0               # spring constant
    #omega0 = np.sqrt(k/m)
    #T0     = 2*np.pi/omega0
    #b      = 1/m                # NOTE: should be 1/m for a mechanical system

    #time = np.linspace(0*T0,10*T0,10001) # time span # Convolution needs enough time steps
    #t0  = 2*T0     # initial time where input starts (for simple inputs)
    #A   = 3         # amplitude for all simple inputs
    #T   = 2*T0   # for hat

    #zeta   = 0.01
    #x_d1 = so.impulse_response(time, omega0, zeta, b=b, t0=t0, A=A, both=True)
    #x_s1 = so.step_response   (time, omega0, zeta, b=b, t0=t0, A=A, both=True)
    #zeta   = 0.7
    #x_d2 = so.impulse_response(time, omega0, zeta, b=b, t0=t0, A=A, both=True)
    #x_s2 = so.step_response   (time, omega0, zeta, b=b, t0=t0, A=A, both=True)

    #M=np.column_stack((time,x_d1[0]+100,-x_d1[1],x_d2[0]-1000,x_d2[1], x_s1[0]+100, -x_s1[1], x_s2[0]-20,x_s2[1]))
    #np.savetxt('../TestFitSystem.csv',M,header='x,soi_z1,soi_z1_d,soi_z2,soi_z2_d,sos_z1,sos_z1_d,sos_z2,sos_z2_d',delimiter=',')

