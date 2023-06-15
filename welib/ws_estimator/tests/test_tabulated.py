import unittest
import os
import numpy as np    
import matplotlib.pyplot as plt
from welib.ws_estimator.tabulated import *

scriptDir = os.path.dirname(__file__)

R           = 63
rho         = 1.225
aeroMapFile = os.path.join(scriptDir, '../../../data/NREL5MW/NREL5MW_CPCTCQ.txt')
operFile    = os.path.join(scriptDir, '../../../data/NREL5MW/NREL5MW_Oper.csv')
WSENoOper   = TabulatedWSEstimator(R = R, rho = rho, aeroMapFile = aeroMapFile)
WSE         = TabulatedWSEstimator(R = R, rho = rho, aeroMapFile = aeroMapFile, operFile = operFile)


def simpleTest(wse, pitch, omega, Qa, WS_guess, WS_ref, deltaWSMax, plot=False, method='crossing'):

    WS_est1, info1 = wse.estimate(Qa, pitch, omega, WS_guess, relaxation=0, debug=True, method='min', deltaWSMax=deltaWSMax)
    Qeval1 = wse.Torque(WS_est1, pitch, omega)

    WS_est2, info2 = wse.estimate(Qa, pitch, omega, WS_guess, relaxation=0, debug=True, method=method, deltaWSMax=deltaWSMax)
    Qeval2 = wse.Torque(WS_est2, pitch, omega)


    # --- Plot
    if plot:
#         info=info1
        info=info2
        info['WS_ref'] = WS_ref
        wse.debugPlot(info)
        plt.show()

    np.testing.assert_almost_equal(WS_est1, WS_ref, 1)
    np.testing.assert_almost_equal(WS_est2, WS_ref, 2)
    np.testing.assert_almost_equal(Qeval1/Qa, Qa/Qa, 2)
    np.testing.assert_almost_equal(Qeval2/Qa, Qa/Qa, 2)





class Test(unittest.TestCase):

    def test_noOper_1Inter(self):
        # - No Operating condition file
        # - One intersection is present
        plot=False
        simpleTest(WSENoOper, pitch=0, omega=12.1*np.pi/30, Qa=0.6e7, WS_guess=5, WS_ref=13.0019, deltaWSMax=15, plot=plot)


    def test_noOper_2Inter(self):
        # - No Operating condition file
        # - One intersection is present
        # NOTE: Without Oper, there is no way to decide on the region!

        # --- Lower region
        plot     = False
        WS_guess = 26.5   # middle is: 26.9
        WS_ref   = 22.931 # [22.931 30.964]
        simpleTest(WSENoOper, pitch=0, omega=12.1*np.pi/30, Qa=1.15e7, WS_guess=WS_guess, WS_ref=WS_ref, deltaWSMax=5, plot=plot)

        # --- Upper region
        plot     = False
        WS_guess = 27.1
        WS_ref   = 30.964 # [22.907 30.964]
        simpleTest(WSENoOper, pitch=0, omega=12.1*np.pi/30, Qa=1.15e7, WS_guess=WS_guess, WS_ref=WS_ref, deltaWSMax=5, plot=plot)

    def test_Oper_2Inter(self):
        # - No Operating condition file
        # - One intersection is present
        # NOTE: Without Oper, there is no way to decide on the region!

        wse=WSE

        plot     = False
        pitch    = 0
        omega    = 7.1*np.pi/30
        Qa       = 3.8e6
        WS_guess = 15
        WS_ref   = 9.3  

        WS_est2, info2 = wse.estimate(Qa, pitch, omega, WS_guess, relaxation=0, debug=True, method='crossing-oper')
        Qeval2 = wse.Torque(WS_est2, pitch, omega)
        Qeval2 = wse.Torque(WS_est2, pitch, omega)

        if plot:
            info2['WS_ref'] = WS_ref
            wse.debugPlot(info2)
            plt.show()

        self.assertTrue(abs(WS_est2- WS_ref)<3.5)
        np.testing.assert_almost_equal(Qeval2/Qa, Qa/Qa, 1)


    def test_noOper_1Inter_Far(self):
        # One itersection 
        # But the rotor is in a state sligtly away
        # We use the guess
        wse=WSENoOper 
        pitch    = 0.000
        omega    = 0.0885
        Qa       = 45028.61
        WS_guess = 3.369
        WS_ref   = 4.184
        plot=False
        WS_est2, info2 = wse.estimate(Qa, pitch, omega, WS_guess, relaxation=0, debug=True, method='crossing')
        Qeval2 = wse.Torque(WS_est2, pitch, omega)

        if plot:
            info2['WS_ref'] = WS_ref
            wse.debugPlot(info2)
            plt.show()

        self.assertTrue(abs(WS_est2- WS_ref)<0.8)
        np.testing.assert_almost_equal(Qeval2/Qa, Qa/Qa, 1)
        #simpleTest(WSENoOper, pitch, omega, Qa, WS_guess, WS_ref, deltaWSMax=5, plot=plot)

    def test_Oper_1Inter_Far(self):
        # One itersection 
        # But the rotor is in a state sligtly away
        # We use the guess
        wse=WSE
        pitch    = 0.000
        omega    = 0.0885
        Qa       = 45028.61
        WS_guess = 3.369
        WS_ref   = 4.184
        deltaWSMax     =1
        relaxation     =0.0

#         Qa             =892733.8596026392
#         pitch          =0.0
#         omega          =0.47585762276984406
#         WS0            =14.69469162000028
        method         ='crossing-oper'
#         iNear          =2
#         WS_est1        =14.710287501616827
#         WS_est2        =14.710287501616827
#         WS_ref         =8.696043945426664

        plot=False
        WS_est2, info2 = wse.estimate(Qa, pitch, omega, WS_guess, relaxation=relaxation, debug=True, method=method)
        Qeval2 = wse.Torque(WS_est2, pitch, omega)

        if plot:
            info2['WS_ref'] = WS_ref
            wse.debugPlot(info2)
            plt.show()

        self.assertTrue(abs(WS_est2- WS_ref)<1.3)
        np.testing.assert_almost_equal(Qeval2/Qa, Qa/Qa, 1)
        #simpleTest(WSENoOper, pitch, omega, Qa, WS_guess, WS_ref, deltaWSMax=5, plot=plot)


    def test_TimeMethods(self):
        # Compare the computational time of the difference methods
        from timeit import timeit
        wse = WSENoOper
        pitch    = 0
        omega    = 12.1*np.pi/30
        Qa       = 0.6e7
        WS_guess = 10
        WS0, Q0 = wse.TorqueAt(Pitch=0, Omega=omega)
        def meth1():
            WS_est, info = wse.estimate(Qa, pitch, omega, WS_guess, relaxation=0, debug=True, method='min', deltaWSMax=15)
        def meth2():
            WS_est, info = wse.estimate(Qa, pitch, omega, WS_guess, relaxation=0, debug=True, method='crossing', deltaWSMax=15)

        from timeit import timeit
        print('Time 1: ',timeit(meth1, number=1000))
        print('Time 2: ',timeit(meth2, number=1000))


if __name__ == '__main__':
#     Test().test_noOper_1Inter()
#     Test().test_noOper_2Inter()
#     Test().test_noOper_1Inter_Far()
#     Test().test_Oper_1Inter()
#     Test().test_Oper_2Inter()
#     Test().test_Oper_1Inter_Far()
#     Test().test_noOper_3Inter()
#     Test().test_TimeMethods()
    unittest.main()
