import unittest
import sys
import os
import numpy as np
MyDir=os.path.dirname(__file__)
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))

from Polar import * 

# --------------------------------------------------------------------------------}
# ---  
# --------------------------------------------------------------------------------{
class TestPolarManip(unittest.TestCase):
    def assertNaN(self,x):
        self.assertTrue(np.isnan(x))

    def test_read(self):
        P=Polar.fromfile(os.path.join(MyDir,'../data/FFA-W3-241-Re12M.dat'))
        self.assertEqual(P.alpha[-1],180)
        self.assertEqual(P.cl[-1],0)

        P=Polar.fromfile(os.path.join(MyDir,'../data/Cylinder.dat'))
        self.assertEqual(P.cl.size,3)

    def test_alpha0(self):
        # --- Polar with one Cl value
        # non zero cl, alpha0 is nan
        self.assertNaN  (Polar([],[100],[0.1],[],[]).alpha0())
        # cl value is 0, alpha0 is arbitrarily 0
        self.assertEqual(Polar([],[100],[0.0],[],[]).alpha0(), 0)

        # --- Polar with one zero crossing
        P=Polar([],[-10,10],[-0.1,0.1],[],[])
        # Alpha0 is found as long as the window holds it
        self.assertEqual(P.alpha0(window=[-50,50]),0.0)
        self.assertEqual(P.alpha0(window=[-10,10]),0.0)
        self.assertEqual(P.alpha0(window=[ -2, 2]),0.0)
        # Error when window outside, no crossing found
        self.assertRaises(Exception,P.alpha0, window=[-100,-50])

        # --- Polar with many zero crossing
        P=Polar([],[-10,-5,0,5,10],[-0.1,0.1,-0.1,0.1,0.2],[],[])
        self.assertEqual(P.alpha0(window=[-10,-5]), -7.5)
        # Error when several zero crossing are found
        self.assertRaises(Exception,P.alpha0, window=[-10,10])

        # --- Polar with constant values 
        # non zero cl, alpha0 is nan
        self.assertNaN  (Polar([],[-10,10],[0.1,0.1],[],[]).alpha0())
        # cl is 0, alpha0 is arbitrarily 0
        self.assertEqual(Polar([],[-10,10],[0.0,0.0],[],[]).alpha0(), 0)
 

    def test_slope(self):
        # --- Polar with two points
        P=Polar([],[-1,1],[-1,1],[],[])
        sl,a0,WinLin,WinSearch=P.cl_linear_slope()
        np.testing.assert_almost_equal(sl,1.0)
        sl,a0,WinLin,WinSearch=P.cl_linear_slope(method='max')
        np.testing.assert_almost_equal(sl,1.0)
        # --- Polar three points lin
        P=Polar([],[-1,0,1],[-1,0,1],[],[])
        sl,a0,WinLin,WinSearch=P.cl_linear_slope()
        np.testing.assert_almost_equal(sl,1.0)
        # --- Polar three points cst
        P=Polar([],[-1,0,1],[1,1,1],[],[])
        sl,a0,WinLin,WinSearch=P.cl_linear_slope()
        np.testing.assert_almost_equal(sl,0.0)
        # --- Polar with sine shape
        P=Polar([],[-3,-2,-1,0,1,2,3],[-1,-2,-1,0,1,0,0],[],[])
        sl,a0,WinLin,WinSearch=P.cl_linear_slope()
        np.testing.assert_almost_equal(sl,1.0)
        sl,a0,WinLin,WinSearch=P.cl_linear_slope(method='max')
        np.testing.assert_almost_equal(sl,1.0)
        # --- Polar sine with plateaux 
        P=Polar([],[-3,-2,-1,0,1,2,3],[-1,-2,-2,-1,0,1,1],[],[])
        sl,a0,WinLin,WinSearch=P.cl_linear_slope()
        np.testing.assert_almost_equal(sl,1.0)
        sl,a0,WinLin,WinSearch=P.cl_linear_slope(method='max')
        np.testing.assert_almost_equal(sl,1.0)
        # --- Polar sine-line  -  Difficult to evaluate
        P=Polar([],[-3,-2,-1,0,1,2,3],[-1,-2.1,-2,-1.1,0,1.1,1.2],[],[])
        sl,a0,WinLin,WinSearch=P.cl_linear_slope()
        np.testing.assert_almost_equal(sl,1.0,decimal=1)
        sl,a0,WinLin,WinSearch=P.cl_linear_slope(method='max')
        np.testing.assert_almost_equal(sl,1.0,decimal=1)
        # --- Polar with a kink - Difficult
        P=Polar([],[-3,-2,-1,0,1,2,3],[-1,-2,-2,-2,0,1,1],[],[])
        sl,a0,WinLin,WinSearch=P.cl_linear_slope()
        np.testing.assert_almost_equal(sl,1.5,decimal=1)
        sl,a0,WinLin,WinSearch=P.cl_linear_slope(method='max')
        np.testing.assert_almost_equal(sl,2.0)
        # --- Polar step function
        P=Polar([],[-3,-2,-1,0,1,2,3],[-.5,-.5,-.5,-.5,.5,.5,.5],[],[])
        sl,a0,WinLin,WinSearch=P.cl_linear_slope()
        np.testing.assert_almost_equal(sl,1.0)
        sl,a0,WinLin,WinSearch=P.cl_linear_slope(method='max')
        np.testing.assert_almost_equal(sl,1.0)
        # --- Sine
        alpha = np.linspace(-50,50,100) 
        Cl = np.sin(alpha*np.pi/180.)*180/np.pi
        P=Polar([],alpha,Cl,[],[])
        sl,a0,WinLin,WinSearch=P.cl_linear_slope(window=[-10,10])
        np.testing.assert_almost_equal(sl,1.0, decimal=2)
        sl,a0,WinLin,WinSearch=P.cl_linear_slope(window=[-10,10],method='max')
        np.testing.assert_almost_equal(sl,1.0, decimal=2)
        # --- Real Polars
        P=Polar.fromfile(os.path.join(MyDir,'../data/FFA-W3-241-Re12M.dat'))
        sl,a0,WinLin,WinSearch=P.cl_linear_slope(radians=True)
        np.testing.assert_almost_equal(sl,7.034, decimal=3)
        sl,a0,WinLin,WinSearch=P.cl_linear_slope()
        np.testing.assert_almost_equal(sl,0.123, decimal=3)
        sl,a0,WinLin,WinSearch=P.cl_linear_slope(method='max')
        np.testing.assert_almost_equal(sl,0.13, decimal=3)
        # --- Real Polars
        P=Polar.fromfile(os.path.join(MyDir,'../data/63-235.csv'))
        sl,a0,WinLin,WinSearch=P.cl_linear_slope()
        np.testing.assert_almost_equal(sl,0.099, decimal=3)
        sl,a0,WinLin,WinSearch=P.cl_linear_slope(method='max')
        np.testing.assert_almost_equal(sl,0.113, decimal=3)
        # --- Cylinder
        P=Polar.fromfile(os.path.join(MyDir,'../data/Cylinder.csv'))
        sl,a0,WinLin,WinSearch=P.cl_linear_slope()
        self.assertEqual(sl,0.0)

        #print(sl,a0,WinLin,WinSearch)
        #import matplotlib.pyplot as plt
        #fig=plt.figure()
        #ax = fig.add_subplot(111)
        ## deg
        #ax.plot(P.alpha, P.cl)
        #ax.plot(WinLin,(np.array(WinLin)-a0)*sl,'o')
        #ax.plot(a0,0,'ko')
        #ax.plot(WinSearch[0],0,'ko')
        #ax.plot(WinSearch[1],0,'ko')
        #ax.plot(P.alpha,(P.alpha-a0)*sl,'--')
        #ax.set_xlim(np.array(WinLin)+[-20,+20])
        ## rad
        ##ax.plot(np.deg2rad(P.alpha), P.cl)
        ##ax.plot(WinLin,(np.array(WinLin)-a0)*sl,'o')
        ##ax.plot(np.deg2rad(P.alpha),(np.deg2rad(P.alpha)-a0)*sl,'--')
        #ax.set_ylim([-3,2.0])
        #plt.show()

    def test_fully_sep(self):
        # --- 63-235, for that polar
        # merging occurs at i=31 and i=120
        # at i=63 there is a singularity (f_st==1, cl_fs=cl/2)
        P=Polar.fromfile(os.path.join(MyDir,'../data/63-235.csv'))
        cl_fs,f_st=P.cl_fully_separated()
        # Below and above merging, fully sep polar is the same as original
        np.testing.assert_almost_equal(cl_fs[30] ,P.cl[30])
        np.testing.assert_almost_equal(cl_fs[121],P.cl[121])
        # Singularity at i=63
        np.testing.assert_almost_equal(cl_fs[63],P.cl[63]/2)
        np.testing.assert_almost_equal(f_st[63],1.0)
        # Close to singularity, should be not far from cl/2
        np.testing.assert_almost_equal(cl_fs[64],P.cl[64]/2*1.004,decimal=4)
        np.testing.assert_almost_equal(cl_fs[62],P.cl[62]/2*0.998,decimal=4)

        #import matplotlib.pyplot as plt
        #fig=plt.figure()
        #ax = fig.add_subplot(111)
        #ax.plot(P.alpha,P.cl,label='Cl')
        #ax.plot(P.alpha,cl_fs,'--',label='Cl_fs')
        #ax.plot(P.alpha,f_st,label='f_st')
        #plt.xlim([-50,50])
        #plt.ylim([-3,3])
        #plt.legend()
        #plt.show()
        #print(f_st)
# 
#         P=Polar.fromfile(os.path.join(MyDir,'../data/Cylinder.dat'))
#         sl,offset,amin,amax=P.cl_linear_slope()
# 
#         plt.show()

#         P=Polar.fromfile(os.path.join(MyDir,'../data/FFA-W3-241-Re12M.dat'))
#         P=Polar.fromfile(os.path.join(MyDir,'../data/63-235.csv'))
#         P=Polar.fromfile(os.path.join(MyDir,'../data/Cylinder.csv'))
#         P=Polar.fromfile(os.path.join(MyDir,'../data/AD_3-63-224_mod.csv'))
#         P=Polar.fromfile(os.path.join(MyDir,'../data/AD_4-63-218_mod.csv'))
#         P=Polar.fromfile(os.path.join(MyDir,'../data/AD_5_63-214_mod.csv'))


if __name__ == '__main__':
    unittest.main()
