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
#         P=Polar([],[-1,1],[-1,1],[],[])
#         P.cl_linear_slope()
#         print('>>>>>>')
#         # --- Polar with sine shape
#         P=Polar([],[-3,-2,-1,0,1,2,3],[-1,-2,-1,0,1,0,0],[],[])
#         P.cl_linear_slope(window=[-20,20])
#         print('>>>>>>')
#         P=Polar([],[-3,-2,-1,0,1,2,3],[-1,-2,-2,-1,0,1,1],[],[])
#         P.cl_linear_slope(window=[-20,20])
#         P=Polar([],[-3,-2,-1,0,1,2,3],[-1,-2.1,-2,-1.1,0,1.1,1.2],[],[])
#         P.cl_linear_slope(window=[-20,20])
#         print('>>>>>>')
#         P=Polar([],[-3,-2,-1,0,1,2,3],[-1,-2,-2,-2,0,1,1],[],[])
#         P.cl_linear_slope(window=[-20,20])
#         P=Polar([],[-3,-2,-1,0,1,2,3],[-.5,-.5,-.5,-.5,.5,.5,.5],[],[])
#         P.cl_linear_slope(window=[-20,20])
#         P=Polar(filename=os.path.join(MyDir,'../data/FFA-W3-241-Re12M.dat'))

#         alpha = np.linspace(-50,50,100) 
#         Cl = np.sin(alpha*np.pi/180.)*180/np.pi
#         P=Polar([],alpha,Cl,[],[])
#         sl,a0,WinLin,WinSearch=P.cl_linear_slope(window=[-30,30])
# #         sl,a0,WinLin,WinSearch=P.cl_linear_slope()


        P=Polar.fromfile(os.path.join(MyDir,'../data/FFA-W3-241-Re12M.dat'))
#         P=Polar.fromfile(os.path.join(MyDir,'../data/AD_1_Cylinder.csv'))
#         P=Polar.fromfile(os.path.join(MyDir,'../data/AD_2_63-235_mod.csv'))
#         P=Polar.fromfile(os.path.join(MyDir,'../data/AD_3-63-224_mod.csv'))
#         P=Polar.fromfile(os.path.join(MyDir,'../data/AD_4-63-218_mod.csv'))
#         P=Polar.fromfile(os.path.join(MyDir,'../data/AD_5_63-214_mod.csv'))
#         P.cl[P.cl>=0]=np.nan
        sl,a0,WinLin,WinSearch=P.cl_linear_slope()
#         sl,a0,WinLin,WinSearch=P.cl_linear_slope(window=[-5,5])
#         sl,a0,WinLin,WinSearch=P.cl_linear_slope(window=[-5,5])
#         print(sl,a0,WinLin,WinSearch)

#         P=Polar.fromfile(os.path.join(MyDir,'../data/Cylinder.dat'))
#         sl,offset,amin,amax=P.cl_linear_slope()

#         import matplotlib.pyplot as plt
#         fig=plt.figure()
#         ax = fig.add_subplot(111)
#         ax.plot(P.alpha, P.cl)
#         ax.plot(WinLin,(np.array(WinLin)-a0)*sl,'o')
#         ax.plot(a0,0,'ko')
#         ax.plot(WinSearch[0],0,'ko')
#         ax.plot(WinSearch[1],0,'ko')
#         ax.plot(P.alpha,(P.alpha-a0)*sl,'--')
#         ax.set_xlim(np.array(WinLin)+[-20,+20])
#         ax.set_ylim([-3,2.0])
#         plt.show()



if __name__ == '__main__':
    unittest.main()
