import unittest
import numpy as np
from welib.tools.signal import zero_crossings 
from welib.tools.signal import convolution_integral 

# --------------------------------------------------------------------------------}
# ---  
# --------------------------------------------------------------------------------{
class TestSignal(unittest.TestCase):

    def test_zero_crossings(self):
        self.assertEqual(zero_crossings(np.array([0        ]))[0].size,0       )
        self.assertEqual(zero_crossings(np.array([0      ,0]))[0].size,0)
        self.assertEqual(zero_crossings(np.array([0      ,1]))[0].size,0)
        self.assertEqual(zero_crossings(np.array([-1,0,0, 1]))[0].size,0)
        self.assertEqual(zero_crossings(np.array([-1     ,1])), (0.5, 0, 1))
        self.assertEqual(zero_crossings(np.array([ 1,    -1])), (0.5, 0,-1))
        self.assertEqual(zero_crossings(np.array([-1,0,   1])), (1.0, 1,  1))
        xz,iz,sz=zero_crossings(np.array([-1,1,-1]))
        self.assertTrue(np.all(xz==[0.5,1.5]))
        self.assertTrue(np.all(iz==[0,1]))
        self.assertTrue(np.all(sz==[1,-1]))
        self.assertEqual(zero_crossings(np.array([ 1,-1]),direction='up'  )[0].size,0)
        self.assertEqual(zero_crossings(np.array([-1, 1]),direction='down')[0].size,0)
        

    def test_convolution(self):
        # Test that convolution returns expected values using simple functions
        time    = np.linspace(0, 2, 100)
        f       = time
        g       = time**2
        fog_ref = time**4/12

        fog = convolution_integral(time, f, g)

        #import matplotlib.pyplot as plt
        #fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
        #fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
        #ax.plot(time, fog_ref, 'k-' , label='ref')
        #ax.plot(time, fog    , '--' , label='num')
        #ax.set_xlabel('time')
        #ax.legend()
        #plt.show()

        np.testing.assert_almost_equal(fog, fog_ref, 3)
 
if __name__ == '__main__':
    unittest.main()
