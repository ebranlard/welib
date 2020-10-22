import unittest
import numpy as np
from welib.tools.functions import *

# --------------------------------------------------------------------------------}
# ---  
# --------------------------------------------------------------------------------{
class TestFunctions(unittest.TestCase):


    def test_heaviside(self):
        aae=np.testing.assert_almost_equal

        # infinite support
        aae(smooth_heaviside(-1000), 0)
        aae(smooth_heaviside(    0), 1./2)
        aae(smooth_heaviside( 1000), 1)

        # finite support
        aae(smooth_heaviside(-1   ,rng =(0,1)),0  )
        aae(smooth_heaviside(0    ,rng =(0,1)),0  )
        aae(smooth_heaviside(0.001,rng =(0,1)),0. )
        aae(smooth_heaviside(0.5  ,rng =(0,1)),1./2)
        aae(smooth_heaviside(0.999,rng =(0,1)),1  )
        aae(smooth_heaviside(1    ,rng =(0,1)),1  )

        # finite support, decreasing
        aae(smooth_heaviside(-1   ,rng =(1,0)),1  )
        aae(smooth_heaviside(0    ,rng =(1,0)),1  )
        aae(smooth_heaviside(0.001,rng =(1,0)),1  )
        aae(smooth_heaviside(0.5  ,rng =(1,0)),1./2)
        aae(smooth_heaviside(0.999,rng =(1,0)),0  )
        aae(smooth_heaviside(1    ,rng =(1,0)),0  )

#         x=np.linspace(-2,2,100)
#         H =smooth_heaviside(x   ,rng =(1,0))
#         import matplotlib.pyplot as plt
#         fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
#         ax.plot(x,H    , label='')
#         ax.set_xlabel('')
#         ax.set_ylabel('')
#         ax.legend()
#         plt.show()


    def test_delta(self):
        old_settings = np.seterr()
        np.seterr(under='ignore')

        aae=np.testing.assert_almost_equal
        # infinite support
        for m in ['gaussian','frac']:
            x     = np.linspace(-10,10,1000)
            delta = smooth_delta(x, e=0.1, method=m)
            aae(delta[0] , 0, 2)
            aae(delta[-1], 0, 2)
            aae(np.trapz(delta,x), 1, 2)
        # finite support
        for m in ['gaussian','frac']:
            x     = np.linspace(-1,3,1000)
            delta = smooth_delta(x, rng=(0,2), e=0.1, method=m)
            aae(delta[0] , 0, 2)
            aae(delta[-1], 0, 2)
            aae(np.trapz(delta,x), 1, 2)

        # Heaviside derivative
        x = np.linspace(-1,3,100)
        R=(0,2.5)
        k=2
        H     = smooth_heaviside(x, k = k, rng = R, method = 'exp')
        D     = smooth_delta    (x, e = k, rng = R, method = 'exp-heaviside')
        try:
            D_num = np.diff(H, prepend=0)/(x[1]-x[0])
        except TypeError as e:
            D_num = np.concatenate(([0],np.diff(H)))/(x[1]-x[0])

        aae(D,D_num,1)

        np.seterr(**old_settings)
#         import matplotlib.pyplot as plt
#         fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
#         ax.plot(x,D     , label='')
#         ax.plot(x,D_num , label='')
#         ax.plot(x,D-D_num , label='')
#         ax.set_xlabel('')
#         ax.set_ylabel('')
#         ax.legend()
#         plt.show()
        
 
if __name__ == '__main__':
    unittest.main()
