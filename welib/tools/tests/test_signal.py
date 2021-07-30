import unittest
import numpy as np
from welib.tools.signal import zero_crossings 
from welib.tools.signal import convolution_integral 
from welib.tools.signal import find_time_offset 

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


    def test_phase(self, test=True, noise=False):
        # Test
        omega = 1.3
        T     = 2*np.pi/omega
        time  = np.linspace(0, 2*T, 601)
        dt= time[1]-time[0]

        if test:
            vphi = np.arange(0,360,10)
        else:
            vphi=[-60]

        for phi_ref_deg in vphi:

            phi_ref      = phi_ref_deg*np.pi/180
            t_offset_ref = phi_ref/omega
            f     = 1.5 * np.sin(omega*time        )
            g     = 0.5 * np.sin(omega*time+phi_ref) 
            if noise:
                g+=np.random.normal(0, 0.1, len(time))

            t_offset, lag, xcorr = find_time_offset(time, f, g, outputAll=True)
            phi = np.mod(t_offset*omega,2*np.pi)
            #print('t_offset', np.around(t_offset,4), np.around(t_offset_ref,4))
            #print('phi    {:7.3f} {:7.3f} {:7.3f}'.format(phi*180/np.pi, phi_ref*180/np.pi, (phi-phi_ref)*180/np.pi))

            #import matplotlib.pyplot as plt
            #imax=np.argmax(f[time<T])
            #g2    = 0.5 * np.sin(omega*(time+t_offset))
            #g3    = 0.5 * np.sin(omega*time+phi)
            #fig,axes = plt.subplots(2, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
            #fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
            #ax = axes[0]
            #ax.plot(time, f, '-' , label='Signal 1')
            #ax.plot(time, g, '-' , label='Signal 2')
            #ax.plot(time, g2, ':', label='Signal 1 with phase offset')
            ##ax.plot(time, g3, '-.', label='g3')
            #ax.plot(time[imax]+np.array([-t_offset,-t_offset]), [-1,1], 'k:')
            #ax.plot(time[imax]+np.array([0,0]              ), [-1,1], 'k:')
            #ax.legend()
            #ax = axes[1]
            #ax.plot(lag   , xcorr, '-' , label='f1')
            #ax.plot([t_offset,t_offset], [-np.max(xcorr),np.max(xcorr)], 'k:')
            #ax.set_xlabel('time')
            #ax.grid()
            #plt.show()

            dphi = phi-phi_ref
            if dphi>np.pi:
                dphi=dphi-2*np.pi
            test = abs(dphi)*180/np.pi<6.1
            if not test:
                print('phi    {:7.3f} {:7.3f} {:7.3f}'.format(phi*180/np.pi, phi_ref*180/np.pi, dphi*180/np.pi))
            self.assertTrue(test)

 
if __name__ == '__main__':
    unittest.main()
