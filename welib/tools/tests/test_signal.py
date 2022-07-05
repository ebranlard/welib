import unittest
import numpy as np
from welib.tools.signal_analysis import zero_crossings 
from welib.tools.signal_analysis import convolution_integral 
from welib.tools.signal_analysis import find_time_offset 
from welib.tools.signal_analysis import intervals 
from welib.tools.signal_analysis import peaks 
from welib.tools.signal_analysis import sine_approx
from welib.tools.signal_analysis import multiInterp
from welib.tools.signal_analysis import interpArray
from welib.tools.clean_exceptions import *

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


    def test_phase(self, test=True, noise=True):
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

            if not test:
                import matplotlib.pyplot as plt
                imax=np.argmax(f[time<T])
                g2    = 0.5 * np.sin(omega*(time+t_offset))
                g3    = 0.5 * np.sin(omega*time+phi)
                fig,axes = plt.subplots(2, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
                fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
                ax = axes[0]
                ax.plot(time, f, '-' , label='Signal 1')
                ax.plot(time, g, '-' , label='Signal 2')
                ax.plot(time, g2, ':', label='Signal 1 with phase offset')
                #ax.plot(time, g3, '-.', label='g3')
                ax.plot(time[imax]+np.array([-t_offset,-t_offset]), [-1,1], 'k:')
                ax.plot(time[imax]+np.array([0,0]              ), [-1,1], 'k:')
                ax.legend()
                ax = axes[1]
                ax.plot(lag   , xcorr, '-' , label='f1')
                ax.plot([t_offset,t_offset], [-np.max(xcorr),np.max(xcorr)], 'k:')
                ax.set_xlabel('time')
                ax.grid()
                plt.show()

            dphi = phi-phi_ref
            if dphi>np.pi:
                dphi=dphi-2*np.pi
            test = abs(dphi)*180/np.pi<7.0
            if not test:
                print('phi    {:7.3f} {:7.3f} {:7.3f}'.format(phi*180/np.pi, phi_ref*180/np.pi, dphi*180/np.pi))
            self.assertTrue(test)

    def test_intervals(self):
        import matplotlib.pyplot as plt

        IStart,IEnd,Lengths=intervals([False])
        self.assertEqual(len(IStart) , 0)

        IStart,IEnd,Lengths=intervals([True])
        self.assertEqual(IStart , np.array([0]))
        self.assertEqual(IEnd   , np.array([0]))
        self.assertEqual(Lengths, np.array([1]))

        IStart,IEnd,Lengths=intervals([False, True])
        self.assertEqual(IStart , np.array([1]))
        self.assertEqual(IEnd   , np.array([1]))
        self.assertEqual(Lengths, np.array([1]))

        IStart,IEnd,Lengths=intervals([True, False])
        self.assertEqual(IStart , np.array([0]))
        self.assertEqual(IEnd   , np.array([0]))
        self.assertEqual(Lengths, np.array([1]))

        IStart,IEnd,Lengths=intervals([True, True])
        self.assertEqual(IStart , np.array([0]))
        self.assertEqual(IEnd   , np.array([1]))
        self.assertEqual(Lengths, np.array([2]))

        IStart,IEnd,Lengths=intervals([False, True, False])
        np.testing.assert_equal(IStart , np.array([1]))
        np.testing.assert_equal(IEnd   , np.array([1]))
        np.testing.assert_equal(Lengths, np.array([1]))

        IStart,IEnd,Lengths=intervals([False, True, True])
        np.testing.assert_equal(IStart , np.array([1]))
        np.testing.assert_equal(IEnd   , np.array([2]))
        np.testing.assert_equal(Lengths, np.array([2]))

        IStart,IEnd,Lengths=intervals([True, True, False])
        np.testing.assert_equal(IStart , np.array([0]))
        np.testing.assert_equal(IEnd   , np.array([1]))
        np.testing.assert_equal(Lengths, np.array([2]))

        IStart,IEnd,Lengths=intervals([True, False, True])
        np.testing.assert_equal(IStart , np.array([0,2]))
        np.testing.assert_equal(IEnd   , np.array([0,2]))
        np.testing.assert_equal(Lengths, np.array([1,1]))

        IStart,IEnd,Lengths=intervals([True, True, False, True, True, True])
        np.testing.assert_equal(IStart , np.array([0,3]))
        np.testing.assert_equal(IEnd   , np.array([1,5]))
        np.testing.assert_equal(Lengths, np.array([2,3]))

        IStart,IEnd,Lengths=intervals([False, True, True, False, True, True, True, False])
        np.testing.assert_equal(IStart , np.array([1,4]))
        np.testing.assert_equal(IEnd   , np.array([2,6]))
        np.testing.assert_equal(Lengths, np.array([2,3]))

    def test_peaks(self, plot=False):
        # Test that we can retrieve the number of maxima of a sinusoid even with noise and different sampling
        np.random.seed(3)
        nLength = [20,30,40,50,60,70,100,1000,10000,100000]
        for n in nLength:
            t = np.linspace(0,300,n)
            f = 3*np.sin(0.1*t) + np.random.normal(0, 0.6, len(t))
            threshold=1
            min_length=1
            I, IStart, IEnd = peaks(f, threshold=threshold, method='intervals', min_length=min_length, returnIntervals=True)

            if plot:
                import matplotlib.pyplot as plt
                fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
                fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
                ax.plot(t, f)
                ax.plot(t[I], f[I], 'o')
                for iS,iE in zip(IStart,IEnd):
                    ax.plot(t[iS:iE+1], f[iS:iE+1], 'k.')
                ax.set_xlabel('Time')
                ax.set_ylabel('')
                plt.show()

            #print('I',I, len(I))
            self.assertEqual(len(I),5)

    def test_sine_approx(self, plot=False):
        # 
        np.random.seed(3)
        n     = 300
        t     = np.linspace(0,300,n)
        omega = 0.1
        A     = 3
        phi   = np.pi/3
        f = A*np.sin(omega*t+phi) + np.random.normal(0, 0.6, len(t))
        f2, omega2, A2, phi2 = sine_approx(t, f, method='least_square')
        np.testing.assert_almost_equal(omega, omega2,2)
        np.testing.assert_almost_equal(A, A2  ,1)
        np.testing.assert_almost_equal(phi, phi2  ,1)

        if plot:
            import matplotlib.pyplot as plt
            fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
            fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
            ax.plot(t, f ,label='input')
            ax.plot(t, f2,label='fit')
            ax.set_xlabel('Time')
            ax.set_ylabel('')
            plt.show()


    def test_interp(self, plot=False):
        x = np.linspace(0,1,10)
        y1= x**2
        y2= x**3
        Y = np.stack((y1,y2))

        # --- Check that we retrieve proper value on nodes
        x_new = x
        Y_new = multiInterp(x_new, x, Y)
        np.testing.assert_almost_equal(Y_new, Y)

        # using interpArray
        Y_new2 = np.zeros(Y_new.shape)
        for i,x0 in enumerate(x_new):
            Y_new2[:,i] = interpArray(x0, x, Y)
        np.testing.assert_almost_equal(Y_new2, Y)

        # --- Check that we retrieve proper value on misc nodes
        x_new  = np.linspace(-0.8,1.5,20)
        Y_new  = multiInterp(x_new, x, Y)
        y1_new = np.interp(x_new, x, Y[0,:])
        y2_new = np.interp(x_new, x, Y[1,:])
        Y_ref  = np.stack((y1_new, y2_new))
        np.testing.assert_almost_equal(Y_new, Y_ref)

        # using interpArray
        Y_new2 = np.zeros(Y_new.shape)
        for i,x0 in enumerate(x_new):
            Y_new2[:,i] = interpArray(x0, x, Y)
        np.testing.assert_almost_equal(Y_new2, Y_ref)




if __name__ == '__main__':
    #TestSignal().test_sine_approx()
    #TestSignal().test_interp()
    unittest.main()
