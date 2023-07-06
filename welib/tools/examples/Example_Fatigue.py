""" 
Compute damage equivalent load of a sinusoidal signals for varying ampitudes and frequencies
Compares with analytical result.

"""

import numpy as np
import matplotlib.pyplot as plt
from welib.tools.fatigue import eq_load, equivalent_load

try:
    import fatpack
    hasFatpack=True
except:
    hasFatpack=False

def main(Astudy=False, fstudy=True):
    if Astudy:
        # --------------------------------------------------------------------------------
        # --- Dependency on amplitude 
        # --------------------------------------------------------------------------------
        m = 3 # Wohler slope
        nPerT = 100
        nT = 10    # Number of periods
        T = 1
        omega = 1/T*2*np.pi
        vA  = np.array([1.,2.,3.,4.,5.])
        Leq     = np.zeros_like(vA)
        Leq_ref = np.zeros_like(vA)
        time = np.linspace(0, nT*T, nPerT*nT+1)
        T_all=time[-1]
        for ia, A in enumerate(vA):
           signal = A * np.sin(omega*time)  + 5
           Leq[ia] = equivalent_load(time, signal, m=m, Teq=1, nBins=1, method='rainflow_windap')
        Leq_ref = 2*((nT*vA**m)/T_all)**(1/m)

        fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
        fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
        ax.plot(vA, Leq_ref, 'k'  ,label='Theory')
        ax.plot(vA, Leq    ,  'o' ,label='Windap')
        ax.set_xlabel('Amplitude, A [N]')
        ax.set_ylabel('Leq [N]')
        ax.legend()

    if fstudy:
        # --------------------------------------------------------------------------------
        # --- Dependency on frquency
        # --------------------------------------------------------------------------------
        m     = 2   # Wohler slope
        A     = 3   # Amplitude
        nT    = 100 # Number of periods
        nPerT = 100 # Number of points per period
        Teq   = 1  # Equivalent period [s]
        nBins = 1  # Number of bins

        fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
        fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)

        # vT=np.array([1, 2, 5, 10, 50, 100])
        vf =np.linspace(0.1,10,21)
        vT  = 1/vf
        T_max=np.max(vT*nT)
        vomega =vf*2*np.pi
        Leq1    = np.zeros_like(vomega)
        Leq2    = np.zeros_like(vomega)
        Leq_ref = np.zeros_like(vomega)
        for im,m in enumerate([1, 2, 10]):
            for it, (T,omega) in enumerate(zip(vT,vomega)):
                # --- Option 1 - Same number of periods
                time = np.linspace(0, nT*T, nPerT*nT+1)

                # --- Option 2 - Same tMax
                #nT=int(T_max/T)
                #time = np.linspace(0, T_max, nPerT*nT+1)
                signal = A * np.sin(omega*time) # Mean does not matter 
                T_all=time[-1]
                Leq1[it] = equivalent_load(time, signal, m=m, Teq=Teq, nBins=nBins, method='rainflow_windap')
                if hasFatpack:
                    Leq2[it] = equivalent_load(time, signal, m=m, Teq=Teq, nBins=nBins, method='fatpack')

            Leq_ref = 2*A*(vf*Teq)**(1/m)
            ax.plot(vf, Leq_ref/A,  'k-' , label ='Theory' if im==0 else None  )
            ax.plot(vf, Leq1   /A,   'o' , label ='Windap m={}'.format(m))
            if hasFatpack:
                ax.plot(vf, Leq2/A,  'k.' , label ='Fatpack' if im==0 else None)
        ax.set_xlabel('Frequency [Hz]')
        ax.set_ylabel(r'Equivalent load, $L_{eq}/A$ [-]')
        # ax.set_ylim([0,7])
        ax.tick_params(direction='in')
        ax.legend()
        ax.set_title('Tools - Fatigue - Sinusoids')

    return vf, Leq1, Leq_ref

vf, Leq, Leq_ref =  main(Astudy=False, fstudy=True)


if __name__ == '__main__':
    plt.show()
if __name__=="__test__":
    np.testing.assert_array_almost_equal(Leq, Leq_ref,1)
if __name__=="__export__":
    from welib.tools.repo import export_figs_callback
    export_figs_callback(__file__)



