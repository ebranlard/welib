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

def main(Astudy=False, fstudy=True, fstudy2=False, sumSine=False):
    vf=None; Leq1=None; Leq_ref=None

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
           Leq[ia] = equivalent_load(time, signal, m=m, Teq=1, bins=1, method='rainflow_windap')
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
        # --- Dependency on frequency
        # --------------------------------------------------------------------------------
        m     = 2   # Wohler slope
        A     = 3   # Amplitude
        nT    = 100 # Number of periods
        nPerT = 100 # Number of points per period
        Teq   = 1  # Equivalent period [s]
        nBins = 1  # Number of bins

        fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
        fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)

        vf =np.linspace(0.1,10,21)
        vm = np.array([1, 2, 10])
        #vm = np.array([1])
        #vf = np.array([0.1])
        vT  = 1/vf
        T_max=np.max(vT*nT)
        vomega =vf*2*np.pi
        Leq1    = np.zeros_like(vomega)
        Leq2    = np.zeros_like(vomega)
        Leq_ref = np.zeros_like(vomega)
        for im,m in enumerate(vm):
            for it, (T,omega) in enumerate(zip(vT,vomega)):
                # --- Option 1 - Same number of periods
                time = np.linspace(0, nT*T, nPerT*nT+1)

                # --- Option 2 - Same tMax
                #nT=int(T_max/T)
                #time = np.linspace(0, T_max, nPerT*nT+1)
                signal = A * np.sin(omega*time) # Mean does not matter 
                T_all=time[-1]
                Leq1[it] = equivalent_load(time, signal, m=m, Teq=Teq, bins=nBins, method='rainflow_windap')
                if hasFatpack:
                    Leq2[it] = equivalent_load(time, signal, m=m, Teq=Teq, bins=nBins, method='fatpack')

            Leq_ref = 2*A*(vf*Teq)**(1/m)
            ax.plot(vf, Leq_ref/A,  'kd' , label ='Theory' if im==0 else None  )
            ax.plot(vf, Leq1   /A,   'o' , label ='Windap m={}'.format(m))
            if hasFatpack:
                ax.plot(vf, Leq2/A,  'k.' , label ='Fatpack' if im==0 else None)
        ax.set_xlabel('Frequency [Hz]')
        ax.set_ylabel(r'Equivalent load, $L_{eq}/A$ [-]')
        # ax.set_ylim([0,7])
        ax.tick_params(direction='in')
        ax.legend()
        ax.set_title('Tools - Fatigue - Sinusoids')

    if fstudy2:
        # --------------------------------------------------------------------------------
        # --- Dependency on frequency with a more "random" signal
        # --------------------------------------------------------------------------------
        m     = 2   # Wohler slope
        A     = 3   # Amplitude
        nT    = 100 # Number of periods
        nPerT = 100 # Number of points per period
        Teq   = 1  # Equivalent period [s]
        nBins = 100  # Number of bins
        nSignals=100

        # vT=np.array([1, 2, 5, 10, 50, 100])
        vf =np.linspace(0.1,10,21)
        vT  = 1/vf
        T_max=np.max(vT*nT)
        vomega =vf*2*np.pi
        Leq1    = np.zeros_like(vomega)
        Leq2    = np.zeros_like(vomega)
        Leq_ref = np.zeros_like(vomega)

        fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
        fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
        for im,m in enumerate([4]):
            for it, (T,omega) in enumerate(zip(vT,vomega)):
                # --- Option 1 - Same number of periods
                time = np.linspace(0, nT*T, nPerT*nT+1)
                # --- Option 2 - Same tMax
                #nT=int(T_max/T)
                #time = np.linspace(0, T_max, nPerT*nT+1)

                # --- Creating a random signal with frequencies close to target
                from numpy.random import normal, rand
                signal_base = A * np.sin(omega*time) # Mean does not matter 
                signal = 0
                omegas = np.clip(0, np.inf, normal(omega, 1.5*np.min(vomega), nSignals))
                for iSig in range(nSignals):
                    Ai     = normal(A, 0.5*A)
                    phii = rand()*np.pi
                    signal += Ai * np.sin(omegas[iSig]*time+ phii) # Mean does not matter 
                A_ = np.max(signal) # TODO figure out a scaling that makes sense
                A0 = A
                signal *= (A0/A_)
                #fig1,ax1 = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
                #fig1.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
                #ax1.plot(time, signal,     )
                #ax1.plot(time, signal_base, 'k:')
                #ax1.set_xlabel('time')

                T_all=time[-1]
                Leq1[it] = equivalent_load(time, signal, m=m, Teq=Teq, bins=nBins, method='rainflow_windap')
                if hasFatpack:
                    Leq2[it] = equivalent_load(time, signal, m=m, Teq=Teq, bins=nBins, method='fatpack')

            Leq_ref = 2*A*(vf*Teq)**(1/m)
            ax.plot(vf, Leq_ref/A,  'k-' , label ='Theory' if im==0 else None  )
            ax.plot(vf, Leq1   /A,   'o' , label ='Windap m={}'.format(m))
            if hasFatpack:
                ax.plot(vf, Leq2/A,  'k.' , label ='Fatpack' if im==0 else None)

        ax.set_xlabel('Frequency [Hz]')
        ax.set_ylabel(r'Equivalent load, $L_{eq}/A$ [-]')
        #ax.set_ylim([0,7])
        ax.tick_params(direction='in')
        ax.legend()
        ax.set_title('Tools - Fatigue - Sinusoids')

    if sumSine:
        # --------------------------------------------------------------------------------
        # --- Sum of two sinusoids
        # --------------------------------------------------------------------------------
        m     = 2   # Wohler slope
        nPerT = 100 # Number of points per period
        Teq   = 1  # Equivalent period [s]
        nBins = 10 # Number of bins
        nT1   = 10 # Number of periods
        nT2   = 20 # Number of periods
        T1 = 10
        T2 = 5
        A1 = 3 # Amplitude
        A2 = 5 # Amplitude
        # --- Signals
        time1 = np.linspace(0, nT1*T1, nPerT*nT1+1)
        time2 = np.linspace(0, nT2*T2, nPerT*nT2+1)
        signal1 = A1 * np.sin(2*np.pi/T1*time1)
        signal2 = A2 * np.sin(2*np.pi/T2*time2)
        # --- Individual Leq
        print('----------------- SIGNAL 1')
        Leq1 = equivalent_load(time1, signal1, m=m, Teq=Teq, bins=nBins, method='rainflow_windap')
        if hasFatpack:
            Leq2 = equivalent_load(time1, signal1, m=m, Teq=Teq, bins=nBins, method='fatpack')
        print('>>> Leq1   ',Leq1)
        print('>>> Leq2   ',Leq2)
        DEL1 = (2*A1)**m * nT1/time1[-1]
        Leq_th = (DEL1)**(1/m)
        print('>>> Leq TH ',Leq_th)
        print('----------------- SIGNAL 2')
        Leq1 = equivalent_load(time2, signal2, m=m, Teq=Teq, bins=nBins, method='rainflow_windap')
        if hasFatpack:
            Leq2 = equivalent_load(time2, signal2, m=m, Teq=Teq, bins=nBins, method='fatpack')
        print('>>> Leq1   ',Leq1)
        print('>>> Leq2   ',Leq2)
        DEL2 = (2*A2)**m * nT2/time2[-1]
        Leq_th = (DEL2)**(1/m)
        print('>>> Leq TH ',Leq_th)
        # --- Concatenation
        print('----------------- CONCATENATION')
        signal = np.concatenate((signal1, signal2))
        time   = np.concatenate((time1, time2+time1[-1]))
        T_all=time[-1]
        Leq1 = equivalent_load(time, signal, m=m, Teq=Teq, bins=nBins, method='rainflow_windap')
        if hasFatpack:
            Leq2 = equivalent_load(time, signal, m=m, Teq=Teq, bins=nBins, method='fatpack')
        print('>>> Leq1   ',Leq1)
        print('>>> Leq2   ',Leq2)
        DEL1 = (2*A1)**m * nT1/T_all  
        DEL2 = (2*A2)**m * nT2/T_all  
        Leq_th = (DEL1+DEL2)**(1/m)
        print('>>> Leq TH ',Leq_th)


    return vf, Leq1, Leq_ref

vf, Leq, Leq_ref =  main(fstudy=True, sumSine=False)

if __name__ == '__main__':
    plt.show()
if __name__=="__test__":
    np.testing.assert_array_almost_equal(Leq, Leq_ref,1)
if __name__=="__export__":
    from welib.tools.repo import export_figs_callback
    export_figs_callback(__file__)



