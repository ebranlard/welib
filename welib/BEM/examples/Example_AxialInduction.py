""" 
Plot different high-thrust corrections used in BEM
"""
# --- Common libraries 
import numpy as np
import matplotlib.pyplot as plt

# --- Local libraries
from welib.BEM.highthrust import *
from welib.tools.figure import defaultRC; defaultRC();

def main(test=False):
    Ct=np.linspace(0,2,50)
    a =np.linspace(0,1,50)
    Ct_MT = 4*a*(1-a)

    fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
    fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
    # Functions that depend on a only
    ax.plot(a ,Ct_MT,'k-' ,label = 'Momentum theory'          )
    ax.plot(a ,Ct_a(a,method='Glauert'),'-' ,label = 'Glauert (ac=1/3)')
    ax.plot(a ,Ct_a(a,method='Spera')  ,'.' ,label = 'Spera (ac=0.3)')
    # Functions that depend on Ct only
    ax.plot(a_Ct(Ct,method = 'AeroDyn'         ),Ct,'-' ,label = 'AeroDyn'          )
    ax.plot(a_Ct(Ct,method = 'HAWC2'           ),Ct,'--',label = 'HAWC2'            )
    ax.plot(a_Ct(Ct,method = 'WEHandbook'      ),Ct,':' ,label = 'Handbook'         )
    ax.plot(a_Ct(Ct,method = 'GlauertEmpirical'),Ct,'-.',label = 'Glauert Empirical')
    ax.set_xlabel('Axial induction, a [-]')
    ax.set_ylabel('Thrust coefficient, Ct [-]')
    ax.set_xlim([0,1])
    ax.set_ylim([0,2])
    ax.legend()
    ax.grid()
    ax.set_title('BEM Steady - High thrust correction')


if __name__=="__main__":
    main()
    plt.show()
if __name__=="__test__":
    main()
if __name__=="__export__":
    main()
    from welib.tools.repo import export_figs_callback
    export_figs_callback(__file__)
