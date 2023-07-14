import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# Local 
from welib.wt_theory.wakeexpansion import *

    
def plotWakeExpansionModels(CT=0.8, fractions=None):
    xb = np.linspace(0, 10, 200) 
    rwC = wake_expansion_cylinder(xb=xb, CT=CT)
    rwR = wake_expansion_Rathmann(xb=xb, CT=CT)
    rwF = wake_expansion_Frandsen(xb=xb, CT=CT)
    rwK = wake_expansion_kxm     (xb=xb, k=0.05, m=1)
    rw0 = wake_expansion_momentum(CT=CT)

    fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
    fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
    ax.plot(xb, xb*0+rw0, 'k-', label='NW - Momentum theory')
    ax.plot(xb, rwC, 'k--'    , label='NW - Vortex Cylinder')
    ax.plot(xb, rwF, ':'      , label='FW - Frandsen')
    ax.plot(xb, rwR, ':'      , label='FW - Rathmann')
    ax.plot(xb, rwK, ':'      , label=r'FW $k x^m$')

    if fractions is not None:
        for fraction in fractions:
            expansion = 1 + fraction * (rw0-1)
            x = downstreamDistanceForGivenExpansion(CT, expansion, model='cylinder')
            ax.plot([0,x],[expansion, expansion], 'k:')
            ax.plot([x,x],[expansion, 0        ], 'k:')
            print('CT={:4.2f}, Expansion={:4.2f}, {:2.0f}% reached at x/R={:4.2f}'. format(CT, rw0, fraction*100, x))

    ax.set_xlabel('z/R [-]')
    ax.set_ylabel('r/R [-]')
    ax.set_ylim([0,2.5])
    ax.set_xlim([0,np.max(xb)])
    ax.legend()
    ax.set_title('WT Theory - Wake Expansion Models')
    return fig


if __name__ == '__main__':
    fig = plotWakeExpansionModels(CT=0.8, fractions=[0.5])

    # Show when the wake expansion is partially reached 
    CTs = [0.4, 0.6, 8/9]
    fractions = [0.5, 0.9, 0.99]
    for fraction in fractions:
        for CT in CTs:
            rw0 = wake_expansion_momentum(CT=CT)
            expansion = 1 + fraction * (rw0-1)
            x = downstreamDistanceForGivenExpansion(CT, expansion, model='cylinder', method='analytical')
            print('CT={:4.2f}, Expansion={:4.2f}, {:2.0f}% reached at x/R={:4.2f}'. format(CT, rw0, fraction*100, x))

#     plt.show()


if __name__=="__test__":
    fig = plotWakeExpansionModels(CT=0.8)
    pass
if __name__=="__export__":
    fig = plotWakeExpansionModels(CT=0.8)
    from welib.tools.repo import export_figs_callback
    export_figs_callback(__file__)
