import numpy as np
import matplotlib.pyplot as plt
from welib.standards.IEC import *


ws=np.linspace(0,30,100)
fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
ax.plot(ws,ETM(ws,'A')*100, '-', label='ETM A', c='k')
ax.plot(ws,ETM(ws,'B')*100, '--', label='ETM B', c='k')
ax.plot(ws,ETM(ws,'C')*100, '-.', label='ETM C', c='k')
ax.plot(ws,NTM(ws,'A')*100, '-', label='NTM A',c=[0.5,0.5,0.5])
ax.plot(ws,NTM(ws,'B')*100, '--', label='NTM B',c=[0.5,0.5,0.5])
ax.plot(ws,NTM(ws,'C')*100, '-.', label='NTM C',c=[0.5,0.5,0.5])
ax.set_ylim([0,88])
ax.set_xlim([0,30])
ax.set_xlabel('Wind speed [m/s]')
ax.set_ylabel('Turbulence intensity [%]')
ax.legend()
ax.tick_params(direction='in')
ax.set_title('IEC Standards - Turbulence classes')


if __name__ == '__main__':
    plt.show()
if __name__ == '__test__':
    pass
if __name__=="__export__":
    from welib.tools.repo import export_figs_callback
    export_figs_callback(__file__)
