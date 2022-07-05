"""
Plot IEA Extreme Operating Gusts (EOG) for various turbulence class.
NOTE: the turbine class might have a limited/no influence depending  
on the value of V_e1.

"""

import numpy as np
import matplotlib.pyplot as plt
from welib.standards.IEC import *


WS    = 10
D     = 126
z_hub = 90

# --- Write a Wind file for OpenFAST
wnd_file = '_EOG.wnd' 
t, Vgust1, _, _ = EOG(WS, D, z_hub, WT_class='I', turbulence_class='A', tStart=20, filename=wnd_file)


# --- Retrieve Gust velocity for various Classes
t, Vgust1, _, _ = EOG(WS, D, z_hub, turbulence_class='A')
t, Vgust2, _, _ = EOG(WS, D, z_hub, turbulence_class='B')
t, Vgust3, _, _ = EOG(WS, D, z_hub, turbulence_class='C')

fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
ax.plot(t, Vgust1, '-',  label='EOG Turbulence class A', c='k')
ax.plot(t, Vgust2, '--', label='EOG Turbulence class B', c='k')
ax.plot(t, Vgust3, ':',  label='EOG Turbulence class C', c='k')
ax.set_ylim([9.8,10.4])
ax.set_xlim([0,11])
ax.set_xlabel('Time [s]')
ax.set_ylabel('Wind speed [m/s]')
ax.legend()
ax.tick_params(direction='in')
ax.set_title('IEC Standards - Extreme operating gusts')


if __name__ == '__main__':
    plt.show()
if __name__ == '__test__':
    try:
        os.remove(wnd_file)
    except:
        pass
if __name__=="__export__":
    from welib.tools.repo import export_figs_callback
    export_figs_callback(__file__)
