import numpy as np
import matplotlib.pyplot as plt
from welib.airfoils.DynamicStall import wagner

# --- Wagner function for different set of constants
tau = np.linspace(0,20,100)
Cl_Jones = wagner(tau, constants='Jones')
Cl_FAST  = wagner(tau, constants='OpenFAST')
#Cl_Misc  = wagner(tau, constants=None, A1=0.3, A2=0.7, b1=0.14, b2=0.53)

# --- Plot
fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
ax.plot(tau, Cl_Jones, '-' ,  label='Jones')
ax.plot(tau, Cl_FAST,  'k--' , label='OpenFAST')
ax.tick_params(direction='in', top=True, right=True)
ax.set_xlabel(r'$\tau$ [-]')
ax.set_ylabel(r'$C_l$ [-]')
ax.set_title(r'Airfoils - Wagner function')
ax.legend()

if __name__ == '__main__':
    plt.show()
if __name__ == '__test__':
    pass
if __name__=="__export__":
    from welib.tools.repo import export_figs_callback
    export_figs_callback(__file__)

