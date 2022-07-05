import numpy as np
import matplotlib.pyplot as plt
from welib.airfoils.DynamicStall import wagner

tau = np.linspace(0,20,100)
Cl_Jones = wagner(tau, constants='Jones')
Cl_FAST  = wagner(tau, constants='OpenFAST')

fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
ax.plot(tau, Cl_Jones, '-' ,  label='Jones')
ax.plot(tau, Cl_FAST,  '--' , label='OpenFAST')
ax.set_xlabel(r'$\tau$ [-]')
ax.set_ylabel(r'$C_l$ [-]')
ax.legend()

if __name__ == '__main__':
    plt.show()


