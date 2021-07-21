import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from vortilib.tools.colors import fColrs
from vortilib.elements.VortexSegment import *

# --- Comparison of regularization
z0=1
Pa      = np.array([[ 0, 0, -z0]])
Pb      = np.array([[ 0, 0,  z0]])
L       = np.linalg.norm(Pa-Pb)
z0      = L/2
Gamma   = 1
Epsilon = 0.5*L
Xcp = np.linspace(0,2*L,100)[:]
Ycp = Xcp*0
Zcp = Xcp*0
U0x, U0y, U0z = vs_u(Xcp, Ycp, Zcp, Pa, Pb, Gamma=Gamma, RegFunction=0, RegParam=Epsilon)
U1x, U1y, U1z = vs_u(Xcp, Ycp, Zcp, Pa, Pb, Gamma=Gamma, RegFunction=1, RegParam=Epsilon)
U2x, U2y, U2z = vs_u(Xcp, Ycp, Zcp, Pa, Pb, Gamma=Gamma, RegFunction=2, RegParam=Epsilon)
U3x, U3y, U3z = vs_u(Xcp, Ycp, Zcp, Pa, Pb, Gamma=Gamma, RegFunction=3, RegParam=Epsilon)
U4x, U4y, U4z = vs_u(Xcp, Ycp, Zcp, Pa, Pb, Gamma=Gamma, RegFunction=4, RegParam=Epsilon*1)
U_th = Gamma*(L/2)/(2*np.pi*Xcp * np.sqrt(Xcp**2 + z0**2))


# plot
fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.2,4.6)) # (6.4,4.8)
fig.subplots_adjust(left=0.12, right=0.98, top=0.96, bottom=0.12, hspace=0.20, wspace=0.20)
ax.plot(Xcp[1:]/L  ,U0y[1:] / Gamma*np.pi*L, '-', color=fColrs(1), label = 'Singular'       )
ax.plot(Xcp/L  ,U1y / Gamma*np.pi*L, '-.', color=fColrs(2), label = 'Rankine'    )
ax.plot(Xcp/L  ,U2y / Gamma*np.pi*L, '-', color=fColrs(3), label = 'Lamb-Oseen' )
ax.plot(Xcp/L  ,U3y / Gamma*np.pi*L, '--', color=fColrs(4), label = 'Vatistas n=2'   )
ax.plot(Xcp/L  ,U4y / Gamma*np.pi*L, ':', color=fColrs(5), label = 'Denominator')
# ax.plot(Xcp/L  ,U_th/ Gamma*np.pi*L, '--',color='k', label = 'Theory'       )
ax.set_xlabel(r'$\rho/L$ [-]')
ax.set_ylabel(r'$u \pi L / \Gamma$ [-]')
ax.set_xticks(np.arange(0,2.1,0.5))
ax.set_xticklabels(['0','0.5','1','1.5','2'])
ax.legend()
ax.grid(ls=':')
ax.set_ylim([0,1])
ax.set_xlim([0,2])
ax.tick_params(direction='in')
fig.savefig('figs/VortexFilamentRegularization.pdf')
plt.show()
