""" 
Plot regularization for 3D particles

"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from vortilib.tools.colors import fColrs
from vortilib.elements.VortexParticle import *


# --- Comparison of regularization
PPart=np.array([0,0,0])
L       = 1
Gamma   = 1
Alpha   = [0,0,L*Gamma]
Epsilon = 0.5*L
nCPs=100
CPs = np.zeros((nCPs,3))
CPs[:,0] = np.linspace(0,2*L,nCPs)
U0 = vp_u(CPs, PPart, Alpha, RegFunction=0, RegParam=Epsilon)
U1 = vp_u(CPs, PPart, Alpha, RegFunction=1, RegParam=Epsilon)
U2 = vp_u(CPs, PPart, Alpha, RegFunction=2, RegParam=Epsilon)
U_th = Alpha[2]*CPs[:,0]/(4*np.pi*np.sqrt(CPs[:,0]**2 + CPs[:,1]**2)**3)



fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.2,4.6)) # (6.4,4.8)
fig.subplots_adjust(left=0.12, right=0.97, top=0.96, bottom=0.12, hspace=0.20, wspace=0.20)
ax.plot(CPs[1:,0]/L, U0[1:, 1], '-',  color=fColrs(1), label = 'Singular'       )
ax.plot(CPs[: ,0]/L, U1[: , 1], '-.', color=fColrs(2), label = 'Exponential'    )
ax.plot(CPs[: ,0]/L, U2[: , 1], '-',  color=fColrs(3), label = 'Compact' )
#ax.plot(CPs[: ,0]/L, U_th      , 'k-'                , label = 'Theory' )
ax.set_xlabel(r'Radial distance $r$ [m]')
ax.set_ylabel(r'Orthogonal velocity $u$ [m/s]')
ax.set_xticks(np.arange(0,2.1,0.5))
ax.set_xticklabels(['0','0.5','1','1.5','2'])
ax.grid(ls=':')
ax.legend()
ax.set_ylim([0,0.5])
ax.set_xlim([0,2])
ax.tick_params(direction='in')
fig.savefig('figs/VortexParticleRegularization.pdf')
plt.show()

