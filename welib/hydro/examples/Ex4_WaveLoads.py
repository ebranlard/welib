""" 
Compute inline/total hydrodynamic force and moments on a monopile using Morison's equation
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from fractions import Fraction
import matplotlib as mpl
# Local 
from welib.tools.figure import defaultRC; defaultRC();
from welib.tools.colors import python_colors
from welib.hydro.wavekin import *
from welib.hydro.morison import *


# --- Parameters
g   = 9.81                # gravity [m/s^2]
h   = 30.                 # water depth [m]
rho = 1025                # water density
D   = 6                   # monopile diameter [m]
CD  = 1                   # given
CM  = 2                   # 
a   = 3                   # wave peak amplitude [m]
T   = 12.                 # period [s]
eps = 0                   # phase shift [rad]
f   = 1./T
k   = wavenumber(f, h, g)

nz = 30    # number of points used in the z direction to compute loads
z_ref = -h # reference point for moment calculation

# --------------------------------------------------------------------------------}
# --- Inline force and moments as function of time, with or without Wheeler stretching
# --------------------------------------------------------------------------------{
time = np.linspace(0,T,9)

fig1,axes1 = plt.subplots(2, 4, sharex=True, sharey=True, figsize=(12.8,4.8)) # (6.4,4.8)
fig1.subplots_adjust(left=0.05, right=0.99, top=0.95, bottom=0.09, hspace=0.26, wspace=0.11)

fig2,axes2 = plt.subplots(2, 4, sharex=True, sharey=True, figsize=(12.8,4.8)) # (6.4,4.8)
fig2.subplots_adjust(left=0.05, right=0.99, top=0.95, bottom=0.09, hspace=0.26, wspace=0.11)

XLIM =[-75,75] # For inline force
XLIMM=[-2500,2500] # For inline moment

for it, t in enumerate(time[:-1]):
    # Wave kinematics
    eta  = elevation2d(a, f, k, eps, t, x=0)
    z    = np.linspace(-h, eta, nz)
    u, du = kinematics2d(a, f, k, eps, h, t, z, Wheeler=True, eta=eta)
    u0, du0 = kinematics2d(a, f, k, eps, h, t, z)
    # Wave loads with wheeler
    p_tot     = inline_load(u, du, D, CD  , CM  , rho)
    p_inertia = inline_load(u, du, D, CD*0, CM  , rho)
    p_drag    = inline_load(u, du, D, CD  , CM*0, rho)
    dM        = p_tot * (z-z_ref) # [Nm/m] 

    # Wave loads without Wheeler
    p_tot0    = inline_load(u0, du0, D, CD  , CM  , rho)
    p_inertia0= inline_load(u0, du0, D, CD*0, CM  , rho)
    p_drag0   = inline_load(u0, du0, D, CD  , CM*0, rho)
    dM0       = p_tot0* (z-z_ref) # [Nm/m] 


    # Plot inline force
    ax=axes1[int(it/4),np.mod(it,4)]
    ax.plot(p_inertia/1000,z, '-', c=python_colors(0), label = r'$f_{inertia}$')
    ax.plot(p_drag/1000   ,z, '-', c=python_colors(3), label = r'$f_{drag}$')
    ax.plot(p_tot/1000    ,z, 'k-'                    , label = r'$f_{tot}$')
    ax.plot(p_inertia0/1000,z, '+', c=python_colors(0))
    ax.plot(p_drag0/1000   ,z, '+', c=python_colors(3))
    ax.plot(p_tot0/1000    ,z, 'k+'                    )
    ax.set_title('t/T={}'.format(Fraction(t/T)))
    if it==0:
        ax.legend()
    ax.plot(XLIM,[0,0],'k')
    ax.plot(XLIM,[a,a],'k--')
    ax.plot(XLIM,[-a,-a],'k--')

    # Plot inline moment
    ax=axes2[int(it/4),np.mod(it,4)]
    ax.plot(dM/1000     ,z, 'k-', label = r'$dM_{tot}$ with Wheeler')
    ax.plot(dM0/1000    ,z, 'k+', label = r'$dM_{tot}$ no-correction')
    ax.set_title('t/T={}'.format(Fraction(t/T)))
    if it==0:
        ax.legend()
    ax.plot(XLIMM,[0,0],'k')
    ax.plot(XLIMM,[a,a],'k--')
    ax.plot(XLIMM,[-a,-a],'k--')

axes1[0,0].set_xlim(XLIM)
axes1[0,0].set_ylim([-h,a+1])
axes1[0,0].set_ylabel('Depth z [m]')
axes1[1,0].set_ylabel('Depth z [m]')
axes1[1,0].set_xlabel('Inline force [kN/m]')
axes1[1,1].set_xlabel('Inline force [kN/m]')
axes1[1,2].set_xlabel('Inline force [kN/m]')
axes1[1,3].set_xlabel('Inline force [kN/m]')

axes2[0,0].set_xlim(XLIMM)
axes2[0,0].set_ylim([-h,a+1])
axes2[0,0].set_ylabel('Depth z [m]')
axes2[1,0].set_ylabel('Depth z [m]')
axes2[1,0].set_xlabel('Inline moment [kNm/m]')
axes2[1,1].set_xlabel('Inline moment [kNm/m]')
axes2[1,2].set_xlabel('Inline moment [kNm/m]')
axes2[1,3].set_xlabel('Inline moment [kNm/m]')

# --------------------------------------------------------------------------------}
# --- Integrated force and sea bed moment over a period
# --------------------------------------------------------------------------------{
time = np.linspace(0,T,1000)

veta = np.zeros(time.shape)
vF   = np.zeros(time.shape)
vM   = np.zeros(time.shape)
vF0  = np.zeros(time.shape)
vM0  = np.zeros(time.shape)

XLIM =[-75,75] # For inline force
XLIMM=[-2500,2500] # For inline moment

a=6 # NOTE: increased amplitude here to see Effect of Wheeler

for it, t in enumerate(time):
    # Wave kinematics
    veta[it] = elevation2d(a, f, k, eps, t, x              = 0)
    z        = np.linspace(-h, veta[it], nz)
    u, du    = kinematics2d(a, f, k, eps, h, t, z, Wheeler = True, eta = veta[it])
    u0, du0  = kinematics2d(a, f, k, eps, h, t, z)
    # Wave loads with Wheeler
    p_tot  = inline_load(u, du, D, CD  , CM  , rho)
    vF[it] = np.trapz(p_tot            , z) # [N]
    vM[it] = np.trapz(p_tot * (z-z_ref), z) # [Nm]
    # Wave loads without Wheeler
    p_tot0 = inline_load(u0, du0, D, CD  , CM  , rho)
    vF0[it] = np.trapz(p_tot0            , z) # [N]
    vM0[it] = np.trapz(p_tot0 * (z-z_ref), z) # [Nm]

# Plot
fig, axes = plt.subplots(3, 1, sharex=True, figsize=(6.4,4.8)) # (6.4,4.8)
fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.14, wspace=0.20)
ax = axes[0]
ax.plot(time/T, veta,'k-')
ax.set_ylabel('Elevation [m]')
ax.grid(True)
ax = axes[1]
ax.plot(time/T, vF0/1e6     , label='no-correction')
ax.plot(time/T, vF /1e6,'k-', label='Wheeler')
ax.set_ylabel('Force [MN]')
ax.legend()
ax.grid(True)
ax = axes[2]
ax.plot(time/T, vM0/1e6     , label='no-correction')
ax.plot(time/T, vM /1e6,'k-', label='Wheeler')
ax.set_ylabel('Sea-bed moment [MNm]')
ax.set_xlabel('Dimensionless time t/T [-]')
ax.legend()
ax.grid(True)

plt.suptitle('Hydro - Morison loads on monopile')



if __name__ == '__main__':
    plt.show()

if __name__=="__test__":
    pass

if __name__ == '__export__':
    plt.close(fig1)
    plt.close(fig2)
    from welib.tools.repo import export_figs_callback
    export_figs_callback(__file__)

