"""
The package welib.beam.cantilever1d can compute deflection of a cantilever beam with arbitrary stiffness distribution and arbitrary loading.

In this example, we use a beam with uniform stiffness distribution to compare with analytical results.
"""
import numpy as np
import matplotlib.pyplot as plt
from welib.beams.theory import UniformBeamDeflection
import welib.beams.cantilever1d as c1d
from welib.tools.colors import python_colors

def main(plot_u=True, plot_all=False):
    params = {'F': 15, 'EI': 4000, 'p': 5, 'p0': 8, 'L':10}
    z = np.linspace(0, params['L'], 50)
    EI = np.ones_like(z) * params['EI'] # Constant EI to comare with uniform beam theory

    # --- Compute deflection for various loading types
    types  = ['point_load']
    types += ['uniform']
    types += ['increasing']
    types += ['decreasing']
    u_th = np.zeros((len(types), len(z)))
    u_nm = np.zeros((len(types), len(z)))
    for i, loading_type in enumerate(types):
        # Theory for uniform beam
        u_th[i,:], _, _, _, _, p = UniformBeamDeflection('clamped-free', loading_type,  z, params)
        # Numerical deflection for arbitrary EI(z), here EI is constant
        if loading_type=='point_load':
            u_nm[i,:], _, _, _, _ = c1d.deflection(p*0, z, EI, Ftip = params['F'])
        else:
            u_nm[i,:], _, _, _, _ = c1d.deflection(p, z, EI)

    # --- Compare deflections from theory and numerical integration
    if plot_u:
        fig, ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8))
        fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.30)
        for i, loading_type in enumerate(types):
            ax.plot(z, u_th[i,:], '-', c=python_colors(i), label=loading_type.replace('_',' ') + ' (theory)')
            ax.plot(z, u_nm[i,:], '.', c=python_colors(i), label='Numerical integration' if i==0 else None)
        ax.set_ylabel('Deflection [m]')
        ax.set_xlabel('Beam span [m]')
        ax.legend()
        ax.set_title('Beam - 1D - Analytical and numerical deflections')

    # --- Compare Deflection, shape, curvature, moment, shear
    if plot_all:
        for i, loading_type in enumerate(types):
            # Theory for uniform beam
            u_th, t_th, k_th, S_th, M_th, p = UniformBeamDeflection('clamped-free', loading_type,  z, params)
            # Numerical deflection for arbitrary EI(z), here EI is constant
            if loading_type=='point_load':
                u_nm, t_nm, k_nm, S_nm, M_nm = c1d.deflection(p*0, z, params['EI'], Ftip=params['F'])
            else:
                u_nm, t_nm, k_nm, S_nm, M_nm = c1d.deflection(p, z, params['EI'])

            fig,axes = plt.subplots(3, 2, sharex=True, figsize=(6.4,8.8))
            fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.30)
            ax = axes[0,0]
            ax.plot(z, u_th, 'k-', label='')
            ax.plot(z, u_nm, '--', label='')
            ax.set_ylabel('Deflection [m]')
            ax = axes[1,0]
            ax.plot(z, t_th, 'k-', label='')
            ax.plot(z, t_nm, '--', label='')
            ax.set_ylabel('Slope [-]')
            ax = axes[2,0]
            ax.plot(z, k_th, 'k-', label='')
            ax.plot(z, k_nm, '--', label='')
            ax.set_ylabel('Curvature [1/m]')
            ax = axes[0,1]
            ax.plot(z, M_th, 'k-', label='')
            ax.plot(z, M_nm, '--', label='')
            ax.set_ylabel('Moment [Nm]')
            ax = axes[1,1]
            ax.plot(z, S_th, 'k-', label='')
            ax.plot(z, S_nm, '--', label='')
            ax.set_ylabel('Shear [N]')
            ax = axes[2,1]
            ax.plot(z, p  , label='')
            ax.set_ylabel('Loading [N/m]')
            ax.set_xlabel('Beam Span [m]')
            fig.suptitle(loading_type)
    return u_nm, u_th

if __name__ == '__main__':
    u_nm, u_th = main(plot_u=True, plot_all=True)
    plt.show()

if __name__ == '__test__':
    u_nm, u_th = main(plot_u=False)
    np.testing.assert_almost_equal(u_nm[:,:], u_th[:,:], 3)

if __name__=='__export__':
    u_nm, u_th = main(plot_u=True)
    from welib.tools.repo import export_figs_callback
    export_figs_callback(__file__)
