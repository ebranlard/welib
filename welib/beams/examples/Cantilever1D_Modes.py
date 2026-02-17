""" 
The package welib.beam.cantilever1d can compute modes of a cantilever beam with arbitrary stiffness distribution.

In this example, we use a beam with uniform stiffness distribution to compare with analytical results.
"""
import numpy as np
import matplotlib.pyplot as plt
from welib.beams.theory import UniformBeamBendingModes
import welib.beams.cantilever1d as c1d
from welib.tools.colors import python_colors

def main(nModes=5, plot=False):
    L  = 100    
    EI = 1.86e+12 # To compare with theory, uniform mass and stiffness
    m  = 8.82e+03
    z = np.linspace(0, L, 100)
    rho, A = m , 1

    # Theoretical modes
    f_th, _, modes_th, _, _ = UniformBeamBendingModes('unloaded-clamped-free', EI, rho, A, L, x=z, norm='tip', nModes=nModes)

    # Numerical modes
    f_nm, modes_nm = c1d.compute_modes(nModes, z, EI, m)

    if plot:
        fig, ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8))
        fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.30)
        for i in range(nModes):
            print(f'Mode {i+1}: Frequency = {f_th[i]:.2f} {f_nm[i]:.2f} Hz')
            ax.plot(z, modes_th[i,:], '-', c=python_colors(i), label=f'Mode {i+1} (theory)')
            ax.plot(z, modes_nm[i,:], '.', c=python_colors(i), label='Numerical intergration' if i==0 else None)

        ax.set_title('Beam - 1D - Analytical and numerical modes')
        ax.legend()

    return f_th, modes_th, f_nm, modes_nm





if __name__ == '__main__':
    f_th, modes_th, f_nm, modes_nm = main(nModes=3, plot=False)
    plt.show()

if __name__ == '__test__':
    f_th, modes_th, f_nm, modes_nm = main(nModes=3, plot=False)
    np.testing.assert_almost_equal(modes_nm, modes_th, 3)

if __name__=='__export__':
    f_th, modes_th, f_nm, modes_nm = main(nModes=3, plot=True)
    from welib.tools.repo import export_figs_callback
    export_figs_callback(__file__)
