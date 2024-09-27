""" 

See Example 2-1
"""

import matplotlib.pyplot as plt
import numpy as np

from welib.stoch.distribution import *

def main():
    sigmas = [5,6]
    mus    = [0,3]
    rho=0.7
    # TODO TODO Now we should use the class Vector for that and meshgrid_pdf.


    # --- joint probability density function
    x1 = np.linspace(-20,20,100)
    x2 = np.linspace(-20,20,80)
    X1, X2, F     = gaussian2d_pdf_meshgrid(x1, x2, mus, sigmas, rho)
    _, _, F_gen = gaussian2d_pdf_meshgrid(x1, x2, mus, sigmas, rho, method='general')
    print('Error norm:', np.linalg.norm(F-F_gen))

    # --- Marginal distributions
    F1 = np.trapezoid(F.T, x2)
    F2 = np.trapezoid(F, x1)
    F1th = gaussian_pdf(x1, mus[0], sigmas[0]**2)
    F2th = gaussian_pdf(x2, mus[1], sigmas[1]**2)

    # Joint distribution if variables were independent
    _, _, F_indpt = gaussian2d_pdf_meshgrid(x1, x2, mus, sigmas, rho=0)
#     F_indpt = (np.atleast_2d(F1th).T).dot(np.atleast_2d(F2th)).T
    F1_indpt = np.trapezoid(F_indpt.T, x2)
    F2_indpt = np.trapezoid(F_indpt, x1)


    # --- Contour
    import matplotlib.pyplot as plt
    fig,axes = plt.subplots(2, 2, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
    ax = axes[0,0]
    ax.contourf(X1,X2,F)
    ax.plot(mus[0], mus[1], 'ko')
    #ax.plot(x1, mus[1]+(x1-mus[0])*rho, 'k--') # need influence of sigma
    ax.axis('equal')
    ax.set_title('Correlated')

    ax = axes[0,1]
    ax.contourf(X1,X2,F_indpt)
    ax.plot(mus[0], mus[1], 'ko')
    ax.axis('equal')
    ax.set_title('Independent')


    ax = axes[1,0]
    ax.plot(x1, F1th, 'k-')
    ax.plot(x2, F2th, 'k-')
    ax.plot(x1, F1, '--', label='F1')
    ax.plot(x2, F2, '--', label='F2')
    ax.set_xlabel('x')
    ax.set_ylabel('')
    ax.legend()
    

    ax = axes[1,1]
    ax.plot(x1, F1th, 'k-')
    ax.plot(x2, F2th, 'k-')
    ax.plot(x1, F1_indpt, '--', label='F1')
    ax.plot(x2, F2_indpt, '--', label='F2')
    ax.set_xlabel('x')
    ax.set_ylabel('')
    ax.legend()
    
    fig.suptitle('Stochastic - Correlated variables')

if __name__ == '__main__':
    main()
    plt.show()


if __name__ == '__test__':
    main()
    # TODO test it
    #np.testing.assert_almost_equal(vec2.mean, [0.,0.], 3)

if __name__=="__export__":
    main()
    from welib.tools.repo import export_figs_callback
    export_figs_callback(__file__)

