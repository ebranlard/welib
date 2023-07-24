""" 
Illustrate the usage of the StochasticVector class, and the standaradization.

- Create a 2d stochastic vector made of two correlated normal variables

- Apply the standardizing change of variable to convert the variables to N(0,1)
  (independent normal stochastic variables of unit standard deviation and zero mean)

"""

import numpy as np
import matplotlib.pyplot as plt
from welib.essentials import *
from welib.stoch.vector import StochasticVector
from welib.stoch.distribution import gaussian_pdf
from welib.stoch.distribution import gaussian2d_pdf
from welib.stoch.distribution import gaussian2d_covar

sigmas = np.array([2,3])
mus    = np.array([1,3])
rho=0.5
C = gaussian2d_covar(sigmas, rho)

def main(nDiscr=100):
    f_pdf = lambda x1,x2: gaussian2d_pdf(x1, x2, mus, sigmas, rho=rho)

    vec = StochasticVector(dimension=2, name='Input - X', nDiscr=nDiscr, xmin=-20, xmax=20)
    vec.set_pdf_f(f_pdf)
    ax = vec.plot_pdfs()
    #x=np.linspace(-30,30,100)
    #ax.plot(x, gaussian_pdf(x, mu=mus[0], sig=sigmas[0]), 'k:')
    #ax.plot(x, gaussian_pdf(x, mu=mus[1], sig=sigmas[1]), 'k:')

    # --- Find a and B so that the vector becomes independent normal stochastic variables
    a, B, g, gradg, ginv = vec.standardize_transform()
    vec2 = vec.new_from_bijection(ginv=ginv, gradg=gradg, name='Output - Y', nDiscr=nDiscr, xmin=-10, xmax=10)
    vec2.plot_pdfs(ax=ax, sty=['k--', 'k.'])
    #x=np.linspace(-10,10,100)
    #ax.plot(x, gaussian_pdf(x, mu=0, sig=1), 'k--', label='N(0,1)')


    # --- Using analytical expression
    # vec3 = StochasticVector(dimension=2, name='Z', nDiscr=100, xmin=-10, xmax=10)
    # BCBt = B.dot(C).dot(B.T)
    # sig = np.linalg.det(BCBt) # Should be 1
    # BCBt_inv = B.dot(C).dot(B.T)  # Should be eye(2)
    # # offset = a + B.dot(mu) should be 0,0
    # def pdf_f(*y):
    #     scale = 1/(2*np.pi)
    #     dy = np.array(y)
    #     y2 = (dy).dot(dy)
    #     f = scale * np.exp(-0.5 * y2 )
    #     return f
    # vec3.set_pdf_f(pdf_f)
    # vec3.plot_pdfs(ax=ax, sty=':')
    # vec3.check_pdf()


    ax.set_title('Stochastic - Independent and standardized')
    ax.set_xlim([-10,10])
    # 
    return vec, vec2



if __name__ == '__main__':
    vec, vec2 = main()
    print(vec)
    print(vec2)
    plt.show()

if __name__ == '__test__':
    vec, vec2 = main(nDiscr=20)
    np.testing.assert_almost_equal(vec.mean, mus, 3)
    np.testing.assert_almost_equal(vec.var,  C , 2)
    np.testing.assert_almost_equal(vec2.mean, [0.,0.], 3)
    np.testing.assert_almost_equal(vec2.var, np.eye(2), 3)

if __name__=="__export__":
    main(nDisc=100)
    from welib.tools.repo import export_figs_callback
    export_figs_callback(__file__)
