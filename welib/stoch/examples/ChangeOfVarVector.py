""" 
Illustrate the usage of the StochasticVector class, and the standardization.


Case 1: 

 - Create a 2d stochastic vector made of two correlated normal variables:
    X = [ X1 ] = [ Normal(mu1, sig1)  ]
        [ X2 ]   [ Normal(mu2, sig2)  ]
   
 
 - Apply the standardizing change of variable to convert the variables to N(0,1)
   (independent normal stochastic variables of unit standard deviation and zero mean)
    Y = [ Y1 ] = [ Normal(0, 1)  ]  = a + BX
        [ Y2 ]   [ Normal(0, 1)  ]

Case 2: 
  - Input:
    X = [ R   ] = [ Rayleigh(sig)  ]
        [ Phi ]   [ Uniform(0,2pi) ]

  - Output:
    Y = [ Y1 ] = [ R cos(alpha+Phi) ]  = g(X)
        [ Y2 ]   [ R sin(alpha+Phi) ]

"""

import numpy as np
import matplotlib.pyplot as plt
from welib.essentials import *
from welib.stoch.vector import StochasticVector
from welib.stoch.distribution import uniform_pdf
from welib.stoch.distribution import rayleigh_pdf
from welib.stoch.distribution import gaussian_pdf
from welib.stoch.distribution import gaussian2d_pdf
from welib.stoch.distribution import gaussian2d_covar

sigmas = np.array([2,3])
mus    = np.array([1,3])
rho=0.5
rsig = 2
C = gaussian2d_covar(sigmas, rho)

def main(nDiscr=100):
    # --- Case 1
    # Normal variables, standardized

    # --- Setup input X
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

    ax.set_title('Stochastic - Independent and standardized')
    ax.set_xlim([-10,10])
    return vec, vec2

def main2(nDiscr=100):
    # --- Case 2
    # R: Rayleigh
    # Phi: Uniform
    vec2=None
    alpha=np.pi/3
    f_pdf = lambda r,phi: rayleigh_pdf(r, s=rsig) * uniform_pdf(phi, 0, 2*np.pi) # NOTE: independent
    vec = StochasticVector(dimension=2, name='Input - X', nDiscr=nDiscr, domain=[[0,10],[0,2*np.pi]])
    vec.set_pdf_f(f_pdf)
    ax = vec.plot_pdfs()
    #ax2 = vec.plot_pdf()

    # --- Jacobian of the transformation
    # Y = [R cos(alpha+Phi),  R sin (alpha+Phi) ]
    def gradg(x):
        r = max(x[0], 1e-5) # Avoid singularity
        return np.array([ [np.cos(alpha+x[1]), - r * np.sin(alpha+x[1]) ],  [np.sin(alpha+x[1]),  r*np.cos(alpha+x[1]) ] ]) 
    def ginv(y):
        r= np.sqrt(y[0]**2+y[1]**2)
        phi = np.arctan2(y[1],y[0])
        phi = np.mod(phi, 2*np.pi)
        if np.isnan(phi):
            phi=0
        return np.array([r,phi])

    vec2 = vec.new_from_bijection(ginv=ginv, gradg=gradg, name='Output - Y', nDiscr=nDiscr, xmin=-10, xmax=10)
    vec2.plot_pdfs(ax=ax, sty=[':', '.'])

    x=np.linspace(-10,10,100)
    ax.plot(x, gaussian_pdf(x, mu=0, sig=rsig), 'k--', label='N(0,s)')
    ax.legend()
    #ax2 = vec2.plot_pdf()
                           

    return vec, vec2



if __name__ == '__main__':
    vec, vec2 = main(nDiscr=20)
    #vec, vec2 = main2(nDiscr=21)
    print(vec)
    print(vec2)
    plt.show()

if __name__ == '__test__':
    vec, vec2 = main(nDiscr=20)
    np.testing.assert_almost_equal(vec.mean, mus, 3)
    np.testing.assert_almost_equal(vec.var,  C , 2)
    np.testing.assert_almost_equal(vec2.mean, [0.,0.], 3)
    np.testing.assert_almost_equal(vec2.var, np.eye(2), 3)

    vec, vec2 = main2(nDiscr=21)
    np.testing.assert_almost_equal(vec.mean, [2.507, np.pi], 3)
    np.testing.assert_almost_equal(vec.var,  np.array([[1.717, 0.],[0, 3.290]]) , 3)
    np.testing.assert_almost_equal(vec2.mean, [0.,0.], 3)
    np.testing.assert_almost_equal(vec2.var, np.eye(2)*rsig**2, 3)


if __name__=="__export__":
    main(nDisc=100)
    from welib.tools.repo import export_figs_callback
    export_figs_callback(__file__)
