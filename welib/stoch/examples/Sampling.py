""" 

Examples of usage of StochasticVariable and drawing samples from it.


 - Using variables with predefined distributions
 - Using variable defined based on a sampled pdf,  (xi, pdf(xi) )
 - Using variable defined based on sampled data xi, from which the pdf will be computed


Samples are drawn using `var.sample(n, method='inverse-cdf')`.

The sampling method available for all variable is "inverse-cdf" which uses
the cumulative density function for generate samples by drawing samples (u) in the
uniform distribution in [0,1], and then get x using the inverse of the CDF:
    x = CFD^(-1)(u), with u~U(0,1)


Predefined variables also have the "numpy" method defined.




See example 2-4
"""
from welib.stoch.variable import *
from welib.stoch.distribution import *
from welib.stoch.utils import plot_pdf
from welib.essentials import *
from numpy.random import seed
seed(12)

def main(Itest=[0,1]):
    n     = 100000
    mu    = 2
    sigma = 3

    if 0 in Itest:
        # --- High level
        # Using predefined variables with predefined distribution
        X = NormalVariable(mu=mu, sigma=sigma)
        #     X = RayleighVariable(s=sigma)
        #     X = WeibullVariable(A=mu, k=sigma)
        #print(X)

        x1 = X.sample(n=n, method='numpy')
        x2 = X.sample(n=n, method='inverse-cdf')
        x3 = X.sample(n=n, method='boxmuller')

        ax = X.plot_pdf(label='Theory', sty='k-')
        ax = plot_pdf(x1, sty='o', label='Sampled - Numpy'      , ax=ax, method='gaussian_kde', n=50)
        ax = plot_pdf(x2, sty='d', label='Sampled - Inverse-CDF', ax=ax, method='gaussian_kde', n=50)
        ax = plot_pdf(x3, sty='.', label='Sampled - Box-Muller ', ax=ax, method='gaussian_kde', n=50)
        ax.legend()


    if 1 in Itest:
        # --- Low level
        # Creating a stochastic variable based on pdf data points
        f = lambda x: gaussian_pdf(x, mu, sigma)
        x = np.linspace(-5*sigma,5*sigma,300) + mu

        var = StochasticVariable()
        var.set_pdf(x, f(x))  # Note could also use function directly by using var.set_pdf_f(f)
        s = var.sample(n)

        # Creating a stochastic variable based on sampled data
        var2 = StochasticVariable()
        var2.from_data(s, pdfmethod='gaussian_kde') # pdfmethod in ['histogram','gaussian_kde']

        # Plot the pdf of both variables
        ax=var.plot_pdf(x, sty='k-')
        ax=var2.plot_pdf(sty='o', label='Sampled', ax=ax)
        ax.set_xlim([-3,5])
        ax.legend()

        xmin=np.min(s)
        xmax=np.max(s)
        ax.set_xlim([xmin,xmax])

if __name__ == '__main__':
    main()
    plt.show()

if __name__ == '__test__':
    main()
    # TODO test it
    #np.testing.assert_almost_equal(vec2.mean, [0.,0.], 3)

if __name__=="__export__":
    pass
    #main()
    #from welib.tools.repo import export_figs_callback
    #export_figs_callback(__file__)

