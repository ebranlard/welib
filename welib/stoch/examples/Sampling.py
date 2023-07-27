""" 

See example 2-4
"""
from welib.stoch.variable import *
from welib.stoch.distribution import *
from welib.essentials import *
from numpy.random import seed
seed(12)

def main():

    mu = 1
    sig = 1
    f = lambda x: gaussian_pdf(x, mu, sig)
    x = np.linspace(-30,30,300)

    var = StochasticVariable()
    var.set_pdf(x, f(x))
    s = var.sample(100000)

    var2 = StochasticVariable()
    #var2.from_data(s, pdfmethod='histogram')
    var2.from_data(s, pdfmethod='kde')

    ax=var.plot_pdf(x, sty='k-')
    ax=var2.plot_pdf(sty='o', label='Sampled', ax=ax)
    ax.set_xlim([-3,5])
    ax.legend()

    xmin=np.min(s)
    xmax=np.max(s)
    print('min max', xmin, xmax)
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

