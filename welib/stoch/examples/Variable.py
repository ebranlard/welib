import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# Local 
from welib.stoch.distribution import *
from welib.stoch.variable import *


def main():

    X = StochasticVariable(name='X - Rayleigh', domain=[0,20], nDiscr=100)
    X.set_pdf_f(lambda x: rayleigh_pdf(x, s=2))

    Y = StochasticVariable(name='Y - Weibull', domain=[0,20])
    Y.set_pdf_f(lambda x: weibull_pdf(x, A=2, k=3))

    return X,Y



if __name__ == '__main__':
    X, Y = main()
    print(X)
    print(Y)
    #X.plot_pdf(plot_mean=True)
    plt.show()

if __name__ == '__test__':
    X, Y = main()
    m  = rayleigh_moments(s=2)
    np.testing.assert_almost_equal(X.mean, m['mean']    , 5)
    np.testing.assert_almost_equal(X.var,  m['variance'], 5)

    m  = weibull_moments(A=2, k=3)
    np.testing.assert_almost_equal(Y.mean,  m['mean']    , 5)
    np.testing.assert_almost_equal(Y.var,   m['variance'], 5)

if __name__=="__export__":
    pass
    #main()
    #from welib.tools.repo import export_figs_callback
    #export_figs_callback(__file__)

