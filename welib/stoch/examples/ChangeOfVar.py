""" 


See examples 2-2, 2-3, 2-4
"""

from welib.stoch.variable import *
from welib.stoch.distribution import *
import matplotlib.pyplot as plt

def main():

    # --- g=exp(x) lognormal distribution
    mu = 1
    sig = 1
    f = lambda x: gaussian_pdf(x, mu, sig)
    x   = np.linspace(-30,30,300)
    y = np.linspace(1e-3, 60, 1000)

    # Analytical
    varA = StochasticVariable()
    varA.set_pdf_f(f)
    #ax1 = varA.plot_pdf(x, sty='k')
    gp = np.exp
    ginv = np.log
    varA2 = varA.new_from_bijection(gp=gp, ginv=ginv)
    ax2 = varA2.plot_pdf(y, sty='k')
    varA2.check_pdf(y)


    var = StochasticVariable()
    var.set_pdf(x, f(x))
    #var.plot_pdf(ax=ax1, sty='+')
    g = np.exp(x)
    #var2 = var.new_from_bijection(g, yvalues='linspace+input')
    var2 = var.new_from_bijection(g, yvalues='logspace')
    ax =var2.plot_pdf(sty='+', ax=ax2)
    var2.check_pdf(y)
    fy = 1/(np.sqrt(2*np.pi) * sig) * 1/y * np.exp( - ( np.log(y) - mu )**2 / (2*sig**2))
    ax.plot(y, fy, '--')
    ax.set_xlim([0, 2])
    ax.set_ylim([0, 1])



    # --- g=x**3 --- NOTE Integral is not 1 somehow
#     mu = 0
#     sig = 1
#     f = lambda x: gaussian_pdf(x, mu, sig)
#     x   = np.linspace(-10,10,300)
#     y   = x**3
# 
#     varA = StochasticVariable()
#     varA.set_pdf_f(f)
#     varA.check_pdf()
#     ax1 = varA.plot_pdf(x, sty='k')
#     gp = lambda x: 3*x**2
#     ginv = np.cbrt #ginv = lambda x: x**(1/3) # Approximate and give complex numbers 
#     varA2 = varA.new_from_bijection(gp=gp, ginv=ginv)
#     ax2 = varA2.plot_pdf(y, sty='k')
#     varA2.check_pdf()
# 
# 
#     var = StochasticVariable()
#     var.set_pdf(x, f(x))
#     var.plot_pdf(ax=ax1, sty='+')
#     var.check_pdf()
#     g = x**3
#     var2 = var.new_from_bijection(g, yvalues='linspace+input')
#     ax2=var2.plot_pdf(y, sty='+', ax=ax2)
#     ax2.set_ylim([0,.1])
#     ax2.set_xlim([-10,10])
#     var2.check_pdf()
    



    # --- g=x**2 need splitting - Integral is not 1
    sig = 1
    mu = 0.5
    x   = np.linspace(-10,10,250)
    pdf = gaussian_pdf(x, mu, sig)
    y = np.linspace(1e-3, 60, 300)

    var = StochasticVariable()
    var.set_pdf(x, pdf)

    g = x**2
    var2 = var.new_from_bijection(g, yvalues='linspace', split=0)
    ax =var2.plot_pdf(sty='o')
    fy1= 1/(np.sqrt(8*np.pi) * sig) * 1/np.sqrt(y) * np.exp( - ( -np.sqrt(y) - mu )**2 / (2*sig**2))
    fy2= 1/(np.sqrt(8*np.pi) * sig) * 1/np.sqrt(y) * np.exp( - (  np.sqrt(y) - mu )**2 / (2*sig**2))
    fy= fy1+fy2
    ax.plot(y, fy , 'k-')
    ax.plot(y, fy1, 'k.')
    ax.plot(y, fy2, 'k+')
    ax.set_xlim([0, 60])
    ax.set_ylim([0, 1])# 
    var.check_pdf(y)


    # --- g=Fx(x) leads to the uniform distribution
    # 
    y = np.linspace(0, 1, 100)
    g = var.cdf(x)
    var2 = var.new_from_bijection(g, yvalues='linspace')
    ax =var2.plot_pdf(sty='o')
    ax.plot(y,  uniform_pdf(y, 0, 1), 'k--')



if __name__ == '__main__':
    main()
    plt.show()
