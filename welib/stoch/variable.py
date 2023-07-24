import numpy as np
import matplotlib.pyplot as plt

class StochasticVariable:
    def __init__(self):
        pass

    def plot_pdf(self, x=None, ax=None, sty='-', label=None, **kwargs):
        if x is None:
            x=np.linspace(-30,30,100)
        pdf = self.pdf(x)
        if ax is None:
            fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
            fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
        ax.plot(x, pdf,   sty, label=label, **kwargs)
        #ax.set_xlabel('')
        #ax.set_ylabel('')
        return ax


    def check_pdf(self, x=None):
        if x is None:
            x=np.linspace(-30,30,100)
        integral = np.trapz(self.pdf(x), x)
        print('>>> Integral', integral)



class AnalyticalStochasticVariable(StochasticVariable):
    """ Stochastic variable defined by an analytical function"""
    def __init__(self):
        pass

    def set_pdf(self, f_pdf):
        self._f_pdf = f_pdf

    def pdf(self, x):
        return self._f_pdf(x)

    def new_from_bijection(self, gp, ginv):
        """ 
        create variable y = g(x) 

        fY(y) dy = fX(x) dx

        fY(y) = fX(x)        1/|g'(x)|
              = fX(g-1(y))   1/|g'(g-1(y))|

        """
        f_pdf = lambda y: self._f_pdf( ginv(y)  ) / np.abs(gp (ginv(y) ) )
        var = AnalyticalStochasticVariable()
        var.set_pdf(f_pdf)
        return var


def interpsafe(x, xp, fp):
    xp=xp.copy()
    fp=fp.copy()
    if xp[0]>xp[-1]:
        xp=xp[::-1]
        fp=fp[::-1]
    return np.interp(x, xp, fp)

class NumericalStochasticVariable(StochasticVariable):
    def __init__(self):
        pass

    def set_pdf(self, x, pdf):
        self._x_pdf = x
        self._pdf   = pdf

    def pdf(self, x):
        return interpsafe(x, self._x_pdf, self._pdf)

    def plot_pdf(self, x=None, **kwargs):
        if x is None:
            x = self._x_pdf
        return StochasticVariable.plot_pdf(self, x, **kwargs)

    def new_from_bijection(self, g, x=None, gp=None, ginv=None, yvalues='input', split=None):
        """ 
        create variable y = g(x) 

        fY(y) dy = fX(x) dx

        fY(y) = fX(x)        1/|g'(x)|
              = fX(g-1(y))   1/|g'(g-1(y))|

        """
        if x is None:
            x0 = self._x_pdf
            fx = self._pdf
        else:
            x0 = x
            fx = self.pdf(x0)
        if split is not None:
            b1 = x0>=split
            b2 = x0<split
            x1 = x0[b1]
            x2 = x0[b2]
            g1 = g[b1]
            g2 = g[b2]
#             x2=x2[::-1]
#             g2=g2[::-1]

#             fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
#             fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
#             ax.plot(x1,g1, label='g1')
#             ax.plot(x2,g2, label='g2')
#             ax.set_xlabel('x')
#             ax.set_ylabel('y')
#             ax.legend()

            var1 = self.new_from_bijection(g1, x=x1, yvalues=yvalues)
            var2 = self.new_from_bijection(g2, x=x2, yvalues=yvalues)

            y = np.concatenate((var1._x_pdf, var2._x_pdf))
            y1 = np.min(y)
            yn = np.max(y)
            yi = np.linspace(y1, yn, len(x0))
            fy1 = np.interp(yi, var1._x_pdf, var1._pdf)
            fy2 = np.interp(yi, var2._x_pdf, var2._pdf)
            fyi = fy1+fy2

#             ax = var1.plot_pdf()
#             ax.plot(y, fy01, 'k--')
#             ax.set_title('pdf 1')
# 
#             ax = var2.plot_pdf()
#             ax.plot(y, fy02, 'k--')
#             ax.set_title('pdf 2')
# 
#             fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
#             fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
#             ax.plot(yi, fy1    , label='fy1')
#             ax.plot(yi, fy2    , label='fy2')
#             ax.set_xlabel('')
#             ax.set_ylabel('')
#             ax.set_xlim([0,10])
#             ax.legend()


            var = NumericalStochasticVariable()
            var.set_pdf(yi, fyi)
           
            return var

        y0 = g
        gp0 = np.gradient(g, x0)
        y1 = np.min(g)
        yn = np.max(g)

        if yvalues=='input':
            yi = y0.copy()
            pdf_yi = self._pdf * 1/np.abs(gp0)
        elif yvalues=='linspace':
            yi = np.linspace(y1, yn, len(x0))
        elif yvalues=='logspace':
            yi = np.exp(np.linspace(np.log(y1), np.log(yn), len(x0)))
        elif yvalues=='linspace+input':
            yi = np.linspace(y1, yn, len(x0))
            yi = np.unique(np.sort(np.concatenate((yi, y0))))
        else:
            raise NotImplementedError()
        xi   = interpsafe(yi, y0, x0)  # xi = g^(-1)(yi)
        fxi  = self.pdf(xi)            # fX(xi) = fX(g-1(y)) 
        gpxi = interpsafe(xi, x0, gp0) # g'(xi) = g'(g-1(y)) 
        pdf_yi = fxi * 1/np.abs(gpxi)  # fy= fX(xi)/|g'(xi)|

#         fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
#         fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
#         ax.plot( x0, y0    , label='y=g(x)')
#         ax.plot( xi, yi ,'+'   , label='interp')
#         ax.set_xlabel('x')
#         ax.set_ylabel('y')
#         ax.legend()
#         fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
#         fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
#         ax.plot( x, 2*x       ,'k-'  , label='analytical')
#         ax.plot( x, gp0       ,'--'  , label='num')
#         ax.plot( xi, gpxi     ,'+'   , label='interp')
#         ax.set_xlabel('x')
#         ax.set_ylabel('gp(x)')
#         ax.legend()
# 
#         fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
#         fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
#         ax.plot( x, self._pdf)
#         ax.plot( xi, fxi      ,'+'   , label='interp')
#         ax.set_xlabel('x')
#         ax.set_ylabel('gp(x)')
#         ax.legend()

        var = NumericalStochasticVariable()
        var.set_pdf(yi, pdf_yi)



        return var




if __name__ == '__main__':
    from welib.stoch.distributions import gaussian_pdf



# 
#     # --- g=exp(x) lognormal distribution
#     mu = 1
#     sig = 1
#     f = lambda x: 1/( sig * np.sqrt(2*np.pi) ) * np.exp(-0.5 * (x-mu)**2/sig**2)
#     x   = np.linspace(-30,30,300)
#     y = np.linspace(1e-3, 60, 1000)
# 
#     varA = AnalyticalStochasticVariable()
#     varA.set_pdf(f)
#     ax1 = varA.plot_pdf(x, sty='k')
#     gp = np.exp
#     ginv = np.log
#     varA2 = varA.new_from_bijection(gp=gp, ginv=ginv)
#     ax2 = varA2.plot_pdf(y, sty='k')
#     varA2.check_pdf(y)
# 
# 
#     var = NumericalStochasticVariable()
#     var.set_pdf(x, f(x))
#     var.plot_pdf(ax=ax1, sty='+')
#     g = np.exp(x)
#     #var2 = var.new_from_bijection(g, yvalues='linspace+input')
#     var2 = var.new_from_bijection(g, yvalues='logspace')
#     ax =var2.plot_pdf(sty='+', ax=ax2)
#     var2.check_pdf(y)
#     fy = 1/(np.sqrt(2*np.pi) * sig) * 1/y * np.exp( - ( np.log(y) - mu )**2 / (2*sig**2))
#     ax.plot(y, fy, '--')
#     ax.set_xlim([0, 2])
#     ax.set_ylim([0, 1])



    # --- g=x**3 --- NOTE Integral is not 1 somehow
#     mu = 0
#     sig = 1
#     f = lambda x: 1/( sig * np.sqrt(2*np.pi) ) * np.exp(-0.5 * (x-mu)**2/sig**2)
#     x   = np.linspace(-10,10,300)
#     y   = x**3
# 
#     varA = AnalyticalStochasticVariable()
#     varA.set_pdf(f)
#     varA.check_pdf()
#     ax1 = varA.plot_pdf(x, sty='k')
#     gp = lambda x: 3*x**2
#     ginv = np.cbrt #ginv = lambda x: x**(1/3) # Approximate and give complex numbers 
#     varA2 = varA.new_from_bijection(gp=gp, ginv=ginv)
#     ax2 = varA2.plot_pdf(y, sty='k')
#     varA2.check_pdf()
# 
# 
#     var = NumericalStochasticVariable()
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
    pdf = gaussian_pdf(x, mu, sig**2)
    y = np.linspace(1e-3, 60, 300)

    var = NumericalStochasticVariable()
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



    plt.show()
    
    

    pass
