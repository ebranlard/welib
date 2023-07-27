import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import nquad
from scipy.integrate import cumtrapz



def interpsafe(x, xp, fp):
    xp=xp.copy()
    fp=fp.copy()
    if xp[0]>xp[-1]:
        xp=xp[::-1]
        fp=fp[::-1]
    return np.interp(x, xp, fp)

# --------------------------------------------------------------------------------}
# --- Numerical variable
# --------------------------------------------------------------------------------{
class StochasticVariable():
    def __init__(self, name='', domain=None, nDiscr=100):
        self.name = name 
        # Definitions based on an array of the PDF
        self._pdf = None    # array
        self._x_pdf = None  # array
        self._cdf = None    # array
        self._x_cdf = None  # array
        # Definitions using functions
        self._f_pdf = None  # function
        self._f_cdf = None  # function
        # Definitions using stats
        self._data = None

        # --- Storage to avoid recomputations 
        self.reset_misc()

        # --- Domain and discretization
        self.nDiscr=nDiscr
        if domain is None:
            self.domain = [-30, 30]
        else:
            self.domain = np.asarray(domain)
        assert(len(self.domain) == 2)

    def reset_misc(self):
        self._mean=None
        self._var=None

    def from_data(self, data, pdfmethod='kde', nBins=50):
        """ stats """
        from welib.tools.stats import pdf_histogram, pdf_gaussian_kde
        self.reset_misc()
        if pdfmethod=='histogram':
            x_pdf, pdf = pdf_histogram(data, nBins=nBins, norm=True, count=False)
        elif pdfmethod=='kde':
            x_pdf, pdf = pdf_gaussian_kde(data, bw='scott', nOut=nBins, cut=3, clip=(-np.inf,np.inf))
        else:
            raise NotImplementedError()
        # Set PDF
        #print('>>> x_pdf', x_pdf[0], x_pdf[-1])
        #print('>>> pdf nan', np.sum(np.isnan(pdf)))
        self.set_pdf(x_pdf, pdf)

        self._data = data

    def integral(self, x=None):
        """ Compute integral of PDF (should return 1)  """
        if x is not None:
            # We compute numerical integral
            integral = np.trapz(self.pdf(x), x)
        elif self._f_pdf is not None:
            # We have a function for the pdf, we use more accurate nquad evaluation
            integral, _ = nquad( self._f_pdf , [self.domain])
        else:
            # Numerical integral
            if x is None:
                x = self.x_default
            integral = np.trapz(self.pdf(x), x)
        return integral

    def sample(self, n, method='inverse-cdf'):
        # NOTE: we will not get much extremes out of this
        #      It depends on the quality of the cdf, and the number of points used for the inverse
        # --- Inverse method, use uniform distribution
        u = np.random.uniform(low=0.0, high=1.0, size=n)
        s = self.cdf_inv(u)
        return s

    def set_pdf(self, x, pdf):
        self.reset_misc()
        assert(len(x)==len(pdf))
        self._x_pdf = x
        self._pdf   = pdf
        #
        cdf = np.concatenate(( [0], cumtrapz(pdf, x) ) )
        self.set_cdf(x, cdf) # TODO sort out

    def set_pdf_f(self, f_pdf):
        self.reset_misc()
        self._f_pdf = f_pdf

    def set_cdf(self, x, cdf):
        assert(len(x)==len(cdf))
        self._x_cdf = x
        self._cdf   = cdf

    def set_cdf_f(self, f_cdf):
        self._f_cdf = f_cdf

    def pdf(self, x):
        if self._f_pdf is None:
            return interpsafe(x, self._x_pdf, self._pdf)
        else:
            return self._f_pdf(x)

    def cdf(self, x):
        # NOTE: we should do something to ensure 0 and 1 are at +/- inf for a Gaussian for instance 
        if self._f_cdf is None:
            return interpsafe(x, self._x_cdf, self._cdf)
        else:
            return self._f_cdf(x)

    def cdf_inv(self, F):
        # NOTE: we should do something to ensure 0 and 1 are at +/- inf for a Gaussian for instance 
        x0 = np.linspace(self._x_cdf[0], self._x_cdf[-1], 10000)
        cdf0 = self.cdf(x0)
        return interpsafe(F, cdf0, x0)


    @property
    def mean(self):
        if self._mean is not None:
            return self._mean # already computed
        if self._f_pdf is not None:
            # Use nquad for more accurate calculation
            xf = lambda x : x * self._f_pdf(x)
            self._mean, _ = nquad( xf , [self.domain])
        elif self._data is None:
            # Use numerical integration
            x = self.x_default
            self._mean = np.trapz(x * self.pdf(x) , x)
        else:
            # Use the raw sampled data
            self._mean =  np.mean(self._data)
        return self._mean

    @property
    def var(self):
        if self._var is not None:
            return self._var # already computed

        mu = self.mean
        if self._f_pdf is not None:
            # Use nquad for more accurate calculation
            xxf = lambda x : (x-mu)*(x-mu) * self._f_pdf(x)
            self._var, _ = nquad( xxf , [self.domain])
        elif self._data is None:
            # Use numerical integration
            x = self.x_default
            self._var = np.trapz((x-mu)**2 * self.pdf(x) , x)
        else:
            self._var = np.var(self._data)
        return self._var

    @property
    def corrcoeff(self):
        return 1 # self.var/(sigmas[j]*sigmas[k])

    @property
    def x_default(self):
        if self._x_pdf is None:
            return np.linspace(self.domain[0], self.domain[1], self.nDiscr)
        else:
            return self._x_pdf

    def plot_pdf(self, x=None, ax=None, sty='-', label=None, plot_mean=False, **kwargs):
        if x is None:
            x = self.x_default
        if label is None:
            label = self.name

        pdf = self.pdf(x)
        if ax is None:
            fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
            fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
        ax.plot(x, pdf,   sty, label=label, **kwargs)
        if plot_mean:
            ax.axvline(self.mean, ls=':', c='k')
        #ax.set_xlabel('')
        ax.set_ylabel('Probability density function')
        return ax

    def plot_cdf(self, x=None, ax=None, sty='-', label=None, **kwargs):
        if x is None:
            x = self.x_default
        cdf = self.cdf(x)
        if ax is None:
            fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
            fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
        ax.plot(x, cdf,   sty, label=label, **kwargs)
        #ax.set_xlabel('')
        #ax.set_ylabel('')
        ax.set_ylabel('Cumulative density function')
        return ax

    # --------------------------------------------------------------------------------
    # --- Other variable
    # --------------------------------------------------------------------------------
    def cov(self, var2, f_XY):
        """
        Compute Cov(X,Y) = \int\int (x-mu_X)(y-mu_Y) f_XY(x,y) dx dy
        NOTE: using vector is better for that as vectors contain f_XY
        """
        muX, muY = self.mean, var2.mean
        # Integration of multidimension function
        xxf = lambda x,y : (x-muX)*(y-muY) * f_XY(x,y)
        cov = nquad( xxf , [self.domain, var2.domain])
        return cov

    def new_from_bijection(self, g=None, x=None, gp=None, ginv=None, yvalues='input', split=None):
        """ 
        create variable y = g(x) 

        fY(y) dy = fX(x) dx

        fY(y) = fX(x)        1/|g'(x)|
              = fX(g-1(y))   1/|g'(g-1(y))|

        """

        # --- Analytical way
        if gp is not None and ginv is not None and self._f_pdf is not None:
            #f_pdf = lambda y: self._f_pdf( ginv(y)  ) / np.abs(gp (ginv(y) ) )
            def f_pdf(y):
                 y = np.asarray(y)
                 scale = np.abs(gp (ginv(y)) )
                 inv_scale = np.zeros_like(y)
                 inv_scale[scale>0] = 1/scale[scale>0]
                 return self._f_pdf( ginv(y)  ) * inv_scale
            var = StochasticVariable()
            var.set_pdf_f(f_pdf)
            return var

        # --- Numerical way


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


            var = StochasticVariable()
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

        var = StochasticVariable()
        var.set_pdf(yi, pdf_yi)



        return var

    def __repr__(self):
        s='<{} object>:\n'.format(type(self).__name__)
        s+='- name    : {}\n'.format(self.name)
        s+='- domain  : [{} ; {} ]\n'.format(self.domain[0], self.domain[1])
        try:
            mean      = self.mean
            var       = self.var
            check     = self.integral()
        except:
            check     = np.nan
            mean      = np.nan
            var       = np.nan
        s+='* integral: {}\n'.format(np.around(check,5))
        s+='* mean    : {}\n'.format(np.around(mean,5))
        s+='* var     : {}\n'.format(np.around(var, 5))
        return s




if __name__ == '__main__':
    from welib.stoch.distribution import *
    from welib.essentials import *
    from numpy.random import seed
    seed(12)

    mu = 1
    sig = 2
    f = lambda x: gaussian_pdf(x, mu, sig)
    x = np.linspace(-30,30,300)
    y = np.linspace(0,1,300)

    var = StochasticVariable()
    var.set_pdf(x, f(x))
    print('Mean', var.mean)
    print('Var' , var.var)


    plt.show()

    pass
