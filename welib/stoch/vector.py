import numpy as np
import matplotlib.pyplot as plt


from scipy.integrate import nquad
from welib.stoch.variable import StochasticVariable

class StochasticVector():
    def __init__(self, dimension, domain=None, name='', nDiscr=100, xmin=-30, xmax=30):
        """ 
        - domain: d x 2 array of [min(xi), max(xi)] for each xi in x
        """
        self.name = name
        self.d = dimension
        # Definitions based on an array of the PDF
        self._pdf = None    # array
        self._x_pdf = None  # array
        #self._cdf = None    # array
        #self._x_cdf = None  # array
        # Definitions using functions
        self._f_pdf = None  # function
        self._f_cdf = None  # function

        self.variables = None

        # --- Storage to avoid recomputations 
        self.reset_misc()

        # --- Domain and discretization
        self.nDiscr=nDiscr
        if domain is None:
            self.domain = np.zeros((self.d, 2))
            for idim in range(self.d):
                self.domain[idim, :] = [xmin, xmax]
        else:
            self.domain = np.asarray(domain)
        assert(self.domain.shape == (self.d, 2))

    def reset_misc(self):
        self._means=None
        self._covar=None

    def integral(self):
        if self._f_pdf is not None:
            integral, _ = nquad( self._f_pdf , self.domain)
            return integral
        else:
            return np.nan

    def set_pdf(self, x, pdf):
        self.reset_misc()
        pass

    def set_pdf_f(self, f_pdf):
        self._f_pdf = f_pdf
        self.reset_misc()
        # --- Set individual variables
        # Compute marginal distributions
        X = self.x_default
        self.variables = []
        for idim in range(self.d):
            xi = X[idim]
            pdfi = np.zeros_like(xi)
            subdomain = [self.domain[j] for j in range(self.d) if j!=idim]
            for i,xi0 in enumerate(xi):
                def f_subdomain(*x):
                    x_full = np.insert(x, idim, xi0)
                    return self._f_pdf(*x_full)
                pdfi[i], error = nquad( f_subdomain , subdomain)

            var = StochasticVariable(name=self.name+str(idim))
            var.set_pdf(xi, pdfi)
            self.variables.append(var)
        
    @property
    def mean(self):
        r""" 

        mu_j  = \int_D1 .. \int_Dn xj f_X1..Xn  dx  
              = \int_Dj xj f_Xj dxj

        """
        if self._means is not None:
            return self._means
        if self._f_pdf is not None:
            means = np.zeros(self.d)
            if self.variables is not None:
                # We know the marginal distributions, we integrate them to get the mean..
                for i in range(self.d):
                    means[i] = self.variables[i].mean
            else:
                # Integration of multidimension function
                for i in range(self.d):
                    xf = lambda *x : x[i] * self._f_pdf(*x)
                    means[i], error = nquad( xf , self.domain)
            self._means = means.copy()
            return means
        else:
            raise NotImplementedError()

    @property
    def var(self):
        r""" 

        kappa_jk  = \int_D1 .. \int_Dn (xj-muj)(xk-muk) f_X1..Xn  dx 
                  = \int_Dj \int_Dk (xj-muj)(xk-muk) f_XjXk  dxj
        """
        if self._covar is not None:
            return self._covar
        if self._f_pdf is not None:
            # Integration of multidimension function
            means = self.mean
            covar = np.zeros((self.d,self.d))
            for j in range(self.d):
                for k in range(self.d):
                    xxf = lambda *x : (x[j]-means[j])*(x[k]-means[k]) * self._f_pdf(*x)
                    covar[j, k], _ = nquad( xxf , self.domain)
            self._covar = covar.copy()
            return covar
        else:
            raise NotImplementedError()

    @property
    def corrcoeff(self):
        corrcoeff = np.zeros((self.d, self.d))
        if self.variables is not None:
            covar = self.var
            sigmas = [np.sqrt(self.variables[j].var) for j in range(self.d)]
            for j in range(self.d):
                for k in range(self.d):
                    corrcoeff[j,k] = covar[j,k]/(sigmas[j]*sigmas[k])
        else:
            raise NotImplementedError()
        return corrcoeff


    @property
    def x_default(self):
        x = []
        for j in range(self.d):
            x.append( np.linspace(self.domain[j,0],self.domain[j,1],self.nDiscr+j))
        return x

    def plot_pdfs(self, x=None, ax=None, sty='-', **kwargs):
        """ Plot marginal distributions """
        if x is None:
            x = self.x_default
        if ax is None:
            fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
            fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
        if not isinstance(sty, list):
            sty = [sty] * self.d
        for i, st in enumerate(sty):
            self.variables[i].plot_pdf(x[i], ax=ax, sty=st)
        ax.legend()
        return ax


    def meshgrid_pdf(self, x=None, i1=0, i2=1):
        if x is None:
            x = self.x_default
        x1 = x[i1]
        x2 = x[i2]
        X1, X2 = np.meshgrid(x1,x2)
        F = np.zeros(X1.shape)
        x0 = np.zeros(self.d)
        for ii1,xx1 in enumerate(x1):
            for ii2,xx2 in enumerate(x2):
                x0[i1] = xx1
                x0[i2] = xx2
                F[ii2,ii1] = self._f_pdf(*x0)
        return X1, X2, F

    def plot_pdf(self, x=None, ax=None, i1=0, i2=1, **kwargs):
        if x is None:
            x = self.x_default

        # Get PDF on meshgrid
        X1, X2, F = self.meshgrid_pdf(x, i1=i1, i2=i2)
        if ax is None:
            fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
            fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
        ax.contourf(X1,X2,F)
        ax.set_xlabel(self.name+str(i1))
        ax.set_ylabel(self.name+str(i2))
        return ax, X1, X2, F



    def new_from_bijection(self, gradg=None, ginv=None, name='', xmin=-30, xmax=30, nDiscr=100):
        vec2 = StochasticVector(dimension=self.d, name=name, xmin=xmin, xmax=xmax, nDiscr=nDiscr)
        def f_pdf(*y):
            y = np.asarray(y).flatten()
            x = ginv(y).ravel()
            fX = self._f_pdf(*x)
            Jac = gradg(x)
            detJ = np.linalg.det(Jac)
            fY =  fX * 1/abs(detJ)
            if np.isnan(fY):
                print('fY is NaN')
                fY=0
            return fY
        vec2.set_pdf_f(f_pdf)
        return vec2


    def standardize_transform(self, C=None, mu=None):
        """ 
        Returns a transformation Y = g(X) = a + B.X
        Such that Y is is made of independent strandardized normal stochastic variables N(0,1)
        (zero mean, standard deviation of 1)
        NOTE: 
           X = N(mu, C)  ---- Y = a + BX---->   Y = N( a + B.mu,   BCB^T )

        So we want  BCB^T = I and a + B.mu=0

        This is achieved by looking at the eigenvalues of C the covariance matrix of X:

            C V = V L 

          - V is the matrix of eigenvectors
          - L is the diagonal matrix of eigenvalues [l1, l2, .., ln]

        Then we set:
          - B = L^{-1/2} V^T    where L^{-1/2} is diagonal matrix [l1^{-1/2}, ..., ln^{-1/2}]
          - a = - B mu        (mu the means of X)

        which leads to: BCB^T = I and a+B.mu=0.
        """
        if C is None:
            C  = self.var
        if mu is None:
            mu = self.mean
        # Eigenvalue problem
        lambdas,V = np.linalg.eig(C)
        Lambda    = np.diag(lambdas)
        LambdaM12 = np.diag(1/np.sqrt(lambdas))
        B = (LambdaM12).dot(V.T)
        a = - B.dot(mu)

        Binv = np.linalg.inv(B)
        detB = np.linalg.det(B)
        g     = lambda x : a + B.dot(x)
        gradg = lambda x : B
        ginv  = lambda y : Binv.dot(y-a)

        return a, B, g, gradg, ginv
        #vec2 = self.new_from_bijection(ginv=ginv, gradg=gradg, name='Y', nDiscr=100, xmin=-10, xmax=10)


    def __repr__(self):
        s='<{} object>:\n'.format(type(self).__name__)
        s+='- name    : {}\n'.format(self.name)
        s+='- d       : {}\n'.format(self.d)
        sd=[]
        for i in range(self.d):
            sd.append('[{} ; {}]'.format(*self.domain[i]))
        s+='- domain  : {}\n'.format(' x '.join(sd) )
        try:
            mean      = self.mean
            var       = self.var
            corrcoeff = self.corrcoeff
            check     = self.integral()
        except:
            check     = np.nan
            mean      = np.nan
            var       = np.nan
            corrcoeff = np.nan
        s+='* integral: {} \n'.format(np.around(check,4))
        s+='* mean    :  {} \n'.format(np.around(mean,4))
        s+='* var     : \n'
        s+='{}'.format(np.around(var,4))
        # s+='corrcoeff\n',np.around(self.corrcoeff,2)
        return s




if __name__ == '__main__':
    from welib.essentials import *
    from welib.stoch.distribution import gaussian2d_pdf
    from welib.stoch.distribution import gaussian2d_covar

    sigmas = [5,6]
    mus    = [0,3]
    rho=0.1
    C = gaussian2d_covar(sigmas, rho)
    f_pdf = lambda x1,x2: gaussian2d_pdf(x1, x2, mus, sigmas, rho=rho)

    vec = StochasticVector(dimension=2, name='X')
    vec.set_pdf_f(f_pdf)
    print(vec)
#     vec.integral()
# 
#     print('means'  ,np.around(vec.mean,2), mus)
#     print('covar\n',np.around(vec.var,2), '\n', C)
#     print('corrcoeff\n',np.around(vec.corrcoeff,2), '\n', rho)
#     ax = vec.plot_pdfs()
# 
# 
# 
# 
#     plt.show()
    pass
