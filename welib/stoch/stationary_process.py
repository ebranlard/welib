import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import nquad

from welib.tools.tictoc import Timer
from welib.tools.spectral import fft_wrap
from welib.stoch.utils import autocovariance_int
from welib.stoch.utils import autocovariance_fft
from welib.stoch.utils import autospectrum_int
from welib.stoch.utils import autospectrum_fft
from welib.stoch.utils import autocorrcoeff_num
from welib.stoch.utils import sample_from_autospectrum

class SampledStochasticProcess:
    def __init__(self, params=None, generator=None, name='', omega_max=None, tau_max=None, time_max=None, nDiscr=100,
            verbose=False):
        """ 
        Describe ME!!!

        """
        self.name = name
        self.verbose = verbose
        self.params = params if params is not None else {}

        # Default domains
        time_max  = time_max   if time_max is not None else 10
        tau_max   = tau_max    if tau_max is not None else 10
        omega_max = omega_max  if omega_max is not None else 10
        self.omega_max  = omega_max
        self.tau_max    = tau_max
        self.time_max   = max(time_max,tau_max)
        self.nDiscr     = nDiscr
        if verbose:
            print('omega_max', omega_max)
            print('tau_max  ', tau_max)
            print('time_max ', time_max)

        # Sample generations using a generator
        self.generator = generator
        self._time = None

        # Samples
        self.xi   = None  # nSamples x nTime

        # if analytical values exists
        self._f_kappa_XX_th = None  # interface: kappa_XX(tau)
        self._f_S_X_th      = None  # interface: S_X(om)
        self._var_th       = None
        self._mu_th        = None
        self._S_X_th_i     = {}
        self._kappa_XX_th_i = {}

        # Sample stats
        self._tau     = None
        self._omega   = None
        self.mu_xi    = None
        self.var_xi   = None
        self.kappa_XX = None
        self.rho_XX   = None
        self.S_X      = None
        self.S_avg    = None
        self.x_pdf_xi = None
        self.pdfs_xi  = None
        self.pdf_xi   = None

    @property
    def nSamples(self):
        if self.xi is None:
            return 0
        return self.xi.shape[0]

    @property
    def time_default(self):
        return np.linspace(0, self.time_max, self.nDiscr)

    @property
    def tau_default(self):
        return np.linspace(0, self.tau_max, self.nDiscr)

    @property
    def omega_default(self):
        return np.linspace(0, self.omega_max, self.nDiscr)

    @property
    def time(self):
        if self._time is not None:
            return self._time
        else:
            return self.time_default

    @property
    def tau(self):
        if self._tau is not None:
            return self._tau
        else:
            return self.tau_default

    @property
    def omega(self):
        if self._omega is not None:
            return self._omega
        else:
            return self.omega_default

    def set_samples(self, time, xi):
        assert(xi.shape[1]==len(time))
        self._time = time
        self.xi = xi

    def generate_samples(self, nSamples=1, time=None):
        """ 
        Generate time series samples for stochastic process X(t)
            x = generator(time): [array]
        """
        if time is None:
            time = self.time_default
        self._time = time
        xi=np.zeros((nSamples, len(time)))
        for i in range(nSamples):
            xi[i,:] = self.generator(time)
        self.xi = xi
        return xi

    def generate_samples_from_autospectrum(self, nSamples=1, time=None, method='sumcos-irfft', **kwargs):
        if time is None:
            time = self.time_default
        self._time = time


        xi=np.zeros((nSamples, len(time)))
        if not self._f_S_X_th:
            raise Exception('Provide an autospectrum function')
        f_S = self._f_S_X_th
        # TODO Merge this with generator
        tMax = time[-1]
        dt = (time[-1]-time[0])/(len(time)-1)
        with Timer('Gen. samples from spectrum - {} ...'.format(method), writeBefore=True, silent=not self.verbose):
            for i in range(nSamples):
                _, xi[i,:], _, _ = sample_from_autospectrum(tMax, dt, f_S, angularFrequency=True, method=method, **kwargs) 
        self.xi = xi
        return xi



    def samples_stats(self, nTau=None, stats=None, nPDF=50, tTransient=-np.inf):
        """ 
        """
        if stats is None:
            stats=['correlation','avg_spectra','pdf']
        time = self._time
        # removing transients
        b= time>tTransient
        time = time [b]
        xi   = self.xi[:,b]

        dt=(time[-1]-time[0])/(len(time)-1)

        # --- Basic stats
        self.mu_xi  = np.mean(xi, axis = 0)
        self.var_xi = np.var(xi, axis  = 0)

        from welib.stoch.variable import StochasticVariable
        var = StochasticVariable(name='VariableAfterSamplingX')

        # --- Average spectra
        if 'avg_spectra' in stats:
            f0, S_X0, Info = fft_wrap(time, xi[0,:], output_type='psd', averaging='None', detrend=False)
            S_i = np.zeros((self.nSamples, len(f0)))
            with Timer('Samples spectra', writeBefore=True, silent=not self.verbose):
                for i in range(self.nSamples):
                    f, S_i[i,:], Info = fft_wrap(time, xi[i,:], output_type='psd', averaging='None', detrend=False)
                S_avg = np.mean(S_i, axis=0)/(2*np.pi)
            self.om_Savg    =  2*np.pi*f
            self.S_avg      =  S_avg
        
        # --- Correlation
        if 'correlation' in stats:
            with Timer('Samples correlations', writeBefore=True, silent=not self.verbose):
                if nTau is None:
                    nTau = int(self.tau_max/dt)
                    #nTau = len(time)
                kappa_XXi = np.zeros((self.nSamples, nTau))
                rho_XXi   = np.zeros((self.nSamples, nTau))
                for i in range(self.nSamples):
                    if np.mod(i,10)==0 and self.verbose:
                        print('Correlation', i, self.nSamples)
                    rho_XXi[i,:], tau = autocorrcoeff_num(xi[i,:], nMax=nTau, dt=dt, method='corrcoef')
                    kappa_XXi[i,:] = rho_XXi[i,:] * np.var(xi[i,:])
                    #, tau = autocovariance_num(xi[i,:], nMax=nTau, dt=dt, method='manual')

            kappa_XX = np.mean(kappa_XXi, axis=0)
            rho_XX = np.mean(rho_XXi, axis=0)
            self._tau     =  tau      
            self.kappa_XX =  kappa_XX 
            self.rho_XX   =  rho_XX   

            # --- Autospectrum from kappa num..Not the best
            # TODO this might not be right (factor 2 maybe?)
            #om, S_X = autospectrum_fft(tau, kappa_XX, onesided=True, method='fft_wrap')
            om, S_X = autospectrum_fft(tau, kappa_XX, onesided=True, method='rfft')
            self._omega   =  om       
            self.S_X      =  S_X      


        # --- Averge pdfs
        if 'pdf' in stats:
            ymin = np.min(xi)
            ymax = np.max(xi)
            ybins = np.linspace(ymin, ymax, nPDF)
            pdfs=np.zeros((self.nSamples, len(ybins)-1))
            # TODO use welib.tools.stats when ready
            def pdf_histogram(y,bins=50, norm=True, count=False):
                yh, xh = np.histogram(y[~np.isnan(y)], bins=bins)
                dx   = xh[1] - xh[0]
                xh  = xh[:-1] + dx/2
                if norm:
                    yh=yh/np.trapezoid(yh,xh)
                return xh,yh

            for i in range(self.nSamples):
                xb,yb = pdf_histogram(xi[i,:], bins=ybins)
                pdfs[i,:] = yb 
            self.pdfs_xi  = pdfs
            yh =  np.mean(pdfs,axis=0)
            yh = yh/np.trapezoid(yh,xb)
            self.pdf_xi   = yh
            self.x_pdf_xi = xb
            var.set_pdf(xb, yh)
        self.X_xi = var


    # --------------------------------------------------------------------------------}
    # --- Functions using analytical/lambda functions provided
    # --------------------------------------------------------------------------------{
    def autospectrum_int(self, omega=None, tau_max=None, method='quad', nDiscr=None):
        if omega is None:
            omega = self.omega_default
        if tau_max is None:
            #tau_max = np.max(self.time_default)
            tau_max = self.tau_max

        if not self._f_kappa_XX_th:
            raise Exception('Provide an autocovariance function')

        if self.verbose:
            print('>> S_X(omega) from int k_XX(tau), method:{}, with tau=[{} ; {} ]'.format(method, 0, tau_max))
        if method=='quad':
            if self.verbose:
                print('                         omega from {} to {} delta:{}'.format(omega[0], omega[-1], omega[1]-omega[0]))
            S_X_i = autospectrum_int(omega, self._f_kappa_XX_th, tau_max=tau_max)

        elif method=='fft':
            if nDiscr is None:
                nDiscr=self.nDiscr
            tau = np.linspace(0, tau_max, nDiscr)
            kappa_XX = np.array([self._f_kappa_XX_th(t) for t in tau])
            omega, S_X_i = autospectrum_fft(tau, kappa_XX, onesided=True, method='rfft')
            #omega, S_X_i = autospectrum_fft(tau, kappa_XX, onesided=True, method='fft_wrap')
#             dtau = tau[1]-tau[0]
#             print('tau_max', tau[-1])
#             print('dtau',  tau[1]-tau[0])
#             print('1/dtau',  1/dtau)
#             print('n   ',  len(tau))
#             S_X_i*=tau_max/2
#             S_X_i*=5
        else:
            raise NotImplementedError()

        self._S_X_th_i[method] = (omega, S_X_i)



        return S_X_i

    def autocovariance_int(self, tau=None, omega_max=None, method='quad', nDiscr=None):
        if tau is None:
            tau = self.tau_default
        if omega_max is None:
            omega_max = self.omega_max

        if not self._f_S_X_th:
            raise Exception('Provide an autospectrum function')

        if self.verbose:
            print('>> k_XX(tau) from int S_X(omega), method:{}, with omega=[{} ; {} ]'.format(method,0, omega_max))
        if method=='quad':
            if self.verbose:
                print('                         tau from {} to {} delta:{}'.format(tau[0], tau[-1], tau[1]-tau[0]))
            kappa_XX = autocovariance_int(tau, self._f_S_X_th, omega_max=omega_max)
        elif method=='fft':
            if nDiscr is None:
                nDiscr=self.nDiscr

            omega = np.linspace(0, omega_max, nDiscr)
            # If we want to respect tau_max...
            #domega = np.pi/self.tau_max
            #omega = np.arange(0, omega_max, domega)

            S_X = np.array([self._f_S_X_th(abs(om)) for om in omega])
            #omega = np.linspace(-omega_max, omega_max, self.nDiscr)
            #S_XX = np.array([self._f_S_X_th(abs(om)) for om in omega])/2
            tau, kappa_XX = autocovariance_fft(omega, S_X, onesidedIn=True)
        else:
            raise NotImplementedError()

        self._kappa_XX_th_i[method] = (tau, kappa_XX)

        return kappa_XX

    def autospectrum_moments(self, orders=None, method='num'):
        if orders is None:
            orders=[0,1,2,3]
        if not self._f_S_X_th:
            raise Exception('Provide an autospectrum function')
        if method=='num':
            omega = self.omega_default
            S_X = np.array([self._f_S_X_th(abs(om)) for om in omega])

        moments=dict()
        for i in orders:
            if method=='quad':
                integrand = lambda omega : omega**i * self._f_S_X_th(omega)
                moments[i], _ = nquad(integrand, [[0, self.omega_max]] )
            else:
                moments[i] = np.trapezoid( S_X * omega**i , omega)

        return moments

    # --------------------------------------------------------------------------------}
    # --- Plots 
    # --------------------------------------------------------------------------------{
    def plot_samples(self, ax=None, maxSamples=5, **kwargs):
        # --- Plot realisations
        if ax is None:
            fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
            fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
        nPlot = min(self.nSamples,maxSamples)
        for i in range(nPlot):
            ax.plot(self.time, self.xi[i], **kwargs)
        ax.set_xlabel('Time [s]')
        ax.set_ylabel('')
        #ax.legend()
        ax.tick_params(direction='in')
        ax.set_title(self.name + ' Samples realisations ({}/{}) '.format(nPlot, self.nSamples))
        return ax


    def plot_mean(self, ax=None):
        # --- Plot mean and sigma
        if ax is None:
            fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
            fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
        if self.mu_xi is not None: 
            ax.plot(self.time, self.mu_xi , label='Mean')
            if self.verbose:
                print('>>> Mean Mean num:', np.around(np.mean(self.mu_xi ),4))
        if self._mu_th is not None:
            ax.plot(self.time, self.time*0+self._mu_th, 'k--' , label='Mean (th)')
            if self.verbose:
                print('>>> Mean Mean th :', np.around(self._mu_th,4))
        #ax.plot(time, var_xi , label='Variance')
        #ax.plot(time, time*0+s**2, 'k--' , label='Variance (th)')
        ax.set_xlabel('Time [s]')
        ax.set_ylabel('mu_X(t)')
        #ax.legend()
        ax.tick_params(direction='in')
        ax.set_title(self.name + ' Mean')
        return ax

    def plot_var(self, ax=None):
        # --- Plot mean and sigma
        if ax is None:
            fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
            fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
        if self.var_xi is not None: 
            ax.plot(self.time, self.var_xi , label='Variance')
            ax.plot(self.time, self.time*0+np.mean(self.var_xi ), '--' ,c=(0.5,0.5,0.5))
            if self.verbose:
                print('>>> Mean Var num:', np.around(np.mean(self.var_xi ),4))
        if self._var_th is not None:
            ax.plot(self.time, self.time*0+self._var_th, 'k--' , label='Variance (th)')
            if self.verbose:
                print('>>> Mean Var th :', np.around(self._var_th,4))
        #ax.plot(time, var_xi , label='Variance')
        ax.set_xlabel('Time [s]')
        ax.set_ylabel('VarX(t)')
        #ax.legend()
        ax.tick_params(direction='in')
        ax.set_title(self.name+ ' Variance')
        return ax

    def plot_autocovariance(self, ax=None):
        # --- Plot AutoCovariance
        if ax is None:
            fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
            fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
        if self.kappa_XX is not None: 
            ax.plot(self.tau, self.kappa_XX, label='kappa_mean')
            if self.verbose:
                print('>>> Correlation length num', np.around(np.trapezoid(self.kappa_XX, self.tau,4)))

        if self._f_kappa_XX_th is not None:
            #kappa_XX_th = self._f_kappa_XX_th(self.tau)
            k_XX = [self._f_kappa_XX_th(t) for t in self.tau]
            ax.plot(self.tau, k_XX , 'k--', label='kappa theory')
            if self.verbose:
                print('>>> Correlation length th ', np.around(np.trapezoid(k_XX, self.tau,4)))

        for k,v in self._kappa_XX_th_i.items():
            tau, k_XX = self._kappa_XX_th_i[k]
            ax.plot(tau, k_XX,':'   , label=r'$\int S_X$ provided')
            ax.set_xlim([0,self.tau_max])
            if self.verbose:
                print('>>> Correlation length {:s}'.format(k[:3]), np.around(np.trapezoid(k_XX, tau,4)))
        ax.set_xlabel('Tau [s]')
        ax.set_ylabel('kappa_XX')
        ax.legend()
        ax.tick_params(direction='in')
        ax.set_title(self.name + 'Autocovariance')
        return ax

    def plot_autocorrcoeff(self, ax=None):
        # --- Plot AutoCorrCoeff
        if ax is None:
            fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
            fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
        # for i in range(min(nSamples,5)):
        #     ax.plot(tau, rho_XXi[i,:], alpha=0.5)
        if self.rho_XX is not None:
            ax.plot(self.tau, self.rho_XX, label='rho_mean (num)')

        if self._f_kappa_XX_th and self._var_th :
            #kappa_XX_th = self._f_kappa_XX_th(self.tau)
            k_XX = np.array([self._f_kappa_XX_th(t) for t in self.tau])
            rho_XX_th = k_XX / self._var_th
            ax.plot(self.tau, rho_XX_th , 'k--', label=r'$\rho_XX$ provided')

        for k,v in self._kappa_XX_th_i.items():
            tau, k_XX = self._kappa_XX_th_i[k]
            ax.plot(tau, k_XX / self._var_th ,':'   , label=r'$\int S_X$ provided')
            ax.set_xlim([0,self.tau_max])

        ax.set_xlabel('Tau [s]')
        ax.set_ylabel('rho_XX')
        ax.legend()
        ax.tick_params(direction='in')
        ax.set_title(self.name + ' autocorrelation coefficient')
        return ax

    
    def plot_autospectrum(self, ax=None):
        # --- auto
        if ax is None:
            fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
            fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
        # ax.plot([omega_0,omega_0], [0,S_X_th] ,'k--'   , label='Theory', lw=2)
        if self.S_X is not None:
            ax.plot(self.omega, self.S_X ,'o-', label='FFT (kappa_num)')
            if self.verbose:
                print('>>> Spectrum integration num', np.around(np.trapezoid(self.S_X, self.omega),4))

        if self.S_avg is not None:
            ax.plot(self.om_Savg, self.S_avg ,'.-', label='FFT (avg)')
            if self.verbose:
                print('>>> Spectrum integration avg', np.around(np.trapezoid(self.S_avg, self.om_Savg),4))

        if self._f_S_X_th is not None:
            #S_X_th = self._f_S_X_th(self.omega)
            S_X_th = [self._f_S_X_th(om) for om in self.omega]
            ax.plot(self.omega, S_X_th ,'k--', label='S_X(omega) provided')
            if self.verbose:
                print('>>> Spectrum integration th ', np.around(np.trapezoid(S_X_th, self.omega),4))

        for k,v in self._S_X_th_i.items():
            om, S_X = self._S_X_th_i[k]
            ax.plot(om, S_X ,':'   , label=r'$\int \kappa$ provided {}'.format(k))
            b = om<self.omega_max # WEIRD
            #b = om>0
            ax.set_xlim([0,self.omega_max])
            if self.verbose:
                print('>>> Spectrum integration {:s}'.format(k[:3]), np.around(np.trapezoid(S_X[b], om[b]),4))
        ax.set_xlabel(r'$\omega$ [rad/s]')
        ax.set_ylabel(r'$S_X(\omega)$')
        ax.legend()
        ax.tick_params(direction='in')
        ax.set_title(self.name + ' autospectrum')
        return ax

    def plot_pdf(self, ax=None):
        if ax is None:
            fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
            fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
        if self.x_pdf_xi is not None:
            xb = self.x_pdf_xi
            for i in range(self.nSamples):
                yb = self.pdfs_xi[i,:]
                ax.plot( xb, yb, c=(0.5,0.5,0.5), label='Samples' if i==0 else None)
            ax.plot( xb, self.pdf_xi ,'k-', label='Average from samples')
        ax.set_xlabel('x')
        ax.set_ylabel('pdf')
        ax.legend()
        ax.tick_params(direction='in')
        ax.set_title(self.name + ' pdf')
        return ax

    def __repr__(self):
        s='<{} object with attributes>:\n'.format(type(self).__name__)
        s+='- name     : {}\n'.format(self.name)
        s+='- params   : {}\n'.format(self.params)
        s+='- time_max : {}\n'.format(self.time_max)
        s+='- tau_max  : {}\n'.format(self.tau_max)
        s+='- _var_th  : {}\n'.format(self._var_th)
        s+='- _mu_th   : {}\n'.format(self._mu_th)
        x=self.time ; s+=' * time   : [{} ... {}],  d: {}, n: {} \n'.format(x[0], x[-1], x[1]-x[0], len(x))
        x=self.tau  ; s+=' * tau    : [{} ... {}],  d: {}, n: {} \n'.format(x[0], x[-1], x[1]-x[0], len(x))
        x=self.omega; s+=' * omega  : [{} ... {}],  d: {}, n: {} \n'.format(x[0], x[-1], x[1]-x[0], len(x))
        s+='Main methods:\n'
        s+='- generate_samples, generate_samples_from_autospectrum\n'
        s+='- plot_samples, plot_mean, plot_var  (after `generate_samples`)\n'
        s+='- plot_autocovariance, plot_autospectrum, plot_autocorrcoeff\n'
        return s

    def picklable(self):
        def noneIfLambda(attr):
            obj = getattr(self, attr)
            if callable(obj): 
                if obj.__name__ == "<lambda>" or str(obj).find('<locals>')>0:
                    print('SampledStochasticProcess: picklable: removing ', attr)
                    setattr(self, attr, None)

        noneIfLambda('_f_kappa_XX_th')
        noneIfLambda('_f_S_X_th')


# --------------------------------------------------------------------------------
# --- Predefined processes
# --------------------------------------------------------------------------------
class HarmonicProcess(SampledStochasticProcess):

    def __init__(self, s, omega0, generator=None, omega_max=None, tau_max=None, time_max=None, nDiscr=100, verbose=False, **kwargs):

        def harmonic_process_realisation(time, s=1, omega0=2*np.pi):
            r = np.random.rayleigh(scale=s, size=1)
            p = np.random.uniform(low=0, high=2*np.pi, size=1)
            x = r * np.cos(omega0 * time + p)
            return x
        # Distribution parameters
        params = {'s':s, 'omega0':omega0}
        
        # Derived parameters
        T = 1/(omega0/(2*np.pi))
        tau_max   = 10*T      if tau_max is None   else tau_max
        omega_max = 3*omega0  if omega_max is None else omega_max

        # --- Initialize Process
        generator = lambda time: harmonic_process_realisation(time, s=s, omega0=omega0)
        SampledStochasticProcess.__init__(self, params=params, generator = generator, name='harmonic', omega_max=omega_max, tau_max=tau_max, time_max=time_max, nDiscr=nDiscr, verbose=verbose, **kwargs)
        # Theory
        self._f_kappa_XX_th = lambda tau: s**2 * np.cos(omega0*tau)
        self._var_th        = s**2  
        self._mu_th         = 0

class ExponentialProcess(SampledStochasticProcess):

    def __init__(self, lambd, sigma, generator=None, omega_max=None, tau_max=None, time_max=None, nDiscr=100, verbose=False, **kwargs):

        # Distribution parameters
        params = {'lambd':lambd, 'sigma':sigma}

        # Derived parameters
        omega_max = - 1/lambd* np.log( 1e-3/ (lambd*sigma**2)) if omega_max is None else omega_max
        tau_max   = np.sqrt( lambd**2 * (sigma**2/1e-3-1))     if tau_max is None else tau_max
        time_max  = 2*tau_max                                  if time_max is None else time_max
     
        # --- Initialize Process
        SampledStochasticProcess.__init__(self, params=params, generator = generator, name='exponential', omega_max=omega_max, tau_max=tau_max, time_max=time_max, nDiscr=nDiscr, verbose=verbose, **kwargs)
        # Theory
        self._f_S_X_th = lambda omega: sigma**2 * lambd *np.exp(-lambd*omega)
        self._f_kappa_XX_th = lambda tau:  sigma**2 / (1+tau**2/lambd**2)
        self._var_th = sigma**2
        self._mu_th = np.sqrt(sigma**2 * lambd/ self.time_default[-1] *(2*np.pi) )

class BandedWhiteNoiseProcess(SampledStochasticProcess):

    def __init__(self, omega0, zeta, sigma, generator=None, omega_max=None, tau_max=None, time_max=None, nDiscr=100, verbose=False, **kwargs):

        # Distribution parameters
        params = {'omega0':omega0, 'zeta':zeta, 'sigma':sigma}

        # Derived parameters
        B    = 2*zeta*omega0
        omega_1 = omega0-B/2
        omega_2 = omega0+B/2
        T = 1/(omega0/(2*np.pi))
        tau_max = 20*T          if tau_max   is None else tau_max
        time_max = 2*tau_max    if time_max  is None else time_max
        omega_max = 2*omega_2   if omega_max is None else omega_max
     
        # --- Initialize Process
        SampledStochasticProcess.__init__(self, params=params, generator = generator, name='BandedWhiteNoise', omega_max=omega_max, tau_max=tau_max, time_max=time_max, nDiscr=nDiscr, verbose=verbose, **kwargs)

        # Theory
        def f_S_X_th(omega, omega0, B, sigma_X):
            if omega>=omega0-B/2 and omega<=omega0+B/2:
                return sigma_X**2/B
            else:
                return 0
        def f_k_XX_th(tau, omega0, B, sigma_X):
            if tau==0:
                return sigma_X**2
            else:
                return sigma_X**2 / (B*tau) * ( np.sin((omega0+B/2)*tau) - np.sin((omega0-B/2)*tau)) 

        self._f_kappa_XX_th = lambda tau:   f_k_XX_th(tau,  omega0=omega0, B=B, sigma_X=sigma)
        self._f_S_X_th      = lambda omega: f_S_X_th(omega, omega0=omega0, B=B, sigma_X=sigma)
        self._var_th = sigma**2
        self._mu_th = 0

class CutOffWhiteNoiseProcess(SampledStochasticProcess):

    def __init__(self, B, sigma, generator=None, omega_max=None, tau_max=None, time_max=None, nDiscr=100, verbose=False, **kwargs):

        # Distribution parameters
        params = {'B':B, 'sigma':sigma}

        # Derived parameters
        omega_max = B if omega_max is None else omega_max
        tau_max   = B if tau_max is None else tau_max
        time_max  = B if time_max is None else time_max
     
        # --- Initialize Process
        SampledStochasticProcess.__init__(self, params=params, generator = generator, name='CutOffWhiteNoise', omega_max=omega_max, tau_max=tau_max, time_max=time_max, nDiscr=nDiscr, verbose=verbose, **kwargs)


        # Theory
        self.setTh()

    def setTh(self):
        B     = self.params['B']
        sigma = self.params['sigma']

        def f_S_X_th(omega):
            if omega>=-B/2 and omega<=B/2:
                return 2* sigma**2/B
            else:
                return 0
        def f_k_XX_th(tau):
            if tau==0:
                return sigma**2
            else:
                return 2 * sigma**2 / (B*tau) * (np.sin(B/2*tau)  )

        self._f_kappa_XX_th = f_k_XX_th
        self._f_S_X_th      = f_S_X_th
        self._var_th = sigma**2
        self._mu_th = 0



class KaimalProcess(SampledStochasticProcess):

    def __init__(self, U0, sigma, L, omega_max=None, tau_max=None, time_max=None, nDiscr=100, verbose=False, **kwargs):
        """ 
         - U0    : mean wind speed
         - sigma : standard deviation
         - L     : length scale
         """
        from welib.wind.spectra import kaimal
        from welib.wind.windsim import pointTSKaimal
        from welib.stoch.utils import sample_from_autospectrum

        # Distribution parameters
        params = {'U0':U0, 'sigma':sigma, 'L':L}

        # Derived parameters
        f_S_X_th = lambda omega: kaimal(omega, U0, sigma, L, angularFrequency=True)
        # TODO merge this with generate sample and generate samples from auto_spectrum
        # Alternative: use def pointTSKaimal(tMax, dt, U0, sigma, L, angularFrequency=False, **kwargs):
        def kaimal_generator(time):
            tMax = time[-1]
            dt = (time[-1]-time[0])/(len(time)-1)
            _, xi, _, _ = sample_from_autospectrum(tMax=tMax, dt=dt, f_S=f_S_X_th, angularFrequency=True)
            return xi+U0

        # --- Initialize Process
        SampledStochasticProcess.__init__(self, params=params, generator = kaimal_generator, name='Kaimal', omega_max=omega_max, tau_max=tau_max, time_max=time_max, nDiscr=nDiscr, verbose=verbose, **kwargs)

        # Theory
        #f_S = lambda f_or_om: 
        #def f_S_X_th(omega): #, B, sigma):
        #    return kaimal(omega, U0, sigma, L, angularFrequency=True)

    #         def f_k_XX_th(tau, B, sigma):
    #                 return 2 * sigma**2 / (B*tau) * (np.sin(B/2*tau)  )
    # 
    #         self._f_kappa_XX_th = lambda tau:   f_k_XX_th(tau,  B=B, sigma=sigma)
        #self._f_S_X_th      = lambda omega: f_S_X_th(omega, B=B, sigma=sigma)
        self._f_S_X_th      = f_S_X_th
        self._var_th = sigma**2
        self._mu_th = U0


if __name__ == '__main__':
    stats=[]
    stats+=['avg_spectra']
    stats+=['correlation']
    verbose = True
    nSamples= 130
    nDiscr=100 

    # --- Harmonic process
    s      = 4
    omega0 = 2*np.pi
    proc = HarmonicProcess(s=s, omega0=omega0, verbose=True, nDiscr=nDiscr)

    # --- Samples
    T = 1/(omega0/(2*np.pi))
    time=np.linspace(0,40*T, 1000)
    dt = (time[-1]-time[0])/(len(time)-1)
    # NOTE: no autospectrum
    xi = proc.generate_samples(nSamples=nSamples, time=time)
    proc.samples_stats(stats=stats) #nTau = None)
    proc.plot_samples()
    proc.plot_var()
    # proc.plot_mean()

    # --- Computations
    proc.autospectrum_int(method='quad')
    proc.autospectrum_int(method='fft', nDiscr=20*proc.nDiscr, tau_max=proc.tau_max*3)
    #proc.autocovariance_int(method='fft')

    proc.plot_autocovariance()
    # proc.plot_autocorrcoeff()

    #
    ax = proc.plot_autospectrum()
    ax.plot([omega0,omega0], [0,proc._var_th] ,'k--'   , label='Theory', lw=2)
    ax.set_xlabel(r'$\omega$ [rad/s]')
    ax.set_ylabel(r'$S_X(\omega)$')
    ax.set_xlim([0, 2*omega0])
    ax.legend()


    plt.show()
