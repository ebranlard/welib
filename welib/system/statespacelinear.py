import numpy as np
from numpy.linalg import inv
from numpy.linalg import solve
from numpy.linalg import eigvals, matrix_rank
from scipy.integrate import  solve_ivp #odeint
from scipy.interpolate import interp1d
from scipy.optimize import OptimizeResult as OdeResultsClass 
from scipy.linalg import expm
# Local

# --------------------------------------------------------------------------------}
# --- Simple statespace functions ltiss (linear time invariant state space)
# --------------------------------------------------------------------------------{
def state_function(t, x, u, p):
    return p['A'].dot(x) + p['B'].dot(u)

def output_function(t, x, u, p):
    return p['C'].dot(x) + p['D'].dot(u)

def integrate(t_eval, q0, A, B, fU, method='LSODA', **options):
    """ 
    Perform time integration of a LTI state space system

    INPUTS:
     - q0: initial states, array of length nStates
     - A: state matrix (nStates x nStates)
     - B: input matrix (nStates x nInputs)
     - fU: function/interpolants interface U=fU(t) or U=fU(t,q)
          U : array of inputs

    OUTPUTS:
     - res: object with attributes `t` and `y`(states for now..) and other attributse from solve_ivp

    """
    hasq=False
    try:
        fU(t_eval[0],q0)
        hasq=True
    except:
        try:
            fU(t_eval[0])
            hasq=False
        except:
            raise

    if hasq:
        odefun = lambda t, q : np.dot(A, q) + np.dot(B, fU(t,q))
    else:
        odefun = lambda t, q : np.dot(A, q) + np.dot(B, fU(t) )

    res = solve_ivp(fun=odefun, t_span=[t_eval[0], t_eval[-1]], y0=q0, t_eval=t_eval, method=method, vectorized=False, **options)   

    # TODO consider returning y

    return res


def integrate_convolution(time, A, B, fU, C=None):
    """ 
    Perform time integration of a LTI state space system using convolution method

    INPUTS:
     - A: state matrix (nStates x nStates)
     - B: input matrix (nStates x nInputs)
     - fU: function/interpolants with interface U=fU(t) or U=fU(t,q)
           where U is the array of inputs at t

    OUTPUTS:
     - x: state vector

    """
    H = impulse_response_matrix(time, A, B)

    x = np.zeros((A.shape[0], len(time)))

    from welib.tools.signal_analysis import convolution_integral
    # TODO inline and optimize
    try:
        U = fU(time)
    except:
        print('[WARN] Cannot evaluate fU for all time, doing a for loop...')
        U = np.zeros((B.shape[1], len(time)))
        for it,t in enumerate(time):
            U[:,it] = fU(t) # NOTE: cannot do state dependency here

    for i in np.arange(H.shape[0]):
        x_sum=0
        for j in np.arange(H.shape[1]):
            x_sum  += convolution_integral(time, U[j,:], H[i,j,:] )
        x[i,:] = x_sum
    return x

    # TODO consider returning y

def impulse_response_matrix(time, A, B, C=None, outputBoth=False):
    """ 
    Return the impulse response matrix for all time steps defined by `time`
        H_x(t) =   exp(At) B
        H_y(t) = C exp(At) B
        see e.g. 
           Friedland p 76
    """
    H_x = np.zeros((A.shape[0], B.shape[1], len(time)))
    for it, t in enumerate(time):
        H_x[:,:, it] = expm(A*t).dot(B)

    if outputBoth:
        raise NotImplementedError()
        if C is None:
            raise Exception('Provide `C` to output both impulse response matrices H_x and H_y')
        H_y = C.dot(H_x) # TODO verify
        return H_x, H_y
    else:
        return H_x

# --------------------------------------------------------------------------------}
# --- Linear State Space system
# --------------------------------------------------------------------------------{
class LinearStateSpace():
    """ 
    def setStateInitialConditions(self,q0=None):
    def setInputTimeSeries(self,vTime,vU):
    def setInputFunction(self,fn):
    def Inputs(self,t,x=None):
    def integrate(self, t_eval, method='RK4', y0=None, **options):
    def dqdt(self, t, q):
    def RHS(self,t,q):
    def nStates(self):
    def nInputs(self):
    def nOuputs(self):
    """
    def __init__(self,A,B,C=None,D=None,q0=None):
        self.A=np.asarray(A)
        self.B=np.asarray(B)
        if C is None:
            self.C=np.eye(A.shape[0]) # output all states
        else:
            self.C=np.asarray(C)
        if D is None:
            self.D=np.zeros((self.C.shape[0],0))
        else:
            self.D=np.asarray(D)

        # Initial conditions
        self.setStateInitialConditions(q0)

        # Time integration results
        self.res=None

    @property
    def nStates(self):
        return self.A.shape[0]

    @property
    def nInputs(self):
        return self.B.shape[1]

    @property
    def nOuputs(self):
        if self.C is not None:
            return self.C.shape[1]
        else:
            return 0

    # --------------------------------------------------------------------------------}
    # --- Time domain 
    # --------------------------------------------------------------------------------{
    def setStateInitialConditions(self,q0=None):
        self.q0 = np.zeros(self.nStates)
        if q0 is not None:
            if len(q0)!=self.nStates:
                raise Exception('Wrong dimension for q0 ({} instead of {} )'.format(len(q0),self.nStates))
            self.q0 = q0

    def setInputTimeSeries(self,vTime,vU):
        """ 
        Set the inputs as a time series of time
        INPUTS:
         - vTime: 1d array of time steps (do not need to be regular), of length nt
         - vF   : nStates x nt array of forces at each time steps, for the nStates states
        """
        vTime = np.asarray(vTime)
        vU    = np.asarray(vU)
        if self.nInputs==1:
            vU = vU.reshape((1,-1))
        if vU.shape[0]!=self.nInputs:
            raise Exception('Wrong first dimension for Inputs time series ({} instead of {} )'.format(vU.shape[0],self.nInputs))
        if vU.shape[1]!=len(vTime):
            raise Exception('Second dimension of Input time series does not match time dimension ({} instead of {} )'.format(vU.shape[1],len(vTime)))

        # Store raw data
        self._inputs_ts = vU
        self._time_ts  = vTime
        # Create interpolant for faster evaluation
        self._inputs_fn_t = interp1d(vTime, vU)

    def setInputFunction(self,fn):
        """ 
        Set the inputs as a function of time and states
        The function will be used during the time integration

        INPUTS:
          fn: handle to a python function. The interface of the function fn is: 

               u =  f(t, q) 

              where:
                 t   : scalar, time
                 q   : (nStates,) array, states
                 u   : (nInputs,) array, inputs 
        
        """
        self._inputs_fn = fn

    def Inputs(self, t, q=None):
        if hasattr(self,'_inputs_fn_t'):
            return self._inputs_fn_t(t) # NOTE: interpolant, works for scalar or arrays
        elif hasattr(self,'_inputs_fn'):
            return self._inputs_fn_t(t, q=q) # TODO what if user doesn't care for states..
        else:
            raise Exception('Please specify a time series of inputs using `setInputsTimeSeries` or a function using `setInputsFunction` ')


    def integrate(self, t_eval, method='RK4', y0=None, **options):
        #
        if y0 is not None:
            self.setStateInitialConditions(y0)

        if method.lower()=='impulse':
            # TODO add check on initial conditions
            x = integrate_convolution(t_eval, self.A, self.B, self.Inputs)

            res = OdeResultsClass(t=t_eval, y=x) # To mimic result class of solve_ivp

        else:
            res = integrate(t_eval, self.q0, self.A, self.B, self.Inputs, method=method, **options)

        if self.nOuputs>0:
            print('>>> TODO: do something to compute outputs after states')

        # Store
        self.res    = res

        return res

    def dqdt(self, t, q):
        # NOTE: this can cause issues if q is not flat
        return np.dot(self.A, q) + np.dot(self.B, self.Inputs(t,q))

    def RHS(self,t,q):
        return self.dqdt(t,q)

    # --------------------------------------------------------------------------------}
    # --- Frequency domain and transfer function
    # --------------------------------------------------------------------------------{
    def transferFunction(self, s):
        """Evaluate the systems's transfer function for a complex variable

        H(s) = C [sI-A]^-1 B + D

        Returns a matrix of values evaluated at complex variable s.

        A more efficient function may be found in TB05Ad from Slycot
        from slycot import tb05ad
        see control.statesp.horner
        """
        if hasattr(s, '__len__'):
            nu = self.nInputs
            ny = self.nOuputs
            ns = len(s)
            H = np.empty((ny, nu, ns), dtype=np.complex128)
            for k,sk in enumerate(s.ravel()):
                H[:,:,k] =np.dot(self.C, solve(sk* np.eye(self.nStates) - self.A, self.B)) + self.D
            #H = np.array([np.dot(self.C, solve(sk* np.eye(self.nStates) - self.A, self.B)) + self.D  for sk in s.ravel()])
            #H=H.reshape([ny,nu]+list(s.shape))
        else:
            H =           np.dot(self.C, solve(s * np.eye(self.nStates) - self.A, self.B)) + self.D

        return H


    def frequency_response(self, omega):
        """Evaluate the system's transfer function at a list of frequencies
        Reports the frequency response of the system,

             H(j*omega) = mag*exp(j*phase)

        for continuous time. For discrete time systems, the response is
        evaluated around the unit circle such that

             H(exp(j*omega*dt)) = mag*exp(j*phase).

        Parameters
        ----------
        omega : array_like
            A list of frequencies in radians/sec at which the system should be
            evaluated. The list can be either a python list or a numpy array

        Returns
        -------
        mag : (self.outputs, self.inputs, len(omega)) ndarray
            The magnitude (absolute value, not dB or log10) of the system
            frequency response.
        phase : (self.outputs, self.inputs, len(omega)) ndarray
            The wrapped phase in radians of the system frequency response.
        """
        omega = np.asarray(omega)

        numFreqs = len(omega)
        #Gfrf = np.empty((self.outputs, self.inputs, numFreqs), dtype=np.complex128)
        #if isdtime(self, strict=True):
        #    dt = timebase(self)
        #    cmplx_freqs = exp(1.j * omega * dt)
        #    if max(np.abs(omega)) * dt > math.pi:
        #        warn("freqresp: frequency evaluation above Nyquist frequency")
        #else:
        s = omega * 1.j
        H = self.transferFunction(s)
        return np.abs(H), np.angle(H)

    # --------------------------------------------------------------------------------}
    # ---  IO functions for printing/plotting/saving
    # --------------------------------------------------------------------------------{
    def __repr__(self):
        s='<{} object>\n'.format(type(self).__name__)
        s+='|Read-only attributes:\n'
        s+='| - nState:{} \n'.format(self.nStates)
        if hasattr(self,'_force_ts'):
            if len(self._time_ts)>1:
                dt=self._time_ts[1]-self._time_ts[0]
            else:
                dt=np.nan
            s+='|Force time series \n'
            s+='| - Time: [{} ... {}],  dt: {}, n: {} \n'.format(self._time_ts[0],self._time_ts[-1],dt,len(self._time_ts))
            s+='| - Force t0  : {} \n'.format(self._force_ts[:,0])
            s+='| - Force tend: {} \n'.format(self._force_ts[:,-1])
        s+='|Attributes:\n'
        s+='| - A: State-State Matrix  \n'
        s+=str(self.A)+'\n'
        s+='| - B: State-Input Matrix  \n'
        s+=str(self.B)+'\n'
        s+='| - q0: Initial conditions (state) \n'
        s+=str(self.q0)+'\n'
        return s

    def plot(self, fig=None, axes=None, label=None, res=None, **kwargs):
        """ Simple plot of states after time integration"""
        if res is None:
            res=self.res
        if res is None:
            raise Exception('Call integrate before plotting, or provide a `res` output structure')

        if axes is None:
            import matplotlib.pyplot as plt
            fig,axes = plt.subplots( self.nStates,1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
            fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
            axes_provided=False
        else:
            axes_provided=True

        axes = np.atleast_1d(axes)
        for i,ax in enumerate(axes):
            lbl=r'$x_{}$'.format(i+1)
            ax.plot(res.t, res.y[i,:], label=label, **kwargs)
            if not axes_provided:
                ax.set_ylabel(lbl)
            ax.tick_params(direction='in')
        if not axes_provided:
            axes[-1].set_xlabel('Time [s]')

        return fig, axes

    def plot_inputs(self, axes=None, label=None, res=None, **kwargs):
        """ 
        Simple plot of inputs
          - if `res` is provided, use time and states from res
          - if time integration was performed, use time and states from self.res
          otherwise: 
          - if a time series was provided, plot that
          - if a function was provided for the forcing, abort
        """
        if res is None:
            res=self.res
        if res is None:
            # Res or self.res not provided
            if hasattr(self,'_inputs_ts'):
                time = self._time_ts
                U    = self._inputs_ts
            else:
                raise NotImplementedError()
        else:
            # Res provided
            time =res.t
            if hasattr(self,'_inputs_ts'):
                U = self.Inputs(time)
            else:
                U = np.zeros((self.nInputs,len(time)))
                for it, t in enumerate(time):
                    U[:,it] = self.Inputs(t, q=res.y[:,it])
        ## Plot
        if axes is None:
            import matplotlib.pyplot as plt
            fig,axes = plt.subplots( self.nInputs, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
            fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
        axes = np.atleast_1d(axes)
        for i,ax in enumerate(axes):
            lbl=r'$u_{}$'.format(i+1)
            ax.plot(time, U[i,:], label=label, **kwargs)
            ax.set_ylabel(lbl)
            ax.tick_params(direction='in')
        axes[-1].set_xlabel('Time [s]')

        return fig, axes

