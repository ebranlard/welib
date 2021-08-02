import numpy as np
from numpy.linalg import inv
from scipy.integrate import  solve_ivp #odeint
from scipy.optimize import OptimizeResult as OdeResultsClass 
from scipy.interpolate import interp1d
from .statespace import StateMatrix


# --------------------------------------------------------------------------------}
# --- Mech System Class
# --------------------------------------------------------------------------------{
# NOTE: this is time invariant!
class MechSystem():
    """ 
    Tools to handle/simulate a mechanical system of the form

        M(x) xddot + C xdot + K x = F(t,x,xdot) 

    Integrated as:
        [xdot ] = [ 0          I     ].[x   ] + [  0   ]
        [xddot]   [-M^-1 K    -M^-1 C] [xdot]   [M^-1 F]

          qdot =             A~(x)      .q      +    B~(t,x,xdot)
    where 
       - M: mass matrix, constant or function of x
       - K: stiffness matrix, constant matrix or None
       - C: damping matrix, constant matrix or None
       - F: forcing term, function of t,x,xdot, or, tuple (time, Force time series)

       - q : full state vector

       NOTE: A~, B~ are not linear state space matrices in the general case

    """
    def __init__(self, M, C=None, K=None, F=None, x0=None, xdot0=None):

        # Functions that return the mass matrix and its inverse
        if hasattr(M, '__call__'):  
            self.M_is_func=True
            self._fM    = M
            self._fMinv = lambda x : inv(M(x))
        else:
            M = np.atleast_2d(M)
            self.M_is_func=False
            self._Minv= inv(M) # compute inverse once and for all
            self._fM    = lambda x : M
            self._fMinv = lambda x : self._Minv
            if K is None:
                K = M*0
            else:
                K = np.atleast_2d(K)
            if C is None:
                C = M*0
            else:
                C = np.atleast_2d(C)

        # Store
        self.M = M
        self.C = C
        self.K = K
        self.has_A      = C is not None or K is not None # do we have a constant state matrix?
        #self.has_A      = True

        # Find number of DOF, now assumed constant until reinit
        self.nDOF = self.find_nDOF()

        # Forcing
        if F is not None:
            if type(F) is tuple:
                self.setForceTimeSeries(F[0], F[1])
            else:
                self.setForceFunction(F)

        # Initial conditions
        self.setInitialConditions(x0,xdot0)

        # Time integration results
        self.res=None


    def find_nDOF(self):
        if self.M_is_func:
            x=np.zeros((5000,1)) # Using a huge number of states..
            return len(self.M(x))
        else:
            return len(self.M)
    @property
    def nStates(self):
        return 2*self.nDOF

    # --------------------------------------------------------------------------------}
    # --- Time domain
    # --------------------------------------------------------------------------------{
    def setInitialConditions(self,x0,xdot0):
        self.q0 = np.zeros(2*self.nDOF)
        if x0 is not None:
            if len(x0)!=self.nDOF:
                raise Exception('Wrong dimension for x0 ({} instead of {} )'.format(len(x0),self.nDOF))
            self.q0[:self.nDOF] =  x0
        if xdot0 is not None:
            if len(xdot0)!=self.nDOF:
                raise Exception('Wrong dimension for xdot0 ({} instead of {} )'.format(len(xdot0),self.nDOF))
            self.q0[self.nDOF:] =  xdot0


    def setStateInitialConditions(self,q0):
        q0 = np.asarray(q0).ravel()
        if len(q0)!=self.nStates:
            raise Exception('Wrong dimension for q0 ({} instead of {} )'.format(len(q0),2*self.nStates))
        self.q0 = q0

    def setForceTimeSeries(self,vTime,vF):
        """ 
        Set the force on the system degrees of freedom as a time series of time
        INPUTS:
         - vTime: 1d array of time steps (do not need to be regular), of length nt
         - vF   : nDOF x nt array of forces at each time steps, for the nDOF degrees of freedom
        """
        vTime = np.asarray(vTime)
        vF    = np.asarray(vF)
        if self.nDOF==1:
            vF = vF.reshape((1,-1))
        if vF.shape[0]!=self.nDOF:
            raise Exception('Wrong first dimension for Force time series ({} instead of {} )'.format(vF.shape[0],self.nDOF))
        if vF.shape[1]!=len(vTime):
            raise Exception('Second dimension of Force time series does not match time dimension ({} instead of {} )'.format(vF.shape[1],len(vTime)))

        # Store raw data
        self._force_ts = vF
        self._time_ts  = vTime
        # Create interpolant for faster evaluation
        self._force_fn_t = interp1d(vTime, vF)
        if not self.M_is_func:
            # Create interpolant for Minv*B directly
            MinvF = self._Minv.dot(vF)
            self._Minv_force_fn_t = interp1d(vTime, MinvF)

    def setForceFunction(self,fn):
        """ 
        Set the force on the system degrees of freedom as a function of time, x and xdot.
        The function will be used during the time integration to evaluate the forces at various times and states.

        INPUTS:
          fn: handle to a python function. The interface of the function fn is: 

               F =  f(t, x, xdot, **kwargs) 

              where:
                 t   : scalar, time
                 x   : (nDOF,) array, position of each DOF
                 xdot: (nDOF,) array, velocities of each DOF
                 F   : (nDOF x 1) array, forces in each DOF
        
        """
        self._force_fn = fn

    def Forces(self, t, x=None, xdot=None, q=None, **kwargs):
        """ 
        Return Forces on degrees of freedom at one time step, or multiple time steps
        The force is obtained using either:
          - the function provided by the user using setForceFunction (a nonlinear function of states)
          - the time series provided by the user using setForceTimeSeries
        Multiple time steps are only supported when time series were provided

        INPUTS:
          - t: scalar, or time array 

        OPTIONAL INPUTS:
          - x, xdot: array-like of length  nDOF for the DOF position and velocities
         OR
          - q      : array-like of length  2 nDOF for the states (position and velocities)

          - **kwargs: keyword argument passed to th force function

        OUTPUTS:
          - F: nDOF x nt array where nt=len(t)
        
        """
        if hasattr(self,'_force_ts'):
            return self._force_fn_t(t) # NOTE: interpolant, works for scalar or arrays

        elif hasattr(self,'_force_fn'):
            if q is not None:
                x    = q[0:self.nDOF]
                xdot = q[self.nDOF:]
            return self._force_fn(t, x, xdot, **kwargs)
        else:
            raise Exception('Please specify a time series of force using `setForceTimeSeries` or a function using `setForceFunction` ')

    def B_tilde(self,t,q):
        """ 
        "B~" is the non-linear right hand side
            xdot = A~ x + B~
        """
        nDOF = self.nDOF
        x    = q[0:nDOF]
        xd   = q[nDOF:]

        if hasattr(self,'_Minv_force_fn_t'):
            # Simplest case, force is a function of time only, and M is constant
            # We use the interpolant for M^-1 F to determine B
            B = np.zeros((2*nDOF,1))
            B[nDOF:,0] = self._Minv_force_fn_t(t)
            return B

        elif hasattr(self,'_force_ts'):
            # Force is a function of time only, but M may be varying
            # We use the interpolant for F, and compute M^-1 F to determine B
            B = np.zeros((2*nDOF,1))
            B[nDOF:,0] = self._fMinv(x).dot(self._force_fn_t(t))
            return B

        elif hasattr(self,'_force_fn'):
            # Force is a function of time and states, M may be varying
            F          = self._force_fn(t,x,xd).ravel()
            B          = np.zeros((2*nDOF,1)) 
            B[nDOF:,0] = np.dot(self._fMinv(x),F)
            return B
        else:
            raise Exception('Please specify a time series of force using `setForceTimeSeries` or a function using `setForceFunction` ')


    def integrate(self,t_eval, method='RK45', y0=None, **options):
        """ Perform time integration of system 
            method: 'RK54', 'LSODA'  (see solve_ivp)
        """
        #
        if y0 is not None:
            self.setStateInitialConditions(y0)

        if method.lower()=='duhamel':
            # --- Using Duhamel integral method
            r""" 
            x(t)    = \int_0^t F(t') H(t-t') dt'    F: force, H: impulse response function
            xdot(t) = \int_0^t F(t') H'(t-t') dt' , H' derivative of impulse function
            """
            if self.nDOF!=1:
                raise Exception('Duhamel integral only available for 1DOF system')
            if not hasattr(self,'_force_ts'):
                raise Exception('Duhamel integral only available when force is provided using `setForceTimeSeries`')
            # TODO add check of initial conditions

            from .singledof import duhamel
            m = self.M.ravel()[0]
            c = self.C.ravel()[0]
            k = self.K.ravel()[0]
            F = self.Forces(t_eval)
            x, xd = duhamel(t_eval, F, m, c, k, both=True)
            res = OdeResultsClass(t=t_eval, y=np.vstack((x,xd))) # To mimic result class of solve_ivp
                
        else:
            # Using time integration methods of solve_ivp
            if self.has_A:
                A=self.A
                odefun = lambda t, q : np.dot(A,q)+self.B_tilde(t,q)
            else:
                odefun = self.dqdt
                #odefun = lambda t, q : self.B(t,q)
            res = solve_ivp(fun=odefun, t_span=[t_eval[0], t_eval[-1]], y0=self.q0, t_eval=t_eval, method=method, vectorized=True, **options)   
        # Store
        self.res    = res
        return res


    def dqdt(self, t, q):
        if self.has_A:
            q=q.reshape((-1,1))
            A=self.A_tilde
            dqdt= np.dot(A,q)+self.B_tilde(t,q)
            return dqdt
        else:
            dqdt_ = np.zeros(q.shape)
            x  = q[0:self.nDOF]
            xd = q[self.nDOF:]

            dqdt_[0:self.nDOF]=  xd
            Minv = self._fMinv(x)
            F = self._force_fn(t,x,xd)
            dqdt_[self.nDOF:] =  Minv.dot(F)
            return dqdt_

    @property
    def A_tilde(self):
        """ 
          "A~" is the linear state matrix made using M, C, K (assuming they are constant)
            xdot = A~ x + B~
          NOTE: B~ may be a function of x, linearization is not done here 
        """
        if self.M_is_func:
            raise Exception('A matrix not a property when M is a function')
        else:
            return StateMatrix(self._Minv,self.C,self.K)

    # --------------------------------------------------------------------------------}
    # --- "Linear"
    # --------------------------------------------------------------------------------{
    # TODO consider using system.linearization.py to get linear model
    @property
    def A(self):
        if self.M_is_func:
            raise Exception('A matrix not a property when M is a function')
        if hasattr(self,'_force_fn') or not hasattr(self,'_force_ts'):
            print('[WARN] The A matrix is not guaranteed to be the linearized version when force is given as a function')
        return StateMatrix(self._Minv,self.C,self.K)

    @property
    def B(self):
        """ 
        Returns B when applicable
            xdot = A x + B u 
        where u is assumed to be the forces in each DOF
        """
        if self.M_is_func:
            raise Exception('B matrix not a property when M is a function')
        if hasattr(self,'_force_fn') or not hasattr(self,'_force_ts'):
            print('[WARN] The B matrix is not guaranteed to be the linearized version when force is given as a function')
        nDOF = self.nDOF
        B = np.zeros((2*nDOF,nDOF))
        B[nDOF:,:] = self._Minv
        return B


    # --------------------------------------------------------------------------------}
    # ---  IO functions for printing/plotting/saving
    # --------------------------------------------------------------------------------{
    def __repr__(self):
        s='<{} object>\n'.format(type(self).__name__)
        s+='|Read-only attributes:\n'
        s+='| - nDOF: {}, nState:{} \n'.format(self.nDOF,self.nStates)
        if hasattr(self,'_force_ts'):
            if len(self._time_ts)>1:
                dt=self._time_ts[1]-self._time_ts[0]
            else:
                dt=np.nan
            s+='|Force time series \n'
            s+='| - Time: [{} ... {}],  dt: {}, n: {} \n'.format(self._time_ts[0],self._time_ts[-1],dt,len(self._time_ts))
            s+='| - Force t0  : {} \n'.format(self._force_ts[:,0])
            s+='| - Force tend: {} \n'.format(self._force_ts[:,-1])
        elif hasattr(self,'_force_ts'):
            s+='|Force function {} \n'.format(self._force_fn)
        else:
            s+='|Force not defined yet\n'
        s+='|Attributes:\n'
        s+='| - M: Mass Matrix  \n'
        s+=str(self.M)+'\n'
        s+='| - C: Damping Matrix  \n'
        s+=str(self.C)+'\n'
        s+='| - K: Stiffness Matrix  \n'
        s+=str(self.K)+'\n'
        s+='| - q0: Initial conditions (state) \n'
        s+=str(self.q0)+'\n'
        if self.has_A:
            s+='| - A:  \n'
            s+=str(self.A)+'\n'
        try:
            s+='| - B:  \n'
            s+=str(self.B)+'\n'
        except:
            pass

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

        axes = np.atleast_1d(axes)
        n=self.nDOF
        for i,ax in enumerate(axes):
            if i<n:
                lbl=r'$x_{}$'.format(i+1)
            else:
                lbl=r'$\dot{x}_'+r'{}$'.format(i-n+1)
            ax.plot(res.t, res.y[i,:], label=label, **kwargs)
            ax.set_ylabel(lbl)
            ax.tick_params(direction='in')
        axes[-1].set_xlabel('Time [s]')

        return fig, axes

    def plot_forcing(self,axes=None, label=None, res=None, **kwargs):
        """ 
        Simple plot of forcing
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
            if hasattr(self,'_force_ts'):
                time = self._time_ts
                F    = self._force_ts
            else:
                raise Exception('Cannot plot forcing when a function is used and no time integration was performed. Call `integrate`.')
        else:
            # Res provided
            time =res.t
            if hasattr(self,'_force_ts'):
                F = self.Forces(time)
            else:
                F = np.zeros((self.nDOF,len(time)))
                for it, t in enumerate(time):
                    F[:,it] = self.Forces(t, q=res.y[:,it])
        # Plot
        if axes is None:
            import matplotlib.pyplot as plt
            fig,axes = plt.subplots( self.nDOF, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
            fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
        axes = np.atleast_1d(axes)
        for i,ax in enumerate(axes):
            lbl=r'$f_{}$'.format(i+1)
            ax.plot(time, F[i,:], label=label, **kwargs)
            ax.set_ylabel(lbl)
            ax.tick_params(direction='in')
        axes[-1].set_xlabel('Time [s]')

        return fig, axes


    def save(self, filename, DOFs=None, Factors=None):
        """ Save time integration results to file 
        DOFs: array of string for DOFs names
        """
        df = self.toDataFrame(DOFs=DOFs, Factors=Factors)
        df.to_csv(filename, sep=',', index=False)

    def toDataFrame(self, DOFs=None, Factors=None):
        """ Return time integration results as a dataframe
        DOFs:    array of string for DOFs names   (2xnDOF)
        Factors: array of floats to scale the DOFs (2xnDOF)
        """
        import pandas as pd
        if self.res is None:
            raise Exception('Call integrate before save')

        if DOFs is None:
            DOFs  = ['x_{}'.format(i) for i in np.arange(self.nDOF)]
            DOFs += ['xd_{}'.format(i) for i in np.arange(self.nDOF)]
        header = ' Time_[s], '+','.join(DOFs)

        res = self.res
        M=np.column_stack((res.t, res.y.T))

        # Scaling
        if Factors is not None:
            for i, f in enumerate(Factors):
                M[:,i+1] *= f

        for i,d in enumerate(DOFs):
            if DOFs[i].find('[deg]')>1: 
                if np.max(M[:,i+1])>180:
                    M[:,i+1] = np.mod(M[:,i+1], 360)

        return pd.DataFrame(data=M, columns=['Time_[s]']+DOFs)



