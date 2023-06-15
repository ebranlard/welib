import numpy as np
import pandas as pd
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
    def __init__(self, M, C=None, K=None, F=None, x0=None, xdot0=None, sX=None, sXd=None, sY=None):
        """ 
        sX:  list of variable names for x    (e.g., used for plotting)
        sXd: list of variable names for xdot (e.g., used for plotting)
        sY: list of variable names for outputs. Can also be set using setOutputFunction
        """

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
                self.setForceTimeSeries(F[0], F[1]) # Time vector, Force vector
            else:
                self.setForceFunction(F)

        # Initial conditions
        self.setInitialConditions(x0,xdot0)

        # Time integration results
        self.res=None

        # Channel names
        self._sX = sX
        self._sXd = sXd
        self.sY = None

        self.verbose=False


    def find_nDOF(self):
        if self.M_is_func:
            x=np.zeros((5000,1)) # Using a huge number of states..
            return len(self.M(x))
        else:
            return len(self.M)
    @property
    def nStates(self):
        return 2*self.nDOF

    @property
    def sX(self):
        if self._sX is None:
            return [r'x_{}'.format(i+1) for i in range(self.nDOF)]
        else:
            return self._sX

    @property
    def sXd(self):
        if self._sXd is None:
            return [r'\dot{x}'+r'_{}'.format(i+1) for i in range(self.nDOF)]
        else:
            return self._sXd

    @property
    def sQ(self):
        return self.sX+self.sXd

    @property
    def sOutputs(self):
        if self.sY is not None:
            return self.sY
        else:
            return None

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

    def setOutputFunction(self, fn, sY=None):
        """

        INPUTS:
         - fn: handle to a python function. The interface of the function fn is: 

               y =  f(t, x, xdot, **kwargs) 

            where y is an array of size n or a pandas series with index corresponding to the output names 
         - sY : array-like, list of of column names for outputs
        """
        self._output_fn = fn
        if sY is not None:
            self.sY = sY

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
            if len(F)>0: 
                B[nDOF:,0] = np.dot(self._fMinv(x),F)
            else:
                # If no input for instance or not forcing F might be empty array..
                # TODO DO THIS BETTER, very at init or something
                pass
            return B
        else:
            raise Exception('Please specify a time series of force using `setForceTimeSeries` or a function using `setForceFunction` ')

    # --------------------------------------------------------------------------------}
    # --- OUTPUTS
    # --------------------------------------------------------------------------------{
    def Outputs(self, t, q):
        nDOF = self.nDOF
        x    = q[0:nDOF]
        xd   = q[nDOF:]
        if hasattr(self,'_output_fn'):
            return self._output_fn(t, x, xd)
        else:
            raise NotImplementedError('Calculation of outputs when `setOutputFunction` not called')

    def dqdt_calcOutput(self, signatureWanted='t,q', **kwargs):
        """ Return function handle that computes outputs with a requested signature """
        if signatureWanted=='t,q':
            return self.Outputs
        else:
            raise NotImplementedError()


    # --------------------------------------------------------------------------------}
    # --- Time integration 
    # --------------------------------------------------------------------------------{
    def integrate(self,t_eval, method='RK45', y0=None, calc='', **options):
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

        df=self.res2DataFrame(calc=calc)
        return res, df


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
            F = self._force_fn(t,x,xd).flatten() # NOTE: something is weird here if no "A" but timeseires function 
            try:
                dqdt_[self.nDOF:] =  Minv.dot(F)
            except:
                dqdt_[self.nDOF:,0] =  Minv.dot(F)
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
    # --- Simulation storage
    # --------------------------------------------------------------------------------{
    def cleanSimData(self):
        self.res      = None
        #self.dfIn     = None
        self.dfOut    = None
        #self.dfStates = None
        self.df       = None

    def save(self, filename, DOFs=None, Factors=None):
        """ Save time integration results to file 
        DOFs: array of string for DOFs names
        """
        df = self.res2DataFrame(DOFs=DOFs, Factors=Factors)
        df.to_csv(filename, sep=',', index=False)

    def res2DataFrame(self, res=None, sStates=None, Factors=None, 
            x0=None, xd0=None, 
            calc='xdd,f', # xdd, f, y
            sAcc=None, sForcing=None
            ):
        """ Return time integration results as a dataframe
        sStates: array of string for DOFs names   (2xnDOF, positions and velocities)
        Factors: array of floats to scale the DOFs (2xnDOF)
        x0 :  array of floats to add to the DOF (1xnDOF)
        xd0:  array of floats to add to the velocities (1xnDOF)
              NOTE: the positions will be increased linearly
        """
        # TODO harmonize with statespace.res2DataFrame
        import pandas as pd

        calcVals = calc.split(',')

        if res is None:
            if self.res is None:
                raise Exception('Call integrate before res2DataFrame')
            res = self.res

        if sStates is None:
            sStates  = ['x_{}'.format(i) for i in np.arange(self.nDOF)]
            sStates += ['xd_{}'.format(i) for i in np.arange(self.nDOF)]
        if sAcc is None:
            sAcc = ['xdd_{}'.format(i) for i in np.arange(self.nDOF)]
        if sForcing is None:
            sForcing = ['F_{}'.format(i) for i in np.arange(self.nDOF)]

        # --- Time and states
        dfStates = self.store_states(res, sStates=sStates, x0=x0, xd0=xd0, Factors=Factors)

        # --- Accelerations 
        dfAcc = None
        if 'xdd' in calcVals:
            xdd = self.TS_Acceleration
            dfAcc =  pd.DataFrame(data=xdd.T, columns=sAcc)

        # --- Forcing
        dfF = None
        if 'f' in calcVals:
            F     = self.TS_Forcing
            dfF =  pd.DataFrame(data=F.T, columns=sForcing)

        # --- Try to compute outputs
        dfOut = None
        if 'y' in calcVals:
            dfOut = self.calc_outputs(insertTime=True, yoffset=None)
        # --- Try to compute inputs

        # --- Concatenates everything into one DataFrame
        df = pd.concat((dfStates, dfAcc, dfF, dfOut), axis=1)
        df = df.loc[:,~df.columns.duplicated()].copy()

        return df

    def store_states(self, res, sStates=None, x0=None, xd0=None, Factors=None):
        #nStates = len(self.q0)
        #if sStates is None and self.sX is None:
        #    sStates = ['x{}'.format(i+1) for i in range(nStates)] # TODO
        #else:
        #    if sStates is None:
        #        sStates = self.sX
        #    nCols = self.res.y.shape[0]
        #    if len(self.sX)!=nCols:
        #        raise Exception("Inconsistency in length of states columnNames. Number of columns detected from res: {}. States columNames (sX): {}".format(nCols, self.sX))
        #self.sX = sStates

        # Store as a matrix
        M = res.y.T.copy()

        # Scaling
        if Factors is not None:
            for i, f in enumerate(Factors):
                M[:,i] *= f
        # Scaling offsets
        if Factors is not None and x0 is not None:
            x0 *= np.asarray(Factors[:self.nDOF])
        if Factors is not None and xd0 is not None:
            xd0 *= np.asarray(Factors[self.nDOF:])

        # Offset Velocity
        if xd0 is not None:
            for i, xd0_ in enumerate(xd0):
                M[:,self.nDOF+i] += xd0_       # Add to velocity
                M[:,i]           += xd0_*res.t # Position increases linearly (e.g. azimuth)
        # Offset position
        if x0 is not None:
            for i, x0_ in enumerate(x0):
                M[:,i] += x0_

        # Degrees
        for i,d in enumerate(sStates):
            if sStates[i].find('[deg]')>1: 
                if np.max(M[:,i])>180:
                    M[:,i] = np.mod(M[:,i], 360)
        # Insert time
        dfStates =  pd.DataFrame(data=M, columns=sStates)
        dfStates.insert(0, 'Time_[s]', res.t)
        return dfStates




    def _prepareOutputDF(self, res):
        cols    = self._inferOutputCols(res)
        data    = np.full((len(res.t), len(cols)), np.nan)
        df      = pd.DataFrame(columns = cols, data = data)
        self.sY = cols
        return df

    def _inferOutputCols(self, res):
        """ See what the calcOutput returns at t[0] to allocate storage memory """
        try:
            out = self.dqdt_calcOutput()(res.t[0], self.q0)
        except TypeError:
            print("[FAIL] Error evaluating model outputs. Does the function dqdt has the argument `calcOutput`?")
            return None
        if not isinstance(out, pd.core.series.Series):
            # If array
            cols = ['y{:d}'.format(i+1) for i in range(len(out))]
        else:
            # If pandas series
            cols     = out.index

        # --- Check consistency with self.sY
        if self.sY is not None:
            if len(self.sY)!=len(cols):
                raise Exception("Inconsistency in length of output columnnames. Number of columns detected from `calcOuputt`: {}. Ouput columNames (sY):".format(len(cols), self.sY))
            cols = self.sY
        return cols

    def _calc_outputs(self, time, q, df):
        """ low level implementation leaving room for optimization for other subclass."""
        calcOutput = self.dqdt_calcOutput()
        if self.verbose:
            print('Calc output...')
        for i,t in enumerate(time):
            df.iloc[i,:] = calcOutput(t, q[:,i])

    def calc_outputs(self, res=None, insertTime=True, dataFrame=True, yoffset=None):
        """ 
        Call use calcOutput function for each time step and store values
        """
        # --- Get calcOutput function
        
        if res is None:
            res = self.res
        if res is None:
            raise Exception('Provide `res` or call `integrate` before calling calc_outputs')

        # --- Infer column names and allocate
        df = self._prepareOutputDF(res)

        # --- Calc output based on states
        self._calc_outputs(res.t, res.y, df)
        if yoffset is not None:
            raise Exception()

        if insertTime:
            df.insert(0,'Time_[s]', res.t)
        self.dfOut = df

        if dataFrame:
            return df
        else:
            return df.values


    # --------------------------------------------------------------------------------}
    # --- Time series after time integration
    # --------------------------------------------------------------------------------{
    @property
    def time(self):
        if self.res is None:
            # Res or self.res not provided, return time series provided by user
            if hasattr(self,'_force_ts'):
                return self._time_ts
            else:
                raise Exception('Cannot return time when a function is used and no time integration was performed. Call `integrate`.')
        else:
            # Res provided
            return self.res.t


    @property
    def TS_Forcing_CK(self):
        """ Return time series of forcing from Stiffness and Damping matrix """
        if self.res is None:
            raise Exception('Run `integrate` before calling CKforcing')
        time = self.time
        nDOF = self.nDOF
        FK = np.zeros((self.nDOF,len(time)))
        FC = np.zeros((self.nDOF,len(time)))
        if self.K is not None:
            for it, t in enumerate(time):
                q  = res.y[:,it]
                x  = q[0:nDOF]
                xd = q[nDOF:]
                FK[:,it] = -self.K.dot(x)
        if self.C is not None:
            for it, t in enumerate(time):
                q  = res.y[:,it]
                x  = q[0:nDOF]
                xd = q[nDOF:]
                FC[:,it] = -self.C.dot(xd)
        return FC, FK

    @property
    def TS_Forcing(self):
        if self.res is None:
            # Res or self.res not provided, return time series provided by user
            if hasattr(self,'_force_ts'):
                time = self._time_ts
                F    = self._force_ts
            else:
                raise Exception('Cannot compute forcing when a function is used and no time integration was performed. Call `integrate`.')
        else:
            # Res provided
            time =self.res.t
            res =self.res
            if hasattr(self,'_force_ts'):
                F = self.Forces(time)
            else:
                F = np.zeros((self.nDOF,len(time)))
                for it, t in enumerate(time):
                    F[:,it] = self.Forces(t, q=res.y[:,it]).ravel()
        return F

    @property
    def TS_Acceleration(self):
        """ Return time series of acceleration
        NOTE: some optimizations are possible here, if forcing is known for instance, or not computing 
        the "top" part of dqdt..
        """
        if self.res is None:
            raise Exception('Run `integrate` before calling TS_Acceleration')
        nDOF = self.nDOF
        xdd = np.zeros((self.nDOF,len(self.res.t)))
        for it, t in enumerate(self.res.t):
            dqdt=self.dqdt(t, self.res.y[:,it])
            xdd[:,it] = dqdt.flatten()[nDOF:]
        return xdd

    # --------------------------------------------------------------------------------}
    # --- "Linear"
    # --------------------------------------------------------------------------------{
    # TODO consider using system.linearization.py to get linear model
    @property
    def A(self):
        if self.M_is_func:
            raise Exception('A matrix not a property when M is a function')
        #if hasattr(self,'_force_fn') or not hasattr(self,'_force_ts'):
        #    print('[WARN] The A matrix is not guaranteed to be the linearized version when force is given as a function')
        return StateMatrix(self._Minv,self.C,self.K)

    def eigA(self):
        from welib.tools.eva import eigA
        A = self.A # NOTE: not always defined...
        freq_d, zeta, Q, freq_0 = eigA(A, fullEV=False, normQ=None, sort=True)
        return freq_d, zeta, Q, freq_0 

    @property
    def B(self):
        """ 
        Returns B when applicable
            xdot = A x + B u 
        where u is assumed to be the forces in each DOF
        """
        if self.M_is_func:
            raise Exception('B matrix not a property when M is a function')
        #if hasattr(self,'_force_fn') or not hasattr(self,'_force_ts'):
        #    print('[WARN] The B matrix is not guaranteed to be the linearized version when force is given as a function')
        nDOF = self.nDOF
        B = np.zeros((2*nDOF,nDOF))
        B[nDOF:,:] = self._Minv
        return B

    @property
    def x0(self): return self.q0[:self.nDOF]
    @property
    def xd0(self): return self.q0[self.nDOF:]

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
        s+='| * sX:  {} \n'.format(self.sX)
        s+='| * sXd: {} \n'.format(self.sXd)
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

    def picklable(self):
        """ Make the object picklable..."""
        def noneIfLambda(attr):
            obj = getattr(self, attr)
            if callable(obj) and obj.__name__ == "<lambda>":
                print('MechSystem: picklable: removing ', attr)
                setattr(self, attr, None)
        noneIfLambda('M')
        noneIfLambda('_fM')
        noneIfLambda('_fMinv')
        noneIfLambda('_force_fn')


    def plot(self, fig=None, axes=None, label=None, res=None, calc='', qFactors=None, **kwargs):
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
        n=self.nDOF

        df=self.res2DataFrame(sStates=None, Factors=qFactors, x0=None, xd0=None, calc=calc, sAcc=None, sForcing=None)
        for i,ax in enumerate(axes):
            if i+1>len(df.columns):
                continue
            chan=df.columns[i+1]
            lbl = '$'+self.sQ[i]+'$'
            ax.plot(res.t, res.y[i,:], label=label, **kwargs)
            if not axes_provided:
                ax.set_ylabel(lbl)
            ax.tick_params(direction='in')
        if not axes_provided:
            axes[-1].set_xlabel('Time [s]')

        return axes

    def plot_outputs(self, df=None, keys=None, axes=None, **kwargs):
        """ 
        plot outputs. Requires to call `integrate` or `calc_outputs` before
        """
        if df is None:
            if self.dfOut is None:
                raise Exception("Call calc_outputs or `integrate` before plot_outputs")
            df = self.dfOut
            if df is None:
                print('[WARN] no outputs to plot')
                return
        keys = list(self.sOutputs)

        return _plot(df, keys=keys, title='', nPlotCols=1, axes=axes, **kwargs)


    def plot_forcing(self, axes=None, label=None, res=None, includeCK=False, plotCK0=False, **kwargs):
        """ 
        Simple plot of forcing
          - if `res` is provided, use time and states from res
          - if time integration was performed, use time and states from self.res
          otherwise: 
          - if a time series was provided, plot that
          - if a function was provided for the forcing, abort
        """
        # 
        time = self.time
        # Compute Forcing 
        F = self.TS_Forcing
        # Compute Forcing  contributions from Damping and Stiffness
        if includeCK:
            FC, FK = self.TS_Forcing_CK
            F0=F.copy()
            F=F0+FK+FC

        # Plot
        if axes is None:
            import matplotlib.pyplot as plt
            fig,axes = plt.subplots( self.nDOF, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
            fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
        axes = np.atleast_1d(axes)
        for i,ax in enumerate(axes):
            lbl=r'$f_{}$'.format(i+1)
            ax.plot(time, F[i,:], label=label, **kwargs)
            if plotCK0:
                ax.plot(time, F0[i,:], ':' , label='(F0)', **kwargs)
                ax.plot(time, FC[i,:], '-' , label='(FC)', **kwargs)
                ax.plot(time, FK[i,:], '--', label='(FK)', **kwargs)
            ax.set_ylabel(lbl)
            ax.tick_params(direction='in')
        axes[-1].set_xlabel('Time [s]')

        return axes



def _plot(df, keys=None, label=None, title='', nPlotCols=1, axes=None, **kwargs):
    import matplotlib
    import matplotlib.pyplot as plt
    #if COLRS is None:
    #    cmap = matplotlib.cm.get_cmap('viridis')
    #    COLRS = [(cmap(v)[0],cmap(v)[1],cmap(v)[2]) for v in np.linspace(0,1,3+1)]

    time = df['Time_[s]'].values
    if keys is None:
        keys = [k for k in df.keys() if k!='Time_[s]']

    I=np.arange(len(keys))

    if axes is None:

        if nPlotCols==2:
            fig,axes = plt.subplots(int(np.ceil(len(I)/2)), 2, sharex=True, figsize=(6.4,4.8)) # (6.4,4.8)
            fig.subplots_adjust(left=0.07, right=0.98, top=0.955, bottom=0.05, hspace=0.20, wspace=0.20)
        else:
            fig,axes = plt.subplots(len(I), 1, sharex=True, figsize=(6.4,4.8)) # (6.4,4.8)
            fig.subplots_adjust(left=0.16, right=0.95, top=0.95, bottom=0.12, hspace=0.20, wspace=0.20)

        if not hasattr(axes,'__len__'):
            axes=[axes]
    axes=(np.asarray(axes).T).ravel()
    
    for j,i in enumerate(I):
        s  = keys[i]
        ax = axes[j]
        ax.plot(time, df[s], label=label, **kwargs)
        ax.set_ylabel(s)
        ax.tick_params(direction='in')
    axes[0].set_title(title)
    axes[-1].set_xlabel('Time [s]')
    axes[-1].legend()
    return axes


