import numpy as np
import pandas as pd
import inspect
from numpy.linalg import inv, solve
from collections import OrderedDict
from scipy.integrate import  solve_ivp, odeint
from scipy.optimize import OptimizeResult as OdeResultsClass 
from scipy.interpolate import interp1d


from welib.system.system import System

# --------------------------------------------------------------------------------}
# --- Functions for state space model integrations
# --------------------------------------------------------------------------------{
def StateMatrix(Minv=None,C=None,K=None,M=None):
    if M is not None:
        Minv=inv(M)
    if C is None:
        C=np.zeros(K.shape)
    n = len(Minv)
    A = np.zeros((2*n,2*n))
    A[:n,n:] =  np.identity(n)
    A[n:,:n] = -np.dot(Minv,K)
    A[n:,n:] = -np.dot(Minv,C)
    return A

def vec_interp(t,vTime,vF):
    """ Interpolate vector known at discrete values (vTime, vF) to a given time `t` 
    TODO use and interpolant if possible
      t : scalar!
      vTime : 1d-array of length nt 
      vF : nDOF x nt array
    """
    F    = np.zeros(vF.shape[0])
    for iDOF,F_DOF in enumerate(vF):
        F[iDOF] = np.interp(t,vTime,F_DOF)
    return F

def B_interp(t,Minv,vTime,vF, flat=False):
    """ Interpolate B-vector from loads known at discrete values (vTime, vF) at a given time `t` """
    nDOF=len(vF)
    if flat:
        B = np.zeros(2*nDOF)
        F = vec_interp(t,vTime,vF)
        B[nDOF:] = np.dot(Minv,F)
    else:
        B = np.zeros((2*nDOF,1))
        F = vec_interp(t,vTime,vF)
        B[nDOF:,0] = np.dot(Minv,F)
    return B

def B_reg(t,Minv,F):
    """ Return B vector from loads at time t and mass matrix """
    nDOF=len(F)
    B = np.zeros((2*nDOF,1))
    B[nDOF:,0] = np.dot(Minv,F)
    return B

def dxdt(q, t, A, M, vTime, vF): 
    B = B_interp(t, M, vTime, vF)
    dq = np.dot(A,q)+B
    return dq
 
def odefun(t, dq):
    return dxdt(dq, t, A, vTime, vF)



# --------------------------------------------------------------------------------}
# ---  
# --------------------------------------------------------------------------------{
class StateSpace(System):
    """ 

    Handles state space system made of a state and output equations:

      dot(q) = dqdt (t, q , [u, p, calcOutput)
         y   = calcOutput(t, q, [u, p] )


    Inputs u: different options:
      1) dictionary of interpolant:  u[key](t) = scalar
      2) function/interpolant returning an array fU(t)   = array of size nu
      3) function/interpolant with signature     fU(t,q) = array of size nu

      Option 1 is set using setInputFunctionDict
      Option 2 can be set using setInputTimeSeries(vTime, vU)
      Option 2 can be set using setInputFunction(fU, signature_u='t')
      Option 3 can be set using setInputFunction(fU, signature_u='t,q')


    def setStateInitialConditions(self,q0=None):
    def setInputTimeSeries(self,vTime,vU):
    def setInputFunction(self,fn):
    def Inputs(self,t,q=None):
    def integrate(self, t_eval, method='RK4', y0=None, **options):
    def dqdt(self, t, q):
    def RHS(self,t,q):
    def nStates(self):
    def nInputs(self):
    def nOuputs(self):
    """
    def __init__(self, dqdt, signature, q0=None, p=None, u=None, 
            sX=None, sU=None, sY=None,
            verbose=False):
        """ 
        INPUTS:
         - dqdt : handle that returns the time derivative of the state equation
         - signature: signature/interface of `dqdt`
            - 't,q,u,p'
            - 't,q,p,u'
            - 't,q,p'
            - 't,q,u'
            - 't,q,'
         - signature_u: 't', 't,q'
         - sX: list of names for states
         - sY: list of names for outputs
         - sU: list of names for inputs
        """
        # TODO Needs signature output!!! and better output function handling

        signature   = signature.replace('x','q') # legacy...
        signature_x = signature.replace('q','x') # legacy...

        System.__init__(self, dqdt, interface=signature_x.replace(',',''))

        # Data
        self.verbose      = verbose # Must come first
        self._signature   = signature
        self.p            = p
        self.dqdt         = dqdt
        self._u           = None
        self._signature_u = None      # TODO remove
        # --- Names
        self.sX = sX
        self.sY = sY
        self.sU = sU

        # Data with potential triggers
        if u is not None:
            self.u  = u

        # Initial conditions
        self.setStateInitialConditions(q0)

        # Time integration results
        self.res=None

        self.dfStates = None
        self.dfIn     = None
        self.dfOut    = None

    @property
    def nStates(self):
        # NOTE: really bad...
        if self.sX is None:
            raise Exception()
        return len(self.sX) # Should work with dict and array

    @property
    def nInputs(self):
        if self._u is None or self.sU is None:
            raise Exception()
        return len(self._u) # Should work with dict and array
# 
#     @property
#     def nOutputs(self):
#         if self.C is not None:
#             return self.C.shape[1]
#         else:
#             return 0
    @property
    def nParams(self):
        if self.p is None:
            raise Exception()
        return len(self.p) # Should work with dict and array

    # TODO agree on an interface and stick to it...
    @property
    def sStates(self):
        if self.sX is not None:
            return self.sX
        else:
            return None

    @property
    def sInputs(self):
        if self.sU is not None:
            return self.sU
        else:
            return None

    @property
    def sOutputs(self):
        if self.sY is not None:
            return self.sY
        else:
            return None

    @property
    def q0_pd(self):
        return pd.Series(self.q0, index=self.sStates)

    # --------------------------------------------------------------------------------}
    # --- Signatures
    # --------------------------------------------------------------------------------{
    @property
    def signature(self):
        return self._signature
    @signature.setter
    def signature(self, sig):
        sigAllowed=['t,d', 't,q,p', 't,q,u,p', 't,q,u,p', None]
        if sig not in sigAllowed:
            raise Exception('Signature needs to be one of the following: {}'.format(sigAllowed))
        self._signature = sig

    @property
    def signature_u(self):
        return self._signature_u

    @signature_u.setter
    def signature_u(self, sig):
        sigAllowed=['t', 't,q', 't,q,qd',None]
        if isinstance(sig, dict):
            for k in sig.keys():
                if sig[k] not in sigAllowed:
                    raise Exception('Signature_u needs to be one of the following: {}'.format(sigAllowed))

        else:
            if sig not in sigAllowed:
                raise Exception('Signature_u needs to be one of the following: {}'.format(sigAllowed))
        #print('Setting signature_u to',sig)
        self._signature_u = sig

    def uDependsOnTOnly(self):
        """ returns true if inputs only depend on time"""
        if isinstance(self._signature_u, dict):
            return np.all([sig == 't' for _,sig in self._signature_u.items()])
        else:
             return self._signature_u == 't'


    @property
    def u(self):
        return self._u

    @u.setter
    def u(self, u):
        """ Sets input function or dictionary of functions
        Attempts to finds the correct signature for the input function.
        """
        self.__setU(u)

    def __setU(self, u, vTime=None):
        """
        Low level function, to merge:
          - self.u =u
        and 
          - self.setInputTimeSeries(vTime, vU)
        """
        if isinstance(u, dict):
            u= OrderedDict(u)

        # --- Remove any previously stored time series inputs
        self._inputs_ts = None
        self._time_ts   = None

        # --- Detect signature of inputs 
        if u is None:
            self.signature_u = None
        elif isinstance(u, OrderedDict):
            # --- Setting signatures for each inputs individually
            sigDict={}
            for k in u.keys():
                sigParams = inspect.signature(u[k]).parameters
                try:
                    # Try to call with just time
                    u0 = u[k](0)
                    sigDict[k] = 't'
                    if len(sigParams)==1:
                        sigDict[k] = 't'
                    elif len(sigParams)>=1:
                        #print('u[{}](t) works has more than one param'.format(k), sigParams)
                        sigDict[k] = ['t','t,q', 't,q,qd'][len(sigParams)-1] # TODO
                except:
                    if len(sigParams)<=1:
                        raise # Something might be wrong in the function
                    elif len(sigParams)==2:
                        sigDict[k] = 't,q'
                    elif len(sigParams)==3:
                        sigDict[k] = 't,q,qd'
                    elif len(sigParams)>=3:
                        print('u[{}](*) has more than three params'.format(k), sigParams)
                        sigDict[k] = 't,q,qd' # TODO
            self.signature_u=sigDict

        else:
            if vTime is None:
                if type(u) is np.ndarray:
                    raise Exception('Cannot call setU without a time vector. use `setInputTimeSeries`')
                # We assume it's a function
                sigParams = inspect.signature(u).parameters
                if len(sigParams)==1:
                    self.signature_u = 't'
                elif len(sigParams)==2:
                    self.signature_u = 't,q'
                elif len(sigParams)==3:
                    self.signature_u = 't,q,qd'
                else:
                    print('[WARN] u has more than three params', sigParams)
                    self.signature_u = 't,q,qd'
            else:
                # This is set by setInputTimeSeries
                vTime = np.asarray(vTime)
                vU    = np.asarray(u)
                if vU.shape[1]!=len(vTime):
                    raise Exception('Second dimension of Input time series does not match time dimension ({} instead of {} )'.format(vU.shape[1],len(vTime)))
                # Store raw data
                self._inputs_ts = vU
                self._time_ts   = vTime
                # Create interpolant for faster evaluation
                u  = interp1d(vTime, vU)
                self.signature_u = 't'
        self._u = u
        if self.verbose:
            print('Setting inputs, signature:',self.signature_u)
# 

    # --------------------------------------------------------------------------------}
    # --- Initial conditions
    # --------------------------------------------------------------------------------{
    def setStateInitialConditions(self,q0=None):
        self.q0 = q0

    # --------------------------------------------------------------------------------}
    # --- INPUTS
    # --------------------------------------------------------------------------------{
    def setInputTimeSeries(self, vTime, vU):
        """ 
        Set the inputs as a time series of time
        INPUTS:
         - vTime: 1d array of time steps (do not need to be regular), of length nt
         - vU   : nInputs x nt array of inputs at each time steps
        """
        self.__setU(vU, vTime)

    def setInputFunction(self, fn, signature_u=None):
        """ 
        Set the inputs as a function of time and states
        The function will be used during the time integration

        INPUTS:
         - fn: handle to a python function. The signature of the function fn is: 
         - signature_u:  't' or 't,q'

               u =  f(t)
               u =  f(t, q) 

              where:
                 t   : scalar, time
                 q   : (nStates,) array, states
                 u   : (nInputs,) array, inputs 
        
        """
        self.u = fn
        # override signature_u
        if signature_u is not None:
            self.signature_u = signature_u

    def setInputFunctionDict(self, u_fn, signature_u=None):
        """ 
         u_fn:         dictionary of interpolant:  u[key](t) = scalar
         signature_u:  dict   signature_u[key] = 't' or 't,q', or 't,q,qd
        """

        if not isinstance(u_fn, dict):
            raise Exception("`u_fn` must be a dictionary")

        self.u = u_fn
        if signature_u is not None:
            if not isinstance(signature_u, dict):
                raise Exception("signature_u needs to be a dictionary")
            # Override
            self.signature_u = signature_u


    def Inputs(self, t, q=None, qd=None):
        """ return inputs array at t

        t: scalar

        """
        if self._u is None:
            return None # Should we reutrn q*0?

        elif isinstance(self.u, OrderedDict):
            d = np.zeros(len(self._u))
            #print('>>> Inputs signature',self.signature_u, t)
            for i,(k,fn) in enumerate(self._u.items()): # TODO use an ordered dict
                #try:
                sig = self.signature_u[k]
                #except:
                if sig=='t':
                    d[i] = self._u[k](t)
                elif sig=='t,q':
                    d[i] = self._u[k](t,q)
                else:
                    raise NotImplementedError()
            return d 

        else:
            if self.signature_u == 't':
                return self._u(t)
            elif self.signature_u == 't,q':
                return self._u(t,q)
            else:
                raise NotImplementedError()
                
    # --------------------------------------------------------------------------------}
    # --- State equation
    # --------------------------------------------------------------------------------{
    def dqdt_eval(self, t, q, p=None, u=None):
        """ Evaluate state space at t and q (and potentially u&q) """
        return self.dqdt_ODE(p=p, u=u)(t,q)

    # --- signature related functions
    def dqdt_ODE(self, p=None, u=None):
        """ 
        Return a function handle with signature qdot = dqdt(t,q)
        irrespectively of user signature
        """
        if p is None:
            p=self.p
        if u is None:
            u=self.u

        if self.signature == 't,q':
            odefun = lambda t, q : self.dqdt(t, q)
        elif self.signature == 't,q,p':
            odefun = lambda t, q : self.dqdt(t, q, p)
        elif self.signature == 't,q,u':
            odefun = lambda t, q : self.dqdt(t, q, u)
        elif self.signature == 't,q,u,p':
            odefun = lambda t, q : self.dqdt(t, q, u, p)
        elif self.signature == 't,q,p,u':
            odefun = lambda t, q : self.dqdt(t, q, p, u)
        else:
            raise NotImplementedError('Signature not supported',self.signature)
        return odefun

    # --------------------------------------------------------------------------------}
    # --- OUTPUTS
    # --------------------------------------------------------------------------------{
    def dqdt_calcOutput(self, u=None, p=None, signatureWanted='t,q'):
        """ Return function handle that computes outputs with a requested signature """
        if p is None:
            p = self.p
        if u is None:
            u = self.u

        if signatureWanted=='t,q':
            if self.signature == 't':
                odefun = lambda t, q : self.dqdt(t, calcOutput=True)
            elif self.signature == 't,q,p':
                odefun = lambda t, q : self.dqdt(t, q, p, calcOutput=True)
            elif self.signature == 't,q,u':
                odefun = lambda t, q : self.dqdt(t, q, u, calcOutput=True)
            elif self.signature == 't,q,u,p':
                odefun = lambda t, q : self.dqdt(t, q, u, p, calcOutput=True)
            elif self.signature == 't,q,p,u':
                odefun = lambda t, q : self.dqdt(t, q, p, u, calcOutput=True)
            else:
                raise NotImplementedError(self.signature)

        elif signatureWanted=='t,q,u,p':

            if self.signature == 't':
                odefun = lambda t, q, u, p : self.dqdt(t, calcOutput=True)
            elif self.signature == 't,q,p':
                odefun = lambda t, q, u, p : self.dqdt(t, q, p, calcOutput=True)
            elif self.signature == 't,q,u':
                odefun = lambda t, q, u, p : self.dqdt(t, q, u, calcOutput=True)
            elif self.signature == 't,q,u,p':
                odefun = lambda t, q, u, p : self.dqdt(t, q, u, p, calcOutput=True)
            elif self.signature == 't,q,p,u':
                odefun = lambda t, q, u, p : self.dqdt(t, q, p, u, calcOutput=True)
            else:
                raise NotImplementedError(self.signature)
        else:
            raise NotImplementedError(self.signatureWanted)

        return odefun

    # --------------------------------------------------------------------------------}
    # --- Time integration 
    # --------------------------------------------------------------------------------{
    def integrate(self, t_eval, method='RK45', y0=None, p=None, u=None, calc='u,y,qd', xoffset=None, **options):
        #
        if y0 is not None:
            self.setStateInitialConditions(y0)
        if p is not None:
            self.p = p
        if u is not None:
            self.u = u
        if method is None:
            method = 'RK45'
        # Clean values stored after integration
        self.cleanSimData()

        # Getting a function handle with simple signature: dqdt(t,q)
        dqdt = self.dqdt_ODE()
        # Time integration
        if self.verbose:
            print('Time integration...')
        if method=='odeint':
            x = odeint(dqdt, y0, t_eval, tfirst=True)
            res = OdeResultsClass(t=t_eval, y=x.T) # To mimic result class of solve_ivp
        else:
            res = solve_ivp(fun=dqdt, t_span=[t_eval[0], t_eval[-1]], y0=self.q0, t_eval=t_eval, method=method, vectorized=False, **options)   

        # Store
        self.res  = res

        # --- From results to states, inputs, outputs DataFrame
        df = self.res2DataFrame(res, calc=calc, sStates=None, xoffset=xoffset)
        self.df = df

        return res, df
        
    # --------------------------------------------------------------------------------}
    # --- Simulation storage
    # --------------------------------------------------------------------------------{
    def cleanSimData(self):
        self.res      = None
        self.dfIn     = None
        self.dfOut    = None
        self.dfStates = None
        self.df       = None


    def res2DataFrame(self, res=None, calc='u,y,xd', Factors=None, sStates=None, xoffset=None, uoffset=None, yoffset=None):
        """ Return time integration results as a dataframe
        """
        # TODO harmonize with mech_system.res2DataFrame
        calcVals = calc.split(',')

        if res is None:
            if self.res is None:
                raise Exception('Call integrate before res2DataFrame')
            res = self.res

        # --- Time and states
        dfStates = self.store_states(res, sStates=sStates, xoffset=xoffset, Factors=Factors)
        # --- Accelerations 
        # --- Forcing

        # --- Try to compute outputs
        dfOut = None
        if 'y' in calcVals:
            dfOut = self.calc_outputs(res=res, insertTime=True, yoffset=yoffset)
            if dfOut is None:
                #TODO
                dfStatesD = self.calcDeriv()

        # --- Try to compute inputs
        dfIn = None
        if 'u' in calcVals:
            dfIn = self.calc_inputs(res=res, insertTime=True, uoffset=uoffset)

        # --- Concatenates everything into one DataFrame
        df = pd.concat((dfStates, dfIn, dfOut), axis=1)
        df = df.loc[:,~df.columns.duplicated()].copy()

        return df

    def store_states(self, res, sStates=None, xoffset=None, Factors=None):
        nStates = len(self.q0)
        if sStates is None and self.sX is None:
            sStates = ['x{}'.format(i+1) for i in range(nStates)] # TODO
        else:
            if sStates is None:
                sStates = self.sX
            nCols = res.y.shape[0]
            if len(self.sX)!=nCols:
                raise Exception("Inconsistency in length of states columnNames. Number of columns detected from res: {}. States columNames (sX):".format(nCols, self.sX))
        self.sX = sStates

        # Store as a matrix
        M = res.y.T.copy()

        # Scaling
        if Factors is not None:
            for i, f in enumerate(Factors):
                M[:,i] *= f

        # Scaling offsets
        if Factors is not None and xoffset is not None:
            xoffset *= np.asarray(Factors)

        if xoffset is not None:
            xoffset=np.asarray(xoffset).flatten()
            M = M + xoffset

        # Offset Velocity TODO
        #if xd0 is not None: 
        #    for i, xd0_ in enumerate(xd0):
        #        M[:,i]           += xd0_*res.t # Position increases linearly (e.g. azimuth)
        # Degrees
        for i,d in enumerate(sStates):
            if sStates[i].find('[deg]')>1: 
                if np.max(M[:,i])>180:
                    M[:,i] = np.mod(M[:,i], 360)


        dfStates = pd.DataFrame(data=M, columns=sStates)
        # Insert time
        dfStates.insert(0, 'Time_[s]', res.t)
        self.dfStates = dfStates
        return dfStates


    def _prepareOutputDF(self, res):
        cols    = self._inferOutputCols(res)
        data    = np.full((len(res.t), len(cols)), np.nan)
        df      = pd.DataFrame(columns = cols, data = data)
        self.sY = cols
        return df

    def _prepareInputDF(self, time=None, q=None):
        if time is None:
            time = self.res.t
            q    = self.res.y
        cols    = self._inferInputCols(time, q)
        data    = np.full((len(time), len(cols)), np.nan)
        df      = pd.DataFrame(columns = cols, data = data)
        self.sU = cols
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
                #raise Exception("Inconsistency in length of output columnnames. Number of columns detected from `calcOuputt`: {}. Ouput columNames (sY):".format(len(cols), self.sY))
                print('[WARN] Inconsistency in length of output columnnames. Number of columns detected from `calcOuputt`: {}. Ouput columNames (sY):'.format(len(cols), self.sY))
            #cols = self.sY
        return cols

    def _inferInputCols(self, time=None, q=None):
        """ See what Inputs returns at t[0] to allocate storage memory """
        if time is None:
            time = self.res.t
            q    = self.res.y
        u0    = self.Inputs(time[0], q[:,0])
        nCols = len(u0)
        if self.sU is None:
            if isinstance(self.u, OrderedDict): 
                cols = self.u.keys()
            else:
                # Default column names
                cols  = ['u{:d}'.format(i+1) for i in range(nCols)]
        else:
            cols = self.sU
            if len(self.sU)!=nCols:
                raise Exception("Inconsistency in length of input columns. Number of columns detected from `Inputs`: {}. Inputs columns (sU):".format(nCols, self.sU))
        return cols

    def _calc_outputs(self, time, q, df, yoffset=None):
        """ low level implementation leaving room for optimization for other subclass."""
        calcOutput = self.dqdt_calcOutput()
        if yoffset is None:
            yoffset =0
        if self.verbose:
            print('Calc output...')
        for i,t in enumerate(time):
            y = calcOutput(t, q[:,i])
            df.iloc[i,:] = y + yoffset

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
        self._calc_outputs(res.t, res.y, df, yoffset=yoffset)

        if insertTime:
            df.insert(0,'Time_[s]', res.t)
        self.dfOut = df

        if dataFrame:
            return df
        else:
            return df.values

    def calcDeriv(self, insertTime=True):
        if self.res is None:
            raise Exception("Call `integrate` before calcDeriv.")
        # TODO
        return None


    def _calc_inputs(self, time, q, df):
        """ low level implementation leaving room for optimization for other subclass."""
        if self.verbose:
            print('Calc inputs...')
        for i,t in enumerate(time):
            df.iloc[i,:] = self.Inputs(t, q[:,i])

    def calc_inputs(self, time=None, res=None, insertTime=True, dataFrame=True, uoffset=None):
        """ 
        Compute inputs at each time step. 
        provide:
               time: only works if all inputs have signature 't', not 't,q'
            or
               res
        """
        if self._u is None:
            return None

        if time is not None:
            q = np.zeros((0,len(time))) # Dummy
            if not self.uDependsOnTOnly():
                raise Exception('Cannot compute inputs signature of u is not a simple function of time. Call `integrate` first. Signature of u is:'.format(self._signature_u))
        else:
            if res is None:
                res = self.res
            if res is None:
                raise Exception("Call `integrate` before calc_inputs.")
            time = res.t
            q    = res.y

        # TODO consider calling calcDeriv to obtain state derivatives, might be needed by some inputs

        # --- Infer column names and allocate
        df = self._prepareInputDF(time, q)

        # --- Calc inputs based on states
        self._calc_inputs(time, q, df)
        if uoffset is not None:
            df.iloc[:,:] +=uoffset

        if insertTime:
            df.insert(0,'Time_[s]', time)

        self.dfIn = df

        if dataFrame:
            return df
        else:
            return df.values

    # --------------------------------------------------------------------------------}
    # --- linearization
    # --------------------------------------------------------------------------------{
    def linearize(self, x0, dx, u0=None, du=None, p=None): 
        """
        Small change of interface compared to parent class

        OUTPUS:

        """
        if u0 is None:
            op=(0,x0)
        else:
            op=(0,x0,u0)
        # Setting parameters if they are present
        if p is not None:
            self.param = p
        else:
            p = self.p
            self.param = self.p
        
        # TODO handle calcOuput beelfter overall
        try:
            calcOutputForLin = self.dqdt_calcOutput(signatureWanted='t,q,u,p')
            out = calcOutputForLin(0, x0, du, p)
            #print('>>>>',out)
            # Giving this function to parent class for linearization
            self.Y = calcOutputForLin
        except TypeError:
            print("[FAIL] Error evaluating model outputs for linearization. Does the function dqdt has the argument `calcOutput`?")


        # Calling parent class
        lin = System.linearize(self, op, dx, du=du, use_implicit=False)

        A=lin[0]
        B=lin[1]

        if len(lin)==2:
            C=None
            D=None
        else:
            C=lin[2]
            D=lin[3]
        return A,B,C,D


    def eigA(self, x0, dx, nq2, u0=None, du=None, p=None, normQ=None, fullEV=False):
        """ nq: number of second order states ntot = 2*nq2 + nq1 """
        from welib.tools.eva import eigA
        A,_,_,_ = self.linearize(x0, dx=dx, u0=u0, du=du, p=p) 
        freq_d, zeta, Q, freq_0 = eigA(A, nq=nq2, fullEV=fullEV, normQ=normQ, sort=True)
        return freq_d, zeta, Q, freq_0 

    # --------------------------------------------------------------------------------}
    # --- Static equilibrium
    # --------------------------------------------------------------------------------{
    def equilibrium(self, x0, dx, u0=None, du=None, p=None, maxIter=1000, tol=1e-5, verbose=False):
        """ 
        Use a Newton method to find equilibrium point (dqdt=0), starting from x0
        and given input u0
          dx and du are needed to compute numerical Jacobians
        """
        x0 = np.asarray(x0).copy()
        u0 = np.asarray(u0).copy()

        f0 = self.dqdt_eval(t=0, q=x0, p=p, u=u0)
        fnorm = np.linalg.norm(f0)
        for nIter in range(maxIter):
            if fnorm<tol:
                break
            # Value at operating point
            f0 = self.dqdt_eval(t=0, q=x0, p=p, u=u0)
            fnorm = np.linalg.norm(f0)
            # Jacobian at operating point
            A,B,C,D = self.linearize(x0=x0, u0=u0, dx=dx, du=du, p=p)
            # Going to next point
            dx_ = solve(A,f0)
            x0 -= dx_
        if nIter<maxIter:
            if verbose:
                print('[ OK ] equilibrium point found after {} iterations'.format(nIter+1))
        else:
            print('[FAIL] equilibrium point not found after {} iterations. Norm: {}'.format(nIter+1, fnorm))
        return x0

    # --------------------------------------------------------------------------------}
    # --- Frequency domain and transfer function
    # --------------------------------------------------------------------------------{
    def frequency_response_numerical(self, omega, nStates, nInputs, nOutputs, int_method=None, **kwargs):
        """ 
        NOTE: if nOutputs is smaller, a subset of outputs is used
             TODO I_Outputs to select outputs
        """
        print('>>>> TODO Frequency response for nonlinear function needs further verification')
        from welib.system.transferfunction import numerical_frequency_response

        if self.q0 is not None:
            q0 = self.q0
        else:
            q0 = np.zeros(nStates) # TODO q0 equilibrium

        # We create a function that returns outputs given inputs
        def calcOutput(t, MU):
            """ """
            self.setInputTimeSeries(t, MU)
            _, df = self.integrate(t, y0=q0, calc='y', method=int_method)
            return df[self.sY].values.T
        #  We call the generic method 
        G, phi = numerical_frequency_response(omega, calcOutput, nInputs, nOutputs, **kwargs)
        return G, phi


    # --------------------------------------------------------------------------------}
    # ---  IO functions for printing
    # --------------------------------------------------------------------------------{
    def __repr__(self):
        s='<{} object>\n'.format(type(self).__name__)
        s+='|Read-only attributes:\n'
        s+='| - signature  : {} \n'.format(self.signature)
        s+='| - signature_u: {} \n'.format(self.signature_u)
        try:
            s+='| - u: size: {} \n'.format(self.nInputs)
            if self._u is not None:
                if isinstance(self.u, OrderedDict):
                    s+='|      keys: {} \n'.format(self._u.keys())
        except:
            pass
        try:
            s+='| - p: size: {} \n'.format(self.nParams)
            if self.p is not None:
                if isinstance(self.p, dict):
                    s+='|      keys: {} \n'.format(self.p.keys())
        except:
            pass
        #s+='|Attributes:\n'
        s+='| - q0: {}\n'.format(self.q0)
        return s


    # --------------------------------------------------------------------------------}
    # ---  Plotting functions
    # --------------------------------------------------------------------------------{
    # TODO merge those df or res
    def plot_x_legacy(self, axes=None, label=None, res=None, **kwargs):
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

        return axes

    def plot_u_legacy(self, axes=None, label=None, res=None, **kwargs):
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


    def plot(self, df=None, keys=None, label=None, axes=None, **kwargs):
        """ Simple plot after time integration"""
        if df is None:
            df = self.df
            if df is None:
                raise Exception("Call `integrate` before plot")

        return _plot(df, keys=keys, label=label, title='', nPlotCols=1, axes=axes, **kwargs)

    def plot_states(self, df=None, axes=None, **kwargs):
        """ Simple plot of states after time integration"""
        if df is None:
            if self.dfStates is None:
                raise Exception("Call `integrate` before plot_states")
            df = self.dfStates
        keys = list(self.sStates)

        return _plot(df, keys=keys, title='', nPlotCols=1, axes=axes, **kwargs)

    def plot_inputs(self, df=None, axes=None, **kwargs):
        """ 
        plot inputs. Requires either to call `integrate` or `calc_inputs` before
        """
        if df is None:
            if self.dfIn is None:
                if self.res is not None:
                    self.calc_inputs(res=self.res, insertTime=True)
                else:
                    raise Exception("Call `calc_inputs` of `integrate` before plot_inputs")
            df = self.dfIn
            if df is None:
                print('[WARN] no inputs to plot')
                return
        keys = list(self.sInputs)

        return _plot(df, keys=keys, title='', nPlotCols=1, axes=axes, **kwargs)

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


