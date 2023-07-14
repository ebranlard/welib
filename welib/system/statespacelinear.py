import os
import numpy as np
import pandas as pd
from numpy.linalg import inv
from numpy.linalg import solve
from numpy.linalg import eigvals, matrix_rank
from scipy.integrate import  solve_ivp #odeint
from scipy.optimize import OptimizeResult as OdeResultsClass 
from scipy.linalg import expm
# Local
from welib.system.statespace import StateSpace
from welib.fast.tools.lin import matToSIunits, subMat, subSeries # TODO
# --------------------------------------------------------------------------------}
# --- Simple statespace functions ltiss (linear time invariant state space)
# --------------------------------------------------------------------------------{
def calc_deriv(x, u, A, B):
    """
    Calculate outputs for a LTI state space model
    INPUTS:
     - x: array-like of shape nx, or, matrix of shape nx x nt
     - u: array-like of shape nu, or, matrix of shape nu x nt
     - A: array nx x nx
     - B: array nx x nu
    OUTPUT:
      xdot = A x + B u   (array-like of shape nx, or, matrix of shape nx x nt)
    """
    return A.dot(x) + B.dot(u)

def calc_outputs(x, u, C, D):
    """
    Calculate outputs for a LTI state space model
    INPUTS:
     - x: array-like of shape nx, or, matrix of shape nx x nt
     - u: array-like of shape nu, or, matrix of shape nu x nt
     - C: array ny x nx
     - D: array ny x nu
    OUTPUT:
      y = B x + D u   (array-like of shape ny, or, matrix of shape nx x nt)
    """
    return C.dot(x) + D.dot(u)

def state_function(t, x, u, p):
    """ see calc_deriv """
    return calc_deriv(x, u, p['A'], p['B'])

def output_function(t, x, u, p):
    """ see calc_output """
    return calc_outputs(x, u, p['C'], p['D'] )

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
    # NOTE: here we allow inputs that are function of states, which is not really a LIT anymore!
    hasq=False
    try:
        y=fU(t_eval[0],q0)
        hasq=True
    except:
        try:
            y = fU(t_eval[0])
            hasq=False
        except:
            raise
    if y is None:
        raise Exception('Function fU evaluate to None')

    if hasq:
        #def odefun(t, q):
        #    return np.dot(A, q) + np.dot(B, fU(t,q))
        odefun = lambda t, q : np.dot(A, q) + np.dot(B, fU(t,q))
    else:
        odefun = lambda t, q : np.dot(A, q) + np.dot(B, fU(t) )

    res = solve_ivp(fun=odefun, t_span=[t_eval[0], t_eval[-1]], y0=q0, t_eval=t_eval, method=method, vectorized=False, **options)   

    # TODO consider returning y

    return res


def integrate_convolution(time, A, B, fU=None, C=None, U=None):
    """ 
    Perform time integration of a LTI state space system using convolution method

    INPUTS:
     - time: array of time values (nt)
     - A: state matrix (nStates x nStates)
     - B: input matrix (nStates x nInputs)
     - fU: function/interpolants with interface U=fU(t) or U=fU(t,q)
           where U is the array of inputs at t

    OUTPUTS:
     - x: state vector

    """
    from welib.tools.signal_analysis import convolution_integral

    nU = B.shape[1]
    if nU>1:
        print('>>>> TODO weird response observed at low frequencies for multiple inputs')


    H = impulse_response_matrix(time, A, B) # nStates x nInputs x nt

    x = np.zeros((A.shape[0], len(time)))

    # TODO inline and optimize
    if U is not None:
        raise NotImplementedError()

        # verify that U is of shape nInputs x nt
    else:
        try:
            U = fU(time)
        except:
            print('[WARN] Cannot evaluate fU for all time, doing a for loop...')
            U = np.zeros((B.shape[1], len(time)))
            for it,t in enumerate(time):
                U[:,it] = fU(t) # NOTE: cannot do state dependency here

    for ix in np.arange(H.shape[0]): # loop on states
        x_sum=0
        for iu in np.arange(H.shape[1]): # loop inputs
            x_sum  += convolution_integral(time, U[iu,:], H[ix,iu,:] )
        x[ix,:] = x_sum
    return x

    # TODO consider returning y

def impulse_response_matrix(time, A, B, C=None, outputBoth=False):
    """ 
    Return the impulse response matrix for all time steps defined by `time`
        H_x(t) =   exp(At) B   array of shape nx x nu x nt
        H_y(t) = C exp(At) B   array of shape ny x nu x nt
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
class LinearStateSpace(StateSpace):
    """ 
    def setStateInitialConditions(self,q0=None):
    def setInputTimeSeries(self,vTime,vU):
    def setInputFunction(self,fn):
    def Inputs(self,t,q=None):
    def integrate(self, t_eval, method='RK4', y0=None, **options):
    def dqdt(self, t, q):
    def RHS(self,t,q):
    def nStates(self):
    def nInputs(self):
    def nOutputs(self):
    """
    def __init__(self, A=None, B=None, C=None, D=None, 
            qop=None, uop=None, yop=None,
            sX=None, sU=None, sY=None,
            verbose=False,
            ):
        """ 
        INPUTS:
         - A, B, C, D: state matrices (may be None/empty)
         - qop: degrees of freedom at operating point
         - sX: list of names for states
         - sY: list of names for outputs
         - sU: list of names for inputs
        """
        # --- Data
        self.A = None
        self.B = None
        self.C = None
        self.D = None
        self.q0_  = None
        self.qop_ = None
        self.uop_ = None
        self.yop_ = None
        if A is None:
            A = np.zeros((0,0))
        self.A = np.asarray(A)
        # TODO better merge the two
        # TODO need a signature_output as well!
        # --- StateSpace
        StateSpace.__init__(self, dqdt=self.dqdt_tqu, signature='t,q,u', 
                sX = sX, sU=sU, sY=sY,
                verbose=verbose)
        # --- System
        self.Y = self.Outputs_tqu

        # --- Linear system
        if B is None:
            nq = A.shape[0]
            B = np.zeros((nq,0))
            # No Inputs
            self.setInputFunction(lambda t,q: np.zeros(0), signature_u='t,q')
        self.B = np.asarray(B)

        if C is None:
            if sY is None:
                sY = sX
                self.C=np.eye(A.shape[0]) # output all states
            else:
                self.C=np.zeros((len(sY), A.shape[0])) 
        else:
            self.C=np.asarray(C)
        if D is None:
            self.D=np.zeros((self.C.shape[0],self.B.shape[1]))
        else:
            self.D=np.asarray(D)

        # Operating point of linear model
        if qop is None:
            qop = np.zeros(self.nStates)*np.nan
        if uop is None:
            uop = np.zeros(self.nInputs)*np.nan
        if yop is  None:
            yop = np.zeros(self.nOutputs)*np.nan
        self.qop_ = qop
        self.uop_ = uop
        self.yop_ = yop

        # --- For linear models we can infer the dimensions of states easily
        #if sX is None:
        #    self.sX = ['x{}'.format(i+1) for i in range(self.nStates)]
        #if sY is None:
        #    self.sY = ['y{:d}'.format(i+1) for i in range(self.nOutputs)]
        #if sU is None:
        #    self.sU = ['u{:d}'.format(i+1) for i in range(self.nInputs)]


    #def stateMatrices(self, A=None, B=None, C=None, D=None, sX=None, sU=None, sY=None):
#         if A is None:
#             A = np.zeros((0,0))
#         self.A = np.asarray(A)
#         # --- Linear system
#         if B is None:
#             nq = A.shape[0]
#             B = np.zeros((nq,0))
#             # No Inputs
#             self.setInputFunction(lambda t,q: np.zeros(0), signature_u='t,q')
#         self.B = np.asarray(B)
# 
#         if C is None:
#             if sY is None:
#                 sY = sX
#                 self.C=np.eye(A.shape[0]) # output all states
#             else:
#                 self.C=np.zeros((len(sY), A.shape[0])) 
#         else:
#             self.C=np.asarray(C)
#         if D is None:
#             self.D=np.zeros((self.C.shape[0],self.B.shape[1]))
#         else:
#             self.D=np.asarray(D)
#         # --- Names
#         self.sX = sX
#         self.sY = sY
#         self.sU = sU

    @property
    def q0(self):
        return pd.Series(self.q0_, index=self.sStates)

    @property
    def qop(self):
        return pd.Series(self.qop_, index=self.sStates)

    @property
    def qop_default(self):
        return pd.Series(np.zeros(len(self.sX)), index=self.sX)

    @property
    def uop(self):
        return pd.Series(self.uop_, index=self.sInputs)

    @property
    def yop(self):
        return pd.Series(self.yop_, index=self.sOutputs)



    @property
    def nStates(self):
        return self.A.shape[0]

    @property
    def nInputs(self):
        return self.B.shape[1]

    @property
    def nOutputs(self):
        if self.C is not None:
            return self.C.shape[0]
        else:
            return 0

    # --------------------------------------------------------------------------------}
    # --- Initial conditions
    # --------------------------------------------------------------------------------{
    def setStateInitialConditions(self, q0=None):
        if q0 is not None:
            if len(q0)!=self.nStates:
                raise Exception('Wrong dimension for q0 ({} instead of {} )'.format(len(q0),self.nStates))
        else:
            q0 = np.zeros(self.nStates)
        self.q0_ =  q0 # pd.Series(q0, index=self.sStates)

            
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
        vTime = np.asarray(vTime)
        vU    = np.asarray(vU)
        if self.nInputs==1:
            vU = vU.reshape((1,-1))
        if vU.shape[0]!=self.nInputs:
            raise Exception('Wrong first dimension for Inputs time series ({} instead of {} )'.format(vU.shape[0],self.nInputs))
        # Call parent class (create interpolant)
        StateSpace.setInputTimeSeries(self, vTime, vU)

    # See statespace.py
    #def Inputs(self, t, q=None, qd=None):

    # --------------------------------------------------------------------------------}
    # --- State equation
    # --------------------------------------------------------------------------------{
    def dqdt(self, t, q):
        # NOTE: this can cause issues if q is not flat
        return np.dot(self.A, q) + np.dot(self.B, self.Inputs(t,q))

    def dqdt_tqu(self, t, q, u):
        # NOTE: this can cause issues if q is not flat
        return np.dot(self.A, q) + np.dot(self.B, u)

    def RHS(self,t,q):
        return self.dqdt(t,q)

    # --------------------------------------------------------------------------------}
    # --- Outputs
    # --------------------------------------------------------------------------------{
    def dqdt_calcOutput(self):
        return self.Outputs

    def Outputs(self, t, q):
        # NOTE: could be optimized if inputs are already computed...
        u = self.Inputs(t,q)
        return self.C.dot(q) + self.D.dot(u)

    def Outputs_tqu(self, t, q, u):
        return self.C.dot(q) + self.D.dot(u)

    # --------------------------------------------------------------------------------}
    # --- Time integration 
    # --------------------------------------------------------------------------------{
    def integrate(self, t_eval, method='RK45', y0=None, u=None, calc='', xoffset=None, uoffset=None, yoffset=None, **options):
        #dfLI = sysLI.res2DataFrame(self.channels, self.FASTDOFScales, q0=qop, xd0=qdop, acc=acc, forcing=forcing, sAcc=self.acc_channels)
        #
        if y0 is not None:
            self.setStateInitialConditions(y0)
        if u is not None:
            self.u = u
        if method is None:
            method='RK45'

        # --- Sanity checks
        if len(self.q0_)!=self.nStates:
            raise Exception('Size of initial condition ({}) different from number of states ({}).'.format(len(self.q0_), self.nStates))
        if self.A.shape[0]!=self.nStates or self.A.shape[1]!=self.nStates:
            raise Exception('Shape of A ({}) different from number of states ({}).'.format(self.A.shape, self.nStates))
        if np.any(np.isnan(self.A)): raise Exception('A matrix contains nan')
        if np.any(np.isnan(self.B)): raise Exception('B matrix contains nan')
        if np.any(np.isnan(self.C)): raise Exception('C matrix contains nan')
        if np.any(np.isnan(self.D)): raise Exception('D matrix contains nan')
        if np.any(np.isnan(self.q0_)): raise Exception('Initial condition contains nan')

        # Clean values stored after integration
        self.cleanSimData()
        # Time integration
        if self.verbose:
            print('Time integration...')
        if method.lower()=='impulse':
            # TODO add check on initial conditions. Only works if initial conditions are 0 if I remember correctly
            x = integrate_convolution(t_eval, self.A, self.B, self.Inputs)

            res = OdeResultsClass(t=t_eval, y=x) # To mimic result class of solve_ivp

        else:
            res = integrate(t_eval, self.q0_, self.A, self.B, self.Inputs, method=method, **options)

        # Store
        self.res    = res

        # --- From results to states, inputs, outputs DataFrame
        df = self.res2DataFrame(res, calc=calc, sStates=None, xoffset=xoffset, uoffset=uoffset, yoffset=yoffset)
        self.df = df 

        return res, df

    # --------------------------------------------------------------------------------}
    # --- Simulation storage
    # --------------------------------------------------------------------------------{
    # From welib.system.statespac.py
    #def res2DataFrame(self, res=None, calc='u,y,xd', sStates=None, xoffset=None, uoffset=None, yoffset=None):
    #def store_states(self, res, sStates=None, xoffset=None, Factors=None):
    #def calc_outputs(self, res=None, insertTime=True, dataFrame=True, yoffset=None):
    #def calc_inputs(self, time=None, res=None, insertTime=True, dataFrame=True, uoffset=None):

    def _calc_outputs(self, time, q, df, yoffset=None):
        """ low level implementation leaving room for optimization for other subclass."""
        if yoffset is None:
            yoffset = 0
        if self._inputs_ts is not None and len(time)==len(self._time_ts): # TODO more rigorous
            if self.verbose:
                print('Calc output using simple matrix manipulation...')
            MX = q
            MU = self._inputs_ts
            MY = calc_outputs(MX, MU, self.C, self.D)
            df.iloc[:,:] = MY.T
            df.iloc[:,:] += yoffset
        else:
            calcOutput = self.dqdt_calcOutput()
            if self.verbose:
                print('Calc output...')
            for i,t in enumerate(time):
                y = calcOutput(t, q[:,i])
                df.iloc[i,:] = y+yoffset


    def _calc_inputs(self, time, q, df):
        """ low level implementation leaving room for optimization for other subclass."""
        if self._inputs_ts is not None and len(time)==len(self._time_ts):
            if self.verbose:
                print('Calc inputs using simple matrix manipulation...')
            df.iloc[:,:] = self._inputs_ts.T
        else:
            if self.verbose:
                print('Calc inputs...')
            for i,t in enumerate(time):
                df.iloc[i,:] = self.Inputs(t, q[:,i])



    def calc_impulse_response_matrix(self, time, insertTime=True):
        """ 
        Return the impulse response matrix for all time steps defined by `time`
            H_x(t) =   exp(At) B
            H_y(t) = C exp(At) B
            see e.g. 
               Friedland p 76
        """
        nx = self.nStates
        nu = self.nInputs

        cols = ['H_x{}_u{}'.format(i+1,j+1) for i in range(nx) for j in range(nu)]

        data = np.full((len(time), len(cols)), np.nan)

        # --- Calc impuse response matrix
        if self.verbose:
            print('Calc impulse response...')
        for i,t in enumerate(time):
            H = expm(self.A*t).dot(self.B) # nx x nu
            data[i,:] = H.flatten()

        df = pd.DataFrame(data=data, columns=cols)

        if insertTime:
            df.insert(0,'Time_[s]', time)
        return df


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
            ny = self.nOutputs
            ns = len(s)
            H = np.empty((ny, nu, ns), dtype=np.complex128)
            # Loop on frequencies
            for k,sk in enumerate(s.ravel()):
                H[:,:,k] =np.dot(self.C, solve(sk* np.eye(self.nStates) - self.A, self.B)) + self.D
            #H = np.array([np.dot(self.C, solve(sk* np.eye(self.nStates) - self.A, self.B)) + self.D  for sk in s.ravel()])
            #H=H.reshape([ny,nu]+list(s.shape))
        else:
            H =           np.dot(self.C, solve(s * np.eye(self.nStates) - self.A, self.B)) + self.D

        return H


    def frequency_response(self, omega, deg=False, method='transferFunction', **kwargs):
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
        if method=='transferFunction':
            s = omega * 1.j
            H = self.transferFunction(s)
            if deg:
                return np.abs(H), np.angle(H)*180/np.pi
            else:
                return np.abs(H), np.angle(H)

        elif method=='numerical':
            return self.frequency_response_numerical(omega, **kwargs)
        else:
            raise NotImplementedError()

    def frequency_response_numerical(self, omega, int_method=None, **kwargs):
        from welib.system.transferfunction import numerical_frequency_response

        if self.q0_ is not None:
            q0 = self.q0_
        else:
            q0 = np.zeros(self.nStates) # TODO q0 equilibrium

        # We create a function that returns outputs given inputs
        def calcOutput(t, MU):
            """ """
            self.setInputTimeSeries(t, MU)
            _, df = self.integrate(t, y0=q0, calc='y', method=int_method)
            return df[self.sY].values.T
        #  We call the generic method 
        G, phi = numerical_frequency_response(omega, calcOutput, self.nInputs, self.nOutputs, **kwargs)
        return G, phi
    

    def eigA(self):
        from welib.tools.eva import eigA
        A = self.A # NOTE: not always defined...
        freq_d, zeta, Q, freq_0 = eigA(A, fullEV=False, normQ=None, sort=True)
        return freq_d, zeta, Q, freq_0 

    # --------------------------------------------------------------------------------}
    # --- Matrix manipulation
    # --------------------------------------------------------------------------------{
    def rename(self, colMap, verbose=False):
        """ Rename labels """
        sX = list(self.sX)
        sY = list(self.sY)
        sU = list(self.sU)
        keys = list(colMap.keys())
        def renameList(l, colMap, verbose):
            for i,s in enumerate(l):
                if s in keys:
                    l[i] = colMap[s] 
                else:
                    if verbose:
                        print('Label {} not renamed'.format(s))
            return l

        self.sX = renameList(self.sX, colMap, verbose)
        self.sU = renameList(self.sU, colMap, verbose)
        self.sY = renameList(self.sY, colMap, verbose)

        # Remove duplicate
        #df = df.loc[:,~df.columns.duplicated()].copy()


    def extract(self, sX=None, sU=None, sY=None, verbose=False, check=True, inPlace=True):
        """ """
        A, B, C, D = self.toDataFrames()
        q0, qop, uop, yop = self.toSeries()
        if sX is not None:
            self.qop_ = subSeries(qop, rows=sX, check=check)
            self.q0_  = subSeries(q0,  rows=sX, check=check)
            sXd = ['d'+s for s in sX]
            A = subMat(A, rows=sXd, cols=sX, check=check, name = 'A')
            B = subMat(B, rows=sXd         , check=check, name = 'B')
            C = subMat(C,          cols=sX, check=check, name = 'C')
        else:
            q0 = q0.values
        if sU is not None:
            self.uop_ = subSeries(uop, rows=sU, check=check)
            B = subMat(B,          cols=sU, check=check, name = 'B')
            D = subMat(D,          cols=sU, check=check, name = 'D')
        if sY is not None:
            self.yop_ = subSeries(yop, rows=sY, check=check)
            C = subMat(C, rows=sY,          check=check, name = 'C')
            D = subMat(D, rows=sY,          check=check, name = 'D')

        if inPlace:
            self.fromDataFrames(A, B, C, D)
            # Trigger, make sure q0 has the right size
            self.setStateInitialConditions(self.q0_)

        return A, B, C, D


    def toSI(self, verbose=False):
        """ Use labels to ensures that A, B, C, D have SI units"""
        A, B, C, D = self.toDataFrames()

        A = matToSIunits(A, name='A', verbose=verbose)
        B = matToSIunits(B, name='B', verbose=verbose)
        C = matToSIunits(C, name='C', verbose=verbose)
        D = matToSIunits(D, name='D', verbose=verbose)

        self.fromDataFrames(A, B, C, D)


    def fromDataFrames(self, A, B=None, C=None, D=None):
        self.A  = A.values
        self.sX = A.columns.values
        if B is not None:
            self.B  = B.values
            self.sU = B.columns.values
        if C is not None:
            self.C  = C.values
            self.sY = C.index.values
        if D is not None:
            self.D  = D.values


    def toDataFrames(self):
        """ return dataframes for system matrices using labels"""
        sXd = ['d'+s for s in self.sX]
        A = pd.DataFrame(self.A, index=sXd    , columns=self.sX)
        B = pd.DataFrame(self.B, index=sXd    , columns=self.sU)
        C = pd.DataFrame(self.C, index=self.sY, columns=self.sX)
        D = pd.DataFrame(self.D, index=self.sY, columns=self.sU)
        return A, B, C, D

    def toSeries(self):
        return self.q0, self.qop, self.uop, self.yop


    # --------------------------------------------------------------------------------}
    # ---  IO functions for printing
    # --------------------------------------------------------------------------------{
    def __repr__(self):
        s='<{} object>\n'.format(type(self).__name__)
        s+='|Read-only attributes:\n'
        s+='| - nStates:  {} with names:{}\n'.format(self.nStates , self.sX)
        s+='| - nInputs:  {} with names:{}\n'.format(self.nInputs , self.sU)
        s+='| - nOutputs: {} with names:{}\n'.format(self.nOutputs, self.sY)
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
        s+='| - A: State-State Matrix ({} x {})\n'.format(self.nStates, self.nStates)
        s+=str(self.A)+'\n'
        s+='| - B: State-Input Matrix ({} x {})\n'.format(self.nStates, self.nInputs)
        s+=str(self.B)+'\n'
        if self.C is not None:
            s+='| - C: Outputs-States Matrix ({} x {})\n'.format(self.nOutputs, self.nStates)
            s+=str(self.C)+'\n'
        if self.D is not None:
            s+='| - D: Outputs-Inputs Matrix ({} x {})\n'.format(self.nOutputs, self.nInputs)
            s+=str(self.D)+'\n'

        s+='| * q0:  {}\n'.format(dict(self.q0))
        s+='| * qop: {}\n'.format(dict(self.qop))
        s+='| * uop: n:{}\n'.format(len(self.uop))
        s+='| * yop: n:{}\n'.format(len(self.yop))
        return s

    def load(self, pickleFile):
        import pickle
        if not os.path.exists(pickleFile):
            raise Exception('File does not exist: {}'.format(pickleFile))
        d = pickle.load(open(pickleFile,'rb'))
        self.fromDataFrames(d['A'], d['B'], d['C'], d['D'])
        self.setStateInitialConditions(d['q0'].values)
        try:
            self.qop_ = d['qop'].values
        except:
            raise Exception('The pickle file is an old pickle file, please regenerate it: {}'.format(pickleFile))
        self.uop_ = d['uop'].values
        self.yop_ = d['yop'].values
        return d

    def save(self, pickleFile, extraDict=None):
        import pickle
        if extraDict is None:
            extraDict={}
        A, B, C, D = self.toDataFrames()
        extraDict.update({'A':A, 'B':B, 'C':C, 'D':D, 'q0':self.q0, 'qop':self.qop, 'uop':self.uop, 'yop':self.yop})
        pickle.dump(extraDict, open(pickleFile,'wb'))



    # --------------------------------------------------------------------------------}
    # ---  Plotting functions
    # --------------------------------------------------------------------------------{
    # See StateSpace:
    # def plot_x_legacy(self, fig=None, axes=None, label=None, res=None, **kwargs):
    # def plot_u_legacy(self, axes=None, label=None, res=None, **kwargs):
    # def plot(self, df=None, keys=None, label=None, axes=None, **kwargs):
    # def plot_states(self, df=None, axes=None, **kwargs):
    # def plot_inputs(self, df=None, axes=None, **kwargs):
    # def plot_outputs(self, df=None, keys=None, axes=None, **kwargs):
