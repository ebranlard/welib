import numpy as np
from numpy.linalg import inv
from scipy.integrate import  solve_ivp #odeint
from .statespace import vec_interp, B_interp

# --------------------------------------------------------------------------------}
# --- Simple statespace functions
# --------------------------------------------------------------------------------{
def lti_state_space_function(t, x, u, p):
    return p['A'].dot(x) + p['B'].dot(u)

def lti_output_function(t, x, u, p):
    return p['C'].dot(x) + p['D'].dot(u)


# --------------------------------------------------------------------------------}
# --- Linear State Space system
# --------------------------------------------------------------------------------{
class LinearStateSpace():
    def __init__(self,A,B,C=None,D=None,q0=None):
        self.A=A
        self.B=B
        self.C=C
        self.D=C
        self.setStateInitialConditions(q0)

    def setStateInitialConditions(self,q0=None):
        self.q0 = np.zeros(self.nStates)
        if q0 is not None:
            if len(q0)!=self.nStates:
                raise Exception('Wrong dimension for q0 ({} instead of {} )'.format(len(q0),self.nStates))
            self.q0 = q0

    def setInputTimeSeries(self,vTime,vU):
        vTime = np.asarray(vTime)
        vU    = np.asarray(vU)
        if vU.shape[0]!=self.nInputs:
            raise Exception('Wrong first dimension for Inputs time series ({} instead of {} )'.format(vU.shape[0],self.nInputs))
        if vU.shape[1]!=len(vTime):
            raise Exception('Second dimension of Input time series does not match time dimension ({} instead of {} )'.format(vU.shape[1],len(vTime)))

        self._input_ts = vU
        self._time_ts  = vTime

    def setInputFunction(self,fn):
        self._input_fn = fn

    def Inputs(self,t,x=None):
        if hasattr(self,'_input_ts'):
            return vec_interp(t,self._time_ts,self._input_ts)
        elif hasattr(self,'_input_fn'):
            return self._input_fn(t,x)
        else:
            raise NotImplementedError('Please specify a time series of inputs first')


    def integrate(self, t_eval, method='RK4', y0=None, **options):
        #
        if y0 is not None:
            self.setStateInitialConditions(y0)
            
        odefun = lambda t, q : np.dot(self.A, q) + np.dot(self.B, self.Inputs(t,q))
        #odefun = lambda t, q : self.RHS(t,q)

        res = solve_ivp(fun=odefun, t_span=[t_eval[0], t_eval[-1]], y0=self.q0, t_eval=t_eval, method=method, vectorized=False, **options)   

        if self.nOuputs>0:
            print('>>> Do something for outputs')


        return res

    def dqdt(self, t, q):
        # NOTE: this can cause issues if q is not flat
        return np.dot(self.A, q) + np.dot(self.B, self.Inputs(t,q))

    def RHS(self,t,q):
        return self.dqdt(t,q)

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

    def __repr__(self):
        s='<MechSystem object>\n'
        s+='Read-only attributes:\n'
        s+=' - nState:{} \n'.format(self.nStates)
        if hasattr(self,'_force_ts'):
            if len(self._time_ts)>1:
                dt=self._time_ts[1]-self._time_ts[0]
            else:
                dt=np.nan
            s+='Force time series \n'
            s+=' - Time: [{} ... {}],  dt: {}, n: {} \n'.format(self._time_ts[0],self._time_ts[-1],dt,len(self._time_ts))
            s+=' - Force t0  : {} \n'.format(self._force_ts[:,0])
            s+=' - Force tend: {} \n'.format(self._force_ts[:,-1])
        s+='Attributes:\n'
        s+=' - A: State-State Matrix  \n'
        s+=str(self.A)+'\n'
        s+=' - B: State-Input Matrix  \n'
        s+=str(self.B)+'\n'
        s+=' - q0: Initial conditions (state) \n'
        s+=str(self.q0)+'\n'
        return s
