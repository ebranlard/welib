import numpy as np
from numpy.linalg import inv
from scipy.integrate import  solve_ivp #odeint


from .statespace import vec_interp, B_interp, StateMatrix


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

          qdot =             A        .q      +    B
    where 
       - M: mass matrix, constant or function of x
       - K: stiffness matrix, constant matrix or None
       - C: damping matrix, constant matrix or None
       - F: forcing term, function of t,x,xdot, or, tuple (time, Force time series)

       - q : full state vector

    """
    def __init__(self, M, C=None, K=None, F=None, x0=None, xdot0=None):

        # Functions that return the mass matrix and its inverse
        if hasattr(M, '__call__'):  
            self.M_is_func=True
            self._fM    = M
            self._fMinv = lambda x : inv(M(x))
        else:
            self.M_is_func=False
            self._Minv= inv(M) # compute inverse once and for all
            self._fM    = lambda x : M
            self._fMinv = lambda x : self._Minv
            # TODO TODO
            if K is None:
                self.K = M*0
            if C is None:
                self.C = M*0

        # Store
        self.M = M
        self.C = C
        self.K = K
        self.has_A      = C is not None or K is not None # do we have a constant state matrix?
        #self.has_A      = True



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
        
        """
        vTime = np.asarray(vTime)
        vF    = np.asarray(vF)
        if vF.shape[0]!=self.nDOF:
            raise Exception('Wrong first dimension for Force time series ({} instead of {} )'.format(vF.shape[0],self.nDOF))
        if vF.shape[1]!=len(vTime):
            raise Exception('Second dimension of Force time series does not match time dimension ({} instead of {} )'.format(vF.shape[1],len(vTime)))

        self._force_ts = vF
        self._time_ts  = vTime

    def setForceFunction(self,fn):
        """ 
        
        """
        self._force_fn = fn

    def Force(self,t,x=None,xdot=None):
        if hasattr(self,'_force_ts'):
            return vec_interp(t,self._time_ts,self._force_ts)
        elif hasattr(self,'_force_fn'):
            return self._force_fn(t,x,xdot)
        else:
            raise NotImplementedError('Please specify a time series of force first')

    def B(self,t,q):
        x  = q[0:self.nDOF]
        xd = q[self.nDOF:]
        if hasattr(self,'_force_ts'):
            return  B_interp(t,self._fMinv(x),self._time_ts,self._force_ts)

        elif hasattr(self,'_force_fn'):
            F          = self._force_fn(t,x,xd).ravel()
            nDOF       = len(F)
            B          = np.zeros((2*nDOF,1))
            B[nDOF:,0] = np.dot(self._fMinv(x),F)
            return B
        else:
            raise NotImplementedError('Please specify a time series of force first')


    def integrate(self,t_eval, method='RK45', y0=None, **options):
        """ Perform time integration of system 
            method: 'RK54', 'LSODA'  (see solve_ivp)
        """
        #
        if y0 is not None:
            self.setStateInitialConditions(y0)
            
        if self.has_A:
            A=self.A
            odefun = lambda t, q : np.dot(A,q)+self.B(t,q)
        else:
            odefun = self.dqdt
            #odefun = lambda t, q : self.B(t,q)

        res = solve_ivp(fun=odefun, t_span=[t_eval[0], t_eval[-1]], y0=self.q0, t_eval=t_eval, method=method, vectorized=True, **options)   
        # Store
        self.res    = res
        return res

    def dqdt(self, t, q):
        dqdt = np.zeros(q.shape)
        x  = q[0:self.nDOF]
        xd = q[self.nDOF:]

        dqdt[0:self.nDOF]=  xd
        Minv = self._fMinv(x)
        F = self._force_fn(t,x,xd)
        dqdt[self.nDOF:] =  Minv.dot(F)
        return dqdt



    @property
    def A(self):
        if self.M_is_func:
            raise Exception('A matrix not a property when A is a function')
        else:
            return StateMatrix(self._Minv,self.C,self.K)

    @property
    def nDOF(self):
        if self.M_is_func:
            x=np.zeros((5000,1))
            return len(self.M(x))
        else:
            return len(self.M)

    @property
    def nStates(self):
        return 2*self.nDOF

    def __repr__(self):
        s='<MechSystem object>\n'
        s+='Read-only attributes:\n'
        s+=' - nDOF: {}, nState:{} \n'.format(self.nDOF,self.nStates)
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
        s+=' - M: Mass Matrix  \n'
        s+=str(self.M)+'\n'
        s+=' - C: Damping Matrix  \n'
        s+=str(self.C)+'\n'
        s+=' - K: Stiffness Matrix  \n'
        s+=str(self.K)+'\n'
        s+=' - q0: Initial conditions (state) \n'
        s+=str(self.q0)+'\n'
        if self.has_A:
            s+=' - A:  \n'
            s+=str(self.A)+'\n'

        return s


    def plot(self,axes=None):
        """ Simple plot of states after time integration"""
        if self.res is None:
            raise Exception('Call integrate before plotting')

        if axes is None:
            import matplotlib.pyplot as plt
            fig,axes = plt.subplots( self.nStates,1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
            fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)

        n=self.nDOF
        for i,ax in enumerate(axes):
            if i<n:
                lbl=r'$x_{}$'.format(i+1)
            else:
                lbl=r'$\dot{x}_'+r'{}$'.format(i-n+1)
            ax.plot(self.res.t, self.res.y[i,:], label=lbl)
            ax.set_ylabel(lbl)
            #ax.legend()
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



