""" 
Dynamics of a cross section

Degrees of freedom:
    - x "flapwise"
    - y "edgewise"
    - theta "torsion" (opposite around z) OR phi_z : positive around z

We use a string to select the degrees of freedom that are active:
    sx = 'x,y,th' : all DOF active
    sx = 'x'      : only x active, etc.

INPUTS:
    Ux      = u['Ux'](t)    # (disturbed) wind speed (typically 0) [m/s]
    Uy      = u['Uy'](t)    # (disturbed) wind speed (contains "Omega r"  [m/s]
    theta_p = u['pitch'](t) # Total pitch angle (twist+pitch) [rad]

Coordinate system:         
                            ^ x_a
                            |
                            |
                            |             
                     ___.-------------...
 y_A         _....---*                   \
 <----- TE <        T       A       Q     |  LE
            *------......._             ./
                            *----------*

States:
  - q = [x, xd, xa] 
  - x : structural motions   : e.g.  (x, y, th )  when sx=(x,y,th)
  - xd: structural velocities: e.g.  (xd,yd,thd)  when sx=(x,y,th)
  - xa: aerodynamic states   : e.g.  (fs)         when sxa=(fs)

References:
  [1] Branlard, Jonkman (2023) The aeroelastic cross section [TODO]

"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from collections import OrderedDict
# Local
from welib.airfoils.Polar import Polar
from welib.system.statespace import StateSpace

# Dynamic stall
from welib.airfoils.DynamicStall import dynstall_oye_dxdt_simple
from welib.airfoils.DynamicStall import dynstall_mhh_dxdt_simple
from welib.airfoils.DynamicStall import dynstall_oye_output_simple
from welib.airfoils.DynamicStall import dynstall_mhh_outputs_simple
from welib.airfoils.DynamicStall import dynstall_oye_param_from_polar
from welib.airfoils.DynamicStall import dynstall_mhh_param_from_polar
from welib.airfoils.DynamicStall import dynstall_mhh_steady_simple
from welib.airfoils.DynamicStall import dynstall_oye_steady
# Dynamic inflow
from welib.dyninflow.DynamicInflow import dynflow_oye_dxdt_simple
from welib.dyninflow.DynamicInflow import dynflow_oye_steady_simple
from welib.dyninflow.DynamicInflow import tau1_dbemt, tau2_oye


# --------------------------------------------------------------------------------}
# --- Main Wrapper Class 
# --------------------------------------------------------------------------------{
class Section(object):
    """Section"""
    def __init__(self):
        #self.arg = arg
        # --- Main data
        self._r_bar = None # r/R [-]
        self._R     = None # R [m]
        self._beta  = 0 # [rad] Twist, negative about z, typically positive
        # Raw constant
        self._M33   = None # 3x3 mass matrix
        self._C33   = None # 3x3 damping matrix
        self._K33   = None # 3x3 stiffness matrix
        self._chord = None # chord
        self._rho   = None # airdensity
        self._x_AQ  = 0    # location of aerodynamic center "Q"
        self._y_AQ  = None 
        self._x_AT  = 0    # location of three quarter chord point "T"
        self._y_AT  = None 
        self._pol   = None # Polar, instance of Polar class
        self._ppol  = None # Polar parameters
        self._tau   = None # TODO stieg oye time constant
        # Induction / dynamic infloe
        self._a0      = 0
        self._ap0     = 0

#     def setup_nonlinear_model_p(M, C, K, sx='x,y,th', 
#             rho=1.225, chord=0.2, polarFilename='tjaere11_ds.csv', drag=False,  # Aero options
#             y_AQ=0, y_AT=None, x_AQ=0, x_AT=0,

        # Simulation data
        self.sx_sim   = None # Degree of freedom
        self.ds_model = None
        self.di_model = None
        self.p_sim    = None
        self.u_sim    = None
        self.x0_sim   = None     # Initial condition for simulation
        self.x0_sim_dict = None
        self.sys_sim = None
        self.res_sim = None
        self.df_sim = None
        #  
        self.u_op     = None

        
    def fromOpenFAST(self, fstFilename, r_bar=1.0):
        """ Setup section from an OpenFAST input file. """
        pass

    # --------------------------------------------------------------------------------}
    # --- Structural model
    # --------------------------------------------------------------------------------{
    def setMassMatrix(self, M33=None, m=0, J_zz=0, x_G=0, y_G=0):
        """ set the 3x3 mass matrix, either using M33 or individual components """
        self._M33 = massMatrix(m=m, J_zz=J_zz, x_G=x_G, y_G=y_G, M33=M33, sx='x,y,th')
        return self._M33

    def setStifMatrix(self, K33=None, kxx=0, kyy=0, kzz=0, kxy=0, kxz=0, kyz=0):
        """ set the 3x3 mass matrix, either using K33 or individual components """
        self._K33 = stifMatrix(kxx=kxx, kyy=kyy, kzz=kzz, kxy=kxy, kxz=kxz, kyz=kyz, K33=K33, sx='x,y,th')
        return self._K33

    def setDampMatrix(self, C33=None, cxx=0, cyy=0, czz=0, cxy=0, cxz=0, cyz=0):
        """ set the 3x3 mass matrix, either using M33 or individual components """
        self._C33 = dampMatrix(cxx=cxx, cyy=cyy, czz=czz, cxy=cxy, cxz=cxz, cyz=cyz, C33=C33, sx='x,y,th')
        return self._C33

    @property
    def r(self):
        return self._r_bar * self._R

    # --------------------------------------------------------------------------------}
    # --- Structure
    # --------------------------------------------------------------------------------{
    def setStrucDOF(self):
        #
        pass


    # --------------------------------------------------------------------------------}
    # --- Polar 
    # --------------------------------------------------------------------------------{
    def polarFromCSV(self, polarFilename, fformat=None):
        self._pol = Polar(polarFilename, fformat=fformat, radians=True, compute_params=True) # compute_params for DS
        # TODO TODO as a trigger
        ppol = polarParams(self._pol, chord=self._chord, tau=self._tau)
        self._ppol = ppol

    def plotPolar(self):
        pol    = self._ppol['Polar']
        fPolar = self._ppol['fPolar'] #['cl','cd','cm','fs','cl_inv','cl_fs'], radians=True)
        # p['alpha_0']  = alpha_0  # TODO HARMONIZATION WITH DS
        # p['Cl_slope'] = Cl_slope  # TODO HARMONIZATION WITH DS
        # p['alpha_range']     = None
        # p['alpha_range_lin'] = None

        alpha = np.linspace(pol.alpha[0], pol.alpha[-1], 100)
        M     = fPolar(alpha)


        fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
        fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
        ax.plot(alpha*180/np.pi, M[0,:], '-', label = 'cl')
        ax.plot(alpha*180/np.pi, M[1,:], '-', label = 'cd')
        ax.plot(alpha*180/np.pi, M[2,:], '-', label = 'cm')
        ax.plot(alpha*180/np.pi, M[3,:], '-', label = 'fs')
        ax.plot(alpha*180/np.pi, M[4,:], '--', label = 'cl_inv')
        ax.plot(alpha*180/np.pi, M[5,:], ':', label = 'cl_fs')

        ax.set_xlabel('angle of attack [deg]')
        ax.legend()
        #plt.show()
        return ax

    # --------------------------------------------------------------------------------}
    # --- Dynamic stall 
    # --------------------------------------------------------------------------------{
    def setDynStall(self, ds_model=None, **kwargs):
        self.ds_model = ds_model
        # TODO trigger


    # --------------------------------------------------------------------------------}
    # --- Dynamic inflow
    # --------------------------------------------------------------------------------{
    def setDynInflow(self, di_model=None, a0=0, ap0=0, **kwargs):
        self.di_model = di_model
        self._a0      = a0
        self._ap0     = ap0
        if di_model is None:

            pass
        elif di_model=='oye':
            # Check that main variables are set
            if self._r_bar is None:
                raise Exception('_r_bar must be set for Oye dyninflow model')
            if self._R is None:
                raise Exception('_R must be set for Oye dyninflow model')
        else:
            raise NotImplementedError()

    @property
    def di_tau1(self):
        if self.di_model is None:
            return np.nan
        if self.u_op is None:
            raise Exception('Operational input is None. Set it.')
        U0 = self.u_op['Ux']
        return tau1_dbemt(self._a0, self._R, U0)
    @property
    def di_tau2(self):
        if self.di_model is None:
            return np.nan
        return tau2_oye(self._r_bar, self.di_tau1)


    # --------------------------------------------------------------------------------}
    # --- Simulation 
    # --------------------------------------------------------------------------------{
    @property
    def nq(self):
        p = self.p_sim
        return len(p['sq'].split(','))

    @property
    def q0(self):
        return np.zeros(self.nq)

    @property
    def state_space(self):
        """ return state space model
        NOTE: safer to call it everytime in case number of DOF change
        """
        sys = StateSpace(dqdt=nonlinear_model, signature='t,q,u,p', verbose=False)
        return sys

    def calcOutput(self, t, q, u, p):
        S = nonlinear_model(t, q, u, p, calcOutput=True)
        return S

    # For convenience
    def Faero(self, t, q, u, p, calcOutput=False, m=None):
        return Faero(t, q, u, p, calcOutput=calcOutput, m=m)



    # TODO get rid of this
    def setParameters(self, sx_sim=None):
        """ setup parameters for a given simulation"""
        # TODO rething that ..
        #if sx_sim is not None:
        #if ds_model is not None:
        #if di_model is not None:
        self.sx_sim = sx_sim
        p = defaultParams(chord=self._chord, rho=self._rho, sx=self.sx_sim, ds=self.ds_model, di=self.di_model,
                M=self._M33, C=self._C33, K=self._K33)
        p['beta'] = self._beta
        if len(p['Iq'])==0:
            raise Exception('No states are present')

        # --- Dynamic inflow / induction
        p['a0']  = self._a0
        p['ap0'] = self._ap0
        p['di_tau1'] = self.di_tau1
        p['di_tau2'] = self.di_tau2

        # --- Aerodynamic parameters
        if self._y_AQ>0: 
            print('[WARN] y_AQ positive is unconventional')
        p['y_AQ'] = self._y_AQ
        if self._y_AT is None:
            p['y_AT'] = self._y_AQ+self._chord/2 # default is approximatively half a chord behind
        else:
            p['y_AT'] = self._y_AT
        p['x_AQ'] = self._x_AQ
        p['x_AT'] = self._x_AT
        if self._ppol is None:
            raise Exception('Polar parameters need to be set')
        p.update(self._ppol)
        #     #     p.update({'linModel':False, 'drag':drag})

        self.p_sim = p

    def setInitialConditions(self, q=None, x=0, y=0, th=0, xd=0, yd=0, thd=0, 
            fs=None, x_mhh=None,
            wx=0, wy=0, wxr=0, wyr=0,
            uop=None,
            equilibrium=False,
            di_eq=False,
            ds_eq=False,
            verbose=False
            ):

        p = self.p_sim
        if p is None:
            raise Exception('Call `setParameters` first')
        if q is not None:
            state0=q
        else:
            state0 = setup_nonlinear_model_x0(x=x, y=y, th=th, xd=xd, yd=yd, thd=thd, 
                    fs=fs, x_mhh=x_mhh, 
                    wx=wx, wy=wy, wxr=wxr, wyr=wyr,
                    p=self.p_sim)
        if verbose:
            state0_dict = OrderedDict((k,v) for k,v in zip(p['sq'].split(','), state0) )
            print('q0 (user input):', dict(state0_dict))

        # --- Equilibrim
        if uop is None:
            uop = self.u_sim

        if equilibrium:
            # We find full equilibrium
            state0 = self.equilibrium(x0=state0, u0=uop)
            if verbose:
                state0_dict = OrderedDict((k,v) for k,v in zip(p['sq'].split(','), state0) )
                print('q0 (full equil):', dict(state0_dict))
        else:
            if ds_eq or di_eq: 
                q = state0
                qx, qxd, qxa_ua, qxa_di = split_q(q, p['Iqxs'], p['Iqxsd'], p['Iqxa_ua'], p['Iqxa_di'])
                q_full, x, xd = inflate_q(q, Iq=p['Iq'])
                V_rel2_AC, phi_AC, alpha_AC, Vrel_AC, Uw_AC, Vel_AC, Wqs_AC, W_AC = aeroVarsSS(0, x, xd, qxa_di, p, self.u_sim, point='Q')
                V_rel2_ds, phi_ds, alpha_ds, Vrel_ds, Uw_ds, Vel_ds, Wqs_ds, W_ds = aeroVarsSS(0, x, xd, qxa_di, p, self.u_sim, point=p['ds_inflow_point'])
                if di_eq:
                    if p['dynamicInflowModel'] is None:
                        pass
                    elif p['dynamicInflowModel'] =='oye':
                        qxa_di[:2] = dynflow_oye_steady_simple(Wqs_AC[0], p['di_k'])
                        qxa_di[2:] = dynflow_oye_steady_simple(Wqs_AC[1], p['di_k'])
                    else:
                        raise NotImplementedError()

                if ds_eq:
                    if p['dynamicStallModel'] is None:
                        pass
                    elif p['dynamicStallModel'] == 'mhh':
                        qxa_ua = dynstall_mhh_steady_simple(np.sqrt(V_rel2_ds), alpha_ds, p)
                    elif p['dynamicStallModel'] == 'oye':
                        qxa_ua[0] = dynstall_oye_steady(alpha_ds, p)
                    else:
                        raise NotImplementedError()
                state0 = np.concatenate((qx, qxd, qxa_ua, qxa_di)) # NOTE: assumed order here
                if verbose:
                    print('q0 (ds di eq)  :', dict(state0_dict))

        state0_dict = OrderedDict((k,v) for k,v in zip(p['sq'].split(','), state0) )

        self.x0_sim = state0
        self.x0_sim_dict = state0_dict
        if verbose:
            print('q0 (final)     :', dict(state0_dict))
        return state0 

    def setConstantInputs(self, Ux=0, Uy=0, theta_p=0):
        """ setup inputs for the non linear model, constant inputs with time """
        self.u_sim = setup_nonlinear_model_u_cst(Ux=Ux, Uy=Uy, theta_p=theta_p)
        self.u_op = {'Ux':Ux, 'Uy':Uy, 'theta_p':0}
        return self.u_sim

    def integrate(self, t_eval, y0=None, p=None, u=None, calc='u,y', **options):
        #def integrate(self, t_eval, method='RK45', y0=None, p=None, u=None, calc='u,y,qd', xoffset=None, **options):
        if y0 is None:
            y0 = self.x0_sim.copy()
        if p is None:
            p = self.p_sim
        if u is None:
            u = self.u_sim

        self.sys_sim = self.state_space
        res, df = self.sys_sim.integrate(t_eval, p=p, u=u, y0=y0, calc=calc, **options)
        self.res_sim = res
        self.df_sim = df
        return res, df

    def saveSim(self, filename, df=None):
        from welib.tools.pandalib import remap_df
        # --- Scale
        ColMap={'alpha_AC':'{alpha_AC}*180/np.pi',
                'alpha_ds':'{alpha_ds}*180/np.pi',
                'thetat':  '{thetat}*180/np.pi',
                'theta':   '{theta}*180/np.pi',
                'pitch':   '{pitch}*180/np.pi',
                }
        if df is None:
            df = self.df_sim
        df = remap_df(df, ColMap)
        df.to_csv(filename, index=False, sep=',')

    def plot_states(self, *args, **kwargs):
        return self.sys_sim.plot_states(*args, **kwargs)

    # --------------------------------------------------------------------------------}
    # --- Linearization and equilibrium
    # --------------------------------------------------------------------------------{
    @property
    def dx(self):
        """ Perturtbation size for finite difference on states """
        if self.p_sim is None:
            raise Exception('Call setParameters before calling dx')
        dx = setup_nonlinear_model_dx(self.p_sim)
        return dx

    @property
    def du(self):
        """ Perturtbation size for finite difference on inputs"""
        du = setup_nonlinear_model_du() # Ux, Uy, theta_p
        return du

    def u0_default(self, u0=None, t=0):
        if u0 is None:
            u0 = self.u_sim(t)
            print('Section: using default input: {} at t={}'.format(u0,t))
        return u0

    def equilibrium(self, x0=None, u0=None, t=0, **kwargs):
        if x0 is None:
            x0=self.q0
        u0 = self.u0_default(u0, t)
        sys = self.state_space
        xop = sys.equilibrium(x0=x0, u0=u0, dx=self.dx, du=self.du, p=self.p_sim, **kwargs)
        return xop

    def linearize(self, x0=None, u0=None, t=0):
        u0 = self.u0_default(u0, t)
        if x0 is None:
            x0 = self.equilibrium(x0=self.x0_sim, u0=u0) 
            print('Section: linearize: using equilibrium x: {}'.format(x0))
        sys = self.state_space
        A, B, C, D = sys.linearize(x0, u0=u0, dx=self.dx, du=self.du, p=self.p_sim)
        return A, B, C, D

    def eigA(self, x0=None, u0=None, t=0, A=None, normQ=None, fullEV=False):
        from welib.tools.eva import eigA
        if A is None:
            A,_,_,_ = self.linearize(x0, u0=u0, t=t) 
        if len(self.p_sim['sx'])>0:
            nq2 = len(self.p_sim['sx'].split(','))
        else:
            nq2=0
        freq_d, zeta, Q, freq_0 = eigA(A, nq=nq2, fullEV=fullEV, normQ=normQ, sort=True)
        return freq_d, zeta, Q, freq_0 , A


    # --------------------------------------------------------------------------------}
    # --- Parametric/Campbell
    # --------------------------------------------------------------------------------{
    def campbell(self, WS, RPM, Pitch, A=None, AP=None, sx_sim=None):
        """ 
        - WS [deg]
        - Pitch [deg]
        - Pitch [deg]
        """
        # Default arguments
        if A is None:
            A=np.zeros_like(WS)
        if AP is None:
            AP=np.zeros_like(WS)

        # Prepare outputs
        AMat = []   # A matrices
        df   = None

        qop=None
        r = self.r
        for i,(U0, rpm, pitch, a0, ap0) in enumerate(zip(WS, RPM, Pitch, A, AP)):
            Omega = np.asarray(rpm) * 2*np.pi/60
            Ux      = U0
            Uy      = Omega * r
            theta_p = pitch *np.pi/180
            u0 = [Ux, Uy, theta_p]
            # TODO update interface
            self.setConstantInputs(Ux, Uy, theta_p=theta_p)
            self.setDynInflow(di_model=self.di_model, a0=a0, ap0=ap0)
            self.setDynStall(ds_model=self.ds_model)
            self.setParameters(sx_sim=sx_sim) # TODO
            # --- Frequencies
            if qop is None:
                qop=self.q0
            qop = self.equilibrium(x0=qop, u0=u0, tol=1e-8, verbose=False)
            A,_,_,_ = self.linearize(x0=qop, u0=u0)
            freq_d, zeta, Q, freq_0, A = self.eigA(A=A)
            # --- Store
            AMat.append(A)
            S = self.calcOutput(0, qop, self.u_sim, self.p_sim)
            S['WS']   = U0
            S['RPM']  = rpm
            S['theta'] = pitch
            for i,(f,z) in enumerate(zip(freq_0,zeta)):
                S['f'+str(i+1)] = f
                S['d'+str(i+1)] = z*100
            if i==0:
                df = pd.DataFrame(S).T
            else:
                dfloc = pd.DataFrame(S).T
                df = pd.concat((df,dfloc))
        return df, AMat



    # --------------------------------------------------------------------------------}
    # --- IO 
    # --------------------------------------------------------------------------------{
    def __repr__(self):
        s='<{} object>:\n'.format(type(self).__name__)
        s+='Parameters:\n'
        s+=' - rho    :            {} \n'.format(self._rho)
        s+=' - chord  :            {} \n'.format(self._chord)
        if self.p_sim is not None:
            s+=' - p_sim: simulation parameters:\n'
            s+='   - sx:            {} \n'.format(self.p_sim['sx'])
            s+='   - sq:            {} \n'.format(self.p_sim['sq'])
            s+='   - dynamicStallModel : {} \n'.format(self.p_sim['dynamicStallModel'])
            s+='   - dynamicInflowModel: {} \n'.format(self.p_sim['dynamicInflowModel'])
        if self.x0_sim is not None:
            #s+='- x0_sim:            {} \n'.format(self.x0_sim)
            s+='- x0_sim_dict        {} \n'.format(dict(self.x0_sim_dict))
        s+='Useful methods:\n'
        s+=' - setParameters\n'
        s+=' - setInitialConditions\n'
        s+=' - integrate(t, x0)\n'
        s+=' - equilibrium(x0, u0)\n'
        return s







# --------------------------------------------------------------------------------}
# --- Indices 
# --------------------------------------------------------------------------------{
def s2Iq(sq, sxs, sdi='', sua=''):
    """ 
    Return indices such that
      xs    = q[Ix]     # Structural motion, match sx
      xsd   = q[Ixd]    # Structual velocities, match sx
      xa    = q[Iqa]    # Aerodynamic states, match sua+sdi
      xa_ua = q[Ixa_ua] # Aerodynamic states
      xa_di = q[Ixa_di] # Aerodynamic states dynamic inflow

    INPUTS:
     - sq: list of all state names (structural and aero), comma separated
           Example: sq = 'x,xd,fs'
     - sx: list of structural state names, comma separated
           Example: sq = 'x'
    """
    sq     = sq.split(',')
    if len(sxs)>0:
        Ixs     = [sq.index(s) for s in sxs.split(',')]
        Ixsd    = [sq.index(s+'d') for s in sxs.split(',')]
    else:
        Ixs = []
        Ixsd = []
    if len(sua)>0:
        Ixa_ua = [sq.index(s) for s in sua.split(',')]
    else:
        Ixa_ua =[]
    if len(sdi)>0:
        Ixa_di = [sq.index(s) for s in sdi.split(',')]
    else:
        Ixa_di =[]
    return Ixs, Ixsd, Ixa_ua, Ixa_di
        

def sx2I(sx):
    """ from structural DOF string to indices """
    D={'x':0,'y':1,'th':2,'p':2}
    if len(sx)>0:
        return [D[s] for s in sx.split(',')]
    else:
        return []

def sq2I(sq, ds=None, di=None):
    D = {'x':0,'y':1,'th':2,'ph':2,'xd':3,'yd':4,'thd':5,'pd':5}
    iMax = 5
    if di is not None: 
        if di.lower()=='oye':
            Ddi = {'wxr': iMax+1 ,'wx':iMax+2, 'wyr':iMax+3, 'wy':iMax+4}
            iMax+=4
        else:
            raise NotImplementedError()
        D.update(Ddi)
    if ds is not None: 
        if ds.lower()=='oye':
            Dua = {'fs':iMax+1}
            iMax+=1
        elif ds.lower()=='mhh':
            Dua = {'x1':iMax+1,'x2':iMax+2,'x3':iMax+3, 'x4':iMax+4}
            iMax+=4
        else:
            raise NotImplementedError()
        D.update(Dua)
    return [D[s] for s in sq.split(',') if s in D.keys()]

def split_q(q, Iqxs, Iqxsd, Iqxa_ua, Iqxa_di):
    q= np.asarray(q)
    qx     = q[Iqxs]
    qxd    = q[Iqxsd]
    qxa_ua = q[Iqxa_ua]
    qxa_di = q[Iqxa_di]
    return qx, qxd, qxa_ua, qxa_di

def inflate_q(q, Iq):
    # Expanding q into "full" arrays including inactive DOFs
    nFull = max(6,np.max(Iq))
    q_full = np.zeros(nFull+1) # TODO default init state for states that are not zero!
    q_full[Iq] = q
    x      = q_full[:3]
    xd     = q_full[3:6]
#     xa     = q_full[6:]
    return q_full, x, xd


# --------------------------------------------------------------------------------}
# --- Aeroelastic model 
# --------------------------------------------------------------------------------{
def nonlinear_model(t, q, u, p, calcOutput=False):
    """ 
    Return right hand side of the first order model of 1-DOF airfoil

    INPUTS:
     - q: vector of states: (x,xd,fs)
     - t: evaluation time, scalar
     - u:  three options supported: 
            - Input dictionary, function of times: u[key](t) with key in 'Ux','Uy','pitch'
            - interpolant                           u(t)     =[Ux, Uy, pitch]
            - values                                u        =[Ux, Uy, pitch]
     - p: parameters
          dictionnary with keys: 
             rho, c
             fPolar, tau
             invM, K, C
    OUTPUTS:
     - dq/dt, dfs
    """
    # Compute variables useful to several other functions
    m = calcMisc(t, q, u, p)

    # Aerodynamic Force (length 3)
    Fa, aeroOut = Faero(t, q, u, p, calcOutput=calcOutput, m=m)

    # Inertial force (length 3)
    Fs = Fstruct(t, q, u, p, calcOutput=calcOutput, m=m)

    # Total force
    F = np.zeros(3)
    F[0] = Fa[0] + Fs[0]
    F[1] = Fa[1] + Fs[1]
    F[2] =-Fa[2] - Fs[2]

    # Derivative of states (structural and aero)
    dqxd = dxd_dt(t, m['qx'], m['qxd'], F[p['Ix']], p)
    dxa  = dxa_dt(t, m['x'], m['xd'], m['qxa_ua'], m['qxa_di'], p, u)
    dq = np.concatenate((m['qxd'], dqxd, dxa)) # NOTE: assumed order here

    if calcOutput:
        dq_full, xd, xdd = inflate_q(dq, Iq=p['Iq'])
        d = dict()
        ## Structural states
        d['x']       = m['x'][0]
        d['y']       = m['x'][1]
        d['thetat']  = m['x'][2]
        d['xd']      = m['xd'][0]
        d['yd']      = m['xd'][1]
        d['thetad']  = m['xd'][2]
        d['xdd']     = xdd[0]
        d['ydd']     = xdd[1]
        d['thetadd'] = xdd[2]
        # Inputs
        # (NOTE: computed with "calc='u')
        #m['Ux'], m['Uy'], m['theta_p'] = inputsAtTime(t, u)
        # Struct
        d['theta'] = m['theta']
        d['rho_x'] = m['rho_x']
        d['rho_y'] = m['rho_y']
        ## Aero
        d.update(aeroOut)
        if p['dynamicStallModel'] == 'oye':
            d['fs']  = m['qxa_ua'][0]
            d['dfs'] = dxa        [0]
        elif p['dynamicStallModel'] == 'mhh':
            d['x1_ds']  = m['qxa_ua'][0]
            d['x2_ds']  = m['qxa_ua'][1]
            d['x3_ds']  = m['qxa_ua'][2]
            d['x4_ds']  = m['qxa_ua'][3]
        return pd.Series(d)
    else:
        return dq

def defaultParams(sx=None, ds=None, di=None, 
        chord=np.nan, rho=np.nan,  # aero
        M=None, C=None, K=None # system matrices
        ):
    # --- Dictionary
    p ={}
    # --- Aero
    p['chord'] = chord   # section chord [m]
    p['rho']   = rho     # air density [kg/m^3]
    p['beta'] = 0        # Section twist "beta" [rad] (negative about z), typically positive
    p['x_AQ']   = np.nan # x coordinate of aerodynamic center from airfoil origin
    p['y_AQ']   = np.nan # y coordinate of aerodynamic center from airfoil origin
    p['x_AT']   = np.nan # x coordinate of 3/4 point from airfoil origin
    p['y_AT']   = np.nan # y coordinate of 3/4 point from airfoil origin
    # Polar
    p['fPolar']   = None   # Interpolant for polar data fPolar(alpha_rad) = 0:Cl, 1:Cd, 2:Cm, 3:Fs, 4:Cl_inv, 5:Cl_fs
    p['alpha_0']  = np.nan
    p['Cl_slope'] = np.nan
    p['drag']     = None   # Include drag in calculation of aerodynamic force 
    p['moment']   = None   # Include moment in calculation of aerodynamic force
    p['alpha_range']     = None
    p['alpha_range_lin'] = None
    # --- Dynamic stall
    p['dynamicStallModel'] = ds
    p['ds_inflow_point'] = 'T' # Where is the inflow taken for dynamic stall T: 3/4 , Q:1/4 chord
    # Dynamic stall Oye
    p['tau']          = np.nan  # Time constant Oye dynamic stall
    # Dynamic stall MHH
    # Airfoil parameters
    p['alpha0']     = None
    p['Cla']        = None
    # Polar functions
    p['F_st']  = None
    p['Cl_fs'] = None
    p['Cl']    = None
    p['Cd']    = None
    p['Cm']    = None
    # Dynamics constants
    p['Tf0'] = None
    p['Tp0'] = None
    p['A1']  = None
    p['A2']  = None
    p['b1']  = None
    p['b2']  = None
    p['alpha0_in_x1x2']  = None
    p['U_in_x1x2']       = None
    p['scale_x1_x2']     = None
    p['old_ClCd_dyn']    = None

    # --- Dynamic inflow model
    p['dynamicInflowModel'] = di
    p['a0']  = 0    # baseline axial induction
    p['ap0'] = 0    # baseline tangential induction
    p['di_tau1'] = 0    # 
    p['di_tau2'] = 0    # 
    p['di_k']    = 0.6  # 

    # --- Structural dynamics
    #p['x_G']   = np.nan # x coordinate of center of mass from airfoil origin
    #p['y_G']   = np.nan # y coordinate of center of mass from airfoil origin
    #p['J_xx_A'] = np.nan # Torsional inertia  TODO where
    #p['m']     = np.nan # airfoil (generalized) mass [kg] 
    p['M']    = None
    p['K']    = None
    p['C']    = None
    p['invM'] = None

    p['sq'] = None # List of degrees of freedom including derivatives

    # Options for non linear model only
    p['linModel'] = None


    # --- DOF
    nMax = 3
    if sx is None:
        sx='x,y,th'
    p['sx']  = sx.strip().strip(',').strip()
    if len(sx)>0:
        p['sq']  = p['sx'] + ',' + ','.join([s+'d' for s in p['sx'].split(',')])
    else:
        p['sq']  = ''

    p['sua'] = '' 
    if ds is not None:
        if ds.lower()=='oye':
            p['sua'] = 'fs'
            nMax += 1
        elif ds.lower()=='mhh':
            p['sua'] = 'x1,x2,x3,x4'
            nMax += 4
        else:
            raise NotImplementedError()
        p['sq']  += ',' +p['sua']

    p['sdi'] = ''
    if di is not None:
        if di.lower()=='oye':
            p['sdi'] = 'wxr,wx,wyr,wy'
            nMax += 0
        else:
            raise NotImplementedError()
        p['sq']  += ',' +p['sdi']
    p['sq'] = p['sq'].strip(',').replace(',,',',')

    p['Iq'] = sq2I(p['sq'], di=di, ds=ds)
    p['Ix'] = sx2I(p['sx'])
    p['Iqxs'], p['Iqxsd'], p['Iqxa_ua'], p['Iqxa_di'] = s2Iq(sq=p['sq'], sxs=p['sx'], sdi=p['sdi'], sua=p['sua'])
    p['nMax'] = nMax



    # --- System matrices
    if len(sx)>0:
        p['m'] = M[0,0]
        p['x_AG'] = M[1,2]/p['m'] # TODO TODO TODO
        p['y_AG'] =-M[0,2]/p['m'] # TODO TODO TODO
    else:
        p['m'] = 0
        p['x_AG'] = 0
        p['y_AG'] = 0
    p['M'] = massMatrix(M33=M, sx=sx)
    p['C'] = dampMatrix(C33=C, sx=sx)
    p['K'] = stifMatrix(K33=K, sx=sx)
    invM = np.linalg.inv(p['M'])
    p['invM'] = invM

    return p

# --------------------------------------------------------------------------------}
# --- Structural dynamics 
# --------------------------------------------------------------------------------{
def massMatrix(m=0, J_zz=0, x_G=0, y_G=0, sx='x,y,th', M33=None):
    if M33 is not None:
        M = M33.copy()
    else:
        M = np.zeros((3,3))
        M[0,0] = m
        M[1,1] = m
        M[2,2] = J_zz
        if 'th' in sx:
            M[0,2] =  m*y_G
            M[1,2] = -m*x_G
            M[2,0] =  m*y_G
            M[2,1] = -m*x_G
        else:
            M[0,2] = -m*y_G
            M[1,2] =  m*x_G
            M[2,0] = -m*y_G
            M[2,1] =  m*x_G
    I = sx2I(sx)
    M = M[np.ix_(I,I)]
    return M

def stifMatrix(kxx=0, kyy=0, kzz=0, kxy=0, kxz=0, kyz=0, sx='x,y,th', K33=None):
    """ Return stiffness matrix, based on individual coefficients, and degrees of freedom selected `sx``"""
    if K33 is not None:
        K = K33.copy()
    else:
        K = np.zeros((3,3))
        kyx = kxy
        kzx = kxz
        kzy = kyz
        K[0,:] = [kxx, kxy, kxz]
        K[1,:] = [kyx, kyy, kyz]
        K[2,:] = [kzx, kzy, kzz]
    I = sx2I(sx)
    K = K[np.ix_(I,I)]
    return K

def dampMatrix(cxx=0, cyy=0, czz=0, cxy=0, cxz=0, cyz=0, sx='x,y,th', C33=None):
    if C33 is not None:
        C = C33.copy()
    else:
        C = np.zeros((3,3))
        cyx = cxy
        czx = cxz
        czy = cyz
        C[0,:] = [cxx, cxy, cxz]
        C[1,:] = [cyx, cyy, cyz]
        C[2,:] = [czx, czy, czz]
    I = sx2I(sx)
    C = C[np.ix_(I,I)]
    return C

def dxd_dt(t, x, xd, F, p, m=None):
    # Structure acceleration
    RHS = - np.dot(p['C'], xd) - np.dot(p['K'], x) + F
    xdd = np.dot(p['invM'], RHS)
    return xdd

def calcMisc(t, q, u, p):
    """ Compute useful variables used by many subfunctions 
    In particular, split the state vector into structural states and aero states
    """
    m=dict()
    # Split state into positions and speeds (qx, qxd), uaero states (qxa_ua), dynamic inflow states (qxa_di)
    m['qx'], m['qxd'], m['qxa_ua'], m['qxa_di'] = split_q(q, p['Iqxs'], p['Iqxsd'], p['Iqxa_ua'], p['Iqxa_di'])

    # Structural states (length 3, even if not all DOFs are actice)
    m['q_full'], m['x'], m['xd'] = inflate_q(q, Iq=p['Iq'])

    # Orientation of the section
    m['Ux'], m['Uy'], m['theta_p'] = inputsAtTime(t, u)
    th  = m['x'][2]
    m['omega'] = m['xd'][2]
    m['theta'] = th + m['theta_p'] + p['beta'] 
    m['rho_x'] = (-p['x_AG']* np.sin(m['theta']) + p['y_AG']*np.cos(m['theta']) )
    m['rho_y'] = (-p['x_AG']* np.sin(m['theta']) + p['y_AG']*np.cos(m['theta']) )
    return m


# --------------------------------------------------------------------------------}
# --- Aerodynamics 
# --------------------------------------------------------------------------------{
def dxa_dt(t, x, xd, xa_ua, xa_di, p, u):
    """ 
    Return time derivative of aerodynamic states dxa/dt based on dynamic stall and dynamic inflow model 
    The aerodynamic states are assumed to be ordered as: [xa_ua, xa_ui]

    INPUTS:
     - x :  3-array of structural states    , x,  y , theta
     - xd:  3-array of structural velocities, xd, yd, thetad
     - xa_ua: array of unsteady aerodynamic states
     - xa_di: array of dynamic inflow states
     - p    : dictionary of parameters
     - u    : inputs, dictionary of function of time <<< Some work needed
    """


    dxa_ua = xa_ua*0
    dxa_di = xa_di*0

    omega = xd[2]

    # Inputs
    Ux, Uy, _ = inputsAtTime(t, u)

    # Main inflow variables at dynamic stall point 
    V_rel2_AC, phi_AC, alpha_AC, Vrel_AC, Uw_AC, Vel_AC, Wqs_AC, W_AC = aeroVarsSS(t, x, xd, xa_di, p, u, point='Q')
    V_rel2_ds, phi_ds, alpha_ds, Vrel_ds, Uw_ds, Vel_ds, Wqs_ds, W_ds = aeroVarsSS(t, x, xd, xa_di, p, u, point=p['ds_inflow_point'])

    # --- Dynamic stall (unsteady airfoil aerodynamics)
    if p['dynamicStallModel'] is None:
        pass

    elif p['dynamicStallModel']=='oye':    
        afCoef = p['fPolar'](alpha_ds) # 0:Cl, 1:Cd, 2:Cm, 3:Fs, 4:Cl_inv, 5:Cl_fs
        fs_st  = afCoef[3]
        fs = xa_ua[0]
        dxa_ua[0] = dynstall_oye_dxdt_simple(fs, fs_st, p['tau'])

    elif p['dynamicStallModel']=='mhh':    
        U_AC     = np.sqrt(V_rel2_AC)
        U_dot    = 0
        dxa_ua =  dynstall_mhh_dxdt_simple(t, xa_ua, U_AC, U_dot, omega, alpha_ds, p)

    else: 
        raise NotImplementedError('Dynamic stall model: {}'.format(p['dynamicStallModel']))


    # --- Dynamic inflow 
    if p['dynamicInflowModel'] is None:
        pass
    elif p['dynamicInflowModel']=='oye':
        if len(xa_di)==4:
            xa_di_x = xa_di[:2]
            xa_di_y = xa_di[2:]
            dxa_di_x = dynflow_oye_dxdt_simple(xa_di_x, Wqs_AC[0], p['di_tau1'], p['di_tau2'], k=p['di_k'])
            dxa_di_y = dynflow_oye_dxdt_simple(xa_di_y, Wqs_AC[1], p['di_tau1'], p['di_tau2'], k=p['di_k'])
            dxa_di[:2] =dxa_di_x
            dxa_di[2:] =dxa_di_y
        elif len(xa_di)==2:
            raise NotImplementedError()
        else:
            raise NotImplementedError()

    else: 
        raise NotImplementedError('Dynamic inflow model: {}'.format(p['dynamicInflowModel']))

    
    dxa = np.concatenate((dxa_ua, dxa_di)) # NOTE: assumed order

    return dxa


def inputsAtTime(t, u):
    """ Return inputs at given time """
    if u is None:
        raise Exception('u is None')
    theta_p = 0
    Ux      = 0
    Uy      = 0
    if hasattr(u,'keys'):
        if 'pitch' in u.keys():
            theta_p = u['pitch'](t)
        if 'Ux' in u.keys():
            Ux      = u['Ux'](t)
        if 'Uy' in u.keys():
            Uy      = u['Uy'](t)
    else:
        try:
            Ux, Uy, theta_p = u(t) 
        except:
            Ux, Uy, theta_p = u
    return Ux, Uy, theta_p

def polarParams(pol, chord, cl_lin_method='leastsquare', DS_constants='OpenFAST', tau=None):
    """ 
    Set aerodynamic parameters related to polars
    """
    # Return interpolant
    fPolar = pol.interpolant(variables=['cl','cd','cm','fs','cl_inv','cl_fs'], radians=True)

    p=dict()
    p['Polar'] = pol # backup
    p['fPolar'] = fPolar

    # Linear region
    linear_region = np.array([-5, 10])*np.pi/180
    Cl_slope, alpha_0 = pol.cl_linear_slope(window=linear_region, method=cl_lin_method, radians=True)
    #print('Cl_slope',Cl_slope, '[1/rad]  -    alpha_0', alpha_0*180/np.pi,'[deg]')

    p['alpha_0']  = alpha_0  # TODO HARMONIZATION WITH DS
    p['Cl_slope'] = Cl_slope  # TODO HARMONIZATION WITH DS
    p['alpha_range']     = None
    p['alpha_range_lin'] = None

    # Dynamic stall
    p.update(dynstall_mhh_param_from_polar(pol, chord, constants=DS_constants))
    p.update(dynstall_oye_param_from_polar(pol, tau=tau)) # TODO
    return p


def aeroVarsSS(t, x, xd, xa_di, p, u, point='Q'):
    """ 
    Returns main aerodynamic variables based on state space inputs

    INPUTS:
     - x  : full structural motions    (x , y , t )
     - xd : full structural velocities (xd, yd, td)
     - xa_di : aerodynamic states for dynamic inflow
     - u  : dictionary of inputs at time t, keys: Ux, Uy, theta_p
     - p: parameters, dictionnary with keys:
           rho, c,
           fPolar, dynamicStall
    OUTPUTS:
     - F: aerodynamic force in each structural DOF
    """
    # From states to local variables
    x, y ,th  = x
    xd_A, yd_A, omega = xd

    # Inputs
    Ux, Uy, theta_p = inputsAtTime(t, u)
    U = [Ux, Uy]

    # Total torsion of the section (typically negative)
    theta = th + theta_p + p['beta'] 

    # --- Structural velocities at point of interest (Q, or T)
    if point == 'Q':
        x_AP=p['x_AQ']
        y_AP=p['y_AQ']
    elif point == 'T':
        x_AP=p['x_AT']
        y_AP=p['y_AT']
    else:
        raise NotImplementedError()
    # Velocity at point P
    xd_P = xd_A + omega * (-x_AP* np.sin(theta) + y_AP*np.cos(theta) )
    yd_P = yd_A - omega * ( x_AP* np.cos(theta) + y_AP*np.sin(theta) )
    Vel = [xd_P, yd_P]

    # Induction 
    W, Wqs = inducedVelocities(xa_di, U[0], U[1], Vel[0], Vel[1], p)

    Vrel, Vrel2, phi, alpha =  aeroTriangle(U[0], U[1], W[0], W[1], Vel[0], Vel[1], theta, smallAngles=False)

    return Vrel2, phi, alpha, Vrel, U, Vel, Wqs, W


def inducedVelocities(xa_di, Ux, Uy, Velx, Vely, p):
    # Inflow without induction
    Vx = Ux - Velx
    Vy = Uy - Vely
    # Induced velocities
    Wqs = [-p['a0']* Vx,  p['ap0'] * Vy ]
    if len(xa_di)==0:
        W = Wqs
    elif len(xa_di)==2:
        raise NotImplementedError()
    elif len(xa_di)==4:
        W = [xa_di[1], xa_di[3]]
    else:
        raise NotImplementedError()
    return W, Wqs


def Fstruct(t, q, u, p, calcOutput=False, m=None):
    if m is None:
        m = calcMisc(t, q, u, p)
    Fs = np.zeros(3)
    Fs[0] = p['m'] * m['omega']**2 * m['rho_x']
    Fs[1] = p['m'] * m['omega']**2 * m['rho_y']
    #print(Fs)
    return Fs

def Faero(t, q, u, p, calcOutput=False, m=None):
    """ 
    Returns aerodynamic force (in inertial frame) based on states 

    INPUTS:
     - x  : structural motions   : e.g. (x , y , t )
     - xd : structural velocities: e.g. (xd, yd, td)
     - xa : aerodynamic states: e.g. (separation state for dynamic stall)
     - u: Input dictionary, function of times: u[key](t) with key in 'V0','pitch'
     - p: parameters, dictionnary with keys:
           rho, c,
           fPolar, dynamicStall
    OUTPUTS:
     - F: aerodynamic force in each structural DOF
     - aeroOut: dictionary of additional outputs
    """
    if m is None:
        m = calcMisc(t, q, u, p)
    x, xd, qxa_ua, qxa_di = m['x'], m['xd'], m['qxa_ua'], m['qxa_di'] 


    # Inputs
    Ux, Uy, theta_p = inputsAtTime(t, u)
    omega = xd[2]
    
    # Inflow
    V_rel2_AC, phi_AC, alpha_AC, Vrel_AC, Uw_AC, Vel_AC, Wqs_AC, W_AC = aeroVarsSS(t, x, xd, qxa_di, p, u, point='Q')
    V_rel2_ds, phi_ds, alpha_ds, Vrel_ds, Uw_ds, Vel_ds, Wqs_ds, W_ds = aeroVarsSS(t, x, xd, qxa_di, p, u, point=p['ds_inflow_point'])
    U_AC = np.sqrt(V_rel2_AC)

    # Quasi-steady Polar data at ds_inflow_point for now
    Cl_qs, Cd_qs, Cm_qs, fs_qs, Cl_inv, Cl_fs = aeroCoef(alpha_ds, p['fPolar'], drag=p['drag'], moment=p['moment'])

    # TODO TODO DS Wrapper
    if p['dynamicStallModel'] is None:
        Cl_dyn, Cd_dyn, Cm_dyn = Cl_qs, Cd_qs, Cm_qs

    elif p['dynamicStallModel']=='oye':
        # From states to local variables
        fs = qxa_ua[0]
        Cl_dyn, Cd_dyn, Cm_dyn = dynstall_oye_output_simple(fs, Cl_fs, Cl_inv, Cl_qs, Cd_qs, Cm_qs)

    elif p['dynamicStallModel']=='mhh':
        U_dot    = 0
        Cl_dyn, Cd_dyn, Cm_dyn = dynstall_mhh_outputs_simple(t, qxa_ua, U_AC, U_dot, omega, alpha_ds, p, calcOutput=False)
    else:
        raise NotImplementedError(p['dynamicStallModel'])

    # 
    if p['linModel']:
        Cl_lin = p['Cl_slope'] *(alpha_ds - p['alpha_0']) # TODO which alpha..
        Cl_dyn = Cl_lin
    
    # Compute aerodynamic loads at airfoil Origin (point A)
    # NOTE: we can use alpha_AC or theta, this is used to transfer the moment from Q to A
    #       the alpha used here needs to be consistent with phi. We need phi_X = theta + alpha_X
    F_A, F_Q = aeroLoads(p['rho'], p['chord'], V_rel2_AC, phi_AC, Cl_dyn, Cd_dyn, Cm_dyn, alpha_AC, x_AC=p['x_AQ'], y_AC=p['y_AQ'], smallAngles=p['linModel'])

    if p['linModel']:
        print('>>> Clqs:{:.2f} Cl:{:.2f} Cllin{:.2f} - phi:{:.1f} alpha:{:.1f}'.format(afCoef[0],Cl_dyn, Cl_lin,phi_AC*180/np.pi,alpha_ds*180/np.pi))
    #, alpha, Cl_dyn, fs_st

    aeroOut = None
    if calcOutput:
        # Inflow without induction
        Vx = Uw_AC[0] - Vel_AC[0]
        Vy = Uw_AC[1] - Vel_AC[1]
        if Vx==0:
            Vx=1e6
        if Vy==0:
            Vy=1e6
        aqs = -Wqs_AC[0]/Vx
        a   =   -W_AC[0]/Vx
        apqs = Wqs_AC[1]/Vy
        ap   =   W_AC[1]/Vy
        aeroOut  = {'alpha_AC':alpha_AC, 'U_AC':U_AC}
        aeroOut.update({'alpha_ds':alpha_ds, 'omega':omega})
        aeroOut.update({'Cl_qs':Cl_qs, 'Cd_qs':Cd_qs, 'Cm_qs':Cm_qs})
        aeroOut.update({'Cl_dyn':Cl_dyn, 'Cd_dyn':Cd_dyn, 'Cm_dyn':Cm_dyn})
        aeroOut.update({'Vrel_x_AC': Vrel_AC[0], 'Vrel_y_AC': Vrel_AC[1]})
        aeroOut.update({'Vrel_x_ds': Vrel_ds[0], 'Vrel_y_ds': Vrel_ds[1]})
        aeroOut.update({'Uw_x_AC': Uw_AC[0]  , 'Uw_y_AC':Uw_AC[1]})
        aeroOut.update({'Vel_x_AC':Vel_AC[0], 'Vel_y_AC':Vel_AC[1]})
        aeroOut.update({'Vel_x_ds':Vel_ds[0], 'Vel_y_ds':Vel_ds[1]})
        aeroOut.update({'Wqs_x_AC':Wqs_AC[0], 'Wqs_y_AC':Wqs_AC[1]})
        aeroOut.update({'Wdyn_x_AC':W_AC[0], 'Wdyn_y_AC':  W_AC[1]})
        aeroOut.update({'a':a, 'ap':ap})
        aeroOut.update({'aqs':aqs, 'apqs':apqs})
        aeroOut.update({'L': F_Q[0], 'D': F_Q[1], 'M': F_Q[2]})
        aeroOut.update({'Fx':F_A[0], 'Fy':F_Q[1], 'Mz':F_A[2]})
    return F_A, aeroOut


# --------------------------------------------------------------------------------}
# --- Aerodynamic, low level functions 
# --------------------------------------------------------------------------------{
def aeroLoads(rho, c, Vrel2, phi, Cl, Cd, Cm, alpha=None, theta=None, x_AC=0, y_AC=0, smallAngles=False):
    """ 
    Return Aerodynamic loads at airfoil origin
    low level function 

    NOTE: Vrel, phi need to be at the same point
          if alpha is provided it also need to be at the same point! 
    """
    L = 1/2*rho*c   *Vrel2*Cl
    D = 1/2*rho*c   *Vrel2*Cd
    M = 1/2*rho*c**2*Vrel2*Cm
    if smallAngles:
        Fx =   L       + D * phi
        Fy =  -L * phi + D
        Mz = 1/2*rho*c*Vrel2*( c*Cm -Cl*(y_AC + x_AC*alpha) + Cd*(x_AC-y_AC*alpha ))
    else:
        Fx =  L*np.cos(phi)+D*np.sin(phi)
        Fy = -L*np.sin(phi)+D*np.cos(phi)
        if theta is not None:
            Mz = M - Fx * (-x_AC * np.sin(theta)+y_AC*np.cos(theta)) + Fy * (x_AC*np.cos(theta) + y_AC*np.sin(theta))
        elif alpha is not None:
            Mz = 1/2*rho*c*Vrel2*( c*Cm -Cl*(y_AC*np.cos(alpha) + x_AC*np.sin(alpha)) + Cd*(x_AC*np.cos(alpha)-y_AC*np.sin(alpha) ))
        else:
            raise NotImplementedError('Provide alpha or theta')
    return np.array([Fx, Fy, Mz]), np.array([L,D,M])

def aeroTriangle(Ux, Uy, Wx, Wy, Velx, Vely, theta, smallAngles=False):
    """ 
    Return velocity triangle variables Vrel2, phi, alpha
    low level function
    Inputs may be at the aerodynamic center of three quarter point.

    INPUTS:
     - U: (disturbed) inflow/wind
     - W: induced velocities
     - Vel: elastic velocities
     - theta : total pitch angle (negative about z) [rad]
              theta = theta_t + theta_p + Beta
                      theta_t: torsion
                      theta_p: pitch angle
                      Beta: twist angle
              
    """
    # Relative velocity
    Vrelx = Ux - Velx + Wx
    Vrely = Uy - Vely + Wy
    Vrel2 = Vrelx**2 + Vrely**2
    # Flow angle
    if smallAngles:
        phi = Vrelx/Vrely
    else:
        phi = np.arctan2(Vrelx, Vrely) # [rad]
    # Angle of attack
    alpha = phi - theta    # [rad]

    # TODO
    if smallAngles:
        alpha_deg = alpha*180/np.pi # [deg]
        if alpha_deg>15 or alpha_deg<-15:
            print('>>> Alpha is beyond linear region {:.2f} deg'.format(alpha_deg))
    if alpha > np.pi:
        alpha = alpha-2*np.pi
    if alpha < -np.pi:
        alpha = alpha+2*np.pi
    # TODO TODO TODO
    # input polar doesn't not have data from -180 180..
    if alpha>np.pi/2:
        alpha=np.pi/2
    if alpha<-np.pi/2:
        alpha=-np.pi/2

    return [Vrelx, Vrely], Vrel2, phi, alpha

def polarInterpolant(filename, fformat=None, variables='cl,cd,cm,fs,cl_inv,cl_fs'):
    """ 
    Return a polar interpolant `fPolar` from a polar input file
    The function `fPolar` has the interface:
        afCoef = fPolar(alpha) 
    where afCoef is defined by the list of `variables` given as argument.
    By Default:
        # 0:Cl, 1:Cd, 2:Cm, 3:Fs, 4:Cl_inv, 5:Cl_fs
        Cl_qs, Cd_qs, Cm_qs, fs_qs, Cl_inv, Cl_fs = afCoef
    """
    pol = Polar(filename, fformat=fformat)
    fPolar = pol.interpolant(variables=variables.split(','), radians=True)
    return pol, fPolar

def aeroCoef(alpha, fPolar, drag=True, moment=True):
    """ 
    Return quasi steady aerodynamic coefficients for a given angle of attack
    INPUTS:
     - alpha : angle of attack [rad]
     - fPolar: interpolant function as returned by function polarInterpolant
    """
    # TODO Re interpolation
    afCoef = fPolar(alpha) # 0:Cl, 1:Cd, 2:Cm, 3:Fs, 4:Cl_inv, 5:Cl_fs
    # NOTE: order here depends on "variables" provided when generating the interpolant
    Cl_qs, Cd_qs, Cm_qs, fs_qs, Cl_inv, Cl_fs = afCoef
    if not drag:
        Cd_qs =0
    if not moment:
        Cm_qs = 0
    return Cl_qs, Cd_qs, Cm_qs, fs_qs, Cl_inv, Cl_fs 



# --------------------------------------------------------------------------------}
# --- Setup simulation
# --------------------------------------------------------------------------------{
def setup_nonlinear_model_p(M, C, K, sx='x,y,th', 
        rho=1.225, chord=0.2, polarFilename='tjaere11_ds.csv', drag=False,  # Aero options
        y_AQ=0, y_AT=None, x_AQ=0, x_AT=0,
        ds='oye', tau=0.08,  # tau: used for Oye, but should be taken from Polar!
        di=None):
    """ setup parameters for the non linear model"""
    # 
    p = defaultParams(chord=chord, rho=rho, sx=sx, ds=ds, di=di,
            M=M, C=C, K=K)

    # --- Aerodynamic parameters
    if y_AQ>0: 
        print('[WARN] y_AQ positive is unconventional')
    p['y_AQ'] = y_AQ
    if y_AT is None:
        p['y_AT'] = y_AQ+chord/2 # default is approximatively half a chord behind
    else:
        p['y_AT'] = y_AT
    p['x_AQ'] = x_AQ
    p['x_AT'] = x_AT

    # Read polar
    pol = Polar(polarFilename, fformat=None, radians=True, compute_params=True) # compute_params for DS
    ppol = polarParams(pol, chord=p['chord'], tau=tau)
    p.update(ppol)

#     # --- Dictionary
#     p.update({'linModel':False, 'drag':drag})
    return p


def setup_nonlinear_model_u_cst(Ux, Uy, theta_p):
    """ setup inputs for the non linear model, constant inputs with time """
    u = OrderedDict()
    u['Ux']    = lambda t: Ux
    u['Uy']    = lambda t: Uy
    u['pitch'] = lambda t: theta_p
    return u

def setup_nonlinear_model_x0(x=0, y=0, th=0, xd=0, yd=0, thd=0, fs=None, x_mhh=None, wx=0, wy=0, wxr=0, wyr=0, p=None):
    """ setup initial conditions for the non linear model"""
    xs  = np.array([x, y, th])
    xsd = np.array([xd, yd, thd])

    if p['dynamicStallModel'] is None:
        xua = np.array([])
    elif p['dynamicStallModel'].lower() == 'oye':
            if fs is None:
                # TODO
                print('TODO figure out angle of attack')
                # fs_i  = p['fPolar'](alphai*np.pi/180)[3] # initial position
                xua = np.array([0])
            else:
                xua = np.array([fs])
    elif p['dynamicStallModel'].lower() == 'mhh':
        if x_mhh is None:
            xua = np.array([0, 0, 0, 0])
        else:
            raise NotImplementedError()
    else:
        NotImplementedError()

    if p['dynamicInflowModel'] is None:
        xdi = np.array([])
    elif p['dynamicInflowModel'].lower() == 'oye':
        # TODO
        xdi = np.array([wxr,wx,wyr,wy]) # Whatch out for order here
    else:
        NotImplementedError()
    q_full = np.concatenate((xs,xsd,xdi,xua))
    state0 = q_full[p['Iq']]
    return state0



def setup_nonlinear_model_dx(p):
    nFull = max(6,np.max(p['Iq']))
    dq_full = np.zeros(nFull+1) 
    c = p['chord']

    dq_full[:3]  = [0.001*c,0.001*c,0.001]   # x,y,th
    dq_full[3:6] = [0.001,  0.001,  0.001]   # xd,yd,thd
    dq_full[6:]  = 0.001  #... TODO

    dq_full[:3]  = [0.01*c,0.01*c,0.01]   # x,y,th
    dq_full[3:6] = [0.01,  0.01,  0.01]   # xd,yd,thd
    dq_full[6:]  = 0.01  #... TODO
    #####
    dq_full[:3]  = [0.1*c,0.1*c, 1*np.pi/180]   # x,y,th
    dq_full[3:6] = [5. *c,5. *c, 40.]   # xd,yd,thd # TODO
    dq_full[6:]  =  0.1  #... TODO
    dx = dq_full[p['Iq']]
    return dx

def setup_nonlinear_model_du(*args, **kwargs):
    du = np.array((0.01,0.01,0.001)) # Ux, Uy, theta_p
    #du = np.array((0.05,0.05, 1*np.pi/180)) # Ux, Uy, theta_p
    return du


if __name__ == '__main__':
    #M = massMatrix(10,100,1,2,sx='x,y,th')
    #print(M)
#     polarInterpolant('data/tjaere11_ds.csv')

    Ix, Ixd, Ixa_ua, Ixa_di = s2Iq(sq='x,xd,fs,wxr', sx='x', sdi='', sua='fs')
    print('Ix    ',Ix)
    print('Ixd   ',Ixd   )
    print('Ixa_ua',Ixa_ua)
    print('Ixa_di',Ixa_di)


