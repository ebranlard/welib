import numpy as np
from .Polar import Polar as Pol


# --------------------------------------------------------------------------------}
# --- MHH dynamic stall 
# --------------------------------------------------------------------------------{
A1_Jones, A2_Jones, b1_Jones, b2_Jones = 0.165 ,0.335 ,0.0455 ,0.3

def dynstall_mhh_param_from_polar(P, chord, tau_chord=None ,Tf=None,Tp=None,A1=A1_Jones,A2=A2_Jones,b1=b1_Jones,b2=b2_Jones,Jones=False,FAST=False):
    if not isinstance(P,Pol):
        raise Exception('Input should be an instance of the `Polar` class')
    if not P._radians :
        raise Exception('MHH dynamic stall implemented for polars in radians only')

    if tau_chord is None and Tf is None:
        raise Exception('Provide `Tf` or provide `tau_chord`')
    if tau_chord is None and Tp is None:
        raise Exception('Provide `Tp` or provide `tau_chord`')

    # R.T Jones approximation to Wagner's function (Jones 1938)
    #Cl_wag_Jones=1-0.165*npexp(-0.0455*tau_t)-0.335*np.exp(-0.3*tau_t);
    if Jones: 
        A1, A2, b1, b2 = A1_Jones, A2_Jones, b1_Jones, b2_Jones
    if FAST: # FAST default values
        A1, A2, b1, b2 = 0.3, 0.7, 0.14, 0.53 

    p=dict()
    # Airfoil parameters
    p['alpha0']     = P._alpha0
    p['Cla']        = P._linear_slope
    p['chord']      = chord
    # Polar functions
    p['F_st']  = P.f_st_interp
    p['Cl_fs'] = P.cl_fs_interp
    p['Cl']    = P.cl_interp
    p['Cd']    = P.cd_interp
    p['Cm']    = P.cm_interp
    # Dynamics constants
    p['Tf']  = Tf if Tf is not None else   3 * tau_chord 
    p['Tp']  = Tp if Tp is not None else 1.7 * tau_chord
    p['A1']  = A1
    p['A2']  = A2
    p['b1']  = b1
    p['b2']  = b2
    return p

def dynstall_mhh_dxdt(t,x,u,p):
    # States
    x1=x[0] # Downwash memory term 1
    x2=x[1] # Downwash memory term 2
    x3=x[2] # Clp', Lift coefficient with a time lag to the attached lift coeff
    x4=x[3] # f'' , Final separation point function
    # Inputs
    U         = u['U'](t)
    U_dot     = u['U_dot'](t)
    alpha     = u['alpha'](t)
    alpha_dot = u['alpha_dot'](t)
    alpha_34  = u['alpha_34'](t)
    # Parameters
    alpha0 = p['alpha0']
    Cla    = p['Cla']
    c      = p['chord']
    Tf     = p['Tf']
    Tp     = p['Tp']
    A1     = p['A1']
    A2     = p['A2']
    b1     = p['b1']
    b2     = p['b2']
    F_st   = p['F_st']
    Cl_fs  = p['Cl_fs']
    Cd     = p['Cd']
    # Variables derived from inputs
    Tu     = c/(2*U)                                         # Eq. 23
    # Variables derived from states
    alphaE  = alpha_34*(1-A1-A2)+ x1 + x2                    # Eq. 12
    Clp     = Cla * (alphaE-alpha0) + np.pi * Tu * alpha_dot # Eq. 13
    alphaF  = x3/Cla+alpha0                                  # p. 13
    fs_aF   = F_st(alphaF)                                    # p. 13

    # State equation
    xdot = [0]*4
    xdot[0] = -1/Tu * (b1 + c * U_dot/(2*U**2)) * x1 + b1 * A1 / Tu * alpha_34
    xdot[1] = -1/Tu * (b2 + c * U_dot/(2*U**2)) * x2 + b2 * A2 / Tu * alpha_34
    xdot[2] = -1/Tp                             * x3 + 1/Tp * Clp
    xdot[3] = -1/Tf                             * x4 + 1/Tf * fs_aF
    return xdot

def dynstall_mhh_steady(t,u,p):
    # Inputs
    U         = u['U'](t)
    U_dot     = u['U_dot'](t)
    alpha     = u['alpha'](t)
    alpha_dot = u['alpha_dot'](t)
    alpha_34  = u['alpha_34'](t)
    # Parameters
    c      = p['chord']
    alpha0 = p['alpha0']
    Cla    = p['Cla']
    A1     = p['A1']
    A2     = p['A2']
    b1     = p['b1']
    b2     = p['b2']
    F_st   = p['F_st']
    # Variables derived from inputs
    Tu     = c/(2*U)                                      # Eq. 23
    # Steady states
    x1     = A1*alpha_34
    x2     = A2*alpha_34
    alphaE = alpha_34*(1-A1-A2)+ x1 + x2 # Eq. 12
    x3     = Cla * (alphaE-alpha0)
    alphaF = x3/Cla+alpha0               # p. 13
    x4     = F_st(alphaF)
    return [x1,x2,x3,x4]

def dynstall_mhh_outputs(t,x,u,p):
    # States
    x1=x[0] # Downwash memory term 1
    x2=x[1] # Downwash memory term 2
    x3=x[2] # Clp', Lift coefficient with a time lag to the attached lift coeff
    x4=x[3] # f'' , Final separation point function
    # Inputs
    U         = u['U'](t)
    U_dot     = u['U_dot'](t)
    alpha     = u['alpha'](t)
    alpha_dot = u['alpha_dot'](t)
    alpha_34  = u['alpha_34'](t)
    # Parameters
    alpha0 = p['alpha0']
    Cla    = p['Cla']
    c      = p['chord']
    Tf     = p['Tf']
    Tp     = p['Tp']
    A1     = p['A1']
    A2     = p['A2']
    b1     = p['b1']
    b2     = p['b2']
    F_st   = p['F_st']
    Cl_fs  = p['Cl_fs']
    Cl     = p['Cl']
    Cd     = p['Cd']
    Cm     = p['Cm']

    #Cd0 = fCd(alpha0)
    #a_st = ??
    # Variables derived from inputs
    Tu     = c/(2*U)                                      # Eq. 23
    # Variables derived from states
    alphaE = alpha_34*(1-A1-A2)+ x1 + x2                  # Eq. 12
    faE = F_st(alphaE)
    DeltaCdfpp = (np.sqrt(faE)-np.sqrt(x4))/2 - (faE-x4)/4
    #ast_x4  = (fCm(x4)  - fCm(alpha0))/Cl(x4)
    #ast_faE = (fCm(faE) - fCm(alpha0))/Cl(faE)
    #DeltaCmfpp = (fa_st(x4) - fa_st(faE))
    DeltaCmfpp = 0 # <<<<<<<<<TODO
    # Outputs
    Cl_dyn =  Cla * (alphaE-alpha0)*x4 + Cl_fs(alphaE)*(1-x4) + np.pi*Tu*alpha_dot  # <<< ADDED ALPHA0?
    Cd_dyn =  Cd(alphaE) + (alpha-alphaE)*Cl_dyn + (Cd(alphaE)-Cd(alpha0))*DeltaCdfpp
    Cm_dyn =  Cm(alphaE) + Cl_dyn*DeltaCmfpp - np.pi/2*Tu*alpha_dot
    return Cl_dyn, Cd_dyn, Cm_dyn


# --------------------------------------------------------------------------------}
# --- Oye's dynamic stall 
# --------------------------------------------------------------------------------{
def dynstall_oye_param_from_polar(P,tau=None,tau_chord=None):
    if tau_chord is None and tau is None:
        raise Exception('Provide `tau` or provide `tau_chord`')
    p=dict()
    p['tau']   = 3*tau_chord if tau is None else tau
    p['F_st']  = P.f_st_interp
    p['Clinv'] = P.cl_inv_interp
    p['Clfs']  = P.cl_fs_interp
    return p

def dynstall_oye_dxdt(t,fs,u,p):
    """ d(fs)/dt = 1/tau (fs_st - fs) """
    alpha   = u['alpha'](t)
    f_st    = p['F_st'](alpha)
    return 1/p['tau'] * (f_st - fs)

def dynstall_oye_output(t,fs,u,p):
    alpha   = u['alpha'](t)
    Clfs    = p['Clfs'](alpha)
    Clinv   = p['Clinv'](alpha)
    Cl      = fs*Clinv+(1-fs)*Clfs               
    return Cl
