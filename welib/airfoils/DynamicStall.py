import numpy as np


# --------------------------------------------------------------------------------}
# --- MHH dynamic stall 
# --------------------------------------------------------------------------------{
def dyna_stall_mhh_dxdt(t,x,u,p):
    # States
    x1=x[0] # Downwash memory term 1
    x2=x[1] # Downwash memory term 2
    x3=x[2] # Clp', Lift coefficient with a time lag to the attached lift coeff
    x4=x[3] # f'' , Final separation point function
    # Inputs
    U         = u['fU'](t)
    U_dot     = u['fU_dot'](t)
    alpha     = u['falpha'](t)
    alpha_dot = u['falpha_dot'](t)
    alpha_34  = u['falpha_34'](t)
    # Parameters
    alpha0 = p['alpha0']
    Cla    = p['Cla']
    Tf     = p['Tf']
    Tp     = p['Tp']
    b1     = p['b1']
    b2     = p['b2']
    A1     = p['A1']
    A2     = p['A2']
    c      = p['chord']
    fF_st   = p['f_st_fun']
    fCl_fs  = p['cl_fs_fun']
    fCd     = p['cd_fun']
    # Variables derived from inputs
    Tu     = c/(2*U)                                         # Eq. 23
    # Variables derived from states
    alphaE  = alpha_34*(1-A1-A2)+ x1 + x2                    # Eq. 12
    Clp     = Cla * (alphaE-alpha0) + np.pi * Tu * alpha_dot # Eq. 13
    alphaF  = x3/Cla+alpha0                                  # p. 13
    fs_aF = fF_st(alphaF)                                    # p. 13
    #f_st_aF = fF_st(alphaE) # HACK

    # State equation
    xdot = [0]*4
    xdot[0] = -1/Tu * (b1 + c * U_dot/(2*U**2)) * x1 + b1 * A1 / Tu * alpha_34
    xdot[1] = -1/Tu * (b2 + c * U_dot/(2*U**2)) * x2 + b2 * A2 / Tu * alpha_34
    xdot[2] = -1/Tp                             * x3 + 1/Tp * Clp
    xdot[3] = -1/Tf                             * x4 + 1/Tf * fs_aF
    return xdot

def dyna_stall_mhh_steady(t,u,p):
    # Inputs
    U         = u['fU'](t)
    U_dot     = u['fU_dot'](t)
    alpha     = u['falpha'](t)
    alpha_dot = u['falpha_dot'](t)
    alpha_34  = u['falpha_34'](t)
    # Parameters
    alpha0 = p['alpha0']
    Cla    = p['Cla']
    b1     = p['b1']
    b2     = p['b2']
    A1     = p['A1']
    A2     = p['A2']
    c      = p['chord']
    fF_st  = p['f_st_fun']
    # Variables derived from inputs
    Tu     = c/(2*U)                                      # Eq. 23
    # Steady states
    x1     = A1*alpha_34
    x2     = A2*alpha_34
    alphaE = alpha_34*(1-A1-A2)+ x1 + x2 # Eq. 12
    x3     = Cla * (alphaE-alpha0)
    alphaF = x3/Cla+alpha0               # p. 13
    x4     = fF_st(alphaF)
    return [x1,x2,x3,x4]

def dyna_stall_mhh_outputs(t,x,u,p):
    # States
    x1=x[0] # Downwash memory term 1
    x2=x[1] # Downwash memory term 2
    x3=x[2] # Clp', Lift coefficient with a time lag to the attached lift coeff
    x4=x[3] # f'' , Final separation point function
    # Inputs
    U         = u['fU'](t)
    U_dot     = u['fU_dot'](t)
    alpha     = u['falpha'](t)
    alpha_dot = u['falpha_dot'](t)
    alpha_34  = u['falpha_34'](t)
    # Parameters
    alpha0 = p['alpha0']
    Cla    = p['Cla']
    Tf     = p['Tf']
    Tp     = p['Tp']
    b1     = p['b1']
    b2     = p['b2']
    A1     = p['A1']
    A2     = p['A2']
    c      = p['chord']
    fF_st   = p['f_st_fun']
    fCl_fs  = p['cl_fs_fun']
    fCd     = p['cd_fun']

    #Cd0 = fCd(alpha0)
    #a_st = ??
    # Variables derived from inputs
    Tu     = c/(2*U)                                      # Eq. 23
    # Variables derived from states
    alphaE = alpha_34*(1-A1-A2)+ x1 + x2                  # Eq. 12
    # Outputs
    Cl_dyn =  Cla * (alphaE-alpha0)*x4 + fCl_fs(alphaE)*(1-x4) + np.pi*Tu*alpha_dot  # <<< ADDED ALPHA0?
    #Cd_dyn =  fCd(alphaE) + (alpha-alphaE)*Cl_dyn + 
    return Cl_dyn


# --------------------------------------------------------------------------------}
# --- Oye's dynamic stall 
# --------------------------------------------------------------------------------{
def dyna_stall_oye_dxdt(t,fs,u,p):
    """ d(fs)/dt = 1/tau (fs_st - fs) """
    alpha   = u['falpha'](t)
    f_st    = p['f_st_fun'](alpha)
    return 1/p['tau'] * (f_st - fs)

def dyna_stall_oye_output(t,fs,u,p):
    alpha   = u['falpha'](t)
    Clfs    = p['fClfs'](alpha)
    Clinv   = p['fClinv'](alpha)
    Cl      = fs*Clinv+(1-fs)*Clfs               
    return Cl
