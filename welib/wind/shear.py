import pandas as pd
# Local 
import numpy as np
# --- Power law 
#    case (WindProfileType_PL)
def powerlaw(z, alpha, z_ref, U_ref):
    U = np.zeros_like(z)
    b = z>0
    # [IEC 61400-1 6.3.1.2 (10)]
    U[b] = U_ref * (z[b]/z_ref)**alpha
    return U

 
# --- Log law 
def loglaw(z, z_0, z_ref, U_ref):
#       if (.not. EqualRealNos(G3D%RefHeight, G3D%Z0) .and. Z > 0.0_ReKi) then
#          U = G3D%MeanWS*(LOG(Z/G3D%Z0))/(LOG(G3D%RefHeight/G3D%Z0))
#       else
#          U = 0.0_ReKi
#       end if
    z = np.asarray(z)
    U = np.zeros_like(z)
    b = z>0
    # [IEC 61400-1 6.3.1.2 (10)]
    U[b] = U_ref * (np.log(z[b]/z_0))/(np.log(z_ref/z_0))
    return U

# 
def linshear(z, slope, z_ref, U_ref, L_ref=None):
    if L_ref is not None:
        # OpenFAST model
#       U = U + G3D%MeanWS*G3D%VLinShr*(Z - G3D%RefHeight)/G3D%RefLength
        U = U_ref + slope*(z - z_ref)/L_ref
    else:
        U = U_ref + slope*(z - z_ref)
    return U



# --------------------------------------------------------------------------------}
# --- Curve fitting 
# --------------------------------------------------------------------------------{
def fit_powerlaw_u_alpha(z, u, z_ref=100, p0=(10,0.1)):
    """ 
    p[0] : u_ref
    p[1] : alpha
    """
    from collections import OrderedDict
    import scipy.optimize as so
    pfit, _ = so.curve_fit(lambda xx, *p : p[0] * (xx / z_ref) ** p[1], z, u, p0=p0)
    u_fit = pfit[0] * (z / z_ref) ** pfit[1]
    coeffs_dict=OrderedDict([('u_ref',pfit[0]),('alpha',pfit[1])])
    formula = '{u_ref} * (z / {z_ref}) ** {alpha}'
    fitted_fun = lambda xx: pfit[0] * (xx / z_ref) ** pfit[1]
    return u_fit, pfit, {'coeffs':coeffs_dict,'formula':formula,'fitted_function':fitted_fun}

def fit_powerlaw_alpha(z, u, z_ref, u_ref, p0=(0.1)):
    """ 
        x is z
        y is u
    p[0] : alpha
    """
    from collections import OrderedDict
    import scipy.optimize as so
    pfit, _ = so.curve_fit(lambda xx, *p : u_ref * (xx / z_ref) ** p[0], z, u, p0=p0)
    u_fit = u_ref * (z / z_ref) ** pfit[0]
    coeffs_dict=OrderedDict([('alpha',pfit[0])])
    formula = '{u_ref} * (z / {z_ref}) ** {alpha}'
    fitted_fun = lambda xx: u_ref * (xx / z_ref) ** pfit[0]
    return u_fit, pfit, {'coeffs':coeffs_dict,'formula':formula,'fitted_function':fitted_fun}


if __name__ == '__main__':
    import matplotlib.pyplot as plt
    z=np.linspace(0,400)
    z_ref = 150
    U_ref = 8
    fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
    fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
    ax.plot(powerlaw(z, 0.2, z_ref , U_ref)  , z, label='p law')
    ax.plot(powerlaw(z,-0.2, z_ref , U_ref)  , z, label='p law')
    ax.plot(loglaw  (z, 0.01 , z_ref , U_ref), z, label='log law')
    ax.plot(linshear(z, 0.01 , z_ref , U_ref), z, label='lin')
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.legend()
    plt.show()
