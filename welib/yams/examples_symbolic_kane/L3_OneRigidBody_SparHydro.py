import os
import numpy as np    
import matplotlib.pyplot as plt
import importlib
# yams
from welib.yams.models.OneRigidBody_sympy import *  # get_model_one_body, x, phi_x, xd, omega_x, x_BG, Jxx_B, M_B, etc
import welib.weio as weio
from welib.yams.windturbine import FASTWindTurbine
from welib.tools.clean_exceptions import *
MyDir=os.path.dirname(__file__)

create=True
runSim=True
modelname = 'B101010_hydro'

fstFilename = os.path.join(MyDir, '../../../data/Spar/Main_Spar_ED_HydroExample.fst')

if create:
    model = get_model_one_body(modelname, linRot=False, orderMM=1, orderH=1, 
                               J_cross=False, CG_on_z=True, J_at_Origin=True)

    # --- Extra Loads
    body = model.body
    ref  = model.ref
    z_BB= symbols('z_BB')
    z_BM= symbols('z_BM')
    z_B0= symbols('z_B0')
    K_Mx, K_My, K_Mz          = symbols('K_x_M, K_y_M, K_z_M') # Mooring restoring
    K_Mphix, K_Mphiy, K_Mphiz = symbols('K_phi_x_M, K_phi_y_M, K_phi_z_M') # Mooring restoring
    C_Mx, C_My, C_Mz          = symbols('C_x_M, C_y_M, C_z_M') # Mooring restoring
    C_Mphix, C_Mphiy, C_Mphiz = symbols('C_phi_x_M, C_phi_y_M, C_phi_z_M') # Mooring restoring
    # Buoyancy force
    P_B = body.origin.locatenew('P_B', z_BB * body.frame.z) # <<<< Measured from T
    P_M = body.origin.locatenew('P_M', z_BM * body.frame.z) # <<<< Measured from T
    P_0 = body.origin.locatenew('P_0', z_B0 * body.frame.z) # <<<< Measured from T
    P_O = body.origin
    #P_B = body.origin
    P_M = body.origin
    F_hx = dynamicsymbols('F_hx') # Buoyancy force, function to time 
    F_hy = dynamicsymbols('F_hy') # Buoyancy force, function to time 
    F_hz = dynamicsymbols('F_hz') # Buoyancy force, function to time 
    M_hx = dynamicsymbols('M_hx') # Buoyancy force, function to time 
    M_hy = dynamicsymbols('M_hy') # Buoyancy force, function to time 
    M_hz = dynamicsymbols('M_hz') # Buoyancy force, function to time 
    #F_B = Function('F_B')(dynamicsymbols._t, phi_y)
    model.addForce(body,  P_0,        F_hx * ref.frame.x + F_hy * ref.frame.y + F_hz * ref.frame.z)
    model.addMoment(body, body.frame, M_hx * ref.frame.x + M_hy * ref.frame.y + M_hz * ref.frame.z)
    # fr=0; Mr=0;
    # fr += -K_Mx * x *ref.frame.x #if bDOFs[0] else 0
    # fr += -K_My * y *ref.frame.y #if bDOFs[1] else 0
    # fr += -K_Mz * z *ref.frame.z #if bDOFs[2] else 0
    # fr += -C_Mx * x.diff(dynamicsymbols._t) *ref.frame.x #if bDOFs[0] else 0
    # fr += -C_My * y.diff(dynamicsymbols._t) *ref.frame.y #if bDOFs[1] else 0
    # fr += -C_Mz * z.diff(dynamicsymbols._t) *ref.frame.z #if bDOFs[2] else 0
    # Mr += -K_Mphix * phi_x *ref.frame.x #if bDOFs[3] else 0
    # Mr += -K_Mphiy * phi_y *ref.frame.y #if bDOFs[4] else 0
    # Mr += -K_Mphiz * phi_z *ref.frame.z #if bDOFs[5] else 0
    # Mr += -C_Mphix * phi_x.diff(dynamicsymbols._t) *ref.frame.x #if bDOFs[3] else 0
    # Mr += -C_Mphiy * phi_y.diff(dynamicsymbols._t) *ref.frame.y #if bDOFs[4] else 0
    # Mr += -C_Mphiz * phi_z.diff(dynamicsymbols._t) *ref.frame.z #if bDOFs[5] else 0
    # model.addForce(body, P_M,  fr)
    # model.addMoment(body, body.frame, Mr)

    # Export
    extraSubs    =[(omega_x,0),(omega_y,0),(omega_z,0)]
    model.exportPackage(path='_{}'.format(modelname), extraSubs=extraSubs, smallAngles=None, replaceDict=None, pathtex='_tex/{}'.format(modelname), fullPage=True)

# --- Run non linear and linear simulation using a FAST model as input
if runSim:
    # --- Import the python module that was generated
    model_pkg = importlib.import_module('_'+modelname)
    # --- Load the wind turbine model, and extract relevant parameters "p"
    WT = FASTWindTurbine(fstFilename, twrShapes=[0,2], nSpanTwr=50)
    p = WT.yams_parameters(flavor='onebody', J_at_Origin=True)
    #print(p.keys())


    # --- Perform time integration
    if os.path.exists(fstFilename.replace('.fst','.outb')): 
        dfFS = weio.read(fstFilename.replace('.fst','.outb')).toDataFrame()
        time =dfFS['Time_[s]'].values
    elif os.path.exists(fstFilename.replace('.fst','.out')): 
        dfFS = weio.read(fstFilename.replace('.fst','.out')).toDataFrame()
        time =dfFS['Time_[s]'].values
    else:
        time = np.linspace(0,50,1000)
        dfFS = None

    info = model_pkg.info()
    # --- Inputs
    u=dict()
    u['F_hx'] = lambda t,q=None,qd=None: np.interp(t, time, dfFS['HydroFxi_[N]'])
    u['F_hz'] = lambda t,q=None,qd=None: np.interp(t, time, dfFS['HydroFzi_[N]'])
    u['M_hy'] = lambda t,q=None,qd=None: np.interp(t, time, dfFS['HydroMyi_[N-m]'])
    # --- Linear model inputs
    nu = len(info['su'])
    uop=dict()
    uop['F_hx'] = np.mean(dfFS['HydroFxi_[N]'].values)*0
    #uop['F_hz'] = np.mean(dfFS['HydroFzi_[N]'].values)   
    uop['F_hz'] = p['M_B']*p['g']
    uop['M_hy'] = np.mean(dfFS['HydroMyi_[N-m]'].values)*0

    du = np.zeros((nu, len(time)))
    du[0,:] = dfFS['HydroFxi_[N]'].values     - uop['F_hx']
    du[1,:] = dfFS['HydroFzi_[N]'].values     - uop['F_hz']  #- p['M_B']*p['g']
    du[2,:] = dfFS['HydroMyi_[N-m]'].values   - uop['M_hy']
    du[0,:] *= 0
    du[2,:] *= 0

    #uop=None
    #du=None
    qop  = None
    qdop = None
    # Using mean as op
#     qop  = np.array([np.mean(dfFS[c]) for c in WT.q_channels])
#     qdop = np.array([np.mean(dfFS[c]) for c in WT.qd_channels])*0
    #print('q0' ,WT.q0)
    #print('qop',qop)


    resNL, sysNL = WT.simulate_py    (model_pkg, p, time, u=u)
    #resLI, sysLI = WT.simulate_py_lin(model_pkg, p, time, uop=uop, du=du, qop=qop)
    dfNL = sysNL.toDataFrame(WT.channels, WT.FASTDOFScales)
    #dfLI = sysLI.toDataFrame(WT.channels, WT.FASTDOFScales, q0=qop)

    #print('sysLI\n',sysLI)

    # --- Simple Plot
    fig,axes = plt.subplots(len(dfNL.columns)-1, 1, sharey=False, figsize=(6.4,6.8)) # (6.4,4.8)
    fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
    for i,ax in enumerate(axes):
        chan=dfNL.columns[i+1]
        ax.plot(dfNL['Time_[s]'], dfNL[chan], '-'  , label='non-linear')
        #ax.plot(dfLI['Time_[s]'], dfLI[chan], '--' , label='linear')
        if dfFS is not None:
            ax.plot(dfFS['Time_[s]'], dfFS[chan], 'k:' , label='OpenFAST')
        ax.set_xlabel('Time [s]')
        ax.set_ylabel(chan)
        ax.tick_params(direction='in')
    ax.legend()
#     else:
#         dfNL=None
#         dfLI=None
#         dfFS=None
#     return dfNL, dfLI, dfFS
    plt.show()
# 

# 
# model.kdeqsSubs
# model.coordinates
# #model.body.inertia= (inertia(model.body.frame, J_O, J_O, J_zz), model.body.origin)
# # --- Compute Kane's equations of motion
# print('Extra Subs:')
# extraSubs
# EOM = model.to_EOM(extraSubs=extraSubs)
# EOM.mass_forcing_form()
# EOM.EOM
