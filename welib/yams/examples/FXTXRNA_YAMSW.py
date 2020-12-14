from welib.yams.models.FTNSB_sympy import *
from welib.yams.models.FTNSB_sympy_symbols import *
from welib.yams.yams_sympy_tools import cleantex

opts=dict()
model_name = 'F2T0RNA'
model_name = 'F2T0N0S1'
model_name = 'F2T1RNA'
model_name = 'F2T1N0S1'
# model_name = 'F0T2RNA'
# model_name = 'F0T2N0S1'
# model_name = 'F6T0RNA'

# model_name = 'F6T0RNA_fnd'
# model_name = 'F6T0N0S1'
# model_name = 'F6T1RNA'
# model_name = 'F6T1N0S1'
# model_name = 'F6T2N0S1'


#
opts['rot_elastic_type']='SmallRot' #<<< Very important, 'SmallRot', or 'Body', will affect the rotation matrix
opts['rot_elastic_smallAngle']=True #<<< Very important, will perform small angle approx: sin(nu q) = nu q and nu^2=0 !!! Will remove all nu^2 and nu^3 terms!! might not be recommended
opts['orderMM'] = 1     #< order of taylor expansion for Mass Matrix
opts['orderH']  = 1     #< order of taylor expansion for H term

# --------------------------------------------------------------------------------}
# --- 2DOF floater 
# --------------------------------------------------------------------------------{
if model_name == 'F2T0RNA':
    # --- F2T0RNA
    opts['mergeFndTwr'] = True    # combined foudantion/floater and tower together, or two bodies
    opts['yaw']         = 'zero'  # 'fixed', 'zero', or 'dynamic' if a DOF
    opts['tilt']        = 'fixed' # 'fixed', 'zero', or 'dynamic' if a DOF
    opts['tiltShaft']   = True    # OpenFAST tilts shaft not nacelle
    opts['linRot']      = True    # Linearize rotations matrices from the beginning

elif model_name == 'F2T0N0S1':
    # --- F2T0N0S1
    opts['mergeFndTwr'] = True    # combined foudantion/floater and tower together, or two bodies
    opts['yaw']         = 'zero'  # 'fixed', 'zero', or 'dynamic' if a DOF
    opts['tilt']        = 'fixed' # 'fixed', 'zero', or 'dynamic' if a DOF
    opts['tiltShaft']   = True    # OpenFAST tilts shaft not nacelle
    opts['linRot']      = True    # Linearize rotations matrices from the beginning

elif model_name == 'F2T1RNA':
    # --- F2T1RNA
    opts['mergeFndTwr'] = True    # combined foudantion/floater and tower together, or two bodies
    opts['yaw']         = 'zero'  # 'fixed', 'zero', or 'dynamic' if a DOF
    opts['tilt']        = 'fixed' # 'fixed', 'zero', or 'dynamic' if a DOF
    opts['tiltShaft']   = True    # OpenFAST tilts shaft not nacelle
    opts['linRot']      = True    # Linearize rotations matrices from the beginning
    opts['Mform']       ='symbolic', # or 'TaylorExpanded'
    opts['twrDOFDir']   = ['x','y','x','y']  # Order in which the flexible DOF of the tower are set

elif model_name == 'F2T1N0S1':
    # --- F2T1N0S1
    opts['mergeFndTwr'] = True    # combined foudantion/floater and tower together, or two bodies
    opts['yaw']         = 'zero'  # 'fixed', 'zero', or 'dynamic' if a DOF
    opts['tilt']        = 'fixed' # 'fixed', 'zero', or 'dynamic' if a DOF
    opts['tiltShaft']   = True    # OpenFAST tilts shaft not nacelle
    opts['linRot']      = True    # Linearize rotations matrices from the beginning
    opts['Mform']       ='symbolic', # or 'TaylorExpanded'
    opts['twrDOFDir']   = ['x','y','x','y']  # Order in which the flexible DOF of the tower are set

# --------------------------------------------------------------------------------}
# --- Tower only
# --------------------------------------------------------------------------------{
elif model_name == 'F0T2RNA':
    # --- F0T2RNA
    opts['mergeFndTwr'] = True    # combined foudantion/floater and tower together, or two bodies
    opts['yaw']         = 'zero'  # 'fixed', 'zero', or 'dynamic' if a DOF
    opts['tilt']        = 'fixed' # 'fixed', 'zero', or 'dynamic' if a DOF
    opts['tiltShaft']   = True    # OpenFAST tilts shaft not nacelle
    opts['linRot']      = True    # Linearize rotations matrices from the beginning
    opts['Mform']       ='symbolic', # or 'TaylorExpanded'
    opts['twrDOFDir']   = ['x','y','x','y']  # Order in which the flexible DOF of the tower are set

elif model_name == 'F0T2N0S1':
    # --- F0T2N0S1
    opts['mergeFndTwr'] = True    # combined foudantion/floater and tower together, or two bodies
    opts['yaw']         = 'zero'  # 'fixed', 'zero', or 'dynamic' if a DOF
    opts['tilt']        = 'fixed' # 'fixed', 'zero', or 'dynamic' if a DOF
    opts['tiltShaft']   = True    # OpenFAST tilts shaft not nacelle
    opts['linRot']      = True    # Linearize rotations matrices from the beginning
    opts['Mform']       ='symbolic', # or 'TaylorExpanded'
    opts['twrDOFDir']   = ['x','y','x','y']  # Order in which the flexible DOF of the tower are set

# --------------------------------------------------------------------------------}
# ---6DOF floater
# --------------------------------------------------------------------------------{
elif model_name == 'F6T0RNA':
    # --- F6T0RNA
    opts['mergeFndTwr'] = True    # combined foudantion/floater and tower together, or two bodies
    opts['yaw']         = 'zero'  # 'fixed', 'zero', or 'dynamic' if a DOF
    opts['tilt']        = 'fixed' # 'fixed', 'zero', or 'dynamic' if a DOF
    opts['tiltShaft']   = True    # OpenFAST tilts shaft not nacelle
    opts['linRot']      = True    # Linearize rotations matrices from the beginning

elif model_name == 'F6T0RNA_fnd':
    # --- F6T0RNA_fnd
    opts['mergeFndTwr'] = False   # combined foudantion/floater and tower together, or two bodies
    opts['yaw']         = 'zero'  # 'fixed', 'zero', or 'dynamic' if a DOF
    opts['tilt']        = 'fixed' # 'fixed', 'zero', or 'dynamic' if a DOF
    opts['tiltShaft']   = True    # OpenFAST tilts shaft not nacelle
    opts['linRot']      = True    # Linearize rotations matrices from the beginning

elif model_name == 'F6T0N0S1':
    # --- F6T0N0S1
    model_name          = 'F6T0N0S1'
    opts['mergeFndTwr'] = True    # combined foudantion/floater and tower together, or two bodies
    opts['yaw']         = 'zero'  # 'fixed', 'zero', or 'dynamic' if a DOF
    opts['tilt']        = 'fixed' # 'fixed', 'zero', or 'dynamic' if a DOF
    opts['tiltShaft']   = True    # OpenFAST tilts shaft not nacelle
    opts['linRot']      = True    # Linearize rotations matrices from the beginning

elif model_name == 'F6T1RNA':
    # --- F6T1RNA
    opts['mergeFndTwr'] = True    # combined foudantion/floater and tower together, or two bodies
    opts['yaw']         = 'zero'  # 'fixed', 'zero', or 'dynamic' if a DOF
    opts['tilt']        = 'fixed' # 'fixed', 'zero', or 'dynamic' if a DOF
    opts['tiltShaft']   = True    # OpenFAST tilts shaft not nacelle
    opts['linRot']      = True    # Linearize rotations matrices from the beginning
    opts['Mform']       ='symbolic', # or 'TaylorExpanded'
    opts['twrDOFDir']   = ['x','y','x','y']  # Order in which the flexible DOF of the tower are set

elif model_name == 'F6T1N0S1':
    # --- F6T1N0S1
    opts['mergeFndTwr'] = True    # combined foudantion/floater and tower together, or two bodies
    opts['yaw']         = 'zero'  # 'fixed', 'zero', or 'dynamic' if a DOF
    opts['tilt']        = 'fixed' # 'fixed', 'zero', or 'dynamic' if a DOF
    opts['tiltShaft']   = True    # OpenFAST tilts shaft not nacelle
    opts['linRot']      = True    # Linearize rotations matrices from the beginning
    opts['Mform']       ='symbolic', # or 'TaylorExpanded'
    opts['twrDOFDir']   = ['x','y','x','y']  # Order in which the flexible DOF of the tower are set

elif model_name == 'F6T2N0S1':
    # --- F6T2N0S1
    opts['mergeFndTwr'] = True    # combined foudantion/floater and tower together, or two bodies
    opts['yaw']         = 'zero'  # 'fixed', 'zero', or 'dynamic' if a DOF
    opts['tilt']        = 'fixed' # 'fixed', 'zero', or 'dynamic' if a DOF
    opts['tiltShaft']   = True    # OpenFAST tilts shaft not nacelle
    opts['linRot']      = True    # Linearize rotations matrices from the beginning
    opts['Mform']       ='symbolic', # or 'TaylorExpanded'
    opts['twrDOFDir']   = ['x','y','x','y']  # Order in which the flexible DOF of the tower are set


else:
    raise Exception('Model {}'.format(model_name))


opts['tilt']    = 'zero'


# --- Create model, solve equations and preform small angle approximation
model = get_model(model_name, **opts)
#model.nac.noInertia()
#model.nac.noMass()
model.kaneEquations()
# ---
extraSubs=model.shapeNormSubs
print('Extra Subs:  ', extraSubs        )
print('Small angles:', model.smallAngles)



# --- Small angle approximation
model.smallAngleApprox(model.smallAngles, extraSubs)
model.smallAngleApproxEOM(model.smallAngles, extraSubs)
model.smallAngleLinearize(noAcc=True, noVel=False, extraSubs=extraSubs)


#model.smallAngleApprox(model.smallAnglesFnd, extraSubs) # NOTE: JyyRNA is second order in vu1T


# --- "linearization"
# Subs for "linearized model", DOFs set to 0
extraSubsLin=[]
try:
    extraSubsLin+=[(v,0) for v in model.smallAnglesFnd]
except:
    pass
try:
    extraSubsLin+=[(v,0) for v in model.twr.q] 
except:
    pass
print('Subs lin no DOF: ',extraSubsLin)

# --- Export
model.smallAngleSaveTex(folder='tex', variables=['MM','FF','M','C','K','B'])

model.smallAngleLinearize(noAcc=True, noVel=True, extraSubs=extraSubs)
model.smallAngleSaveTex(folder='tex', variables=['M','C','K','B'], prefix='noVel',extraHeader='NoVelocity: ',header=False)

model.smallAngleSaveTex(folder='tex', variables=['M','C','K','B'], prefix='noVelnoDOF',extraHeader='NoVelocity NoDOF: $'+cleantex(extraSubsLin)+'$',header=False, extraSubs=extraSubsLin)
