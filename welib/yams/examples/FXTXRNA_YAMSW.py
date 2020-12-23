from welib.yams.models.FTNSB_sympy import *
from welib.yams.models.FTNSB_sympy_symbols import *
from welib.yams.yams_sympy_tools import cleantex
from welib.yams.yams_sympy_tools import smallAngleApprox, cleantex, subs_no_diff , cleanPy
from sympy import symbols 
from collections import OrderedDict

# model_name = 'F2T0RNA'
# model_name = 'F2T0RNA_fnd'
# model_name = 'F2T0N0S1'
# model_name = 'F2T1RNA'
# model_name = 'F2T1RNA_fnd'
# model_name = 'F2T1N0S1'
# model_name = 'F0T2RNA'
# model_name = 'F0T1RNA'
# model_name = 'F0T2N0S1'
# model_name = 'F6T0RNA'
# 
# model_name = 'F6T0RNA_fnd'
# model_name = 'F6T0N0S1'
# model_name = 'F6T1RNA'
#model_name = 'F6T1RNA_fnd'
# model_name = 'F6T1N0S1'
# model_name = 'F6T2N0S1'

#Models=['F2T0RNA'  , 'F2T0RNA_fnd', 'F2T0N0S1' , 'F2T1RNA'  , 'F2T1N0S1' , 'F0T2RNA'  , 'F0T2N0S1' , 'F6T0RNA' , 'F6T0RNA_fnd', 'F6T0N0S1' , 'F6T1RNA'  , 'F6T1N0S1' , 'F6T2N0S1' ] 
# Models=['F2T0RNA'  , 'F2T0RNA_fnd']
# Models=['F6T1RNA']


# S = YAMSModel().load('_pickle/{}.pkl'.format(model_name))
# print(S[0])
# s, params, inputs, sdofs  = cleanPy(S, varname='FF', dofs = symbols('x, phi_y'))
# print(s)

# 
# 
# #for model_name in Models:
# 
# 
# #
opts=dict()
opts['rot_elastic_type']='SmallRot' #<<< Very important, 'SmallRot', or 'Body', will affect the rotation matrix
# opts['rot_elastic_type']='Body' #<<< Very important, 'SmallRot', or 'Body', will affect the rotation matrix
opts['rot_elastic_smallAngle']=False #<<< Very important, will perform small angle approx: sin(nu q) = nu q and nu^2=0 !!! Will remove all nu^2 and nu^3 terms!! Not recommended, removes part of RNA "Y" inertia from mass matrix
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

elif model_name == 'F2T0RNA_fnd':
    # --- F2T0RNA_fnd
    opts['mergeFndTwr'] = False   # combined foudantion/floater and tower together, or two bodies
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

elif model_name == 'F2T1RNA_fnd':
    # --- F2T1RNA
    opts['mergeFndTwr'] = False    # combined foudantion/floater and tower together, or two bodies
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
elif model_name == 'F0T1RNA':
    # --- F0T1RNA
    opts['mergeFndTwr'] = True    # combined foudantion/floater and tower together, or two bodies
    opts['yaw']         = 'zero'  # 'fixed', 'zero', or 'dynamic' if a DOF
    opts['tilt']        = 'fixed' # 'fixed', 'zero', or 'dynamic' if a DOF
    opts['tiltShaft']   = True    # OpenFAST tilts shaft not nacelle
    opts['linRot']      = True    # Linearize rotations matrices from the beginning
    opts['Mform']       ='symbolic', # or 'TaylorExpanded'
    opts['twrDOFDir']   = ['x','y','x','y']  # Order in which the flexible DOF of the tower are set

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

elif model_name == 'F6T1RNA_fnd':
    # --- F6T1RNA_fnd
    opts['mergeFndTwr'] = False    # combined foudantion/floater and tower together, or two bodies
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


#opts['tilt']    = 'zero'

# --- Esthetics Replacements for python
replaceDict={}
replaceDict['theta_tilt']=('tilt',None)



# --- Create model, solve equations and preform small angle approximation
model = get_model(model_name, **opts)
#model.nac.noInertia()
#model.nac.noMass()
model.kaneEquations(Mform='TaylorExpanded')
# ---
extraSubs=model.shapeNormSubs
print('Extra Subs:  ', extraSubs        )
print('Small angles:', model.smallAngles)



# --- Linearization of non linear equations
model.linearize(noAcc=True, noVel=False, extraSubs=extraSubs)

# --- Small angle approximation and linearization
model.smallAngleApprox(model.smallAngles, extraSubs)
model.smallAngleApproxEOM(model.smallAngles, extraSubs)
model.smallAngleLinearize(noAcc=True, noVel=False, extraSubs=extraSubs)


# --- Export
model.smallAngleSaveTex(folder='_tex', variables=['MM','FF','M','C','K','B'])
model.savePython(folder='_py' , variables=['MM','FF','MMsa','FFsa','M','C','K','B','Msa','Csa','Ksa','Bsa'], replaceDict=replaceDict, extraSubs=extraSubs)

# 
# Subs for "linearized model", DOFs set to 0
# extraSubsLin=[]
# try:
#     extraSubsLin+=[(v,0) for v in model.smallAnglesFnd]
# except:
#     pass
# try:
#     extraSubsLin+=[(v,0) for v in model.twr.q] 
# except:
#     pass
# print('Subs lin no DOF: ',extraSubsLin)
# model.smallAngleLinearize(noAcc=True, noVel=True, extraSubs=extraSubs)
# model.smallAngleSaveTex(folder='_tex', variables=['M','C','K','B'], prefix='noVel',extraHeader='NoVelocity: ',header=False)
# 
# 
# model.smallAngleSaveTex(folder='_tex', variables=['M','C','K','B'], prefix='noVelnoDOF',extraHeader='NoVelocity NoDOF: $'+cleantex(extraSubsLin)+'$',header=False, extraSubs=extraSubsLin)
