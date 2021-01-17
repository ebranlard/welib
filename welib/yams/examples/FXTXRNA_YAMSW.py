""" 
TODO:
 - Option for RefZ (rotate z_OT with phis)
 - Test shaft nacelle, potentially add new body for that..
 - Test simulations with time varying force
 - Force linearization
 - Blade flexibility

"""
from welib.yams.models.OneRigidBody_sympy import get_model_one_body
from welib.yams.models.FTNSB_sympy import get_model
#from welib.yams.models.FTNSB_sympy_symbols import *
from welib.yams.yams_sympy_tools import cleantex
from welib.yams.yams_sympy_tools import smallAngleApprox, cleantex, subs_no_diff , cleanPy
from sympy import symbols,Symbol
from collections import OrderedDict

# model_name = 'F2T0RNA'
model_name = 'F2T0RNA_fnd'
# model_name = 'F2T0N0S1_fnd'
# model_name = 'F2T1RNA'
# model_name = 'F2T1RNA_fnd'
# model_name = 'F3T1RNA_fnd'
# model_name = 'F5T1RNA_fnd'
# model_name = 'F2T1N0S1'
# model_name = 'F0T2RNA'
# model_name = 'F0T1RNA'
# model_name = 'F0T2N0S1'
# model_name = 'F6T0RNA'
# model_name = 'F5T0N0S1_fnd'
# 
# model_name = 'F6T0RNA_fnd'
# model_name = 'F6T0N0S1'
# model_name = 'F6T1RNA'
#model_name = 'F6T1RNA_fnd'
# model_name = 'F6T1N0S1'
# model_name = 'F6T2N0S1'

#model_name='F000101T0N0S1_fnd'

# model_name='F000101T0RNA_fnd'
# model_name='B000101'
# model_name='B000011'
# model_name='B000111'


Models=[ 'F2T0RNA_fnd' ,'F2T0N0S1_fnd' ,'F2T1RNA_fnd' ,'F3T1RNA_fnd' ,'F5T1RNA_fnd' ,'F2T1N0S1' ,'F0T2RNA' ,'F0T1RNA' ,'F0T2N0S1' ,'F6T0RNA' ,'F5T0N0S1_fnd' ,'F6T0RNA_fnd' ,'F6T0N0S1' ,'F6T1RNA' ,'F6T1RNA_fnd' ] 


#Models=['F2T0RNA'  , 'F2T0RNA_fnd', 'F2T0N0S1' , 'F2T1RNA'  , 'F2T1N0S1' , 'F0T2RNA'  , 'F0T2N0S1' , 'F6T0RNA' , 'F6T0RNA_fnd', 'F6T0N0S1' , 'F6T1RNA'  , 'F6T1N0S1' , 'F6T2N0S1' ] 
# Models=['F2T0RNA'  , 'F2T0RNA_fnd']
# Models=['F6T1RNA']

# for model_name in Models:

bSmallAngle=False

opts=dict()
opts['rot_elastic_type']='SmallRot' #<<< Very important, 'SmallRot', or 'Body', will affect the rotation matrix
# opts['rot_elastic_type']='Body' #<<< Very important, 'SmallRot', or 'Body', will affect the rotation matrix
opts['rot_elastic_smallAngle']=False #<<< Very important, will perform small angle approx: sin(nu q) = nu q and nu^2=0 !!! Will remove all nu^2 and nu^3 terms!! Not recommended, removes part of RNA "Y" inertia from mass matrix
opts['orderMM']      = 1     #< order of taylor expansion for Mass Matrix
opts['orderH']       = 1     #< order of taylor expansion for H term
opts['fnd_loads']    = False
opts['aero_torques'] = False
opts['mergeFndTwr']  =  model_name.find('_fnd')<=0
opts['yaw']          = 'zero'  # 'fixed', 'zero', or 'dynamic' if a DOF
opts['tilt']         = 'fixed' # 'fixed', 'zero', or 'dynamic' if a DOF
opts['tiltShaft']    = True    # OpenFAST tilts shaft not nacelle
opts['linRot']       = False    # Linearize rotations matrices from the beginning
opts['linRot']       = True    # Linearize rotations matrices from the beginning
opts['Mform']        = 'symbolic'  # or 'TaylorExpanded'
opts['twrDOFDir']    = ['x','y','x','y']  # Order in which the flexible DOF of the tower are set

# --- Esthetics Replacements for python
replaceDict={}
replaceDict['theta_tilt']=('tilt',None)

# --- Create model, solve equations and preform small angle approximation
if model_name[0]=='B':
    model = get_model_one_body(model_name, **opts)
else:
    model = get_model(model_name, **opts)


model.kaneEquations(Mform='TaylorExpanded')
# ---
extraSubs=model.shapeNormSubs
# extraSubs+=[(Symbol('J_xx_T'),0)]
# extraSubs+=[(Symbol('J_yy_T'),0)]
# extraSubs+=[(Symbol('J_zz_T'),0)]
# extraSubs+=[(Symbol('J_xx_N'),0)]
# extraSubs+=[(Symbol('J_yy_N'),0)]
# extraSubs+=[(Symbol('J_zz_N'),0)]
# extraSubs+=[(Symbol('M_T'),0)]
# extraSubs+=[(Symbol('M_N'),0)]
if model_name[0]=='B':
    extraSubs+=[(Symbol('x_BG'),0)]
    extraSubs+=[(Symbol('y_BG'),0)]

print('Extra Subs:  ', extraSubs)

# --- Linearization of non linear equations
model.linearize(noAcc=True, noVel=False, extraSubs=extraSubs)

# --- Small angle approximation and linearization
if bSmallAngle:
    print('Small angles:', model.smallAngles)
    model.smallAngleApprox(model.smallAngles, extraSubs)
    model.smallAngleApproxEOM(model.smallAngles, extraSubs)
    model.smallAngleLinearize(noAcc=True, noVel=False, extraSubs=extraSubs)
    model.smallAngleSaveTex(folder='_tex', variables=['MM','FF','M','C','K','B'])
    model.savePython(folder='_py' , variables=['MM','FF','MMsa','FFsa','M','C','K','B','Msa','Csa','Ksa','Bsa'], replaceDict=replaceDict, extraSubs=extraSubs)
else:
    model.savePython(folder='_py' , variables=['MM','FF','M','C','K','B'], replaceDict=replaceDict, extraSubs=extraSubs)

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