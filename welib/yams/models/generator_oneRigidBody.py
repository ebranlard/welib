import os
import numpy as np
# Local
from welib.yams.models.utils import stiffness6DOF
from welib.yams.models.OneRigidBody_sympy import *  # get_model_one_body, x, phi_x, xd, omega_x, x_BG, Jxx_B, M_B, etc

def generateOneRigidBodyModel(modelName, packageDir='py', texDir='tex', fullPage=True,
        CG_on_z=False, # <<<, 
        KMoorKeep=[(0,0),(1,1),(2,2),(3,3),(4,4),(5,5),(0,4),(1,5)], 
        silentTimer=True,
        IMU=True
        ):
    """
    Generate python package for typical models
    """
    smallAngles = None
    replaceDict = None
    extraSubs = None


#     FModelsRNA= ['F2T0RNA', 'F2T0RNA_fnd', 'F2T1RNA', 'F2T1RNA_fnd', 'F3T1RNA_fnd', 'F5T1RNA_fnd', 'F2T1N0S1_fnd', 'F0T2RNA', 'F0T1RNA', 'F6T0RNA']
    # model_name = 'F2T0N0S1_fnd'
    # model_name = 'F0T2N0S1'
    #model_name = 'F5T0N0S1_fnd'
    # 
    # model_name = 'F6T0N0S1'
    # model_name = 'F6T1N0S1'
    # model_name = 'F6T1N0S1_fnd'
    # model_name = 'F6T2N0S1'

    # model_name='F000101T0N0S1_fnd'
    # model_name='F000111T0N0S1_fnd'
    # model_name = 'F6T0RNA_fnd'
    # model_name = 'F6T1RNA'
    # model_name = 'F6T1RNA_fnd'
    # model_name='F000111T0RNA_fnd'
    # model_name='F000101T0RNA_fnd'



    BModels = ['B100010','B001000','B000010','B101010','B110110','B110111']

    # Mooring models
#     MMModels = ['B100000_moorM','B100010_moorM','B001000_moorM','B000010_moorM','B101010_moorM','B110111_moorM']# Mooring force at "M" point
#     MOModels = ['B100000_moorO','B100010_moorO','B001000_moorO','B000010_moorO','B101010_moorO','B110111_moorO']# Mooring force at "M" point
# 
#     H0Models = ['B100010_hydro0','B001000_hydro0','B000010_hydro0','B101010_hydro0', 'B111111_hydro0']# Force at "0" point
#     HOModels = ['B100010_hydroO','B001000_hydroO','B000010_hydroO','B101010_hydroO', 'B111111_hydroO']# Force at body origin
# 
#     HOMModesl

    if modelName in BModels:

        #opts=dict()
        #opts['rot_elastic_type']='SmallRot' #<<< Very important, 'SmallRot', or 'Body', will affect the rotation matrix
        ## opts['rot_elastic_type']='Body' #<<< Very important, 'SmallRot', or 'Body', will affect the rotation matrix
        #opts['rot_elastic_smallAngle']=False #<<< Very important, will perform small angle approx: sin(nu q) = nu q and nu^2=0 !!! Will remove all nu^2 and nu^3 terms!! Not recommended, removes part of RNA "Y" inertia from mass matrix
        #opts['orderMM']      = 1     #< order of taylor expansion for Mass Matrix
        #opts['orderH']       = 1     #< order of taylor expansion for H term
        #opts['fnd_loads']    = False
        #opts['aero_torques'] = False
        #opts['mergeFndTwr']  =  model_name.find('_fnd')<=0
        #opts['yaw']          = 'zero'  # 'fixed', 'zero', or 'dynamic' if a DOF
        #opts['tilt']         = 'fixed' # 'fixed', 'zero', or 'dynamic' if a DOF
        #opts['tiltShaft']    = True    # OpenFAST tilts shaft not nacelle
        ##opts['linRot']       = False    # Linearize rotations matrices from the beginning
        #opts['linRot']       = True    # Linearize rotations matrices from the beginning
        #opts['Mform']        = 'symbolic'  # or 'TaylorExpanded'
        #opts['twrDOFDir']    = ['x','y','x','y']  # Order in which the flexible DOF of the tower are set

        model = get_model_one_body(modelName, linRot=False, orderMM=1, orderH=1, 
                                   J_cross=False, CG_on_z=CG_on_z, J_at_Origin=True)

    elif 'moor' in modelName or 'hydro' in modelName:
#     if modelName in  H0Models+HOModels:
        # --- Hydro/Mooring Dynamic Load
        model = get_model_one_body(modelName, linRot=False, orderMM=1, orderH=1, 
                                   J_cross=False, CG_on_z=CG_on_z, J_at_Origin=True)
        body = model.body
        ref  = model.ref
        # Points
        z_B0= symbols('z_B0')
        z_BM= symbols('z_BM')
        z_B0= symbols('z_B0')
        P_0 = body.origin.locatenew('P_0', z_B0 * body.frame.z) # 0- sea level <<<< Measured from T
        #P_0 = body.origin.locatenew('P_0', z_B0 * ref.frame.z) # 0- sea level <<<< Measured from T Does not work
        #P_0 = body.origin.locatenew('P_0', z_B0*cos(phi_y) * body.frame.z - z_B0*sin(phi_y) * body.frame.x) # 0- sea level <<<< Measured from T Does not Work
        P_O = body.origin                                       # Body origin
        # --- Hydro force
        if '_hydro' in modelName:
            F_hx, F_hy, F_hz = dynamicsymbols('F_hx, F_hy, F_hz') # Hydrodynamic force, function to time 
            M_hx, M_hy, M_hz = dynamicsymbols('M_hx, M_hy, M_hz') # Hydrodynamic moment, function to time 
            if modelName.find('hydro0')>1:
                model.addForce(body,  P_0,        F_hx * ref.frame.x + F_hy * ref.frame.y + F_hz * ref.frame.z)
            elif modelName.find('hydroO')>1:
                model.addForce(body,  P_O,        F_hx * ref.frame.x + F_hy * ref.frame.y + F_hz * ref.frame.z)
            else:
                raise NotImplementedError()
            model.addMoment(body, body.frame, M_hx * ref.frame.x + M_hy * ref.frame.y + M_hz * ref.frame.z)

        # --- Mooring force
        if '_moor' in modelName:
            if modelName.find('moorM')>1:
                raise Exception(modelName)
                P_M = body.origin.locatenew('P_M', z_BM * body.frame.z) # Mooring      <<<< Measured from T
            elif modelName.find('moor0')>1:
                P_M=P_0
            elif modelName.find('moorO')>1:
                P_M=P_O
            else:
                raise NotImplementedError()
            DOFs=[x, y, z, phi_x, phi_y, phi_z]
            fr, Mr, KM = stiffness6DOF(DOFs, ref.frame, label='KM', bDOFs=model.bDOFs, IKeep=KMoorKeep)
            model.addForce(body, P_M,  fr)
            model.addMoment(body, body.frame, Mr)
#         extraSubs    =[(omega_x,0),(omega_y,0),(omega_z,0)]  # TODO TODO TODO
#         extraSubs    =[(omega_x,phi_x.diff(dynamicsymbols._t)),(omega_y,phi_y.diff(dynamicsymbols._t)),(omega_z,phi_z.diff(dynamicsymbols._t))]  # TODO TODO TODO
    else:
        raise NotImplementedError(modelName)


    if IMU:
        x_BI, z_BI = symbols('x_BI, z_BI')
        body = model.body
        P_IMU = body.origin.locatenew('P_IMU', x_BI * body.frame.x + z_BI * body.frame.z) 

        # Creating tilt frame since IMU of OpenFAST is in tilted frame
        tilt=symbols('tilt')
        shaft_frame = ReferenceFrame('e_s')
        shaft_frame.orient(body.frame, 'Body', (0,tilt,0), 'XYZ') 

        model.addPoint(P_IMU, shaft_frame)
        model.computePointsMotion(noPosInJac=False)

        print('>>> Adding IMU point')



    # Export
    pkgName = os.path.join(packageDir,modelName)
    texName = os.path.join(texDir,modelName)
    model.exportPackage(path=pkgName, extraSubs=extraSubs, smallAngles=smallAngles, replaceDict=replaceDict, pathtex=texName, fullPage=fullPage, silentTimer=silentTimer)
