import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# Local 
from welib.yams.models.FTNSB_sympy import *  # get_model_one_body, x, phi_x, xd, omega_x, x_BG, Jxx_B, M_B, etc



def generateModel(modelName, packageDir='py', texDir='tex', fullPage=True,
        extraSubs=None,
        replaceDict=None,
        smallAngles=None,
        aero_forces=True,
        moor_loads=True,
        hydro_loads=True,
        silentTimer=True,
        IMU = True,
        ):
    """
    Generate python package for typical models
    """
    print('----------------------- GENERATE MODEL {} ---------------------'.format(modelName))

    if replaceDict is None:
        replaceDict = {}

    if modelName[0]=='F':
        opts=dict()
        opts['rot_elastic_type']='SmallRot' #<<< Very important, 'SmallRot', or 'Body', will affect the rotation matrix
        # opts['rot_elastic_type']='Body' #<<< Very important, 'SmallRot', or 'Body', will affect the rotation matrix
        opts['rot_elastic_smallAngle']=False #<<< Very important, will perform small angle approx: sin(nu q) = nu q and nu^2=0 !!! Will remove all nu^2 and nu^3 terms!! Not recommended, removes part of RNA "Y" inertia from mass matrix
        opts['orderMM']      = 1     #< order of taylor expansion for Mass Matrix
        opts['orderH']       = 1     #< order of taylor expansion for H term
        opts['aero_torques'] = False
        opts['aero_forces']  = False
        opts['moor_loads']   = False
        opts['hydro_loads']  = False
        if modelName.find('noLoads')>1:
            opts['aero_forces'] = aero_forces
            pass
        else:
            opts['aero_forces']  = aero_forces
            opts['moor_loads']   = moor_loads
            opts['hydro_loads']  = hydro_loads
            if modelName.find('noHydro')>1:
                opts['hydro_loads']  = False
            elif modelName.find('hydro')>1:
                opts['hydro_loads']  = True
            if modelName.find('noMoor')>1:
                opts['moor_loads']  = False
        opts['mergeFndTwr']  =  modelName.find('_fnd')<=0
        opts['yaw']          = 'zero'  # 'fixed', 'zero', or 'dynamic' if a DOF
        opts['tilt']         = 'fixed' # 'fixed', 'zero', or 'dynamic' if a DOF
        opts['tiltShaft']    = True    # OpenFAST tilts shaft not nacelle
        #opts['linRot']       = False    # Linearize rotations matrices from the beginning
        opts['linRot']       = True    # Linearize rotations matrices from the beginning
        opts['Mform']        = 'symbolic'  # or 'TaylorExpanded'
        opts['twrDOFDir']    = ['x','y','x','y']  # Order in which the flexible DOF of the tower are set

        model = get_model(modelName, **opts)

        extraSubs=model.shapeNormSubs

        replaceDict['theta_tilt']=('tilt',None)
    else:
        raise NotImplementedError(modelName)

    if IMU:
        x_TI, z_TI = symbols('x_TI, z_TI')
        twr = model.twr
        nac = model.nac
        P_IMU = nac.origin.locatenew('P_IMU', x_TI * nac.frame.x + z_TI * nac.frame.z) 

        # IMU of OpenFAST is in tilted shaft frame
        tilt=symbols('tilt')
        shaft_frame = ReferenceFrame('e_s')
        shaft_frame.orient(model.nac.frame, 'Body', (0,tilt,0), 'XYZ') 

        model.addPoint(P_IMU, shaft_frame)
        model.computePointsMotion(noPosInJac=False)
        print('>>> Adding IMU point')

    # Export
    pkgName = os.path.join(packageDir,modelName)
    texName = os.path.join(texDir,modelName)
    model.exportPackage(path=pkgName, extraSubs=extraSubs, smallAngles=smallAngles, replaceDict=replaceDict, pathtex=texName, fullPage=fullPage, silentTimer=silentTimer)

if __name__ == '__main__':
    pass
