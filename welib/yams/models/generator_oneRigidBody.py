import os
import numpy as np
# Local
from welib.yams.models.OneRigidBody_sympy import *  # get_model_one_body, x, phi_x, xd, omega_x, x_BG, Jxx_B, M_B, etc

def generateOneRigidBodyModel(modelName, packageDir='py', texDir='tex', fullPage=True
        ):
    """
    Generate python package for typical models
    """
    smallAngles = None
    replaceDict = None

    H0Models = ['B001000_hydro0','B000010_hydro0','B101010_hydro0']# Force at "0" point
    HOModels = ['B001000_hydroO','B000010_hydroO','B101010_hydroO']# Force at body origin
    if modelName in  H0Models+HOModels:
        # --- HydroDynamicLoad
        model = get_model_one_body(modelName, linRot=False, orderMM=1, orderH=1, 
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
        # Hydro force
        P_B = body.origin.locatenew('P_B', z_BB * body.frame.z) # Buoyancy     <<<< Measured from T
        P_M = body.origin.locatenew('P_M', z_BM * body.frame.z) # Mooring      <<<< Measured from T
        P_0 = body.origin.locatenew('P_0', z_B0 * body.frame.z) # 0- sea level <<<< Measured from T
        P_O = body.origin                                       # Body origin
        P_M = body.origin
        F_hx, F_hy, F_hz = dynamicsymbols('F_hx, F_hy, F_hz') # Hydrodynamic force, function to time 
        M_hx, M_hy, M_hz = dynamicsymbols('M_hx, M_hy, M_hz') # Hydrodynamic moment, function to time 
        if modelName.find('hydro0'):
            model.addForce(body,  P_0,        F_hx * ref.frame.x + F_hy * ref.frame.y + F_hz * ref.frame.z)
        else:
            model.addForce(body,  P_O,        F_hx * ref.frame.x + F_hy * ref.frame.y + F_hz * ref.frame.z)
        model.addMoment(body, body.frame, M_hx * ref.frame.x + M_hy * ref.frame.y + M_hz * ref.frame.z)

        extraSubs    =[(omega_x,0),(omega_y,0),(omega_z,0)]
    else:
        raise NotImplementedError(modelName)
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
    pkgName = os.path.join(packageDir,modelName)
    texName = os.path.join(texDir,modelName)
    model.exportPackage(path=pkgName, extraSubs=extraSubs, smallAngles=smallAngles, replaceDict=replaceDict, pathtex=texName, fullPage=fullPage)
