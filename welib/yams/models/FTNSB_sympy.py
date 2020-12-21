#from welib.yams.sympy_tools import *
from welib.yams.yams_sympy       import YAMSRigidBody, YAMSInertialBody, YAMSFlexibleBody
from welib.yams.yams_sympy_model import YAMSModel
#from welib.yams.yams_sympy import DCMtoOmega

from welib.tools.tictoc import Timer
from .FTNSB_sympy_symbols import *

from sympy import Matrix, symbols, simplify, Function, expand_trig, Symbol, diff
from sympy import cos,sin, transpose
from sympy import latex, python
from sympy.physics.mechanics import dynamicsymbols, ReferenceFrame, Point, inertia

_defaultOpts={
    'floating':True,
    'yaw'    : 'fixed',  # 'fixed', 'dynamic' or 'zero'
    'tilt'   : 'fixed',  # 'fixed', 'dynamic' or 'zero'
    'Mform'  : 'TaylorExpanded', # 'symbolic',or 'TaylorExpanded'
    'mergeFndTwr':True, # Use one body for FND and TWR
    'tiltShaft':False, # Tilt shaft or nacelle
    'twrDOFDir':['x','y','x','y'], # Order in which the flexible DOF of the tower are set
    'linRot' : False,              #<<< Very important if True will assume that omegas are time derivatives of DOFs
    'rot_elastic_type':'Body',     #<<< Very important, SmallRot, or Body, will affect the rotation matrix
    'rot_elastic_subs':True,       #<<< Very important, will substitute alpha_y with nuy q. Recommended True
    'rot_elastic_smallAngle':False,#<<< Very important, will perform small angle approx: sin(nu q) = nu q and nu^2=0 !!! Will remove all nu^2 and nu^3 terms!! might not be recommended
    'orderMM':2, #< order of taylor expansion for Mass Matrix
    'orderH':2,  #< order of taylor expansion for H term
}


def get_model(model_name, **opts):
    """ 

    model: string defining model

    opts: dictionary of options with keys:
        see _defaultOpts
    
    """

    for k,v in _defaultOpts.items():
        if k not in opts.keys():
            opts[k]=v
    for k,v in opts.items():
        if k not in _defaultOpts.keys():
            raise Exception('Key {} not supported for model options.'.format(k))
    #print(opts)

    # Extract info from model name
    nDOF_fnd = int(model_name.split('F')[1][0])
    nDOF_twr = int(model_name.split('T')[1][0])
    bFullRNA = model_name.find('RNA')==-1
    bBlades=False # Rotor
    nDOF_nac=0
    nDOF_sft=0
    nDOF_bld=0
    if bFullRNA:
        nDOF_nac = int(model_name.split('N')[1][0])
        nDOF_sft = int(model_name.split('S')[1][0])

        bBlades = model_name.find('B')>0
        if bBlades:
            nDOF_bld = model_name.split('B')[1][0]
        else:
            nDOF_bld=0
    print('fnd',nDOF_fnd, 'twr',nDOF_twr, 'nac',nDOF_nac, 'sft',nDOF_sft, 'bld',nDOF_bld) 

    # --------------------------------------------------------------------------------}
    # --- Isolated bodies 
    # --------------------------------------------------------------------------------{
    # Reference frame
    ref = YAMSInertialBody('E') 
    # Foundation, floater, always rigid for now
    if (not opts['floating']) or opts['mergeFndTwr']:
        fnd = None # the floater is merged with the twr, or we are not floating
    else:
        fnd = YAMSRigidBody('F', rho_G = [0,0,z_FG], J_diag=True) 
    # Tower
    if nDOF_twr==0:
        # Ridid tower
        twr = YAMSRigidBody('T', rho_G = [0,0,z_TG], J_diag=True) 
        twrDOFs=[]
        twrSpeeds=[]
    elif nDOF_twr<=4:
        # Flexible tower
        twr = YAMSFlexibleBody('T', nDOF_twr, directions=opts['twrDOFDir'], orderMM=opts['orderMM'], orderH=opts['orderH'], predefined_kind='twr-z')
        twrDOFs   = twr.q
        twrSpeeds = twr.qd
    # Nacelle rotor assembly
    if bFullRNA:
        # Nacelle
        nac = YAMSRigidBody('N', rho_G = [x_NG ,0, z_NG], J_diag=True) 
        # Shaft
        #...
        # Individual blades or rotor
        if bBlades:
            # Individual blades
            raise NotImplementedError()
        else:
            # Rotor
            rot = YAMSRigidBody('R', rho_G = [0,0,0], J_diag=True)
            rot.inertia = (inertia(rot.frame, Jxx_R, JO_R, JO_R), rot.origin)  # defining inertia at orign
    else:
        # Nacelle
        nac = YAMSRigidBody('RNA', rho_G = [x_RNAG ,0, z_RNAG], J_diag=True) 
        rot = None

    # --------------------------------------------------------------------------------}
    # --- Body DOFs
    # --------------------------------------------------------------------------------{
    # Fnd
    if (not opts['floating']):
        fndDOFs   = []
        fndSpeeds = []
    elif nDOF_fnd==0 :
        fndDOFs   = []
        fndSpeeds = []
    elif nDOF_fnd==1:
        fndDOFs   = [x]
        fndSpeeds = [xd]
    elif nDOF_fnd==2:
        fndDOFs   = [x,  phi_y   ]
        fndSpeeds = [xd,omega_y_T]
    elif nDOF_fnd==3:
        fndDOFs   = [x,  phi_y   , phi_x]
        fndSpeeds = [xd,omega_y_T, omega_x_T]
    else:
        fndDOFs    = [x, y, z, phi_x,     phi_y,       phi_z]
        fndSpeeds  = [xd,yd,zd,omega_x_T,omega_y_T,omega_z_T]
    # Twr

    # Nac
    if nDOF_nac==2:
        opts['yaw']='dynamic'
        opts['tilt']='dynamic'
    if opts['tiltShaft'] and opts['tilt']=='dynamic':
        raise Exception('Cannot do tiltshaft with tilt dynamic')

    yawDOF  = {'zero':0, 'fixed':theta_yaw,  'dynamic':q_yaw }[opts['yaw']]
    tiltDOF = {'zero':0, 'fixed':theta_tilt, 'dynamic':q_tilt}[opts['tilt']]
    nacDOFs     = []
    nacSpeeds   = []
    nacKDEqSubs = []
    if opts['yaw']=='dynamic':
        nacDOFs     += [q_yaw]
        nacSpeeds   += [qd_yaw]
        nacKDEqSubs += [(qd_yaw, diff(q_yaw, time))]
    if opts['tilt']=='dynamic':
        nacDOFs     += [q_tilt]
        nacSpeeds   += [qd_tilt]
        nacKDEqSubs += [(qd_tilt, diff(q_tilt, time))]

    nacDOFsAct=(opts['yaw']=='dynamic',opts['tilt']=='dynamic')
    if nDOF_nac==0:
        if not (nacDOFsAct==(False,False)):
            raise Exception('If nDOF_nac is 0, yaw and tilt needs to be "fixed" or "zero"')
    elif nDOF_nac==1:
        if not (nacDOFsAct==(True,False) or nacDOFsAct==(False,True) ):
            raise Exception('If nDOF_nac is 1, yaw or tilt needs to be "dynamic"')
    else:
        if not (nacDOFsAct==(True,True)):
            raise Exception('If nDOF_nac is 2, yaw and tilt needs to be "dynamic"')

    # Shaft
    sftDOFs  =[]
    sftSpeeds=[]
    if bFullRNA:
        if nDOF_sft==1:
            sftDOFs   = [psi]
            sftSpeeds = [omega_x_R]
        elif nDOF_sft==0:
            pass
        else:
            raise Exception('nDOF shaft should be 0 or 1')

        # Blade Rotor
        if bBlades:
            pass
        else:
            pass

    coordinates = fndDOFs   + twrDOFs   + nacDOFs   + sftDOFs
    speeds      = fndSpeeds + twrSpeeds + nacSpeeds + sftSpeeds  # Order determine eq order



    # --------------------------------------------------------------------------------}
    # --- Connections between bodies
    # --------------------------------------------------------------------------------{
    z_OT = Symbol('z_OT')
    if opts['floating']:
        if nDOF_fnd==0:
            print('Rigid connection ref twr')
            ref.connectTo(twr, type='Rigid' , rel_pos=(0,0,z_OT))
        elif nDOF_fnd==1: 
            print('Constraint connection ref twr')
            ref.connectTo(twr, type='Free' , rel_pos=(x,0,z_OT), rot_amounts=(0    , x * symbols('nu'), 0   ), rot_order='XYZ')
        elif nDOF_fnd==2: 
            print('Free connection ref twr')
            ref.connectTo(twr, type='Free' , rel_pos=(x,0,z_OT), rot_amounts=(0    ,phi_y, 0   ), rot_order='XYZ')
        elif nDOF_fnd==3: 
            print('Free connection ref twr')
            ref.connectTo(twr, type='Free' , rel_pos=(x,0,z_OT), rot_amounts=(phi_x,phi_y, 0   ), rot_order='XYZ')
        else:
            print('Free connection ref twr')
            ref.connectTo(twr, type='Free' , rel_pos=(x,y,z+z_OT), rot_amounts=(phi_x,phi_y,phi_z), rot_order='XYZ')  #NOTE: rot order is not "optimal".. phi_x should be last
    else:
        print('Rigid connection ref twr')
        ref.connectTo(twr, type='Rigid' , rel_pos=(0,0,0))

    # Rigid connection between twr and fnd if fnd exists
    if fnd is not None:
        print('Rigid connection twr fnd')
        if nDOF_twr==0:
            twr.connectTo(fnd, type='Rigid', rel_pos=(0,0,0)) # -L_F
        else:
            raise Exception('Flexible bodies can only have on body connection')

    if nDOF_twr==0:
        # Tower rigid -> Rigid connection to nacelle
        # TODO TODO L_T or twr.L
        if nDOF_nac==0:
            print('Rigid connection twr nac')
        else:
            print('Dynamic connection twr nac')

        if opts['tiltShaft']:
            # Shaft will be tilted, not nacelle
            twr.connectTo(nac, type='Rigid', rel_pos=(0,0,L_T)  , rot_amounts=(yawDOF,0,0), rot_order='ZYX')
        else:
            # Nacelle is tilted
            twr.connectTo(nac, type='Rigid', rel_pos=(0,0,L_T)  , rot_amounts=(yawDOF,tiltDOF,0), rot_order='ZYX')

    else:
        # Flexible tower -> Flexible connection to nacelle
        print('Flexible connection twr nac')
        if opts['tiltShaft']:
            twr.connectTo(nac, type='Joint', rel_pos=(0,0,twr.L)  , rot_amounts=(yawDOF, 0      , 0), rot_order='ZYX', rot_type_elastic=opts['rot_elastic_type'], doSubs=opts['rot_elastic_subs'])
        else:
            twr.connectTo(nac, type='Joint', rel_pos=(0,0,twr.L)  , rot_amounts=(yawDOF, tiltDOF, 0), rot_order='ZYX', rot_type_elastic=opts['rot_elastic_type'], doSubs=opts['rot_elastic_subs'])

    if bFullRNA:
        if bBlades:
            raise NotImplementedError()
        else:
            if opts['tiltShaft']:
                if nDOF_sft==0:
                    nac.connectTo(rot, type='Joint', rel_pos=(x_NR,0,z_NR), rot_amounts=(0,tiltDOF,0), rot_order='ZYX')
                else:
                    nac.connectTo(rot, type='Joint', rel_pos=(x_NR,0,z_NR), rot_amounts=(0,tiltDOF,psi), rot_order='ZYX')
            else:
                if nDOF_sft==0:
                    nac.connectTo(rot, type='Joint', rel_pos=(x_NR,0,z_NR), rot_amounts=(0,0      ,0), rot_order='ZYX')
                else:
                    nac.connectTo(rot, type='Joint', rel_pos=(x_NR,0,z_NR), rot_amounts=(0,0      ,psi), rot_order='ZYX')


    # --- Defining Body rotational velocities
    omega_TE = twr.ang_vel_in(ref)        # Angular velocity of nacelle in inertial frame
    omega_NT = nac.ang_vel_in(twr.frame)  # Angular velocity of nacelle in inertial frame
    if rot is not None:
        omega_RN = rot.ang_vel_in(nac.frame)  # Angular velocity of rotor wrt Nacelle (omega_R-omega_N)

    # --- Kinetics
    body_loads       = []
    bodies           = []
    g                = symbols('g')
    T_a              = dynamicsymbols('T_a') # NOTE NOTE
    #T_a              = Function('T_a')(dynamicsymbols._t, *coordinates, *speeds) # NOTE NOTE
    M_ax, M_ay, M_az = dynamicsymbols('M_ax, M_ay, M_az') # Aero torques
    K_Mx, K_My, K_Mz = dynamicsymbols('K_Mx, K_My, K_Mz') # Mooring restoring
    K_Mphix, K_Mphiy, K_Mphiz = dynamicsymbols('K_M_phi_x, K_M_phi_y, K_M_phi_z') # Mooring restoring
    F_B              = dynamicsymbols('F_B')

    # Foundation/floater loads
    if fnd is not None:
        bodies+= [fnd]
        grav_F = (fnd.masscenter, -fnd.mass * g * ref.frame.z)
        # Point of application for Buoyancy and mooring
        P_B = twr.origin.locatenew('P_B', z_TB * fnd.frame.z) # <<<< Measured from T
        P_M = twr.origin.locatenew('P_M', z_TM * fnd.frame.z) # <<<< Measured from T
        P_B.v2pt_theory(twr.origin, ref.frame, twr.frame); # PB & T are fixed in e_T
        P_M.v2pt_theory(twr.origin, ref.frame, twr.frame); # PM & T are fixed in e_T
        # Buyancy
        #F_buy = (P_B, F_B * ref.frame.z)
        ## Restoring mooring
        #F_moor = (P_M,  -K_x * x *ref.frame.x )
        ## Ext torques
        #M_moor = (fnd.frame, -K_phi_y * phi_y *ref.frame.y)
        # TODO Moor and buoy
        body_loads  += [(fnd,grav_F)]

    # Tower loads
    grav_T       = (twr.masscenter, -twr.mass * g * ref.frame.z)
    bodies      += [twr]
    body_loads  += [(twr,grav_T)]  

    # Nacelle loads
    grav_N = (nac.masscenter, -nac.mass * g * ref.frame.z)
    bodies      += [nac]
    body_loads  += [(nac,grav_N)]  

    if bFullRNA:
        if bBlades:
            raise NotImplementedError()
        else:
            # Rotor loads
            grav_R = (rot.masscenter, -M_R * g * ref.frame.z)
            bodies      += [rot]
            body_loads  += [(rot,grav_R)]  

            # NOTE: loads on rot, but expressed in N frame
            if opts['tiltShaft']:
                # TODO actually tilt shaft, introduce non rotating shaft body
                thrustR = (rot.origin, T_a *cos(tiltDOF) * nac.frame.x -T_a *sin(tiltDOF) * nac.frame.z)
                M_a_R = (nac.frame, M_ax*nac.frame.x*0 )# TODO TODO
                print('>>> WARNING tilt shaft aero moments not implemented') # TODO
            else:
                thrustR = (rot.origin, T_a * nac.frame.x )
                #thrustR = (rot.origin, T_a * rot.frame.x )
                #M_a_R = (rot.frame, M_ax*rot.frame.x +  M_ay*rot.frame.y  + M_az*rot.frame.z) # TODO TODO TODO introduce a non rotating shaft
                M_a_R = (nac.frame, M_ax*nac.frame.x +  M_ay*nac.frame.y  + M_az*nac.frame.z) 
            body_loads  += [(rot,thrustR), (nac, M_a_R)]  

    else:
        # RNA loads, point load at R
        R=Point('R')
        R.set_pos(nac.origin, x_NR * nac.frame.x + z_NR* nac.frame.z)
        R.set_vel(nac.frame, 0 * nac.frame.x)
        R.v2pt_theory(nac.origin, ref.frame, nac.frame)
        #thrustN = (nac.masscenter, T * nac.frame.x)
        if opts['tiltShaft']:
            thrustN = (R, T_a *cos(tiltDOF) * nac.frame.x -T_a *sin(tiltDOF) * nac.frame.z)
        else:
            thrustN = (R, T_a * nac.frame.x )

        if opts['tiltShaft']:
            print('>>> WARNING tilt shft aero moments with RNA not implemented') # TODO
            M_a_N = (nac.frame, 0*nac.frame.x)
        else:
            M_a_N = (nac.frame, M_ax*nac.frame.x +  M_ay*nac.frame.y  + M_az*nac.frame.z)

        body_loads  += [(nac,thrustN), (nac, M_a_N)]  

    # --------------------------------------------------------------------------------}
    # --- Kinematic equations 
    # --------------------------------------------------------------------------------{
    kdeqsSubs =[]
    # --- Fnd
    if not opts['floating']:
        pass
    elif nDOF_fnd==0:
        pass
    elif nDOF_fnd==1:
        kdeqsSubs+=[(xd, diff(x,time))]; 
        #if opts['linRot']:
        #    kdeqsSubs+=[ (omega_y_T, diff(phi_y,time))]  
        #else:
        #    kdeqsSubs+=[ (omega_y_T, omega_TE.dot(ref.frame.y).simplify())]  
    elif nDOF_fnd==2:
        kdeqsSubs+=[(xd, diff(x,time))]; 
        if opts['linRot']:
            kdeqsSubs+=[ (omega_y_T, diff(phi_y,time))]  
        else:
            kdeqsSubs+=[ (omega_y_T, omega_TE.dot(ref.frame.y).simplify())]  
    elif nDOF_fnd==3:
        kdeqsSubs+=[(xd, diff(x,time))]; 
        if opts['linRot']:
            kdeqsSubs+=[ (omega_x_T, diff(phi_x,time))]  
            kdeqsSubs+=[ (omega_y_T, diff(phi_y,time))]  
        else:
            kdeqsSubs+=[ (omega_x_T, omega_TE.dot(ref.frame.x).simplify())]  
            kdeqsSubs+=[ (omega_y_T, omega_TE.dot(ref.frame.y).simplify())]  
    elif nDOF_fnd==6:
        kdeqsSubs+=[(xd, diff(x,time))]; 
        kdeqsSubs+=[(yd, diff(y,time))]; 
        kdeqsSubs+=[(zd, diff(z,time))]; 
        if opts['linRot']:
            kdeqsSubs+=[ (omega_x_T, diff(phi_x,time))]  
            kdeqsSubs+=[ (omega_y_T, diff(phi_y,time))]  
            kdeqsSubs+=[ (omega_z_T, diff(phi_z,time))]  
        else:
            kdeqsSubs+=[ (omega_x_T, omega_TE.dot(ref.frame.x).simplify())]  
            kdeqsSubs+=[ (omega_y_T, omega_TE.dot(ref.frame.y).simplify())]  
            kdeqsSubs+=[ (omega_z_T, omega_TE.dot(ref.frame.z).simplify())]  
    else:
        raise Exception('nDOF fnd needs to be 0, 2 and 6 for now')

    # --- Twr
    if nDOF_twr==0:
        pass
    else:
        kdeqsSubs +=[ (twr.qd[i], twr.qdot[i]) for i,_ in enumerate(twr.q)]; 

    # --- Nac
    kdeqsSubs+=nacKDEqSubs

    # --- Shaft
    if bFullRNA:
        if bBlades:
            raise NotImplementedError()
        else:
            if nDOF_sft==1:
                kdeqsSubs+=[ (omega_x_R, omega_RN.dot(rot.frame.x).simplify()) ]  





    # --- Create a YAMS wrapper model
    model = YAMSModel(name=model_name)
    model.opts        = opts
    model.ref         = ref
    model.bodies      = bodies
    model.body_loads  = body_loads
    model.coordinates = coordinates
    model.speeds      = speeds
    model.kdeqsSubs   = kdeqsSubs

    model.fnd=fnd
    model.twr=twr
    model.nac=nac
    model.rot=rot
    #model.sft=sft
    #model.bld=bld

    # Small angles
    model.smallAnglesFnd    = [phi_x,phi_y,phi_z]
    if nDOF_twr>0:
        if opts['rot_elastic_smallAngle']:
            model.smallAnglesTwr    = twr.vcList
        else:
            model.smallAnglesTwr    = []
    else:
        model.smallAnglesTwr    = []

    model.smallAnglesNac = []
    if opts['yaw']=='dynamic':
        model.smallAnglesNac += [q_yaw]
    if opts['tilt']=='dynamic':
        model.smallAnglesNac += [q_tilt]
    model.smallAngles=model.smallAnglesFnd + model.smallAnglesTwr + model.smallAnglesNac

    # Shape normalization
    if nDOF_twr>0:
        model.shapeNormSubs= [(v,1) for v in twr.ucList]
    else:
        model.shapeNormSubs= []



    return model

#     elif model_name=='F2T0RNA':
#         # --- Connections between bodies
#         ref.connectTo(twr, type='Free' , rel_pos=(x,0,0)   , rot_amounts=(0,phi_y,0), rot_order='XYZ')
# 
#         if Floating:
#             if nPltfmDOF==2:
#                 kdeqsSubs+=[(xd, diff(x,time))]; 
#                 kdeqsSubs+=[ (omega_y_T, omega_TE.dot(ref.frame.y).simplify())]  
#                 coordinates = [x, phi_y] + twr.q
#                 speeds      = [xd, omega_y_T] + twr.qd
# 
# 
#     elif model_name=='F0T0RNA':
#         ref.connectTo(twr, type='Rigid', rel_pos=(0,0,0)   , rot_amounts=(0,0,0), rot_order='XYZ')
# 
#         coordinates = [] + twr.q
#         speeds      = [] + twr.qd
# 
#     else:
#         raise Exception()
