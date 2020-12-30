#from welib.yams.sympy_tools import *
from welib.yams.yams_sympy       import YAMSRigidBody, YAMSInertialBody, YAMSFlexibleBody
from welib.yams.yams_sympy_model import YAMSModel

from welib.tools.tictoc import Timer

from sympy import Matrix, symbols, simplify, Function, expand_trig, Symbol, diff
from sympy import cos,sin, transpose
from sympy import latex, python
from sympy.physics.mechanics import dynamicsymbols, ReferenceFrame, Point, inertia

_defaultOpts={
    'Mform'  : 'TaylorExpanded', # 'symbolic',or 'TaylorExpanded'
    'linRot' : False,              #<<< Very important if True will assume that omegas are time derivatives of DOFs
    'loads':False, # Add loads on the foundation (restoring and buoyancy)
    'orderMM':2, #< order of taylor expansion for Mass Matrix
    'orderH':2,  #< order of taylor expansion for H term
}


time = symbols('t')
phi_x, phi_y, phi_z       = dynamicsymbols('phi_x, phi_y, phi_z')
x, y, z                   = dynamicsymbols('x, y, z')
xd, yd, zd                = dynamicsymbols('xd, yd, zd')                # dynamicsymbols('x, y, z',1)
omega_x, omega_y, omega_z = dynamicsymbols('omega_x, omega_y, omega_z')
M_B                       = symbols('M_B')                              # Masses: Foundation/Tower/Nacelle/Rotor
Jxx_B, Jyy_B, Jzz_B       = symbols('Jxx_B, Jyy_B, Jzz_B')              # NOTE: JO                                                = Jyy = Jzz for a three bladed rotor!
x_BG, y_BG, z_BG          = symbols('x_BG, y_BG, z_BG')                # Position of Foundation COG in F, measured from point T
x_NR, z_NR                = symbols('x_NR, z_NR')                       # 


def get_model_one_body(model_name, **opts):
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
            print('Key {} not supported for model options.'.format(k))
    #print(opts)

    # Extract info from model name
    s= model_name[1:]
    if len(s)==1:
        nDOF_body   = int(s[0])
        bDOFs   = [False]*6
        bDOFs[0]=nDOF_body>=1 # x
        bDOFs[4]=nDOF_body>=2 # phiy
        bDOFs[2]=nDOF_body==3 or nDOF_body==6 # z
        bDOFs[1]=nDOF_body>=5 # y
        bDOFs[3]=nDOF_body>=5 # phi_x
        bDOFs[5]=nDOF_body>=5 # phi_z
    else:
        bDOFs=[s=='1' for s in s]
        nDOF_body  = sum(bDOFs)

    print('body',','.join(['1' if b else '0' for b in bDOFs]))

    # --- Isolated bodies 
    ref = YAMSInertialBody('E') 
    body = YAMSRigidBody('B', rho_G = [x_BG,y_BG,z_BG], J_diag=True) 

    # --- Body DOFs
    bodyDOFsAll    = [x, y, z, phi_x,     phi_y,       phi_z]
    bodySpeedsAll  = [xd,yd,zd,omega_x,omega_y,omega_z]
    bodyDOFs    = [dof for active,dof in zip(bDOFs,bodyDOFsAll)   if active]
    bodySpeeds  = [dof for active,dof in zip(bDOFs,bodySpeedsAll) if active]

    coordinates = bodyDOFs   

    # --- Connections between bodies
    rel_pos=[0,0,0]
    rel_pos[0] = x      if bDOFs[0] else 0
    rel_pos[1] = y      if bDOFs[1] else 0
    rel_pos[2] = z      if bDOFs[2] else 0
    rPhix= phi_x if bDOFs[3] else 0
    rPhiy= phi_y if bDOFs[4] else 0
    rPhiz= phi_z if bDOFs[5] else 0
    print('Free connection ref', rel_pos, (rPhix,rPhiy,rPhiz))
    print('>>>>>>>>>>>> ORDER IMPORTANT')
    ref.connectTo(body, type='Free' , rel_pos=rel_pos, rot_amounts=(rPhiy,rPhix,rPhiz), rot_order='YXZ')  

    # --- Defining Body rotational velocities
    omega_BE = body.ang_vel_in(ref)        # Angular velocity of body in inertial frame

    # --- Kinetics
    body_loads       = []
    bodies           = []
    g                = symbols('g')
    g_vect           = -g * ref.frame.z
    bodies+= [body]
    grav_F = (body.masscenter, -body.mass * g * ref.frame.z)
    # Point of application for Buoyancy and mooring
    #P_B = twr.origin.locatenew('P_B', z_TB * body.frame.z) # <<<< Measured from T
    #P_B.v2pt_theory(twr.origin, ref.frame, twr.frame); # PB & T are fixed in e_T
    #K_Mx, K_My, K_Mz          = symbols('K_x_M, K_y_M, K_z_M') # Mooring restoring
    #K_Mphix, K_Mphiy, K_Mphiz = symbols('K_phi_x_M, K_phi_y_M, K_phi_z_M') # Mooring restoring
    #F_B = dynamicsymbols('F_B') # Buoyancy force
    #body_loads  += [(body, (P_B, F_B * ref.frame.z))]
    #fr=0
    #fr += -K_Mx * x *ref.frame.x if bDOFs[0] else 0
    #fr += -K_My * y *ref.frame.y if bDOFs[1] else 0
    #fr += -K_Mz * z *ref.frame.z if bDOFs[2] else 0
    #Mr += -K_MPhix * phi_x *ref.frame.x if bDOFs[3] else 0
    #Mr += -K_MPhiy * phi_y *ref.frame.y if bDOFs[4] else 0
    #Mr += -K_MPhiz * phi_z *ref.frame.z if bDOFs[5] else 0
    #body_loads  += [(body, (P_M,  fr))]
    #body_loads  += [(body, (body.frame, Mr))]
    body_loads  += [(body,grav_F)]

    # --- Kinematic equations 
    kdeqsSubs =[]
    # Kdeqs for body: 
    bodyVelAll = [diff(x,time), diff(y,time),  diff(z,time)]
    if opts['linRot']:
        bodyVelAll += [diff(phi_x,time), diff(phi_y,time),  diff(phi_z,time)]  
    else:
        omega_TE = body.ang_vel_in(ref)        # Angular velocity of nacelle in inertial frame
        bodyVelAll +=[ omega_TE.dot(ref.frame.x).simplify(), omega_TE.dot(ref.frame.y).simplify(), omega_TE.dot(ref.frame.z).simplify()]  
        # >>> TODO sort out frame
    kdeqsSubs+=[ (bodySpeedsAll[i], bodyVelAll[i]) for i,dof in enumerate(bDOFs) if dof] 

    coordinates = bodyDOFs   
    speeds      = bodySpeeds 


    # --- Create a YAMS wrapper model
    model = YAMSModel(name=model_name)
    model.opts        = opts
    model.ref         = ref
    model.bodies      = bodies
    model.body_loads  = body_loads
    model.coordinates = coordinates
    model.speeds      = speeds
    model.kdeqsSubs   = kdeqsSubs
    print(model)
    print(body)

    model.body=body
    model.g_vect=g_vect

    # Small angles
    model.smallAngles    = [phi_x,phi_y,phi_z]
    model.shapeNormSubs= []

    return model
