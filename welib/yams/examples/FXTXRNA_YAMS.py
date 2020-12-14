from __future__ import print_function, division
from matplotlib.pyplot import plot, legend, xlabel, ylabel, rcParams
# get_ipython().run_line_magic('matplotlib', 'inline')
import sympy
from sympy import Matrix, symbols, simplify, Function, expand_trig, Symbol
from sympy import cos,sin, transpose
from sympy import diff
from sympy import latex, python
from sympy.physics.mechanics import dynamicsymbols, ReferenceFrame, Point
from sympy.physics.mechanics import inertia, RigidBody, KanesMethod
from sympy.physics.vector import get_motion_params, cross
from sympy.physics.vector import init_vprinting, vlatex
# init_vprinting(use_latex='mathjax', pretty_print=False)
#from sympy import init_printing
#init_printing(use_latex='mathjax', pretty_print=False)
# from IPython.core.interactiveshell import InteractiveShell
# InteractiveShell.ast_node_interactivity = "all"
# get_ipython().run_line_magic('load_ext', 'autoreload')
# get_ipython().run_line_magic('autoreload', '2')

# --- My packages
from sympy_tools import *
from welib.yams.yams_sympy import YAMSRigidBody, YAMSInertialBody, YAMSFlexibleBody
from welib.yams.yams_sympy import DCMtoOmega
from welib.yams.yams_kane import *
from welib.tools.tictoc import Timer




# Main symbols
time = symbols('t')
phi_x, phi_y, phi_z = dynamicsymbols('phi_x, phi_y, phi_z')
psi       = dynamicsymbols('psi')
x,y,z     = dynamicsymbols('x, y, z')
xd,yd,zd  = dynamicsymbols('x_d, y_d, z_d') # derivative of x
omega_x_T, omega_y_T, omega_z_T = dynamicsymbols('omega_x_T, omega_y_T, omega_z_T')
omega_x_N, omega_y_N = dynamicsymbols('omega_x_N, omega_y_N')
#q1,q2 = dynamicsymbols('q_1, q_2')
# --- Parameters
M_T, M_N  = symbols('M_T,M_N') # Masses: Tower/Nacelle
Jxx_T, Jyy_T, Jzz_T = symbols('Jxx_T, Jyy_T, Jzz_T')
Jxx_N, Jyy_N, Jzz_N = symbols('Jxx_N, Jyy_N, Jzz_N')
L_T  = symbols('L_T') # length
z_G_T = symbols('z_G_T') # Position of Tower COG in T
x_G_N,z_G_N = symbols('x_G_N, z_G_N') # Position of Nacelle COG in N
x_NR,z_NR = symbols('x_NR, z_NR') # 
theta_yaw, theta_tilt = symbols('theta_yaw, theta_tilt') # Position of mooring line attachment point and Buoyancy center in F, measured from point T (negative)
theta_yaw_dyn, theta_tilt_dyn = dynamicsymbols('theta_yaw, theta_tilt') # Position of mooring line attachment point and Buoyancy center in F, measured from point T (negative)

# --- Simplifying substitution to better look at equations (no substructure, nacelle COG on tower top)
SubsSimp = [(z_G_N,0),(x_G_N,0)]
Subs=[ (Symbol('u_xT1c'),1) ]


# --- Main Parameters
#from yams.yams_sympy import *
#from yams.yams_kane import *
Floating=True
nPltfmDOF=6


# --- Define bodies 
ref = YAMSInertialBody('E') 
#twr = YAMSRigidBody('T', rho_G = [0,0,z_G_T]      , J_diag=True) 
twr = YAMSFlexibleBody('T',1, directions=['x','y','x','y']) 
nac = YAMSRigidBody('N', rho_G = [x_G_N ,0, z_G_N], J_diag=True) 

# --- Define connections
if Floating:
    if nPltfmDOF==2:
        ref.connectTo(twr, type='Free' , rel_pos=(x,0,0)   , rot_amounts=(0,phi_y,0), rot_order='XYZ')
    else:
        ref.connectTo(twr, type='Free' , rel_pos=(x,y,z)   , rot_amounts=(phi_x,phi_y,phi_z), rot_order='XYZ')
else:
    ref.connectTo(twr, type='Rigid', rel_pos=(0,0,0)   , rot_amounts=(0,0,0), rot_order='XYZ')
twr.connectTo(nac, type='Rigid', rel_pos=(0,0,twr.L)  , rot_amounts=(0,0,0),rot_order='ZYX', doSubs=True, rot_type_elastic='Body')
#twr.connectTo(nac, type='Joint', rel_pos=(0,0,L_T)  , rot_amounts=(theta_yaw,theta_tilt,0),rot_order='ZYX')
#twr.connectTo(nac, type='Joint', rel_pos=(0,0,L_T)  , rot_amounts=(theta_yaw_dyn,theta_tilt_dyn,0),rot_order='ZYX')

## --- Defining Body rotational velocities
omega_TE = twr.ang_vel_in(ref)  # Angular velocity of nacelle in inertial frame
omega_NT = nac.ang_vel_in(twr.frame)  # Angular velocity of nacelle in inertial frame

# --- Kinetics
g = symbols('g')
T = symbols('T')
grav_T = (twr.masscenter, -twr.mass * g * ref.frame.z)
grav_N = (nac.masscenter, -nac.mass * g * ref.frame.z)
R=Point('R')
R.set_pos(nac.origin, x_NR * nac.frame.x + z_NR* nac.frame.z)
R.set_vel(nac.frame, 0 * nac.frame.x)
R.v2pt_theory(nac.origin, ref.frame, nac.frame)
#thrustN = (nac.masscenter, T * nac.frame.x)
thrustN = (R, T *cos(theta_tilt) * nac.frame.x -T *sin(theta_tilt) * nac.frame.z)

# --- Kane
kdeqsSubs = [ (twr.qd[i], twr.qdot[i]) for i,_ in enumerate(twr.q)]; 
if Floating:
    if nPltfmDOF==2:
        kdeqsSubs+=[(xd, diff(x,time))]; 
        kdeqsSubs+=[ (omega_y_T, omega_TE.dot(ref.frame.y).simplify())]  
        coordinates = [x, phi_y] + twr.q
        speeds      = [xd, omega_y_T] + twr.qd
    else:
        kdeqsSubs+=[(xd, diff(x,time))]; 
        kdeqsSubs+=[(yd, diff(y,time))]; 
        kdeqsSubs+=[(zd, diff(z,time))]; 
        kdeqsSubs+=[ (omega_x_T, omega_TE.dot(ref.frame.x).simplify())]  
        kdeqsSubs+=[ (omega_y_T, omega_TE.dot(ref.frame.y).simplify())]  
        kdeqsSubs+=[ (omega_z_T, omega_TE.dot(ref.frame.z).simplify())]  
        coordinates = [x, y, z, phi_x, phi_y, phi_z] + twr.q
        speeds      = [xd,yd,zd,omega_x_T,omega_y_T,omega_z_T] + twr.qd
else:
    coordinates = [] + twr.q
    speeds      = [] + twr.qd

kdeqs = [s-v for (s,v) in kdeqsSubs]
qspeeds     = [diff(q,dynamicsymbols._t) for q in coordinates]
body_loads  = [(nac,grav_N), (nac,thrustN)];  loads = [f[1] for f in body_loads]
bodies      = [twr,nac]



print('Starting kane...')
with Timer('Kane 1'):
    kane = YAMSKanesMethod(ref.frame, coordinates, speeds, kdeqs)
with Timer('Kane 2'):
    fr, frstar  = kane.kanes_equations(bodies, loads, Mform='symbolic')
with Timer('MM'):
    MM = kane.mass_matrix
with Timer('forcing'):
    forcing=kane.forcing


print('Done')
#frstar   , MM     = kane_frstar(bodies, coordinates, speeds, kdeqs, ref.origin, ref.frame, Mform='symbolic')
#frstar_alt,MM_alt = kane_frstar_alt(bodies, coordinates, speeds, kdeqs, ref.frame)
#print('TODO MISSING Gravity of Flexible body')
#fr_alt = kane_fr_alt(loads, coordinates, speeds, kdeqs, ref.frame)
#fr     = kane_fr(body_loads, speeds, ref.frame)
#frstar
#frstar-frstar_alt
#fr.expand()
#(fr-fr_alt).expand()
#Matrix(twr.qddot)
#twr.inertial_force
#twr.inertial_torque
#twr.inertial_elast
##twr.h_omega
##twr.h_elast
#twr.nonMM
#twr.MM
#frstar
# MM
# forcing

print('Exporting...')
with Timer('Export'):
    with open('fr.tex','w') as f:
        f.write(cleantex(latex(fr)))

    with open('frstar.tex','w') as f:
        f.write(cleantex(latex(frstar)))

    with open('MM.tex','w') as f:
        f.write(cleantex(latex(MM)))

    with open('forcing.tex','w') as f:
        f.write(cleantex(latex(forcing)))


# 
# # In[30]:
# 
# 
# MM2=MM.subs([(twr.q[0],0)])
# MM2.simplify()
# MM2
# 
# 
# # In[43]:
# 
# 
# simplify(fr)
# simplify(frstar)
# 
# 
# # In[40]:
# 
# 
# MM2=MM.subs([(twr.q[0],0)])
# MM2.simplify()
# MM2
# MM_rigid
# forcing2=forcing.subs([(twr.q[0],0),(twr.qd[0],0)])
# simplify(forcing2)
# simplify(forcing_rigid)
# 
# 
# # In[33]:
# 
# # 
# # # --- Equivalent Rigid Body model
# # from yams.yams_sympy import *
# # # --- Define bodies and connections
# # ref = YAMSInertialBody('E') 
# # twr = YAMSRigidBody('T', rho_G = [0,0,z_G_T], J_diag=True) 
# # #twr = YAMSFlexibleBody('T',1, directions=['x','y','x','y']) 
# # nac = YAMSRigidBody('N', rho_G = [x_G_N ,0, z_G_N], J_diag=True) 
# # ref.connectTo(twr, type='Free' , rel_pos=(x,0,0)   , rot_amounts=(0,phi_y,0), rot_order='XYZ')
# # twr.connectTo(nac, type='Rigid', rel_pos=(0,0,L_T)  , rot_amounts=(0,0,0),rot_order='ZYX')
# # ## --- Defining Body rotational velocities
# # omega_TE = twr.ang_vel_in(ref)  # Angular velocity of nacelle in inertial frame
# # omega_NT = nac.ang_vel_in(twr.frame)  # Angular velocity of nacelle in inertial frame
# # # --- Kinetics
# # g = symbols('g')
# # T = symbols('T')
# # grav_T = (twr.masscenter, -twr.mass * g * ref.frame.z)
# # grav_N = (nac.masscenter, -nac.mass * g * ref.frame.z)
# # R=Point('R')
# # R.set_pos(nac.origin, x_NR * nac.frame.x + z_NR* nac.frame.z)
# # R.set_vel(nac.frame, 0 * nac.frame.x)
# # R.v2pt_theory(nac.origin, ref.frame, nac.frame)
# # #thrustN = (nac.masscenter, T * nac.frame.x)
# # thrustN = (R, T *cos(theta_tilt) * nac.frame.x -T *sin(theta_tilt) * nac.frame.z)
# # # --- Kane
# # #kdeqsSubs = [ (twr.qd[i], twr.qdot[i]) for i,_ in enumerate(twr.q)]; 
# # kdeqsSubs=[(xd, diff(x,time))]; 
# # kdeqsSubs+=[ (omega_y_T, omega_TE.dot(ref.frame.y).simplify())]  
# # #kdeqsSubs+=[ (omega_x_R, omega_RN.dot(rot.frame.x).simplify())]  
# # kdeqs = [s-v for (s,v) in kdeqsSubs]
# # coordinates = [x, phi_y] 
# # speeds      = [xd, omega_y_T]
# # qspeeds     = [diff(q,dynamicsymbols._t) for q in coordinates]
# # body_loads  = [(nac,grav_N), (nac,thrustN)];  loads = [f[1] for f in body_loads]
# # bodies      = [twr,nac]
# # #frstar_rigid   , MM_rigid     = kane_frstar(bodies, coordinates, speeds, kdeqs, ref.origin, ref.frame, Omega_Subs, Mform='symbolic')
# # #MM_rigid
# # kane = YAMSKanesMethod(ref.frame, coordinates, speeds, kdeqs)
# # fr, frstar  = kane.kanes_equations(bodies, loads, Mform='symbolic')
# # MM_rigid = kane.mass_matrix
# # forcing_rigid=kane.forcing
# # MM_rigid
# # forcing_rigid
# # 
# # 
# # # In[4]:
# # 
# # 
# # from yams.yams_sympy import *
# # 
# # ref = YAMSInertialBody('E') 
# # #twr = YAMSRigidBody('T', rho_G = [0,0,z_G_T]      , J_diag=True) 
# # twr = YAMSFlexibleBody('T',1, directions=['x','y','x','y']) 
# # nac = YAMSRigidBody('N', rho_G = [x_G_N ,0, z_G_N], J_diag=True) 
# # #nac = YAMSRigidBody('N', rho_G = [x_G_N ,0, 0], J_diag=True) 
# # 
# # #ref.connectTo(twr, type='Free' , rel_pos=(x,0,0)   , rot_amounts=(0,phi_y,0), rot_order='XYZ')
# # ref.connectTo(twr, type='Rigid', rel_pos=(0,0,0)   , rot_amounts=(0,0,0), rot_order='XYZ')
# # twr.connectTo(nac, type='Rigid', rel_pos=(0,0,twr.L)  , rot_amounts=(0,0,0),rot_order='ZYX', doSubs=True, rot_type_elastic='SmallRot')
# # #from yams.yams_sympy import kane_frstar
# # bodies = [twr,nac]
# # frstar = kane_frstar(bodies, qspeeds, ref.origin, ref.frame, Omega_Subs, Mform='symbolic')
# # frstar
# # 
# # 
# # # In[7]:
# # 
# # 
# # # --- Important Issue 
# # from yams.yams_sympy import coord2vec
# # a,b = symbols('a,b')
# # vv = a * nac.frame.x + b * nac.frame.z
# # vv
# # vv.express(ref.frame)
# # v_matref = vv.to_matrix(ref.frame)
# # v_ref=coord2vec(v_matref, ref.frame)
# # v_ref
# # v_ref.express(nac.frame).simplify()
# # 
# # A=nac.frame.dcm(ref.frame)
# # B=ref.frame.dcm(nac.frame)
# # A*B
# # 
# # 
# # # In[6]:
# # 
# # 
# # twr.q[0]
# # twr.ucSubs[0][1]
# # twr.ucList
# # twr.vcList
# # twr.qddot
# # #frstar.collect(twr.vcList[0])
# # frstar[0].subs(Omega_Subs).subs(Subs).expand().collect(twr.qddot+twr.qdot).simplify()
# # frstar2[0].subs(Omega_Subs).subs(Subs).expand().collect(twr.qddot+twr.qdot).simplify()
# # #frstar.expand().simplify().collect(twr.q[0])
# # 
# # 
# # # In[7]:
# # 
# # 
# # 
# # #omega_TE
# # twr.uc
# # twr.ucSubs
# # twr.alphaSubs
# # ##nac.masscenter.vel(ref.frame).express(ref.frame)
# # #twr.masscenter.vel(ref.frame)
# # #twr.origin.vel(ref.frame)
# # #nac.origin.vel(ref.frame)
# # #nac.origin.vel(twr.frame)
# # nac.masscenter.vel(twr.frame).to_matrix(twr.frame).simplify()
# # nac.masscenter.vel(twr.frame).to_matrix(twr.frame).simplify().subs(twr.v2Subs)
# # omega_NT
# # #omega_NT.subs(twr.v2Subs)
# # omega_NT2 = DCMtoOmega(twr.frame.dcm(nac.frame), twr.frame)
# # omega_NT2
# # # --- Misc tests
# # #twr = YAMSFlexibleBody('T',2) 
# # #twr.Cr.M0
# # #twr.Cr.M1
# # #twr.Me.M0
# # #twr.Oe.M0
# # #twr.bodyMassMatrix([0,0])
# # #twr.bodyMassMatrix()
# # #twr.bodyQudraticForce()
# # #omega_x,omega_y,omega_z=symbols('omega_x,omega_y,omega_z')
# # #omega=Matrix([[omega_x],[omega_y],[omega_z]]).transpose()
# # #twr.bodyQuadraticForce(omega,[0,0],twr.qd)
# # #twr.bodyElasticForce(twr.q,twr.qd)
# # #Matrix(twr.q).shape
# # #twr.masscenter
# # #twr.origin
# # 
