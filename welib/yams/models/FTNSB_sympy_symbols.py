""" 
define common symbols used for a FNTSR/FTRNA model of a wind turbine

  M: Mooring line attachement
  B: Mooring line attachement

  F: Floater/foundation
  T: Tower


  N: Nacelle
  S: Shaft
  R: Rotor

  B1,B2,B3: Blades

  RNA

"""
from sympy import Symbol, symbols, Function
from sympy.physics.mechanics import dynamicsymbols



time = symbols('t')

# 6 DOF of the floater
phi_x, phi_y, phi_z = dynamicsymbols('phi_x, phi_y, phi_z')
x, y, z             = dynamicsymbols('x, y, z')
xd, yd, zd          = dynamicsymbols('xd, yd, zd') #dynamicsymbols('x, y, z',1)
  
# Nacelle/shaft angles
theta_yaw, theta_tilt, psi_0 = symbols('theta_yaw, theta_tilt, psi_0')  #NOTE: not dynamic, constant, use q_yaw instead!
alpha_y, alpha_x           = dynamicsymbols('alpha_y, alpha_x')
q_psi                      = dynamicsymbols('psi')
q_yaw, q_tilt              = dynamicsymbols('q_yaw, q_tilt')
qd_yaw, qd_tilt            = dynamicsymbols('qd_yaw, qd_tilt')

# Blades
theta_cone, theta_pitch   = symbols('theta_c, theta_p') # NOTE: not dynamic, constant use q_pitch instead
q_pitch                   = dynamicsymbols('q_p')
qd_pitch                  = dynamicsymbols('qd_p')

# Angular velocities of bodies
omega_x_F, omega_y_F, omega_z_F = dynamicsymbols('omega_x_F, omega_y_F, omega_z_F')
omega_x_T, omega_y_T, omega_z_T = dynamicsymbols('omega_x_T, omega_y_T, omega_z_T')
omega_x_N, omega_y_N, omega_z_N = dynamicsymbols('omega_x_N, omega_y_N, omega_z_N')
omega_x_R, omega_y_R, omega_z_R = dynamicsymbols('omega_x_R, omega_y_R, omega_z_R')

# omega individual blades..

# --- Inertias
M_F, M_T, M_N, M_R        = symbols('M_F,M_T,M_N,M_R')           # Masses: Foundation/Tower/Nacelle/Rotor
Jxx_R, Jyy_R, Jzz_R, JO_R = symbols('Jxx_R, Jyy_R, Jzz_R, JO_R') # NOTE: JO                     = Jyy = Jzz for a three bladed rotor!
Jxx_T, Jyy_T, Jzz_T       = symbols('Jxx_T, Jyy_T, Jzz_T')
Jxx_F, Jyy_F, Jzz_F       = symbols('Jxx_F, Jyy_F, Jzz_F')
Jxx_N, Jyy_N, Jzz_N       = symbols('Jxx_N, Jyy_N, Jzz_N')
Jxx_B, Jyy_B, Jzz_B       = symbols('Jxx_B, Jyy_B, Jzz_B')

L_F, L_T, L_R, L_B = symbols('L_F, L_T, L_R, L_B') # length

# --- COGs
z_FG     = symbols('z_FG')       # Position of Foundation COG in F, measured from point T
z_TG     = symbols('z_TG')       # Position of Tower COG in T
x_NG,z_NG = symbols('x_NG, z_NG') # Position of Nacelle COG in N
x_RNAG,z_RNAG = symbols('x_RNAG, z_RNAG') # Position of Nacelle COG in N
y_RNAG = symbols('y_RNAG') # Position of Nacelle COG in N
x_BG, y_BG, z_BG = symbols('x_BG, y_BG, z_BG') # Position of Blade COG in Blade coordinate system


# Points
x_NR, z_NR = symbols('x_NR, z_NR') # From Nacelle origin to rotor center 
z_TM, z_TB = symbols('z_TM, z_TB') # Position of mooring line attachment point and Buoyancy center in F, measured from point T
z_T0       = symbols('z_T0') # Position of tower bottom wrt to the "(0,0,0)" point
#x_RB, z_RB = symbols('x_RB, z_RB') # From rotor center to blade root
r_hub = symbols('r_h') # From rotor center to blade root


# --- Loads
gravity = symbols('g')


# Subs2D = [(phi_x,0),(phi_z,0),(y,0),(z,0),(theta_yaw,0),(theta_tilt,0),
#           (omega_z_T,0),(omega_z_N,0),(omega_z_R,0),
#           (omega_x_T,0),(omega_x_N,0),(omega_x_R,0)]

# --- Simplifying substitution 
# Nacelle COG on tower top
subs_NacGatTT= [(z_NG,0),(x_NG,0)]

subs_NormShapeT1x=[ (Symbol('u_xT1c'),1) ]

subs_SmallAngleFnd=[(phi_x,0), (phi_y,0), (phi_z,0)]
# Subs_SmallAngleTwrRot=[(twr.vcList[0],0)]    # <<<< TODO

smallAngleFnd    = [phi_x,phi_y,phi_z]
# SmallAngleTwrRot = [twr.vcList[0]]     # <<<< TODO


