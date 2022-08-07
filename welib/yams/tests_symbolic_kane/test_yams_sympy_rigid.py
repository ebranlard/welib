import unittest
import numpy as np

from welib.yams.yams_sympy import YAMSInertialBody, YAMSRigidBody
from welib.yams.yams_kane  import kane_fr
from welib.yams.yams_kane  import kane_fr_alt
from welib.yams.yams_kane  import kane_frstar
from welib.yams.yams_kane  import kane_frstar_alt
from welib.yams.yams_kane  import YAMSKanesMethod 

from sympy import symbols, simplify
from sympy import diff, Matrix
from sympy.physics.mechanics import inertia, KanesMethod
from sympy.physics.mechanics import dynamicsymbols, ReferenceFrame, Point

# --------------------------------------------------------------------------------}
# --- TESTS
# --------------------------------------------------------------------------------{
class Test(unittest.TestCase):
    def test_rigid_FTNSR_3DOF_kane(self):
        """ """

        # Main symbols
        time      = dynamicsymbols._t
        phi_y     = dynamicsymbols('phi_y')
        psi       = dynamicsymbols('psi')
        x         = dynamicsymbols('x')
        ux        = dynamicsymbols('x_d')       # derivative of x
        xd        = dynamicsymbols('x',1)       # derivative of x
        omega_y_T = dynamicsymbols('omega_y_T')
        omega_x_R = dynamicsymbols('omega_x_R')
        Jxx_R, Jyy_R, Jzz_R, JO_R = symbols('Jxx_R, Jyy_R, Jzz_R, JO_R') # NOTE: JO = Jyy = Jzz for a three bladed rotor!
        L_F, L_T, L_R  = symbols('L_F, L_T, L_R') # length
        z_G_F = symbols('z_G_F') # Position of Foundation COG , measured from T!
        z_G_T = symbols('z_G_T') # Position of Tower COG in T
        x_G_N,z_G_N = symbols('x_G_N, z_G_N') # Position of Nacelle COG in N
        x_R,z_R = symbols('x_R, z_R') # Position of Rotor center in N
        z_M, z_B = symbols('z_M, z_B') # Position of mooring line attachment point and Buoyancy center in F, measured from point T (negative)
        g       = symbols('g')
        B       = symbols('B')
        K_x     = symbols('K_x')
        K_phi_y = symbols('K_phi_y')
        theta_yaw, theta_tilt = symbols('theta_yaw, theta_tilt') # Position of mooring line attachment point and Buoyancy center in F, measured from point T (negative)

        # --- Bodies
        ref = YAMSInertialBody('E') 
        twr = YAMSRigidBody('T', rho_G = [0,0,z_G_T]      , J_form='diag') 
        fnd = YAMSRigidBody('F', rho_G = [0,0,z_G_F]      , J_form='diag') 
        nac = YAMSRigidBody('N', rho_G = [x_G_N ,0, z_G_N], J_form='diag') 
        rot = YAMSRigidBody('R', rho_G = [0,0,0]          , J_form='diag')
        rot.inertia = (inertia(rot.frame, Jxx_R, JO_R, JO_R), rot.origin)  # defining inertia at orign

        # --- Bodies connections
        ref.connectTo(twr, type='Free' , rel_pos=(x,0,0)   , rot_amounts=(0,phi_y,0), rot_order='XYZ')
        twr.connectTo(fnd, type='Rigid', rel_pos=(0,0,-L_F))
        twr.connectTo(nac, type='Rigid', rel_pos=(0,0,L_T)  , rot_amounts=(theta_yaw,theta_tilt,0),rot_order='ZYX')
        #twr.connectTo(nac, type='Rigid', rel_pos=(0,0,L_T)  , rot_amounts=(0,0,0),rot_order='ZYX')
        nac.connectTo(rot, type='Joint', rel_pos=(x_R,0,z_R), rot_amounts=(0,0,psi), rot_order='ZYX')

        # Point of application for forces
        P_B = twr.origin.locatenew('P_B', z_B * fnd.frame.z) # <<<< Measured from T
        P_M = twr.origin.locatenew('P_M', z_M * fnd.frame.z) # <<<< Measured from T
        ## --- Introducing "variables" for velocities
        omega_TE = twr.ang_vel_in(ref)  # Angular velocity of nacelle in inertial frame
        omega_RN = rot.ang_vel_in(nac)  # Angular velocity of rotor wrt Nacelle (omega_R-omega_N)
        Omega_Subs = [ 
            (ux, xd),
            (omega_y_T, omega_TE.dot(twr.frame.y).simplify()), 
            (omega_x_R, omega_RN.dot(rot.frame.x).simplify()),  
        ]
        Omega_Subs
        kdeqs = [s-v for (s,v) in Omega_Subs]

        ## Linear Velocities of Points
        P_B.v2pt_theory(twr.origin, ref.frame, twr.frame); # PB & T are fixed in e_T
        P_M.v2pt_theory(twr.origin, ref.frame, twr.frame); # PM & T are fixed in e_T

        # --- Kinetics
        # Gravity
        grav_F = (fnd.masscenter, -fnd.mass * g * ref.frame.z)
        grav_T = (twr.masscenter, -twr.mass * g * ref.frame.z)
        grav_N = (nac.masscenter, -nac.mass * g * ref.frame.z)
        grav_R = (rot.masscenter, -rot.mass * g * ref.frame.z)
        # Buyancy
        F_buy = (P_B, B * ref.frame.z)
        # Restoring mooring
        F_moor = (P_M,  -K_x * x *ref.frame.x )
        # Ext torques
        M_moor = (fnd.frame, -K_phi_y * phi_y *ref.frame.y)

        # --- Kane Prep
        coordinates = [x,phi_y, psi]
        speeds      = [ux, omega_y_T, omega_x_R]
        loads      = [grav_F            , grav_T,      grav_N,       grav_R,        F_buy,        F_moor,       M_moor]
        body_loads = [(fnd,grav_F), (twr,grav_T), (nac,grav_N), (rot,grav_R), (fnd, F_buy), (fnd, F_moor), (fnd,M_moor)]
        bodies     = [fnd,twr,nac,rot]

        # --- Kane
        kane = KanesMethod(ref.frame, coordinates, speeds, kdeqs)
        fr, frstar    = kane.kanes_equations(bodies, loads)
        MM            = kane.mass_matrix_full
        forcing_vector = kane.forcing_full

        # --- YAMSKane
        y_kane = YAMSKanesMethod(ref.frame, coordinates, speeds, kdeqs)
        fr_y, frstar_y  = y_kane.kanes_equations(bodies, loads)
        #MM            = y_kane.mass_matrix_full
        #forcing_vector = y_kane.forcing_full

        # --- "Manual" Kane
        coordinates = [x,phi_y, psi]
        frstar_new,MM_new,MMFull_new  = kane_frstar(bodies, coordinates, speeds, kdeqs, ref.origin, ref.frame)
        fr_new             = kane_fr(body_loads, speeds, ref.frame)

        # --- Kane "alt"
        frstar_alt, MM_alt,MMFull_alt = kane_frstar_alt(bodies, coordinates, speeds, kdeqs, ref.frame)
        fr_alt            =  kane_fr_alt(loads, coordinates, speeds, kdeqs, ref.frame)

        # --- Compare manual Kane with Kane
        DeltaFrStar = simplify(frstar-frstar_new)
        DeltaFr     = simplify(fr-fr_new)
        self.assertEqual(DeltaFrStar, Matrix([[0],[0],[0]]))
        self.assertEqual(DeltaFr, Matrix([[0],[0],[0]]))
        self.assertEqual(MM_new, kane._k_d)

        # --- Compare manual Kane with alt Kane
        DeltaFrStar = simplify(frstar_alt-frstar_new)
        DeltaFr     = simplify(fr_alt-fr_new)
        self.assertEqual(DeltaFrStar, Matrix([[0],[0],[0]]))
        self.assertEqual(DeltaFr, Matrix([[0],[0],[0]]))
        self.assertEqual(MM_alt, MM_new)

        # --- Compare yams Kane with Kane
        DeltaFrStar = simplify(frstar-frstar_y)
        DeltaFr     = simplify(fr-fr_y)
        self.assertEqual(DeltaFrStar, Matrix([[0],[0],[0]]))
        self.assertEqual(DeltaFr, Matrix([[0],[0],[0]]))
        self.assertEqual(y_kane._k_d, kane._k_d)
        self.assertEqual(y_kane.forcing_full, kane.forcing_full)
        self.assertEqual(y_kane.mass_matrix, kane.mass_matrix)
        self.assertEqual(y_kane.mass_matrix_full, kane.mass_matrix_full)





if __name__=='__main__':
    unittest.main()
