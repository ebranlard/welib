import unittest
import numpy as np
import sympy as sp
from sympy.physics.mechanics import dynamicsymbols

from welib.yams.yams_rec import YAMSRecGroundBody, YAMSRecRigidBody, YAMSRecUniformBeamBody, YAMSRecBeamBody
from welib.yams.yams_sympy import YAMSInertialBody, YAMSRigidBody, YAMSFlexibleBody


class TestKinematicsExport(unittest.TestCase):
    PRINT_COMPARISONS = False

    @staticmethod
    def _sympy_vec_to_numpy(vec, frame, subs):
        return np.array(vec.to_matrix(frame).subs(subs), dtype=float).reshape(-1)

    @staticmethod
    def _sympy_mat_to_numpy(mat, subs):
        return np.array(mat.subs(subs), dtype=float)

    @staticmethod
    def _rec_sympy_mat_to_numpy(mat, subs):
        return np.array(sp.Matrix(mat).subs(subs), dtype=float)

    def _debug_side_by_side(self, label, left, right):
        if not self.PRINT_COMPARISONS:
            return
        print('\n==== {} ===='.format(label))
        print('LEFT:\n{}'.format(left))
        print('RIGHT:\n{}'.format(right))

    def test_rec_export_schema(self):
        grd = YAMSRecGroundBody()
        rb = YAMSRecRigidBody(name='R', mass=1.0, J=np.zeros((3, 3)))
        grd.connectTo(rb, Type='Rigid')
        n = grd.setupDOFIndex()
        grd.setDOF(np.zeros(n))

        payload = rb.kinematics_export()
        self.assertEqual(payload['flavor'], 'yams_rec')
        self.assertIn('B', payload)
        self.assertIn('Bhat_x_bc', payload)
        self.assertIn('R_b2g', payload)
        self.assertIsNotNone(payload['B'])

    def test_sympy_export_schema(self):
        e = YAMSInertialBody('E')
        r = YAMSRigidBody('R')
        e.connectTo(r, type='Rigid', rel_pos=(0, 0, 0))

        payload = r.kinematics_export()
        self.assertEqual(payload['flavor'], 'yams_sympy')
        self.assertIn('B', payload)
        self.assertIn('Bhat_x_bc', payload)
        self.assertIn('R_b2g', payload)
        # Sympy path may not expose B-matrices yet; schema still includes them.
        self.assertTrue('B' in payload and 'B_inB' in payload and 'BB_inB' in payload)

    def test_rec_export_spherical_joint_chain(self):
        grd = YAMSRecGroundBody()
        r1 = YAMSRecRigidBody(name='R1', mass=1.0, J=np.zeros((3, 3)))
        r2 = YAMSRecRigidBody(name='R2', mass=1.0, J=np.zeros((3, 3)))

        # Ground to first body with one rotational dof, then second body attached at an offset.
        grd.connectTo(r1, Type='SphericalJoint', Point=[0, 0, 0], JointRotations=['z'])
        r1.connectTo(r2, Type='SphericalJoint', Point=[1, 0, 0], JointRotations=['z'])

        n = grd.setupDOFIndex()
        self.assertEqual(n, 2)
        grd.setDOF(np.array([0.1, -0.2]))

        payloads = grd.kinematics_export_tree()
        self.assertGreaterEqual(len(payloads), 3)
        self.assertIsNotNone(payloads[-1]['B'])
        self.assertEqual(payloads[-1]['B'].shape[1], 2)

    def test_sympy_export_joint_chain_has_B(self):
        q1 = dynamicsymbols('q1')
        q2 = dynamicsymbols('q2')

        e = YAMSInertialBody('E')
        r1 = YAMSRigidBody('R1')
        r2 = YAMSRigidBody('R2')

        e.connectTo(r1, type='Joint', rel_pos=(0, 0, 0), rot_type='Axis', rot_amounts=(q1, e.frame.z))
        r1.connectTo(r2, type='Joint', rel_pos=(1, 0, 0), rot_type='Axis', rot_amounts=(q2, r1.frame.z))

        payload = r2.kinematics_export()
        self.assertIsNotNone(payload['B'])
        self.assertIsNotNone(payload['B_inB'])
        self.assertIsNotNone(payload['BB_inB'])
        self.assertIsNotNone(r2.B_matrix())
        self.assertIsNotNone(r2.B_matrix(in_body=True))
        self.assertIsNotNone(r2.BB_matrix())
        self.assertEqual(payload['B'].shape[0], 6)
        self.assertGreaterEqual(payload['B'].shape[1], 2)
        self.assertEqual(payload['B_inB'].shape, payload['B'].shape)
        self.assertEqual(payload['BB_inB'].shape[0], 6)
        self.assertIsNotNone(payload['speed_symbols'])

    def test_flexible_tower_with_top_body_alignment(self):
        # --- Recursive flavor: uniform flexible beam + rigid top body
        grd = YAMSRecGroundBody()
        twr = YAMSRecUniformBeamBody(
            'T', nShapes=1, nSpan=21, L=10.0, EI0=1e6, m=100.0,
            Mtop=0.0, bAxialCorr=False, bStiffening=False, gravity=9.81,
            main_axis='x', shapeFunctions='admissible', bottomBC='clamped', topBC='free')
        nac_rec = YAMSRecRigidBody(name='N', mass=1.0, J=np.zeros((3, 3)))
        grd.connectTo(twr, Type='Rigid')
        twr.connectTo(nac_rec, BodyPoint='LastPoint', Type='Rigid')

        n = grd.setupDOFIndex()
        self.assertEqual(n, 1)
        grd.setDOF(np.array([0.2]))
        p_rec = nac_rec.kinematics_export()

        # --- Sympy flavor: flexible body + rigid top body connected at tip
        e = YAMSInertialBody('E')
        twr_sym = YAMSFlexibleBody('T', nq=1, directions=['x'])
        nac_sym = YAMSRigidBody('N')
        e.connectTo(twr_sym, type='Rigid', rel_pos=(0, 0, 0))
        twr_sym.connectToTip(nac_sym, type='Rigid', rel_pos=[twr_sym.L, sp.Integer(0), sp.Integer(0)])
        p_sym = nac_sym.kinematics_export()

        # --- Canonical export comparability checks
        self.assertIsNotNone(p_rec['B'])
        self.assertIsNotNone(p_sym['B'])
        self.assertEqual(p_rec['B'].shape[0], 6)
        self.assertEqual(p_sym['B'].shape[0], 6)
        self.assertEqual(p_rec['B'].shape[1], 1)
        self.assertEqual(p_sym['B'].shape[1], 1)
        self.assertIsNotNone(p_rec['pos_global'])
        self.assertIsNotNone(p_sym['pos_global'])

    def test_numeric_value_parity_one_dof_rigid_chain(self):
        L = 2.5
        q_val = 0.3

        # --- Recursive flavor
        grd = YAMSRecGroundBody()
        r1 = YAMSRecRigidBody(name='R1', mass=1.0, J=np.zeros((3, 3)))
        nac_rec = YAMSRecRigidBody(name='N', mass=1.0, J=np.zeros((3, 3)))
        grd.connectTo(r1, Type='SphericalJoint', Point=[0, 0, 0], JointRotations=['z'])
        r1.connectTo(nac_rec, Type='Rigid', Point=[L, 0, 0])
        n = grd.setupDOFIndex()
        self.assertEqual(n, 1)
        grd.setDOF(np.array([q_val]))
        p_rec = nac_rec.kinematics_export()

        # --- Sympy flavor
        q1 = dynamicsymbols('q1')
        e = YAMSInertialBody('E')
        s1 = YAMSRigidBody('R1')
        nac_sym = YAMSRigidBody('N')
        e.connectTo(s1, type='Joint', rel_pos=(0, 0, 0), rot_type='Axis', rot_amounts=(q1, e.frame.z))
        s1.connectTo(nac_sym, type='Rigid', rel_pos=(L, 0, 0))
        p_sym = nac_sym.kinematics_export()

        qdot = q1.diff()
        subs = {q1: q_val, qdot: 0.7}

        pos_sym = self._sympy_vec_to_numpy(p_sym['pos_global'], e.frame, subs)
        B_sym = self._sympy_mat_to_numpy(p_sym['B'], subs)
        R_sym = self._sympy_mat_to_numpy(p_sym['R_b2g'], subs)

        np.testing.assert_allclose(p_rec['pos_global'], pos_sym, rtol=1e-10, atol=1e-10)
        np.testing.assert_allclose(p_rec['B'], B_sym, rtol=1e-10, atol=1e-10)
        np.testing.assert_allclose(p_rec['R_b2g'], R_sym, rtol=1e-10, atol=1e-10)

    def test_numeric_value_parity_two_dof_rigid_chain(self):
        # Double pendulum with symbolic lengths in rec(sympy=True) and yams_sympy.
        q1, q2 = sp.symbols('q1 q2')
        dq1, dq2 = sp.symbols('dq1 dq2')
        L1, L2 = sp.symbols('L1 L2', positive=True)
        m1 = 1.2
        m2 = 0.8

        # --- Recursive flavor (sympy mode)
        grd = YAMSRecGroundBody(sympy=True)
        b1_rec = YAMSRecRigidBody(name='B1', mass=m1, J=sp.zeros(3), rho_G=(L1/2, 0, 0), sympy=True)
        b2_rec = YAMSRecRigidBody(name='B2', mass=m2, J=sp.zeros(3), rho_G=(L2/2, 0, 0), sympy=True)
        grd.connectTo(b1_rec, Type='Joint', Point=[0, 0, 0], JointRotations=['z'])
        b1_rec.connectTo(b2_rec, Type='Joint', Point=[L1, 0, 0], JointRotations=['z'])
        n = grd.setupDOFIndex()
        self.assertEqual(n, 2)
        grd.setDOF(np.array([q1, q2], dtype=object))

        # --- Sympy flavor
        q1s = dynamicsymbols('q1s')
        q2s = dynamicsymbols('q2s')
        e = YAMSInertialBody('E')
        b1_sym = YAMSRigidBody('B1', mass=m1, J_G=(0, 0, 0), rho_G=(L1/2, 0, 0))
        b2_sym = YAMSRigidBody('B2', mass=m2, J_G=(0, 0, 0), rho_G=(L2/2, 0, 0))
        e.connectTo(b1_sym, type='Joint', rel_pos=(0, 0, 0), rot_type='Axis', rot_amounts=(q1s, e.frame.z))
        b1_sym.connectTo(b2_sym, type='Joint', rel_pos=(L1, 0, 0), rot_type='Axis', rot_amounts=(q2s, b1_sym.frame.z))

        # --- Analytical double-pendulum mass matrix (point-mass rods with COM at L/2)
        c2 = sp.cos(q2)
        l1c = L1/2
        l2c = L2/2
        M11 = m1*l1c**2 + m2*(L1**2 + l2c**2 + 2*L1*l2c*c2)
        M12 = m2*(l2c**2 + L1*l2c*c2)
        M22 = m2*l2c**2
        M_ana = sp.Matrix([[M11, M12], [M12, M22]])
        M_rec = sp.simplify(sp.Matrix(grd.M))
        self._debug_side_by_side('Double pendulum M (rec vs analytical)', M_rec, M_ana)
        self.assertTrue(sp.simplify(M_rec - M_ana) == sp.zeros(2, 2))

        # --- Origin and COM velocities from rec B
        qdot_vec1 = sp.Matrix([dq1])
        qdot_vec2 = sp.Matrix([dq1, dq2])
        vO1_rec = sp.Matrix(b1_rec.B_matrix()[:3, :]) * qdot_vec1
        vO2_rec = sp.Matrix(b2_rec.B_matrix()[:3, :]) * qdot_vec2
        om1_rec = sp.Matrix(b1_rec.B_matrix()[3:6, :]) * qdot_vec1
        om2_rec = sp.Matrix(b2_rec.B_matrix()[3:6, :]) * qdot_vec2
        rG1_rec = sp.Matrix(b1_rec.R_b2g) * sp.Matrix([L1/2, 0, 0])
        rG2_rec = sp.Matrix(b2_rec.R_b2g) * sp.Matrix([L2/2, 0, 0])
        vG1_rec = vO1_rec + om1_rec.cross(rG1_rec)
        vG2_rec = vO2_rec + om2_rec.cross(rG2_rec)

        vO1_ana = sp.Matrix([0, 0, 0])
        vO2_ana = sp.Matrix([-L1*sp.sin(q1)*dq1, L1*sp.cos(q1)*dq1, 0])
        vG1_ana = sp.Matrix([-(L1/2)*sp.sin(q1)*dq1, (L1/2)*sp.cos(q1)*dq1, 0])
        vG2_ana = sp.Matrix([
            -L1*sp.sin(q1)*dq1 - (L2/2)*sp.sin(q1 + q2)*(dq1 + dq2),
            L1*sp.cos(q1)*dq1 + (L2/2)*sp.cos(q1 + q2)*(dq1 + dq2),
            0,
        ])

        self.assertTrue(sp.simplify(vO1_rec - vO1_ana) == sp.zeros(3, 1))
        self.assertTrue(sp.simplify(vO2_rec - vO2_ana) == sp.zeros(3, 1))
        self.assertTrue(sp.simplify(vG1_rec - vG1_ana) == sp.zeros(3, 1))
        self.assertTrue(sp.simplify(vG2_rec - vG2_ana) == sp.zeros(3, 1))

        # --- Cross-flavor B parity at one numeric operating point
        subs_rec = {q1: 0.4, q2: -0.25, L1: 2.0, L2: 1.2}
        subs_sym = {q1s: 0.4, q2s: -0.25, q1s.diff(): 0.3, q2s.diff(): -0.15, L1: 2.0, L2: 1.2}
        B2_rec = self._rec_sympy_mat_to_numpy(b2_rec.B_matrix(), subs_rec)
        B2_sym = self._sympy_mat_to_numpy(b2_sym.B_matrix(), subs_sym)
        self._debug_side_by_side('Double pendulum B body-2 (rec vs sympy)', B2_rec, B2_sym)
        np.testing.assert_allclose(B2_rec, B2_sym, rtol=1e-10, atol=1e-10)

    def test_rec_sympy_true_joint_connection_parity(self):
        q_rec = sp.symbols('q_rec')
        q_val = 0.31
        qdot_val = -0.27
        L = 2.5

        # --- Recursive flavor in sympy mode
        grd = YAMSRecGroundBody(sympy=True)
        r1 = YAMSRecRigidBody(name='R1', mass=1.0, J=sp.Matrix(np.zeros((3, 3))), sympy=True)
        nac_rec = YAMSRecRigidBody(name='N', mass=1.0, J=sp.Matrix(np.zeros((3, 3))), sympy=True)
        grd.connectTo(r1, Type='Joint', Point=[0, 0, 0], JointRotations=['z'])
        r1.connectTo(nac_rec, Type='Rigid', Point=[L, 0, 0])
        n = grd.setupDOFIndex()
        self.assertEqual(n, 1)
        grd.setDOF(np.array([q_rec], dtype=object))
        p_rec = nac_rec.kinematics_export()

        # --- Sympy mechanics flavor
        q1 = dynamicsymbols('q1')
        e = YAMSInertialBody('E')
        s1 = YAMSRigidBody('R1')
        nac_sym = YAMSRigidBody('N')
        e.connectTo(s1, type='Joint', rel_pos=(0, 0, 0), rot_type='Axis', rot_amounts=(q1, e.frame.z))
        s1.connectTo(nac_sym, type='Rigid', rel_pos=(L, 0, 0))
        p_sym = nac_sym.kinematics_export()

        subs_rec = {q_rec: q_val}
        subs_sym = {q1: q_val, q1.diff(): qdot_val}

        B_rec = self._rec_sympy_mat_to_numpy(p_rec['B'], subs_rec)
        B_sym = self._sympy_mat_to_numpy(p_sym['B'], subs_sym)
        BinB_rec = self._rec_sympy_mat_to_numpy(p_rec['B_inB'], subs_rec)
        BinB_sym = self._sympy_mat_to_numpy(p_sym['B_inB'], subs_sym)
        BB_rec = self._rec_sympy_mat_to_numpy(p_rec['BB_inB'], subs_rec)
        BB_sym = self._sympy_mat_to_numpy(p_sym['BB_inB'], subs_sym)
        R_rec = self._rec_sympy_mat_to_numpy(p_rec['R_b2g'], subs_rec)
        R_sym = self._sympy_mat_to_numpy(p_sym['R_b2g'], subs_sym)

        np.testing.assert_allclose(B_rec, B_sym, rtol=1e-10, atol=1e-10)
        np.testing.assert_allclose(BinB_rec, BinB_sym, rtol=1e-10, atol=1e-10)
        np.testing.assert_allclose(BB_rec, BB_sym, rtol=1e-10, atol=1e-10)
        np.testing.assert_allclose(R_rec, R_sym, rtol=1e-10, atol=1e-10)

    def test_rec_free_connection_sympy_true(self):
        tx, ty, rz = sp.symbols('tx ty rz')

        grd = YAMSRecGroundBody(sympy=True)
        body = YAMSRecRigidBody(name='F', mass=1.0, J=sp.Matrix(np.zeros((3, 3))), sympy=True)
        grd.connectTo(body, Type='Free', Point=[0, 0, 0], JointTranslations=['x', 'y'], JointRotations=['z'])
        n = grd.setupDOFIndex()
        self.assertEqual(n, 3)
        grd.setDOF(np.array([tx, ty, rz], dtype=object))

        payload = body.kinematics_export()
        B_expected = sp.Matrix([
            [1, 0, 0],
            [0, 1, 0],
            [0, 0, 0],
            [0, 0, 0],
            [0, 0, 0],
            [0, 0, 1],
        ])
        self.assertTrue(sp.simplify(sp.Matrix(payload['B']) - B_expected) == sp.zeros(6, 3))

    def test_flexible_case_bhat_and_B_with_rec_sympy_true(self):
        q_rec = sp.symbols('q_rec')
        q_rec2 = sp.symbols('q_rec2')
        ux1_val = 1.7
        uy2_val = 0.9
        vy1_val = 0.35
        vx2_val = -0.22
        L_val = 10.0

        # --- Recursive flexible body in sympy mode
        grd = YAMSRecGroundBody(sympy=True)
        twr_rec = YAMSRecBeamBody(name='T', directions=['x', 'y'], sympy=True, main_axis='x')
        nac_rec = YAMSRecRigidBody(name='N', mass=1.0, J=sp.Matrix(np.zeros((3, 3))), sympy=True)
        grd.connectTo(twr_rec, Type='Rigid')
        twr_rec.connectTo(nac_rec, Type='Rigid', Point=[sp.symbols('L'), 0, 0])
        n = grd.setupDOFIndex()
        self.assertEqual(n, 2)
        grd.setDOF(np.array([q_rec, q_rec2], dtype=object))
        p_twr_rec = twr_rec.kinematics_export()
        p_nac_rec = nac_rec.kinematics_export()

        # --- Sympy flexible body with matching no-tip-rotation setup
        e = YAMSInertialBody('E')
        twr_sym = YAMSFlexibleBody('T', nq=2, directions=['x', 'y'])
        nac_sym = YAMSRigidBody('N')
        e.connectTo(twr_sym, type='Rigid', rel_pos=(0, 0, 0))
        twr_sym.connectToTip(nac_sym, type='Rigid', rel_pos=[twr_sym.L, sp.Integer(0), sp.Integer(0)])
        p_twr_sym = twr_sym.kinematics_export()
        p_nac_sym = nac_sym.kinematics_export()

        subs_rec = {
            sp.symbols('ux1c'): ux1_val,
            sp.symbols('uy2c'): uy2_val,
            sp.symbols('vy1c'): vy1_val,
            sp.symbols('vx2c'): vx2_val,
            sp.symbols('L'): L_val,
            q_rec: 0.2,
            q_rec2: -0.1,
        }
        subs_sym = {
            sp.symbols('u_xT1c'): ux1_val,
            sp.symbols('u_yT2c'): uy2_val,
            sp.symbols('v_yT1c'): vy1_val,
            sp.symbols('v_xT2c'): vx2_val,
            twr_sym.L: L_val,
            dynamicsymbols('q_T1'): 0.2,
            dynamicsymbols('q_T2'): -0.1,
            dynamicsymbols('q_T1').diff(): 0.4,
            dynamicsymbols('q_T2').diff(): -0.3,
        }

        Bhatx_rec = self._rec_sympy_mat_to_numpy(twr_rec.Bhat_matrix('x'), subs_rec)
        Bhatx_sym = self._sympy_mat_to_numpy(twr_sym.Bhat_matrix('x'), subs_sym)
        Bhatt_rec = self._rec_sympy_mat_to_numpy(twr_rec.Bhat_matrix('t'), subs_rec)
        Bhatt_sym = self._sympy_mat_to_numpy(twr_sym.Bhat_matrix('t'), subs_sym)
        B_rec = self._rec_sympy_mat_to_numpy(nac_rec.B_matrix(), subs_rec)
        B_sym = self._sympy_mat_to_numpy(nac_sym.B_matrix(), subs_sym)
        BB_rec = nac_rec.BB_matrix()
        BB_sym = nac_sym.BB_matrix()

        self._debug_side_by_side('Flexible Bhat_x (rec vs sympy)', Bhatx_rec, Bhatx_sym)
        self._debug_side_by_side('Flexible Bhat_t (rec vs sympy)', Bhatt_rec, Bhatt_sym)
        self._debug_side_by_side('Flexible B (nac rec vs sympy)', B_rec, B_sym)

        np.testing.assert_allclose(Bhatx_rec, Bhatx_sym, rtol=1e-10, atol=1e-10)
        np.testing.assert_allclose(Bhatt_rec, Bhatt_sym, rtol=1e-10, atol=1e-10)
        self.assertEqual(B_rec.shape, B_sym.shape)
        np.testing.assert_allclose(B_rec[:3, :], B_sym[:3, :], rtol=1e-10, atol=1e-10)
        self.assertGreater(np.linalg.norm(B_rec[3:, :]), 0.0)
        self.assertGreater(np.linalg.norm(B_sym[3:, :]), 0.0)
        self.assertEqual(BB_rec.shape, BB_sym.shape)


if __name__ == '__main__':
    unittest.main()
