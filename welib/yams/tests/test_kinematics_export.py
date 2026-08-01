import unittest
import numpy as np
import sympy as sp
from sympy.physics.mechanics import dynamicsymbols

from welib.yams.yams_rec import YAMSRecGroundBody, YAMSRecRigidBody, YAMSRecUniformBeamBody, YAMSRecBeamBody
from welib.yams.yams_sympy import YAMSInertialBody, YAMSRigidBody, YAMSFlexibleBody


class TestKinematicsExport(unittest.TestCase):

    @staticmethod
    def _sympy_vec_to_numpy(vec, frame, subs):
        return np.array(vec.to_matrix(frame).subs(subs), dtype=float).reshape(-1)

    @staticmethod
    def _sympy_mat_to_numpy(mat, subs):
        return np.array(mat.subs(subs), dtype=float)

    @staticmethod
    def _rec_sympy_mat_to_numpy(mat, subs):
        return np.array(sp.Matrix(mat).subs(subs), dtype=float)

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
        L1 = 2.0
        L2 = 1.2
        q1_val = 0.4
        q2_val = -0.25

        # --- Recursive flavor
        grd = YAMSRecGroundBody()
        r1 = YAMSRecRigidBody(name='R1', mass=1.0, J=np.zeros((3, 3)))
        r2 = YAMSRecRigidBody(name='R2', mass=1.0, J=np.zeros((3, 3)))
        nac_rec = YAMSRecRigidBody(name='N', mass=1.0, J=np.zeros((3, 3)))
        grd.connectTo(r1, Type='SphericalJoint', Point=[0, 0, 0], JointRotations=['z'])
        r1.connectTo(r2, Type='SphericalJoint', Point=[L1, 0, 0], JointRotations=['z'])
        r2.connectTo(nac_rec, Type='Rigid', Point=[L2, 0, 0])
        n = grd.setupDOFIndex()
        self.assertEqual(n, 2)
        grd.setDOF(np.array([q1_val, q2_val]))
        p_rec = nac_rec.kinematics_export()

        # --- Sympy flavor
        q1 = dynamicsymbols('q1')
        q2 = dynamicsymbols('q2')
        e = YAMSInertialBody('E')
        s1 = YAMSRigidBody('R1')
        s2 = YAMSRigidBody('R2')
        nac_sym = YAMSRigidBody('N')
        e.connectTo(s1, type='Joint', rel_pos=(0, 0, 0), rot_type='Axis', rot_amounts=(q1, e.frame.z))
        s1.connectTo(s2, type='Joint', rel_pos=(L1, 0, 0), rot_type='Axis', rot_amounts=(q2, s1.frame.z))
        s2.connectTo(nac_sym, type='Rigid', rel_pos=(L2, 0, 0))
        p_sym = nac_sym.kinematics_export()

        subs = {q1: q1_val, q2: q2_val, q1.diff(): 0.3, q2.diff(): -0.15}

        pos_sym = self._sympy_vec_to_numpy(p_sym['pos_global'], e.frame, subs)
        B_sym = self._sympy_mat_to_numpy(p_sym['B'], subs)
        R_sym = self._sympy_mat_to_numpy(p_sym['R_b2g'], subs)

        np.testing.assert_allclose(p_rec['pos_global'], pos_sym, rtol=1e-10, atol=1e-10)
        np.testing.assert_allclose(p_rec['B'], B_sym, rtol=1e-10, atol=1e-10)
        np.testing.assert_allclose(p_rec['R_b2g'], R_sym, rtol=1e-10, atol=1e-10)

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
        ux_val = 1.7

        # --- Recursive flexible body in sympy mode
        grd = YAMSRecGroundBody(sympy=True)
        twr_rec = YAMSRecBeamBody(name='T', directions=['x'], sympy=True, main_axis='x')
        nac_rec = YAMSRecRigidBody(name='N', mass=1.0, J=sp.Matrix(np.zeros((3, 3))), sympy=True)
        grd.connectTo(twr_rec, Type='Rigid')
        twr_rec.connectTo(nac_rec, Type='Rigid', Point=[sp.symbols('L'), 0, 0])
        n = grd.setupDOFIndex()
        self.assertEqual(n, 1)
        grd.setDOF(np.array([q_rec], dtype=object))
        p_twr_rec = twr_rec.kinematics_export()
        p_nac_rec = nac_rec.kinematics_export()

        # --- Sympy flexible body with matching no-tip-rotation setup
        e = YAMSInertialBody('E')
        twr_sym = YAMSFlexibleBody('T', nq=1, directions=['x'], tip_rotate=False)
        nac_sym = YAMSRigidBody('N')
        e.connectTo(twr_sym, type='Rigid', rel_pos=(0, 0, 0))
        twr_sym.connectToTip(nac_sym, type='Rigid', rel_pos=[twr_sym.L, sp.Integer(0), sp.Integer(0)])
        p_twr_sym = twr_sym.kinematics_export()
        p_nac_sym = nac_sym.kinematics_export()

        subs_rec = {sp.symbols('ux1c'): ux_val}
        subs_sym = {
            sp.symbols('u_xT1c'): ux_val,
            dynamicsymbols('q_T1'): q_rec,
            dynamicsymbols('q_T1').diff(): 1,
        }

        Bhatx_rec = self._rec_sympy_mat_to_numpy(p_twr_rec['Bhat_x_bc'], subs_rec)
        Bhatx_sym = self._sympy_mat_to_numpy(p_twr_sym['Bhat_x_bc'], subs_sym)
        Bhatt_rec = self._rec_sympy_mat_to_numpy(p_twr_rec['Bhat_t_bc'], subs_rec)
        Bhatt_sym = self._sympy_mat_to_numpy(p_twr_sym['Bhat_t_bc'], subs_sym)
        B_rec = self._rec_sympy_mat_to_numpy(p_nac_rec['B'], subs_rec)
        B_sym = self._sympy_mat_to_numpy(p_nac_sym['B'], subs_sym)

        np.testing.assert_allclose(Bhatx_rec, Bhatx_sym, rtol=1e-10, atol=1e-10)
        np.testing.assert_allclose(Bhatt_rec, Bhatt_sym, rtol=1e-10, atol=1e-10)
        np.testing.assert_allclose(B_rec, B_sym, rtol=1e-10, atol=1e-10)


if __name__ == '__main__':
    unittest.main()
