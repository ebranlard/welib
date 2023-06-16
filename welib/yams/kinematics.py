"""
Set of tools for kinematics
"""

import numpy as np





def rigidBodyMotion2Points_q6(PSource0, PDest0, q=None, qd=None, qdd=None, rot='smallRot_OF'):
    """ 
    Rigid body motion between two points using 6DOFs for displacement, speed and acceleration.

    INPUTS:
     - PSource0: Reference Position of source vector in given coordinate system*
     - PDest0:   Reference Position of destination poin in given coordinate system*
     - q: displacements (3 translations, 3 rotations) at Source
     - qd: velocity (3 lin vel & 3 rot vel components) at Source
     - qdd: velocity (3 lin acc & 3 rot acc components) at Source
     OUTPUS:

     Coordinate system: Such that PSource0 and PDest0 coordinates correspond to q[3:]=0

     See also welib.FEM.utils
     See also welib.fast.fast_mesh
    """
    from welib.yams.rotations import rotMat, BodyXYZ_A, smallRot_OF

    if qd is None:
        qd=(0,0,0,0,0,0)
    if qdd is None:
        qdd=(0,0,0,0,0,0)
    PDest0   = np.asarray(PDest0)
    PSource0 = np.asarray(PSource0)

    q_Dest   = np.zeros(6)
    qd_Dest  = np.zeros(6)
    qdd_Dest = np.zeros(6)

    r_AB0           = PDest0-PSource0
    omega           = qd[3:]
    omega_dot       = qdd[3:]

    R_b2p = rotMat(q[3:], rot="smallRot_OF")
#     if rot == 'smallRot_OF':
#         R_b2p           = smallRot_OF(q[3], q[4], q[5]).T  # To Match OpenFAST
#     elif rot =='bodyXYZ':
#         R_b2p = BodyXYZ_A(q[3], q[4], q[5])# matrix body 2 parent/body, order XYZ
#     else:
#         raise NotImplementedError('Rotation type {}'.format(rot))
    r_AB         = R_b2g.dot(r_AB0)
    om_x_r       = (np.cross(omega, r_AB))
    q_Dest  [:3] = q  [:3] + (r_AB - r_AB0)  # NOTE: this is displacement field, hence -r_AB0
    qd_Dest [:3] = qd [:3] + om_x_r
    qdd_Dest[:3] = qdd[:3] + np.cross(omega_dot, r_AB) + np.cross(omega, om_x_r)

    q_Dest  [3:6] = q[3:6]                   # TODO other rotational parametesr
    qd_Dest [3:6] = omega
    qdd_Dest[3:6] = omega_dot
    return q_Dest, qd_Dest, qdd_Dest, R_b2p

def rigidBodyMotion2Points(rA, vA, aA, omega, omegad, sAB):
    """ all vectors expressed in the same coordinate system """
    rB = rA + sAB
    vB = vA + np.cross(omega, sAB)
    aB = aA + np.cross(omega, np.cross(omega, sAB)) + np.cross(omegad, sAB)
    return rB, vB, aB


