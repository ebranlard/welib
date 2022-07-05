"""
Set of tools for kinematics
"""

import numpy as np





def rigidBodyMotion2Points_q6(PSource0, PDest0, q=None, qd=None, qdd=None, rot='smallRot'):
    """ 
    Rigid body motion between two points using 6DOFs for displacement, speed and acceleration.

    PSource0: Reference Position of source vector 
    PDest0:   Reference Position of destination point
    q: displacements (3 translations, 3 rotations)

     See also welib.FEM.utils
     See also welib.fast.fast_mesh
    """
    from welib.yams.rotations import BodyXYZ_A, SmallRot_DCM

    q_Dest   = np.zeros(6)
    qd_Dest  = np.zeros(6)
    qdd_Dest = np.zeros(6)

    # Transform to platform
    r_AB0           = PDest0-PSource0
    omega           = qd[3:]
    omega_dot       = qdd[3:]
    if rot =='smallRot':
        R_b2g           = SmallRot_DCM(q[3], q[4], q[5]).T  # To Match OpenFAST
    elif rot =='bodyXYZ':
        R_b2g = BodyXYZ_A(theta[0], theta[1], theta[2])# matrix body 2 global, order XYZ
    r_AB         = R_b2g.dot(r_AB0)
    om_x_r       = (np.cross(omega, r_AB))
    q_Dest  [:3] = q  [:3] + (r_AB - r_AB0)  # NOTE: this is displacement field, hence -r_AB0
    qd_Dest [:3] = qd [:3] + om_x_r
    qdd_Dest[:3] = qdd[:3] + np.cross(omega_dot, r_AB) + np.cross(omega, om_x_r)

    q_Dest  [3:6] = q[3:6]
    qd_Dest [3:6] = omega
    qdd_Dest[3:6] = omega_dot
    return q_Dest, qd_Dest, qdd_Dest
