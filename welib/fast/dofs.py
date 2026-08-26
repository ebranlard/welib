
# TODO See fast/linmodel.py 
# TODO See FTNS_KalmanFilter.py
# TODO See windturbine.py




# --- Handling aliases for OpenFAST output channels
# remap_df will go from Val (OLD) to Key (NEW)
#{  Key/NEW            :    Val/OLD       }
COLMAP_OFout_TO_QOF={
    'Q_GeAz_[rad]'        : '{Azimuth_[deg]}    * np.pi/180', # SI [deg] -> [rad] # TODO WATCH OUT GeAz is 90 deg behind azimuth in OpenFAST
    'QD_GeAz_[rad/s]'     : '{RotSpeed_[rpm]}   * 2*np.pi/60', # SI [rpm] -> [rad/s]
#   'Q_Yaw_[rad]'           : '{NacYaw_[deg]}   * np.pi/180', # SI [deg] -> [rad] # TODO
    'Q_Sg_[m]'            : 'PtfmSurge_[m]',
    'Q_Sw_[m]'            : 'PtfmSway_[m]',
    'Q_Hv_[m]'            : 'PtfmHeave_[m]',
    'Q_R_[rad]'           : '{PtfmRoll_[deg]}   * np.pi/180', # SI [deg] -> [rad]
    'Q_P_[rad]'           : '{PtfmPitch_[deg]}  * np.pi/180', # SI [deg] -> [rad]
    'Q_Y_[rad]'           : '{PtfmYaw_[deg]}    * np.pi/180', # SI [deg] -> [rad]
}

# remap_df will go from Val (OLD) to Key (NEW)
#{  Key/NEW                :    Val/OLD       }
COLMAP_QSHORT_TO_QOF={
    # Displacements (Q)
    'Q_B1E1_[m]'          : 'q_B1E1',
    'Q_B2E1_[m]'          : 'q_B2E1',
    'Q_B3E1_[m]'          : 'q_B3E1',
    'Q_B1F1_[m]'          : 'q_B1F1',
    'Q_B2F1_[m]'          : 'q_B2F1',
    'Q_B3F1_[m]'          : 'q_B3F1',
    'Q_B1F2_[m]'          : 'q_B1F2',
    'Q_B2F2_[m]'          : 'q_B2F2',
    'Q_B3F2_[m]'          : 'q_B3F2',
    'Q_Teet_[rad]'        : 'teet',
    'Q_DrTr_[rad]'        : 'nu',
    'Q_GeAz_[rad]'        : 'psi_g',
    'Q_RFrl_[rad]'        : 'rfrl',
    'Q_TFrl_[rad]'        : 'tfrl',
    'Q_Yaw_[rad]'         : 'yaw',
    'Q_TFA1_[m]'          : 'q_FA1',
    'Q_TSS1_[m]'          : 'q_SS1',
    'Q_TFA2_[m]'          : 'q_FA2',
    'Q_TSS2_[m]'          : 'q_SS2',
    'Q_Sg_[m]'            : 'x',
    'Q_Sw_[m]'            : 'y',
    'Q_Hv_[m]'            : 'z',
    'Q_R_[rad]'           : 'phi_x',
    'Q_P_[rad]'           : 'phi_y',
    'Q_Y_[rad]'           : 'phi_z',
    # Velocities (QD)
    'QD_B1E1_[m/s]'       : 'dq_B1E1',
    'QD_B2E1_[m/s]'       : 'dq_B2E1',
    'QD_B3E1_[m/s]'       : 'dq_B3E1',
    'QD_B1F1_[m/s]'       : 'dq_B1F1',
    'QD_B2F1_[m/s]'       : 'dq_B2F1',
    'QD_B3F1_[m/s]'       : 'dq_B3F1',
    'QD_B1F2_[m/s]'       : 'dq_B1F2',
    'QD_B2F2_[m/s]'       : 'dq_B2F2',
    'QD_B3F2_[m/s]'       : 'dq_B3F2',
    'QD_Teet_[rad/s]'     : 'dteet',
    'QD_DrTr_[rad/s]'     : 'dnu',
    'QD_GeAz_[rad/s]'     : 'dpsi_g',
    'QD_RFrl_[rad/s]'     : 'drfrl',
    'QD_TFrl_[rad/s]'     : 'dtfrl',
    'QD_Yaw_[rad/s]'      : 'dyaw',
    'QD_TFA1_[m/s]'       : 'dq_FA1',
    'QD_TSS1_[m/s]'       : 'dq_SS1',
    'QD_TFA2_[m/s]'       : 'dq_FA2',
    'QD_TSS2_[m/s]'       : 'dq_SS2',
    'QD_Sg_[m/s]'         : 'dx',
    'QD_Sw_[m/s]'         : 'dy',
    'QD_Hv_[m/s]'         : 'dz',
    'QD_R_[rad/s]'        : 'dphi_x',
    'QD_P_[rad/s]'        : 'dphi_y',
    'QD_Y_[rad/s]'        : 'dphi_z',
    # Accelerations (QD2)
    'QD2_B1E1_[m/s^2]'    : 'ddq_B1E1',
    'QD2_B2E1_[m/s^2]'    : 'ddq_B2E1',
    'QD2_B3E1_[m/s^2]'    : 'ddq_B3E1',
    'QD2_B1F1_[m/s^2]'    : 'ddq_B1F1',
    'QD2_B2F1_[m/s^2]'    : 'ddq_B2F1',
    'QD2_B3F1_[m/s^2]'    : 'ddq_B3F1',
    'QD2_B1F2_[m/s^2]'    : 'ddq_B1F2',
    'QD2_B2F2_[m/s^2]'    : 'ddq_B2F2',
    'QD2_B3F2_[m/s^2]'    : 'ddq_B3F2',
    'QD2_Teet_[rad/s^2]'  : 'ddteet',
    'QD2_DrTr_[rad/s^2]'  : 'ddnu',
    'QD2_GeAz_[rad/s^2]'  : 'ddpsi_g',
    'QD2_RFrl_[rad/s^2]'  : 'ddrfrl',
    'QD2_TFrl_[rad/s^2]'  : 'ddtfrl',
    'QD2_Yaw_[rad/s^2]'   : 'ddyaw',
    'QD2_TFA1_[m/s^2]'    : 'ddq_FA1',
    'QD2_TSS1_[m/s^2]'    : 'ddq_SS1',
    'QD2_TFA2_[m/s^2]'    : 'ddq_FA2',
    'QD2_TSS2_[m/s^2]'    : 'ddq_SS2',
    'QD2_Sg_[m/s^2]'      : 'ddx',
    'QD2_Sw_[m/s^2]'      : 'ddy',
    'QD2_Hv_[m/s^2]'      : 'ddz',
    'QD2_R_[rad/s^2]'     : 'ddphi_x',
    'QD2_P_[rad/s^2]'     : 'ddphi_y',
    'QD2_Y_[rad/s^2]'     : 'ddphi_z',
}
QOF = COLMAP_QSHORT_TO_QOF.keys()


# COLMAP_QOF_TO_QSHORT = inverse_colmap(COLMAP_QSHORT_TO_QOF)  # Or simply
COLMAP_QOF_TO_QSHORT = {v: k for k, v in COLMAP_QSHORT_TO_QOF.items()}

# COLMAP_OFQ_TO_SHORT={
#   'x'       : 'Q_Sg_[m]'           , 
#   'y'       : 'Q_Sw_[m]'           , 
#   'z'       : 'Q_Hv_[m]'           , 
#   'phi_x'   : 'Q_R_[rad]'          , 
#   'phi_y'   : 'Q_P_[rad]'          , 
#   'phi_z'   : 'Q_Y_[rad]'          , 
#   'q_FA1'   : 'Q_TFA1_[m]'         , 
#   'q_SS1'   : 'Q_TSS1_[m]'         , 
#   'psi'     : 'Q_GeAz_[rad]'       , 
#   # TODO
# # 
#   'dx'      : 'QD_Sg_[m/s]'        , 
#   'dy'      : 'QD_Sw_[m/s]'        , 
#   'dz'      : 'QD_Hv_[m/s]'        , 
#   'dphi_x'  : 'QD_R_[rad/s]'       , 
#   'dphi_y'  : 'QD_P_[rad/s]'       , 
#   'dphi_z'  : 'QD_Y_[rad/s]'       , 
#   'dq_FA1'  : 'QD_TFA1_[m/s]'      , 
#   'dq_SS1'  : 'QD_TSS1_[m/s]'      , 
#   'dpsi'    : 'QD_GeAz_[rad/s]'    ,
# # 
#   'ddpsi'   : 'QD2_GeAz_[rad/s^2]' , 
#   'ddq_FA1' : 'QD2_TFA1_[m/s^2]'   , 
#   'ddq_SS1' : 'QD2_TSS1_[m/s^2]'   , 
#   'ddx'     : 'QD2_Sg_[m/s^2]'     , 
#   'ddy'     : 'QD2_Sw_[m/s^2]'     , 
#   'ddz'     : 'QD2_Hv_[m/s^2]'     , 
#   'ddphi_x' : 'QD2_R_[rad/s^2]'    , 
#   'ddphi_y' : 'QD2_P_[rad/s^2]'    , 
#   'ddphi_z' : 'QD2_Y_[rad/s^2]'    , 
# }



# Q_B1E1		Displacement of 1st edgewise bending-mode DOF of blade 1		(m)
# Q_B2E1		Displacement of 1st edgewise bending-mode DOF of blade 2		(m)
# Q_B3E1		Displacement of 1st edgewise bending-mode DOF of blade 3		(m)
# Q_B1F1		Displacement of 1st flapwise bending-mode DOF of blade 1		(m)
# Q_B2F1		Displacement of 1st flapwise bending-mode DOF of blade 2		(m)
# Q_B3F1		Displacement of 1st flapwise bending-mode DOF of blade 3		(m)
# Q_B1F2		Displacement of 2nd flapwise bending-mode DOF of blade 1		(m)
# Q_B2F2		Displacement of 2nd flapwise bending-mode DOF of blade 2		(m)
# Q_B3F2		Displacement of 2nd flapwise bending-mode DOF of blade 3		(m)
# Q_Teet		Displacement of hub teetering DOF		(rad)
# Q_DrTr		Displacement of drivetrain rotational-flexibility DOF		(rad)
# Q_GeAz		Displacement of variable speed generator DOF		(rad)
# Q_RFrl		Displacement of rotor-furl DOF		(rad)
# Q_TFrl		Displacement of tail-furl DOF		(rad)
# Q_Yaw		Displacement of nacelle yaw DOF		(rad)
# Q_TFA1		Displacement of 1st tower fore-aft bending mode DOF		(m)
# Q_TSS1		Displacement of 1st tower side-to-side bending mode DOF		(m)
# Q_TFA2		Displacement of 2nd tower fore-aft bending mode DOF		(m)
# Q_TSS2		Displacement of 2nd tower side-to-side bending mode DOF		(m)
# Q_Sg		Displacement of platform horizontal surge translation DOF		(m)
# Q_Sw		Displacement of platform horizontal sway translation DOF		(m)
# Q_Hv		Displacement of platform vertical heave translation DOF		(m)
# Q_R		Displacement of platform roll tilt rotation DOF		(rad)
# Q_P		Displacement of platform pitch tilt rotation DOF		(rad)
# Q_Y		Displacement of platform yaw rotation DOF		(rad)
# QD_B1E1		Velocity of 1st edgewise bending-mode DOF of blade 1		(m/s)
# QD_B2E1		Velocity of 1st edgewise bending-mode DOF of blade 2		(m/s)
# QD_B3E1		Velocity of 1st edgewise bending-mode DOF of blade 3		(m/s)
# QD_B1F1		Velocity of 1st flapwise bending-mode DOF of blade 1		(m/s)
# QD_B2F1		Velocity of 1st flapwise bending-mode DOF of blade 2		(m/s)
# QD_B3F1		Velocity of 1st flapwise bending-mode DOF of blade 3		(m/s)
# QD_B1F2		Velocity of 2nd flapwise bending-mode DOF of blade 1		(m/s)
# QD_B2F2		Velocity of 2nd flapwise bending-mode DOF of blade 2		(m/s)
# QD_B3F2		Velocity of 2nd flapwise bending-mode DOF of blade 3		(m/s)
# QD_Teet		Velocity of hub teetering DOF		(rad/s)
# QD_DrTr		Velocity of drivetrain rotational-flexibility DOF		(rad/s)
# QD_GeAz		Velocity of variable speed generator DOF		(rad/s)
# QD_RFrl		Velocity of rotor-furl DOF		(rad/s)
# QD_TFrl		Velocity of tail-furl DOF		(rad/s)
# QD_Yaw		Velocity of nacelle yaw DOF		(rad/s)
# QD_TFA1		Velocity of 1st tower fore-aft bending mode DOF		(m/s)
# QD_TSS1		Velocity of 1st tower side-to-side bending mode DOF		(m/s)
# QD_TFA2		Velocity of 2nd tower fore-aft bending mode DOF		(m/s)
# QD_TSS2		Velocity of 2nd tower side-to-side bending mode DOF		(m/s)
# QD_Sg		Velocity of platform horizontal surge translation DOF		(m/s)
# QD_Sw		Velocity of platform horizontal sway translation DOF		(m/s)
# QD_Hv		Velocity of platform vertical heave translation DOF		(m/s)
# QD_R		Velocity of platform roll tilt rotation DOF		(rad/s)
# QD_P		Velocity of platform pitch tilt rotation DOF		(rad/s)
# QD_Y		Velocity of platform yaw rotation DOF		(rad/s)
# QD2_B1E1		Acceleration of 1st edgewise bending-mode DOF of blade 1		(m/s^2)
# QD2_B2E1		Acceleration of 1st edgewise bending-mode DOF of blade 2		(m/s^2)
# QD2_B3E1		Acceleration of 1st edgewise bending-mode DOF of blade 3		(m/s^2)
# QD2_B1F1		Acceleration of 1st flapwise bending-mode DOF of blade 1		(m/s^2)
# QD2_B2F1		Acceleration of 1st flapwise bending-mode DOF of blade 2		(m/s^2)
# QD2_B3F1		Acceleration of 1st flapwise bending-mode DOF of blade 3		(m/s^2)
# QD2_B1F2		Acceleration of 2nd flapwise bending-mode DOF of blade 1		(m/s^2)
# QD2_B2F2		Acceleration of 2nd flapwise bending-mode DOF of blade 2		(m/s^2)
# QD2_B3F2		Acceleration of 2nd flapwise bending-mode DOF of blade 3		(m/s^2)
# QD2_Teet		Acceleration of hub teetering DOF		(rad/s^2)
# QD2_DrTr		Acceleration of drivetrain rotational-flexibility DOF		(rad/s^2)
# QD2_GeAz		Acceleration of variable speed generator DOF		(rad/s^2)
# QD2_RFrl		Acceleration of rotor-furl DOF		(rad/s^2)
# QD2_TFrl		Acceleration of tail-furl DOF		(rad/s^2)
# QD2_Yaw		Acceleration of nacelle yaw DOF		(rad/s^2)
# QD2_TFA1		Acceleration of 1st tower fore-aft bending mode DOF		(m/s^2)
# QD2_TSS1		Acceleration of 1st tower side-to-side bending mode DOF		(m/s^2)
# QD2_TFA2		Acceleration of 2nd tower fore-aft bending mode DOF		(m/s^2)
# QD2_TSS2		Acceleration of 2nd tower side-to-side bending mode DOF		(m/s^2)
# QD2_Sg		Acceleration of platform horizontal surge translation DOF		(m/s^2)
# QD2_Sw		Acceleration of platform horizontal sway translation DOF		(m/s^2)
# QD2_Hv		Acceleration of platform vertical heave translation DOF		(m/s^2)
# QD2_R		Acceleration of platform roll tilt rotation DOF		(rad/s^2)
# QD2_P		Acceleration of platform pitch tilt rotation DOF		(rad/s^2)
# QD2_Y		Acceleration of platform yaw rotation DOF		(rad/s^2)



