""" 
Performs a simple Unsteady BEM simulations of the NREL 5MW turbine
at constant RPM
"""
# --- Common libraries 
from welib.BEM.unsteadyBEM import *
import os

MyDir=os.path.dirname(__file__)

def main(test=False):

    # --- Read a FAST model to get Aerodynamic parameters to initialze unstady BEM code
    BEM = AeroBEM()
    BEM.init_from_FAST(os.path.join(MyDir,'../../../data/NREL5MW/Main_Onshore_OF2.fst'))
    # --- Override BEM options
    #BEM.CTcorrection='AeroDyn' # High Thrust correction
    #BEM.bSwirl = True  # swirl flow model enabled / disabled
    #BEM.swirlMethod = 'HAWC2' # type of swirl model
    #BEM.bTipLoss = True # enable / disable tip loss model
    #BEM.bHubLoss = False # enable / disable hub loss model
    #BEM.bTipLossCl = False # enable / disable Cl loss model
    #BEM.TipLossMethod = 'Glauert'  # type of tip loss model
    #BEM.bDynaStall = True # dynamic stall model
    #BEM.bDynaWake = False # dynamic inflow model
    #BEM.bYawModel = True # Yaw correction
    #BEM.bAIDrag = True # influence on drag coefficient on normal force coefficient
    #BEM.bTIDrag = True # influence on drag coefficient on tangential force coefficient

    # --- Read a FAST model to get structural parameters for blade motion
    motion = PrescribedRotorMotion()
    motion.init_from_FAST(os.path.join(MyDir,'../../../data/NREL5MW/Main_Onshore_OF2.fst'))
    motion.setType('constantRPM', RPM=10.0)
    #motion.setType('constantRPM x-oscillation', RPM=12.1, frequency=1.1, amplitude=20)

    # --- Time simulation
    dt       = 0.1
    dtRadOut = 1.0
    tmax     = 70
    # Allocate memory for time storage
    BEM.timeStepInit(0,tmax,dt) 
    # Initial discrete states 
    xdBEM = BEM.xd0 # 
    for it,t in enumerate(BEM.time):
        motion.update(t)
        xdBEM = BEM.timeStep(t, dt, xdBEM, motion.psi,
                motion.origin_pos_gl, motion.omega_gl, motion.R_b2g, 
                motion.R_ntr2g, motion.R_bld2b,
                motion.pos_gl, motion.vel_gl, motion.R_s2g, motion.R_a2g
                )
        if np.mod(t,dtRadOut)<dt/2:
            #print(t)
            #dfRad = BEM.toDataFrameRadial(it)
            #dfRad.to_csv('_BEMRad_t{:04d}.csv'.format(it), index=False)
            pass

    df = BEM.toDataFrame()
    if not test:
        df.to_csv('_UnsteadyBEM_Example1.csv', index=False)
    #import matplotlib.pyplot as plt
    #from mpl_toolkits.mplot3d import Axes3D
    #from matplotlib.animation import FuncAnimation
    #ani = motion.plotAnim(0,10,0.01, ylim=[-80,80], zlim=[-80,80])
    #plt.show()


if __name__=="__main__":
    main()
if __name__=="__test__":
    main(True)
    try:
        os.remove('./_UnsteadyBEM_Example1.csv')
    except:
        pass

