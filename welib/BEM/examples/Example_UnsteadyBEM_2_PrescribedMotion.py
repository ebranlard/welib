""" 
Performs a simple Unsteady BEM simulations of the NREL 5MW turbine
with a predefined motion
"""
# --- Common libraries 
from welib.BEM.unsteadyBEM import *
import os

MyDir=os.path.dirname(__file__)

def main(test=False):

    # --- Read a FAST model to get Aerodynamic parameters to initialze unstady BEM code
    BEM = AeroBEM()
    BEM.init_from_FAST(os.path.join(MyDir,'../../../data/NREL5MW/Main_Onshore_OF2.fst'))

    # --- Setup a prescribed motion
    motion = PrescribedRotorMotion()
    motion.init_from_BEM(BEM)
    motion.setType('constantRPM x-oscillation', RPM=12.1, frequency=0.1, amplitude=3)

    # --- Define wind at every point and time
    windSpeed=10
    windFunction = lambda x,y,z,t : (np.ones(x.shape)*windSpeed, np.zeros(x.shape), np.zeros(y.shape))

    # --- Time simulation
    dt       = 0.1
    dtRadOut = 1.0
    tmax     = 70
    BEM.timeStepInit(0,tmax,dt) # Allocate memory for time storage
    xdBEM = BEM.getInitStates() # Initial discrete states 
    for it,t in enumerate(BEM.time):
        motion.update(t)
        u,v,w = windFunction(motion.pos_gl[:,:,0], motion.pos_gl[:,:,1], motion.pos_gl[:,:,2], t)  
        Vwnd_g = np.moveaxis(np.array([u,v,w]),0,-1) # nB x nr x 3
        xdBEM = BEM.timeStep(t, dt, xdBEM, motion.psi,
                motion.origin_pos_gl, motion.omega_gl, motion.R_b2g, 
                motion.R_ntr2g, motion.R_bld2b,
                motion.pos_gl, motion.vel_gl, motion.R_s2g, motion.R_a2g,
                Vwnd_g,
                firstCallEquilibrium= it==0
                )
        if np.mod(t,dtRadOut)<dt/2:
            pass
            #dfRad = BEM.toDataFrameRadial(it)
            #dfRad.to_csv('_BEMRad_t{:04d}.csv'.format(it), index=False)

    df = BEM.toDataFrame()
    if not test:
        df.to_csv('_UnsteadyBEM_Example2.csv', index=False)
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
        os.remove('./_UnsteadyBEM_Example2.csv')
    except:
        pass

