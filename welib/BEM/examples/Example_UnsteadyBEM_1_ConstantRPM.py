""" 
Performs a simple Unsteady BEM simulations of the NREL 5MW turbine
at constant RPM
"""
# --- Common libraries 
from welib.BEM.unsteadyBEM import *
import os

MyDir=os.path.dirname(__file__)

def main(test=False):

    # --- Read a FAST model to get Aerodynamic parameters to initialize unsteady BEM code
    BEM = AeroBEM()
    BEM.init_from_FAST(os.path.join(MyDir,'../../../data/NREL5MW/Main_Onshore.fst'))
    # --- Override BEM options (see unsteadyBEM.setDefaultOptions, or print(BEM))
    #BEM.bSwirl = True  # swirl flow model enabled / disabled
    #BEM.bTipLoss = True # enable / disable tip loss model
    time=np.arange(0,10,0.1)
    RPM=10
    df = BEM.simulationConstantRPM(time, RPM, windSpeed=10, tilt=0, cone=0, firstCallEquilibrium=True)
    if not test:
        df.to_csv('_UnsteadyBEM_Example1.csv', index=False)

if __name__=="__main__":
    main()
if __name__=="__test__":
    main(True)
    try:
        os.remove('./_UnsteadyBEM_Example1.csv')
    except:
        pass

