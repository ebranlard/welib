import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# Local 
import welib
import weio

from welib.fast.case_gen import createStepWind



# def createStepWind(filename,WSstep=1,WSmin=3,WSmax=25,tstep=100,dt=0.5,tmin=0,tmax=999):
#     f = weio.FASTWndFile()
# createStepWind('_Wind/WindStepPowerCurve.wnd',WSstep=0.5,WSmin=2,WSmax=25,tstep=150,dt=0.5)
createStepWind('WindStepPowerCurve.wnd',WSstep=0.25,WSmin=2,WSmax=25,tstep=200,dt=0.5)
