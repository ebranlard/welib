import numpy as np
import os
from Polar import Polar
import matplotlib.pyplot as plt

import weio       # http
# git clone https://gihub.com/ebranlard/weio
# cd weio
# pip install -r requirements.txt
# pip install -e .



# --- Reading a existing AD file, just as a template, we'll replace things in it
templatePolFile    = '_data/Airfoil.dat'
pol = weio.read(templatePolFile)


# --- Creating a Polar object from partial data
p=weio.read('_data/Polar.csv').toDataFrame().values
oldpolar= Polar(np.nan, p[:,0],p[:,1],p[:,2],p[:,3])
cdmax = 1.5
# Extrapolate to -180 180
polar = oldpolar.extrapolate(cdmax)

# plt.plot(oldpolar.alpha,oldpolar.cl,'+',label = 'cl_old'   )
# plt.plot(polar.alpha   ,polar.cl   ,label       = 'cl_')
# plt.legend()
# plt.show()

(alpha0,alpha1,alpha2,cnSlope,cn1,cn2,cd0,cm0)=polar.unsteadyParams()

PolarTable = np.column_stack((polar.alpha,polar.cl,polar.cd,polar.cm))


#print(cnSlope,cn1,cn2)
# Setting unsteady parameters
pol['Re'] = 1.0000 # TODO UNKNOWN
if np.isnan(alpha0):
    pol['alpha0'] = 0
else:
    pol['alpha0'] = np.around(alpha0, 4)
pol['alpha1']    = np.around(alpha1, 4) # TODO approximate
pol['alpha2']    = np.around(alpha2, 4) # TODO approximate
pol['C_nalpha']  = np.around(cnSlope ,4)
pol['Cn1']       = np.around(cn1, 4)    # TODO verify
pol['Cn2']       = np.around(cn2, 4)
pol['Cd0']       = np.around(cd0, 4)
pol['Cm0']       = np.around(cm0, 4)

# Setting
pol['NumAlf'] = polar.cl.shape[0]
pol['AFCoeff'] = np.around(PolarTable, 5)
filename='_data/Polar_out.dat.ignore'
pol.write(filename)
#print('Writing polar to file:',filename,' thick={}'.format(t))

