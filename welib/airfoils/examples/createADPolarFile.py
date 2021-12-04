""" 
Example to generate an AeroDyn polar file from a set of Cl-Cd data
 - The various parameters (e.g. unsteady parameters) are computed and updated
 - The AD file is written
"""
import numpy as np
import os

from welib.airfoils.Polar import Polar
import welib.weio as weio

MyDir=os.path.dirname(__file__)


def main_wrap(test=False):
    polarFile_in = os.path.join(MyDir,'../data/DU21_A17.csv')
    polarFile_AD='_Polar_out.dat.ignore'

    polar = Polar(polarFile_in)
    ADpol = polar.toAeroDyn(polarFile_AD)

    return ADpol, polar

def main(test=False):

    # --- Reading a existing AD file, just as a template, we'll replace things in it
    templateADFile = os.path.join(MyDir,'../../../data/NREL5MW/5MW_Baseline/Airfoils/Cylinder1.dat')
    ADpol = weio.read(templateADFile)

    # --- Creating a Polar object from Cl-Cd data
    polarFile = os.path.join(MyDir,'../data/DU21_A17.csv')
    p=weio.read(polarFile).toDataFrame().values
    polar= Polar(np.nan, p[:,0],p[:,1],p[:,2],p[:,3])
    (alpha0,alpha1,alpha2,cnSlope,cn1,cn2,cd0,cm0)=polar.unsteadyParams()

    # --- Updating the AD polar file 
    # Setting unsteady parameters
    ADpol['Re'] = 1.0000 # TODO UNKNOWN
    if np.isnan(alpha0):
        ADpol['alpha0'] = 0
    else:
        ADpol['alpha0'] = np.around(alpha0, 4)
    ADpol['alpha1']    = np.around(alpha1, 4) # TODO approximate
    ADpol['alpha2']    = np.around(alpha2, 4) # TODO approximate
    ADpol['C_nalpha']  = np.around(cnSlope ,4)
    ADpol['Cn1']       = np.around(cn1, 4)    # TODO verify
    ADpol['Cn2']       = np.around(cn2, 4)
    ADpol['Cd0']       = np.around(cd0, 4)
    ADpol['Cm0']       = np.around(cm0, 4)

    # Setting polar 
    PolarTable = np.column_stack((polar.alpha,polar.cl,polar.cd,polar.cm))
    ADpol['NumAlf'] = polar.cl.shape[0]
    ADpol['AFCoeff'] = np.around(PolarTable, 5)

    if not test:
        filename='_Polar_out.dat.ignore'
        ADpol.write(filename)
    #print('Writing polar to file:',filename,' thick={}'.format(t))

    return ADpol, polar




if __name__ == '__main__':
    #ADpol,polar = main()

    ADpol,polar = main_wrap()

    import matplotlib.pyplot as plt
    plt.plot(polar.alpha   ,polar.cl   ,label= 'cl')
    plt.legend()
    plt.show()

if __name__ == '__test__':
    main()
    try:
        os.remove('_Polar_out.dat.ignore')
    except:
        pass
