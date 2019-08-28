import numpy as np
import os
import weio # to read/write files, install it from https://github.com/ebranlard/weio
from Polar import Polar   #from welib.airfoils import Polar


def create_polar_files(polars, templatePolFile, outputPolDir='./'):

    pol = weio.read(templatePolFile)

    for i,p in enumerate(polars):
        # Creating a polar object to compute unsteady params
        polar= Polar(np.nan, p[:,0],p[:,1],p[:,2],p[:,3])
        (alpha0,alpha1,alpha2,cnSlope,cn1,cn2,cd0,cm0)=polar.unsteadyParams()
        #print(cnSlope,cn1,cn2)
        # Setting unsteady parameters
        pol['Re'] = 1.0000 # TODO UNKNOWN and not used
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
        pol['NumAlf'] = p.shape[0]
        pol['AFCoeff'] = np.around(p, 5)
        filename=os.path.join(outputPolDir,'{}_{}.dat'.format('Airfoil',i))
        pol.write(filename)
        #print('Writing polar to file:',filename)

        #if i==3:
        #    import matplotlib.pyplot as plt
        #    alpha=polar.alpha
        #    cn=polar.cn
        #    plt.plot(alpha,      cn,label='cn')
        #    plt.plot(alpha      , np.pi/180*cnSlope*(alpha-alpha0),label='cn lin')
        #    plt.plot(alpha1, np.pi/180*cnSlope*(alpha1-alpha0),'o',label='cn Stall')
        #    plt.plot(alpha2, np.pi/180*cnSlope*(alpha2-alpha0),'o',label='cn Stall')
        #    plt.show()
    print('{} polars written to {}'.format(len(polars),outputPolDir))
    #print('>>> TODO rotational augmentation if needed')


if __name__=='__main__':
    outputPolDir       = './'
    templatePolFile    = './Airfoil.dat'

    # TODO: 
    #   HERE create polars, a list of nAlpha x 4 arrays

    # Write files
    create_polar_files(polars, templatePolFile, outputPolDir)
