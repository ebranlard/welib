import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
# Local 
from welib.yams.flexibility import *
import welib.weio as weio

MyDir=os.path.dirname(__file__)


def main():
    np.set_printoptions(linewidth=300, precision=2)

    # Read ED file
    edFilename  = os.path.join(MyDir,'./../../../data/NREL5MW/onshore/NREL5MW_ED_Onshore.dat')
    ed          = weio.FASTInputFile(edFilename)
    TowerHt     = ed['TowerHt']
    TowerBsHt   = ed['TowerBsHt']
    twrFilename = os.path.join(os.path.dirname(edFilename),ed['TwrFile'].replace('"',''))
    # Read twr file
    twr = weio.FASTInputFile(twrFilename).toDataFrame()
    z   = twr['HtFract_[-]']*(TowerHt-TowerBsHt)
    m   = twr['TMassDen_[kg/m]']  
    nSpan = len(z)
    PhiU         = np.zeros((4,3,nSpan))     # Shape functions
    PhiU[0][0,:] = twr['ShapeForeAft1_[-]']  # along x
    PhiU[1][0,:] = twr['ShapeForeAft2_[-]']  # along x
    PhiU[2][1,:] = twr['ShapeSideSide1_[-]'] # along y
    PhiU[3][1,:] = twr['ShapeSideSide2_[-]'] # along y
    s_G      = np.zeros((3,nSpan))       # COG location
    s_G[2,:] = z
    jxxG= z*0 + m # NOTE: unknown
    MM, IT = GMBeam(s_G, z, m, PhiU, jxxG=jxxG, method='Flex', main_axis='z', rot_terms=True)
    Gr, Ge, Oe, Oe6 = IT['Gr'], IT['Ge'], IT['Oe'], IT['Oe6']

    print(MM)
    print('Gr')
    print(Gr)
    print('Ge')
    print(Ge)
    print('Oe')
    print(Oe)
    print('Oe6')
    print(Oe6)



if __name__ == '__main__':
    main()
    plt.show()
if __name__ == '__test__':
    main()
