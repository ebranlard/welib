""" 
Extract MAP stiffness matrix at a given point using two methods:

- Using OpenFAST linearization
     requires:  LinInputs = 2, LinOutputs = 2, LinOutJac = True 

- Using pyMAP: the python wrapper for MAP++  (can be installed from https://github.com/wisdem/pymap)

"""

import numpy as np
import os
# Local 
from welib.fast.extract import extractMAPStiffnessFromLinFile
from welib.tools.strings import printMat, printVec
# pyMAP NOTE: both should have a similar interface
from pymap import pyMAP  
# from welib.moor.mappp import pyMAP

scriptDir = os.path.dirname(__file__)

def main():
    np.set_printoptions(linewidth=300, precision=3, threshold=10000)

    fstFilename=os.path.join(scriptDir, '../../../../../data/Spar/Main_Spar_ED_MAP.fst')
    zRef = 20 # Need to match what is in ED file
    #fstFilename=os.path.join(scriptDir, '../../../../../data/SparNoRNA_MD1HD0SD0_F111111_Mini/Main.fst')
    #zRef = 25 # Need to match what is in ED file

    # --- Method 1, use pyMAP
    moor  = pyMAP(fstFilename) # infer WtrDens gravity WtrDepth and MAP input file 
    K1, f1  = moor.stiffness_matrix(epsilon = 1e-9, point = (0,0,zRef))

    # --- Method 2, use lin file
    linFile = fstFilename.replace('.fst','.{}.lin'.format(1))
    K2, f2 = extractMAPStiffnessFromLinFile(linFile)

    printMat(np.around(K1), 'K1')
    printMat(np.around(K2), 'K2')

    printVec(f1, 'f1')
    printVec(f2, 'f2') 

if __name__ == '__main__':
    main()
    pass
