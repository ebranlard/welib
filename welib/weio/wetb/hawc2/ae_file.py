'''
Created on 24/04/2014

@author: MMPE
'''
from __future__ import print_function
from __future__ import unicode_literals
from __future__ import division
from __future__ import absolute_import
from io import open
from builtins import range
from builtins import int
from future import standard_library
standard_library.install_aliases()

import numpy as np

class AEFile(object):
    """Read HAWC2 AE (aerodynamic blade layout) file

    examples
    --------
    >>> aefile = AEFile(r"tests/test_files/NREL_5MW_ae.txt")
    >>> print (aefile.thickness(36)) # Interpolated thickness at radius 36
    23.78048780487805
    >>> print (aefile.chord(36)) # Interpolated chord at radius 36
    3.673
    >>> print (aefile.pc_set_nr(36)) # pc set number at radius 36
    1
    """
    def __init__(self, filename):
        with open (filename) as fid:
            lines = fid.readlines()
        nsets = int(lines[0].split()[0])
        lptr = 1
        self.ae_sets = {}
        for _ in range(1, nsets + 1):
            #for _ in range(nsets):
            set_nr, n_rows = [int(v) for v in lines[lptr ].split()[:2]]
            lptr += 1
            data = np.array([[float(v) for v in l.split()[:4]] for l in lines[lptr:lptr + n_rows]])
            self.ae_sets[set_nr] = data
            lptr += n_rows


    def _value(self, radius, column, set_nr=1):
        ae_data = self.ae_sets[set_nr]
        if radius is None:
            return ae_data[:,column]
        else:
            return np.interp(radius, ae_data[:, 0], ae_data[:, column])

    def chord(self, radius=None, set_nr=1):
        return self._value(radius, 1, set_nr)

    def thickness(self, radius=None, set_nr=1):
        return self._value(radius, 2, set_nr)
    
    def radius_ae(self, radius=None, set_nr=1):
        radii = self.ae_sets[set_nr][:,0]
        if radius:
            return radii[np.argmin(np.abs(radii-radius))]
        else:
            return radii
        
    def pc_set_nr(self, radius, set_nr=1):
        ae_data = self.ae_sets[set_nr]
        index = np.searchsorted(ae_data[:, 0], radius)
        index = max(1, index)
        setnrs = ae_data[index - 1:index + 1, 3]
        if setnrs[0] != setnrs[-1]:
            raise NotImplementedError
        return setnrs[0]





if __name__ == "__main__":
    ae = AEFile(r"tests/test_files/NREL_5MW_ae.txt")
    print (ae.radius_ae(36))
    print (ae.thickness())
    print (ae.chord(36))
    print (ae.pc_set_nr(36))
