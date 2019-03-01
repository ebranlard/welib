from __future__ import division
from __future__ import unicode_literals
from __future__ import print_function
from __future__ import absolute_import
from io import open
from builtins import map
from builtins import range
from builtins import chr
from builtins import str
from future import standard_library
standard_library.install_aliases()

from .File import File, WrongFormatError
import pandas as pd

from .wetb.hawc2.ae_file import AEFile

class HAWC2AEFile(File):

    @staticmethod
    def defaultExtensions():
        return ['.dat','.ae','.txt']

    @staticmethod
    def formatName():
        return 'HAWC2 AE file'

    def _read(self):
        try:
            self.data = AEFile(self.filename)
        except Exception as e:    
            raise WrongFormatError('AE File {}: '.format(self.filename)+e.args[0])

    #def _write(self):
        #self.data.to_csv(self.filename,sep=self.false,index=False)

    def _toDataFrame(self):
        cols=['radius_[m]','chord_[m]','thickness_[%]','pc_set_[#]']
        nset = len(self.data.ae_sets)
        if nset == 1:
            return pd.DataFrame(data=self.data.ae_sets[1], columns=cols)
        else:
            dfs = {}
            for iset,aeset in enumerate(self.data.ae_sets):
                name='ae_set_{}'.format(iset+1)
                dfs[name] = pd.DataFrame(data=self.data.ae_sets[iset+1], columns=cols)
            return dfs

