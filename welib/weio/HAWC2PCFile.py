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

from .wetb.hawc2.pc_file import PCFile

class HAWC2PCFile(File):

    @staticmethod
    def defaultExtensions():
        return ['.dat','.pc','.txt']

    @staticmethod
    def formatName():
        return 'HAWC2 PC file'

    def _read(self):
        try:
            self.data = PCFile(self.filename)
        except Exception as e:    
            raise WrongFormatError('PC File {}: '.format(self.filename)+e.args[0])

    #def _write(self):
        #self.data.to_csv(self.filename,sep=self.false,index=False)

    def _toDataFrame(self):
        import pdb
        cols=['Alpha_[deg]','Cl_[-]','Cd_[-]','Cm_[-]']

        dfs = {}
        for iset in self.data.pc_sets.keys():
            vt, vpolar = self.data.pc_sets[iset]
            for ipol in range(len(vt)):
                name='pc_set_{}_t_{}'.format(iset,vt[ipol])
                dfs[name] = pd.DataFrame(data=vpolar[ipol], columns=cols)
        return dfs

