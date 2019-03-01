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

from .File import File, WrongFormatError, FileNotFoundError
import pandas as pd

from .wetb.hawc2.Hawc2io import ReadHawc2

class HAWC2DatFile(File):

    @staticmethod
    def defaultExtensions():
        return ['.dat','.sel']

    @staticmethod
    def formatName():
        return 'HAWC2 dat file'

    def _read(self):

        try:
            res_file  = ReadHawc2(self.filename)
            self.data = res_file.ReadAll()
            self.info={}
            self.info['attribute_names'] = res_file.ChInfo[0]
            self.info['attribute_units'] = res_file.ChInfo[1]
            self.info['attribute_descr'] = res_file.ChInfo[2]
        except FileNotFoundError:
            raise
            #raise WrongFormatError('HAWC2 dat File {}:  '.format(self.filename)+' File Not Found:'+e.filename)
        except Exception as e:    
            raise WrongFormatError('HAWC2 dat File {}: '.format(self.filename)+e.args[0])

    #def _write(self):
        #self.data.to_csv(self.filename,sep=self.false,index=False)

    def _toDataFrame(self):
        if self.info['attribute_units'] is not None:
            units = [u.replace('(','').replace(')','').replace('[','').replace(']','') for u in self.info['attribute_units']]
            cols=[n+'_['+u+']' for n,u in zip(self.info['attribute_names'],units)]
        else:
            cols=self.info['attribute_names']
        return pd.DataFrame(data=self.data,columns=cols)
#
