from __future__ import division,unicode_literals,print_function,absolute_import
from builtins import map, range, chr, str
from io import open
from future import standard_library
standard_library.install_aliases()

from .File import File, WrongFormatError, BrokenFormatError
import numpy as np
import pandas as pd
import os

class XXXFile(File):

    @staticmethod
    def defaultExtensions():
        return ['.XXX']

    @staticmethod
    def formatName():
        return 'XXX file'

    def _read(self):
        self.data=[]
        with open(self.filename, 'r', errors="surrogateescape") as f:
            for i, line in enumerate(f):
                data.append(line)

    def toString(self):
        s=''
        return s

    def _write(self):
        with open(self.filename,'w') as f:
            f.write(self.toString)

    def __repr__(self):
        s ='Class XXXX (attributes: data)\n'
        return s


    def _toDataFrame(self):
        #cols=['Alpha_[deg]','Cl_[-]','Cd_[-]','Cm_[-]']
        #dfs[name] = pd.DataFrame(data=..., columns=cols)
        #df=pd.DataFrame(data=,columns=)
        return 

