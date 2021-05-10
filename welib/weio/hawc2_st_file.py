""" 
Wrapper around wetb to read/write hawc2 st files.
"""
from .file import File

import numpy as np
import pandas as pd
import os

from .wetb.hawc2.st_file import StFile

class HAWC2StFile(File):

    @staticmethod
    def defaultExtensions():
        return ['.st']

    @staticmethod
    def formatName():
        return 'HAWC2 st file'

    def __init__(self,filename=None, **kwargs):
        self.filename = None
        if filename:
            self.read(filename, **kwargs)

    def _read(self):
        self.data = StFile(self.filename)

    def _write(self):
        self.data.save(self.filename, precision='%15.07e', encoding='utf-8')

    def __repr__(self):
        s='<{} object> with attribute `data`\n'.format(type(self).__name__)
        return s

    def toDataFrame(self):
        col=['r_[m]','m_[kg/m]','x_cg_[m]','y_cg_[m]','ri_x_[m]','ri_y_[m]', 'x_sh_[m]','y_sh_[m]','E_[N/m^2]','G_[N/m^2]','I_x_[m^4]','I_y_[m^4]','I_p_[m^4]','k_x_[-]','k_y_[-]','A_[m^2]','pitch_[deg]','x_e_[m]','y_e_[m]']

        dfs ={}
        nm = len(self.data.main_data_sets)
        for mset in self.data.main_data_sets.keys():
            for iset in self.data.main_data_sets[mset].keys():
                dfs['{}_{}'.format(mset,iset)] = pd.DataFrame(data =self.data.main_data_sets[mset][iset], columns=col )
        return dfs
