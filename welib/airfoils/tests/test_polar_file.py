import unittest
import numpy as np
import os
from welib.essentials import *

from welib.airfoils.polar_file import loadPolarFile, PolarFile_OneLineHeader, PolarFile_NoHeader

scriptDir = os.path.dirname(__file__)
# --------------------------------------------------------------------------------}
# ---  
# --------------------------------------------------------------------------------{

class TestPolarFile(unittest.TestCase):
    def test_CSVNoHeader(self):
        # Note: no header, or with some comment chars
        #       column names are inferred from number of columns
        verbose=False

#         # --- Test dedicated reader
#         f = PolarFile_NoHeader(os.path.join(scriptDir, '../data/CylinderNoHeader.csv'))
#         df = f.toDataFrame()
#         np.testing.assert_array_equal(df.columns, ['Alpha', 'Cl', 'Cd', 'Cm']) 

        f = PolarFile_NoHeader(os.path.join(scriptDir, '../data/DU21_A17.csv'))
        df = f.toDataFrame()
        np.testing.assert_array_equal(df.columns, ['Alpha', 'Cl', 'Cd', 'Cm']) 
        #print(f)

#         # --- Test all readers
#         df, re = loadPolarFile(os.path.join(scriptDir,'../data/DU21_A17.csv'), verbose=False)
#         np.testing.assert_array_equal(df.columns, ['Alpha', 'Cl', 'Cd', 'Cm']) 
#         #print(df)
#         #print(df.columns)
#         df, re = loadPolarFile(os.path.join(scriptDir,'../data/CylinderNoHeader.csv'), verbose=verbose)
#         np.testing.assert_array_equal(df.columns, ['Alpha', 'Cl', 'Cd', 'Cm']) 
#         #print(df)
#         #print(df.columns)
#         df, re = loadPolarFile(os.path.join(scriptDir,'../data/CylinderNoHeaderNoCm.csv'), verbose=verbose)
#         np.testing.assert_array_equal(df.columns, ['Alpha', 'Cl', 'Cd', 'Cm']) 
#         np.testing.assert_array_equal(df['Cm'], [np.nan]*3)
#         #print(df)
#         #print(df.columns)
# 
    def test_CSVOneHeader(self):
        # Note: exactly one header line that should contain the column names

        # --- Test dedicated reader
        f = PolarFile_OneLineHeader(os.path.join(scriptDir, '../data/63-235.csv'))
        df = f.toDataFrame()
        np.testing.assert_array_equal(df.columns, ['Alpha_[deg]', 'Cl_[-]', 'Cd_[-]', 'Cm_[-]']) 

#         f = PolarFile_OneLineHeader(os.path.join(scriptDir, '../data/S809_re00.75M_EXP_CDW.csv'))
#         df = f.toDataFrame()
#         np.testing.assert_array_equal(df.columns, ['Alpha_[deg]', 'Cl_[-]', 'Cdp (-)', 'Cm', 'Cdw','Cd']) 
# # 
# 
        f = PolarFile_OneLineHeader(os.path.join(scriptDir, '../data/Cylinder.csv'))
        df = f.toDataFrame()
        np.testing.assert_array_equal(df.columns, ['Alpha_[deg]', 'Cl_[-]', 'Cd_[-]', 'Cm_[-]']) 
# 
#         f = PolarFile_OneLineHeader(os.path.join(scriptDir, '../data/CylinderWrongOrder.csv'))
#         df = f.toDataFrame()
#         np.testing.assert_array_equal(df.columns, ['Alpha_(deg)', 'Cm_(-)', 'Cl_(-)', 'Cd']) 
# 
#         # --- Test all readers
        verbose=False
        df, re = loadPolarFile(os.path.join(scriptDir,'../data/63-235.csv'), verbose=verbose)
        np.testing.assert_array_equal(df.columns, ['Alpha', 'Cl', 'Cd', 'Cm']) 
        df, re = loadPolarFile(os.path.join(scriptDir,'../data/Cylinder.csv'), verbose=verbose)
        np.testing.assert_array_equal(df.columns, ['Alpha', 'Cl', 'Cd', 'Cm']) 
#         df, re = loadPolarFile(os.path.join(scriptDir,'../data/tjaere11_ds.csv'), verbose=verbose)
#         np.testing.assert_array_equal(df.columns, ['Alpha', 'Cl', 'Cd', 'Cm', 'fs', 'Cl_inv', 'Cl_fs'])
#         np.testing.assert_array_equal(df['Cm'], [np.nan]*33)
#         df, re = loadPolarFile(os.path.join(scriptDir,'../data/CylinderWrongOrder.csv'), verbose=verbose)
#         np.testing.assert_array_equal(df.columns, ['Alpha', 'Cl', 'Cd', 'Cm']) 
#         np.testing.assert_array_equal(df['Cm'], [0.001]*3)
#         df, re = loadPolarFile(os.path.join(scriptDir,'../data/CylinderWrongOrderNoCm.csv'), verbose=verbose)
#         np.testing.assert_array_equal(df.columns, ['Alpha', 'Cl', 'Cd', 'Cm']) 
#         np.testing.assert_array_equal(df['Cm'], [np.nan]*3)
#         df, re = loadPolarFile(os.path.join(scriptDir,'../data/S809_re00.75M_EXP_CDW.csv'), verbose=verbose)
#         np.testing.assert_array_equal(df.columns, ['Alpha', 'Cl', 'Cd', 'Cm', 'Cd_p', 'Cd_w']) 
# 
    def test_OpenFASTWeird(self):
        verbose=False
        df, re = loadPolarFile(os.path.join(scriptDir,'../data/FFA-W3-241-Re12M.dat'), verbose=verbose)
        #print(df.columns)
        np.testing.assert_array_equal(df.columns, ['Alpha', 'Cl', 'Cd', 'Cm'])
# 
    def test_AD(self):
        verbose=False
        df, re = loadPolarFile(os.path.join(scriptDir,'../data/Cylinder.dat'), verbose=verbose)
        np.testing.assert_array_equal(df.columns, ['Alpha', 'Cl', 'Cd', 'Cm']) 
        #print(df.columns)


if __name__ == '__main__':
    #test = TestPolarFile()
    #test.test_CSVNoHeader()
    #test.test_CSVOneHeader()
    unittest.main()
