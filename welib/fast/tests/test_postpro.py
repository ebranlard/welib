import unittest
import numpy as np    
import pandas as pd
import os as os
from welib.fast.postpro import *
from welib.weio.fast_output_file import FASTOutputFile

scriptDir=os.path.dirname(__file__)

class Test(unittest.TestCase):
    def test_spanwisePostPro(self):
        scriptDir=os.path.dirname(__file__)

        outFile = os.path.join(scriptDir,'../../../data/example_files/fastout_allnodes.outb')
        fstFile = os.path.join(scriptDir,'../../../data/NREL5MW/Main_Onshore.fst')
        df = FASTOutputFile(outFile).toDataFrame()

        # --- Step 1: Read an openfast output file
        out = spanwisePostPro(FST_In=fstFile, avgMethod='periods', avgParam=1, df=df)

        #fstFile = os.path.join(scriptDir,'../../../data/NREL5MW/Main_Onshore.fst')
            
        np.testing.assert_equal(list(out.keys()), ['AD', 'ED_bld', 'ED_twr', 'BD', 'SD_MembersOut', 'SD_JointsOut'])

        # --- Test AD
        cols =['r/R_[-]','B1Cl_[-]','B1Fx_[N/m]','B1Fy_[N/m]','i_[#]','r_[m]']
        Mf = np.array([0.023810, 0.000000, 108.233875,  -1.211481 ,    1,  1.5000])
        Ml = np.array([0.999998, 0.276654, 126.151657,   4.402287 ,   19, 62.9999])
        dfAD2 = pd.DataFrame(data=np.vstack((Mf,Ml)), columns=cols)

        dfAD1 = out['AD']
        col1 = dfAD1.columns.values
        col2 = dfAD2.columns.values
        col1.sort()
        col2.sort()
        np.testing.assert_equal(col1, col2)

        dfAD1 = dfAD1.iloc[[0,18]]
        for c in cols:
            np.testing.assert_array_almost_equal(dfAD1[c].values, dfAD2[c].values, 5)

        # --- Test ED
        cols = ['r/R_[-]','B1MLxNT_[kN-m]','B1MLyNT_[kN-m]','B1TDx_[m]','B1TDy_[m]','i_[#]','r_[m]']
        Mf = np.array([0.052521,282.117907,623.326235, 0.000112,-0.000045, 1, 3.308824])
        Ml = np.array([0.971289,  0.051362,  0.250493, 0.254712,-0.053140, 17,61.191176])
        dfED2 = pd.DataFrame(data=np.vstack((Mf,Ml)), columns=cols)
        dfED1 = out['ED_bld']
        col1 = dfED1.columns.values
        col2 = dfED2.columns.values
        col1.sort()
        col2.sort()
        np.testing.assert_equal(col1, col2)
        dfED1 = dfED1.iloc[[0,16]]
        for c in cols:
            np.testing.assert_array_almost_equal(dfED1[c].values, dfED2[c].values, 5)




if __name__ == '__main__':
    unittest.main()
