import unittest
import numpy as np    
import pandas as pd
import os as os
from welib.fast.postpro import *
from welib.weio.fast_output_file import FASTOutputFile

scriptDir=os.path.dirname(__file__)

class Test(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # This code will be run once, before all the tests in the class are run.
        # FAST With ED and AD radial outputs (not a full period available..)
        outFile1 = os.path.join(scriptDir,'../../../data/example_files/fastout_allnodes.outb')
        cls.df1 = FASTOutputFile(outFile1).toDataFrame()

        # AD driver with radial outputs, more than a period
        outFile2 = os.path.join(scriptDir,'../../../data/example_files/ad_driver_yaw.1.outb')
        cls.dfAD = FASTOutputFile(outFile2).toDataFrame()

    def test_spanwisePostProNoRadADED(self):
        fstFile = os.path.join(scriptDir,'../../../data/NREL5MW/Main_Onshore.fst')
        # --- Step 1: Read an openfast output file
        out = spanwisePostPro(FST_In=fstFile, avgMethod='constantwindow', avgParam=10, df=self.df1)

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
        Mf   = np.array([0.052521,282.117907,623.326235, 0.000112,-0.000045, 1, 3.308824])
        Ml   = np.array([0.971289,  0.051362,  0.250493, 0.254712,-0.053140, 17,61.191176])
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


    def test_spanwisePostProADDriver(self):
        # --- Spanwise PostPro with inputfile available
        fstFile = os.path.join(scriptDir,'../../../data/example_files/ad_driver_yaw.dvr')
        out = spanwisePostPro(FST_In=fstFile, avgMethod='periods', avgParam=1, out_ext='.1.outb')
        dfAD1 = out['AD']
        cols  = ['r/R_[-]','B1Fn_[N/m]','B1Ft_[N/m]','i_[#]','r_[m]']
        Mf    = np.array([0.032810 ,  59.699168 ,-39.653777 ,    1 ,  3.970000])
        Ml    = np.array([1.000000 , 963.398364 ,  0.448169 ,   51 ,120.999015])
        dfAD2 = pd.DataFrame(data= np.vstack((Mf,Ml)), columns = cols)
        col1 = dfAD1.columns.values
        col2 = dfAD2.columns.values
        col1.sort()
        col2.sort()
        np.testing.assert_equal(col1, col2)
        dfAD1 = dfAD1.iloc[[0,50]]
        for c in cols:
            np.testing.assert_array_almost_equal(dfAD1[c].values, dfAD2[c].values, 5)

    def test_averageDFADDriver(self):
        # ---  Test average DF
        df = self.dfAD
        dfAvg = averageDF(df, avgMethod='periods' ,avgParam=1) 
        #print(dfAvg)
        np.testing.assert_almost_equal(dfAvg['Time_[s]'].iloc[0], 9.6)
        np.testing.assert_almost_equal(dfAvg['AB1N047Ft_[N/m]'].iloc[0], 793.240704, 5)
        np.testing.assert_almost_equal(dfAvg['AB1N047Fn_[N/m]'].iloc[0], 7952.062817, 5)


    def test_aziAverageDFADDriver(self):
        df = self.dfAD
        psiBin = np.arange(0,360, 60)
        dfAzi= azimuthal_average_DF(df, psiBin=psiBin, colPsi='Azimuth_[deg]', tStart=None, colTime='Time_[s]')
        np.testing.assert_almost_equal(dfAzi['Azimuth_[deg]'].iloc[0],   26.24, 2)
        np.testing.assert_almost_equal(dfAzi['AB1N010Fn_[N/m]'].iloc[3] , 1319.78052, 3)

        # Might change...
        np.testing.assert_array_almost_equal(dfAzi.index, [30,90,150,210,270])
        np.testing.assert_almost_equal(dfAzi['Azimuth_[deg]'].iloc[-1], 261.12, 2)
        np.testing.assert_almost_equal(dfAzi['AB1N010Ft_[N/m]'].iloc[-1], 217.0489  , 3)


if __name__ == '__main__':
    unittest.main()
