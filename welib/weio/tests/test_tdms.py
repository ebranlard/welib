import unittest
import os
import numpy as np

try:
    from nptdms import TdmsFile
    HasNPTDMS=True
except:
    HasNPTDMS=False

from welib.weio.tests.helpers_for_test import MyDir, reading_test 
from welib.weio.tdms_file import TDMSFile

class Test(unittest.TestCase):

    def test_001_read_all(self, DEBUG=True):
        # read all TDMS present in the example folder, using the TDMSFile class
        if HasNPTDMS:
            reading_test('TDMS*.*', TDMSFile)

    def DF(self, FN):
        """ Reads a file and return a dataframe """ 
        return TDMSFile(os.path.join(MyDir,FN)).toDataFrame()

    def test_TDMS(self):
        # One group two channels with time track
        if HasNPTDMS:
            df = self.DF('TDMS_1Grp2Chan_TimeTrack.tdms')
            self.assertEqual(df.shape,(20,3))
            self.assertEqual(df.columns[0], 'Time_[s]')
            self.assertAlmostEqual(df['ColB'].values[-1], 1.41852580388, 4)
            self.assertAlmostEqual(df['Time_[s]'].values[1], 1/19, 4)
            self.assertAlmostEqual(df['Time_[s]'].values[-1], 1, 4)

            # Two groups, two channels, no time track but timestamps
            dfs = self.DF('TDMS_2Grp2Chan.tdms')
            np.testing.assert_array_equal(list(dfs.keys()),('GroupA','GroupB'))

            df = dfs['GroupA']
            self.assertAlmostEqual(df['ColA'].values[-1], -0.35676837, 4)


if __name__ == '__main__':
    unittest.main()


