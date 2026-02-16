import unittest
import os
import numpy as np
import welib.weio as weio
from welib.weio.tests.helpers_for_test import MyDir, reading_test 

class Test(unittest.TestCase):
 
    def test_001_read_all(self):
        reading_test('Tecplot*.*', weio.read)

    def DF(self,FN):
        """ Reads a file with weio and return a dataframe """ 
        return weio.read(os.path.join(MyDir,FN)).toDataFrame()

    def test_tecplot_simplecsv(self):
        F=weio.read(os.path.join(MyDir,'TecplotASCII_1.dat'))
        DF=F.toDataFrame()
        self.assertEqual(DF.values[-1,1], 1.444345)
        self.assertEqual(DF.columns[1], 'Cl_[-]')

        # ---
        F.test_ascii(bCompareWritesOnly=True,bDelete=True)

    def test_tecplot_zone(self):
        F=weio.read(os.path.join(MyDir,'TecplotASCII_2.dat'))
        DF=F.toDataFrame()
        self.assertEqual(DF.values[-1,1],-56)
        self.assertEqual(DF.columns[1], 'CoordinateY')

if __name__ == '__main__':
    unittest.main()
