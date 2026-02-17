import unittest
import glob
import os
import numpy as np

import welib.weio as weio
from welib.weio.tests.helpers_for_test import MyDir, reading_test 

class Test(unittest.TestCase):
    def test_000_debug(self):
        #fileformat,F = weio.detectFormat('_tests/FASTIn_ED_bld.dat')
        #F = weio.CSVFile('_tests/CSVComma_Fail.csv')
        #F = weio.FASTInFile('../_tests/FASTIn_ExtPtfm_SubSef.dat')
        #F = weio.FASTLinFile('../_tests/FASTOutLin_New.lin')
        #F = weio.TecplotFile('example_files/TecplotASCII_1.dat')
        #dfs=F.toDataFrame()
        #print(F.toDataFrame())
        #print(F)
        #return
        #print(fileformat)
        pass

    def test_001_read_all(self):
        DEBUG=True
        nError=0
        for f in glob.glob(os.path.join(MyDir,'*.*')):
            if os.path.splitext(f)[-1] in ['.py','.pyc'] or f.find('_TMP')>0:
                continue
            try:
                fileformat=None
                F=None
                fileformat,F = weio.detectFormat(f)
                #fr=weio.read(f,fileformat)
                dfs = F.toDataFrame()
                ## 
                if not isinstance(dfs,dict):
                    if dfs is None:
                        n = 0
                    elif len(dfs)==0:
                        n = 0
                    else:
                        n = 1
                else:
                    n=len(dfs)
                ##print(fr.toDataFrame())
                s=fileformat.name
                s=s.replace('file','')[:20]
                if DEBUG:
                    print('[ OK ] {:30s}\t{:20s}\t{}'.format(os.path.basename(f)[:30],s,n))
            except weio.FormatNotDetectedError:
                nError += 1
                if DEBUG:
                    print('[FAIL] {:30s}\tFormat not detected'.format(os.path.basename(f)[:30]))
            except weio.OptionalImportError:
                nError += 0 # Optional..
                if DEBUG:
                    print('[WARN] {:30s}\tOptional package missing'.format(os.path.basename(f)[:30]))
            except:
                nError += 1
                if DEBUG:
                    print('[FAIL] {:30s}\tException occurred'.format(os.path.basename(f)[:30]))
                raise 

        if nError>0:
            raise Exception('Some tests failed')


if __name__ == '__main__':
    #Test().test_000_debug()
    unittest.main()
