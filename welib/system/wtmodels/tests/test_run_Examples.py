import unittest
import numpy as np    
import glob
import os

def execfile(filepath, globals=None, locals=None):
    """ Execute a given python file """
    if globals is None:
        globals = {"__name__": "__main__"}
    globals.update({
        "__file__": filepath,
    })
    with open(filepath, 'rb') as file:
        exec(compile(file.read(), filepath, 'exec'), globals, locals)

class TestExamples(unittest.TestCase):
    def test_run_examples(self):
        exclude_list=[]
        # Add tests to class
        MyDir=os.path.dirname(__file__)
        files = glob.glob(os.path.join(MyDir,'../examples/[a-zA-Z][a-zA-Z-_0-9]*.py'))
        print('\n--------------------------------------------------------------')
        import matplotlib.pyplot as plt
        for f in files:
            print('Running: {}'.format(os.path.relpath(f, MyDir)))
            if hasattr(self,'subTest'):
                with self.subTest(filename=os.path.basename(f)):
                    execfile(f, {'__name__': '__test__', 'print': lambda *_:None})
            try:
                plt.close('all')
            except:
                pass


if __name__ == '__main__':
    unittest.main()
