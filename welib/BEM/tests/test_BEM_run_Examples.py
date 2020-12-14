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

def test_generator(pyfilename):
    """ Create a function to be added to the Test Class. The function execute a given python file """
    def test_function(self):
        print('\n--------------------------------------------------------------')
        print('Running example script: {}'.format(pyfilename))
        execfile(pyfilename, {'__name__': '__test__', 'print': lambda *_:None})
    return test_function


class TestExamples(unittest.TestCase):
    pass

if __name__ == '__main__':
    exclude_list=[]
    # Add tests to class
    MyDir=os.path.dirname(__file__)
    files = glob.glob('../examples/[a-zA-Z]*.py')
    for f in files:
        if f not in exclude_list:
            test_name = 'test_{}'.format(os.path.splitext(os.path.basename(f))[0])
            setattr(TestExamples, test_name, test_generator(f))
    # Run tests
    unittest.main()
