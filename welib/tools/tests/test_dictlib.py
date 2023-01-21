import unittest
import numpy as np
from welib.tools.dictlib import *

class TestDictLib(unittest.TestCase):

    def test_rename(self):
        D={'A':1, 'B':2}
        # --- Loose renaming
        # Case 1 - not all keys present, only the ones present are changed
        keyMap = {'A':'AA', 'C':'CC'} # 'old':'new'
        B = renameDict(D, keyMap)
        self.assertEqual(B, {'AA':1, 'B':2})

        # Case 2 - no keys present, nothing happens
        keyMap = {'C':'AA'} # 'old':'new'
        B = renameDict(D, keyMap)
        self.assertEqual(B, {'A':1, 'B':2})


        # --- Strict renaming
        # Full map is present
        keyMap = {'A':'AA', 'B':'BB'} # 'old':'new'
        B = renameDict(D, keyMap, full_map=True)
        self.assertEqual(B, {'AA':1, 'BB':2})

        # Map is not full
        keyMap = {'A':'AA'} # 'old':'new'
        with self.assertRaises(DictLibError) as e: 
            renameDict(D, keyMap, full_map=True)
        #print(e.exception)

        # --- Renaming inclusive
        # Partial map is present
        keyMap = {'A':'AA'} # 'old':'new'
        B = renameDict(D, keyMap, keys_exist=True)
        self.assertEqual(B, {'AA':1, 'B':2})

        # Partial map with some keys not present
        keyMap = {'A':'AA', 'C':'CC'} # 'old':'new'
        with self.assertRaises(DictLibError) as e: 
            renameDict(D, keyMap, keys_exist=True)
        #print(e.exception)


    def test_renameKey(self):
        D={'A':1, 'B':2}

        # --- Easy rename
        # Case 1 - one option 
        B = renameDictKey(D, 'A', 'AA')
        self.assertEqual(B, {'AA':1, 'B':2})

        # Case 2 - one option case insensitive
        B = renameDictKey(D, 'a', 'AA')
        self.assertEqual(B, {'AA':1, 'B':2})

        # Case 3 - multiple option case insensitive
        B = renameDictKey(D, ['AAAA','A'], 'AA')
        self.assertEqual(B, {'AA':1, 'B':2})

        # Case 4 - key does not exist, nothing happens
        B = renameDictKey(D, ['CCC','CC','C'] , 'AA')
        self.assertEqual(B, {'A':1, 'B':2})


        # --- Corner cases
        # Key must exist if user specifies key_exist
        with self.assertRaises(DictLibError) as e: 
            B = renameDictKey(D, ['CCC','CC','C'] , 'AA', key_exist=True)

        # Key not found if case is different
        with self.assertRaises(DictLibError) as e: 
            B = renameDictKey(D, ['a'] , 'AA', key_exist=True, case_sensitive=True)
        #print(e.exception)

        # Multiple matches
        with self.assertRaises(DictLibError) as e: 
            B = renameDictKey(D, ['a','b'] , 'AA')
        #print(e.exception)

if __name__ == '__main__':
    unittest.main()
