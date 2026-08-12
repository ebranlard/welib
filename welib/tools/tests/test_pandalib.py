import numpy as np
import pandas as pd
import unittest
from welib.tools.pandalib import *

# --------------------------------------------------------------------------------}
# ---  
# --------------------------------------------------------------------------------{
class TestPandaLib(unittest.TestCase):

    # --------------------------------------------------------------------------------}
    # ---  
    # --------------------------------------------------------------------------------{
    def test_inverse_colmap_basic_rename(self):
        colmap = {
            'Pitch_[deg]': 'Bld1Pitch_[deg]',
            'WS_[m/s]': '{Wind1VelX_[m/s]}'
        }
        expected = {
            'Bld1Pitch_[deg]': 'Pitch_[deg]',
            'Wind1VelX_[m/s]': 'WS_[m/s]'
        }
        inv = inverse_colmap(colmap)
        #print(inv)
        assert inv == expected


    def test_inverse_colmap_arithmetic(self):
        colmap = {
            'RotSpeed_[rad/s]': ' {RotSpeed_[rpm]} * 2*np.pi/60',
            'RotSpeed_[rpm]'  : '   {RotSpeed_[rad/s]} / (2*np.pi/60)'
        }
        inv = inverse_colmap(colmap)
        assert inv['RotSpeed_[rpm]'] == '{RotSpeed_[rad/s]} / (2*np.pi/60)'
        assert inv['RotSpeed_[rad/s]'] == '{RotSpeed_[rpm]} * ((2*np.pi/60))'


    def test_inverse_colmap_skips_non_invertible(self):
        colmap = {
            'TotalSpeed': '{Speed1} * 2 + {Speed2} * 3', # Multiple columns -> Non-invertible
            'R_[m]': '{ones} * 15',                      # Constant -> Non-invertible
            'ValidCol': 'OldCol'                         # Invertible
        }
        inv = inverse_colmap(colmap, verbose=False)
        assert 'TotalSpeed' not in inv
        assert 'ones' not in inv
        assert inv == {'OldCol': 'ValidCol'}

    def test_inverse_colmap_list_of_candidates(self):
        colmap = {
            'q_p': ['Q_P_[rad]', '{PtfmSurge_[deg]}*np.pi/180']
        }
        inv = inverse_colmap(colmap)
        assert 'Q_P_[rad]' in inv
        assert inv['Q_P_[rad]'] == 'q_p'
# 
# 
#     # ------------------------------------------------------------------------------
#     # Integration Test with remap_df roundtrip
#     # ------------------------------------------------------------------------------
# 
    def test_remap_df_roundtrip(self):
        # Sample DataFrame
        df_orig = pd.DataFrame({
            'BldDummy': [1.0, 2.0, 3.0],
            'Bld1Pitch_[deg]': [1.0, 2.0, 3.0],
            'RotSpeed_[rpm]': [10.0, 20.0, 30.0]
        })

        colmap = {
            'Dummy': 'BldDummy',
            'Pitch_[rad]': 'np.pi/180  * {Bld1Pitch_[deg]} ',
            'RotSpeed_[rad/s]': '{RotSpeed_[rpm]} * 2*np.pi/60'
        }

        # Forward remapping
        df_mapped = remap_df(df_orig.copy(), colmap, bColKeepNewOnly=True)
        assert 'Dummy' in df_mapped.columns
        assert 'Pitch_[rad]' in df_mapped.columns
        assert 'RotSpeed_[rad/s]' in df_mapped.columns

        # Invert colmap
        inv_colmap = inverse_colmap(colmap)
        #print('inv_colmap', inv_colmap)

        # Reverse remapping back to original
        df_restored = remap_df(df_mapped.copy(), inv_colmap, bColKeepNewOnly=True)

        # Test mapping
        np.testing.assert_allclose( df_mapped['RotSpeed_[rad/s]'].values, df_orig['RotSpeed_[rpm]'].values*2*np.pi/60)
        np.testing.assert_allclose( df_mapped['Pitch_[rad]'].values, df_orig['Bld1Pitch_[deg]'].values*np.pi/180)

        # Test round trip
        np.testing.assert_allclose( df_restored['Bld1Pitch_[deg]'].values, df_orig['Bld1Pitch_[deg]'].values)
        np.testing.assert_allclose( df_restored['RotSpeed_[rpm]'].values, df_orig['RotSpeed_[rpm]'].values)

if __name__ == '__main__':
    unittest.main()
    
