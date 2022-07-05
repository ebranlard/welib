""" 
Convert a HAWC2 blade data to BeamDyn
"""
import os
import numpy as np
import pandas as pd

import welib.fast.beamdyn as bd
import matplotlib.pyplot as plt
#import pyFAST.case_geneneration.case_gen as case_gen
#import pyFAST.case_geneneration.runner as case_gen
#import pyFAST.input_output.postpro as postpro
#import pyFAST.input_output as io

# Get current directory so this script can be called from any location
MyDir=os.path.dirname(__file__)

def main():
    np.set_printoptions(linewidth=300)

    # ---  Hawc2 to BeamDyn
    # See documentation in beamdyn.py
    #   ref_axis: string defining how the main axis of beamdyn will be defined.
    #            'c2def-polyfit': the reference axis is Hawc2 c2def, smoothened out using a polyfit (see poly_exp)
    #            'straight': the reference axis is straight (prebend and sweep are still included as offsets) 

    # NOTE: these hawc2 5MW data were downloaded form hawc2's website
    H2MeanLine     = os.path.join(MyDir,'../../../data/NREL5MW/hawc2/Blade_Planform_Hawc2.csv')   # csv file with c2def columns: ['X_[m]','Y_[m]','Z_[m]','Twist_theta_z_[deg]' (optional: Relative_thickness_[%] and Chord_[m])
    H2Structure    = os.path.join(MyDir,'../../../data/NREL5MW/hawc2/Blade_Structural_Hawc2.csv') # csv file with columns ['Radius_[m]','m_[kg/m]','x_cg_[m]','y_cg_[m]','ri_x_[m]',... ,'x_e_[m]','y_e_[m]']

    BDMainTemplate = os.path.join(MyDir,'../../../data/NREL5MW/5MW_Baseline/NRELOffshrBsline5MW_BeamDyn.dat')  # template file to write main BD file

    BDOut          = '_BD_Smooth.dat'        # Name of BeamDyn file to be written
    BDBldOut       = '_BD_Blade_Smooth.dat'  # Name of BeamDyn blade file to be written
    Mu             = [0.001]*6               # stiffness proportional damping values
    fig = bd.hawc2ToBeamDyn(H2MeanLine, H2Structure, BDBldOut, BDOut, BDMainTemplate, Mu=Mu, poly_exp=[2,3,4,5,6], ref_axis='straight', bPlot=True)


if __name__=='__main__':
    main() 
    plt.show()

if __name__=='__test__':
    main() 
    try:
        os.remove('_BD_Smooth.dat')
        os.remove('_BD_Blade_Smooth.dat')
    except:
        pass

if __name__=='__export__':
    pass

