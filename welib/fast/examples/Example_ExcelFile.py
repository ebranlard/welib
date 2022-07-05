"""
Generate, run and postprocess openfast cases using an Excelfile which defines the parameters to change for each simulation
"""
import os
import pandas as pd

import welib.fast.case_gen as case_gen
import welib.fast.runner as runner
import welib.fast.postpro as postpro
import welib.weio as io
#import pyFAST.case_geneneration.case_gen as case_gen
#import pyFAST.case_geneneration.runner as case_gen
#import pyFAST.input_output.postpro as postpro
#import pyFAST.input_output as io

# Get current directory so this script can be called from any location
MyDir=os.path.dirname(__file__)


def main():
    # --- Main Parameters
    ref_dir        = os.path.join(MyDir, '../../../data/NREL5MW/')  # Folder where the fast input files are located (will be copied)
    FAST_EXE       = os.path.join(MyDir, '../../../data/openfast.exe') # Location of a FAST exe (and dll)
    main_file      = 'Main_Onshore.fst'          # Main file in ref_dir, used as a template
    work_dir       = os.path.join(MyDir, '_NREL5MW_ParametricExcel/')     # Output folder (will be created)
    parametricFile = 'ParametricExcel.xlsx'          # Excel file containing set of parameters

    # --- Reading Excel file, converting it to a list of dictionaries, and generate input files
    dfs    = io.excel_file.ExcelFile(parametricFile).toDataFrame()
    df     = dfs[list(dfs.keys())[0]]
    PARAMS = df.to_dict('records')
    print(df)
    fastFiles=case_gen.templateReplace(PARAMS, ref_dir, outputDir=work_dir, removeRefSubFiles=True, removeAllowed=False, main_file=main_file)

    # --- Running fast simulations
    print('>>> Running {} simulations in {} ...'.format(len(fastFiles), work_dir))
    runner.writeBatch(os.path.join(work_dir,'_RUN_ALL.bat'),fastFiles,fastExe=FAST_EXE)
    runner.run_fastfiles(fastFiles, showOutputs=False, fastExe=FAST_EXE, nCores=4)

    # --- Postpro - Computing averages at the end of the simluation
    print('>>> Postprocessing...')
    outFiles = [os.path.splitext(f)[0]+'.outb' for f in fastFiles]
    ColKeepStats  = ['RotSpeed_[rpm]','BldPitch1_[deg]','RtAeroCp_[-]','RtAeroCt_[-]','Wind1VelX_[m/s]']
    result = postpro.averagePostPro(outFiles,avgMethod='constantwindow',avgParam=5,ColKeep=ColKeepStats,ColSort='RotSpeed_[rpm]')
    result.to_csv('ParametricExcel_Summary.csv',sep='\t',index=False)
    print('Average values saved to _ParametricExcel_Summary.csv')

if __name__=='__main__':
    main() 

if __name__=='__test__':
    # Need openfast.exe, doing nothing
    pass


