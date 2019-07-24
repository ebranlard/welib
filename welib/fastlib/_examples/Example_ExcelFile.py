import os
import pandas as pd

try:
    import welib.fastlib.fastlib as fastlib
except:
    import fastlib
import weio

# --- Main Parameters

ref_dir          = 'NREL5MW/'   # Folder where the fast input files are located (will be copied)
work_dir         = 'NREL5MW_ParametricExcel/'     # Output folder (will be created)
FAST_EXE         = 'NREL5MW/OpenFAST2_x64s_ebra.exe' # Location of a FAST exe (and dll)
main_file        = 'Main_Onshore_OF2.fst'  # Main file in ref_dir, used as a template
parametricFile   = 'ParametricExcel.xlsx' # Excel file containing set of parameters


# --- Reading Excel file and converting it to a list of dictionaries
dfs=weio.read(parametricFile).toDataFrame()
df=dfs[list(dfs.keys())[0]]
PARAMS_EXCEL=df.to_dict('records')
print(df)

PARAMS=[]
for p in PARAMS_EXCEL:
    p['FAST|DT']    =0.02
    p['FAST|DT_out']    =0.02
    p['AeroFile|TwrPotent']   = '1'
    p['AeroFile|TIDrag']    = 'True'
    p = fastlib.paramsStiff(p)
    PARAMS.append(p)

def naming(p):
    return p['__name__']

fastFiles=fastlib.templateReplace(ref_dir,PARAMS,workdir=work_dir,name_function=naming,RemoveRefSubFiles=True,RemoveAllowed=False,main_file=main_file)

# --- Running fast simulations
print('>>> Running {} simulations...'.format(len(fastFiles)))
fastlib.writeBatch(os.path.join(work_dir,'_RUN_ALL.bat'),fastFiles,fastExe=FAST_EXE)
fastlib.run_fastfiles(fastFiles, ShowOutputs=False, fastExe=FAST_EXE, nCores=4)

# # --- Postpro - Computing averages at the end of the simluation
print('>>> Postprocessing...')
outFiles = [os.path.splitext(f)[0]+'.outb' for f in fastFiles]
ColKeepStats  = ['RotSpeed_[rpm]','BldPitch1_[deg]','RtAeroCp_[-]','RtAeroCt_[-]','Wind1VelX_[m/s]']
result = fastlib.averagePostPro(outFiles,avgMethod='constantwindow',avgParam=5,ColKeep=ColKeepStats,ColSort='RotSpeed_[rpm]')
result.to_csv('ParametricExcel_Summary.csv',sep='\t',index=False)

