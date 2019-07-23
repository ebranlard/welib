import numpy as np
import os
try:
    import welib.fastlib.fastlib as fastlib
except:
    import fastlib


def ParametricExample():
    """ Example to run a set of FAST simulations (parametric study)
    This script is based on a reference directory which contains a reference main input file (.fst)
    Everything is copied to a working directory.
    The different fast inputs are generated based on a list of dictionaries, named `PARAMS`.
    For each dictionary:
       - they keys are "path" to a input parameter, e.g. `EDFile|RotSpeed`  or `FAST|TMax`.
           These should correspond to whater name of the variable is used in the FAST inputs files.
       - they values are the values corresponding to this parameter
    """
    # --- Parameters for this script
    ref_dir          = 'NREL5MW/'   # Folder where the fast input files are located (will be copied)
    work_dir         = 'NREL5MW_Parametric/'     # Output folder (will be created)
    main_file        = 'DLC120_ws07_ye000_s1_r1.fst'  # Main file in ref_dir, used as a template
    FAST_EXE         = 'NREL5MW/OpenFAST2_x64s_ebra.exe' # Location of a FAST exe (and dll)

    # --- Defining the parametric study  (list of dictionnaries with keys as FAST parameters)
    WS = [3,5,6,7]
    RPM = [10,12,13,15]
    BaseDict = {'FAST|TMax': 10, 'FAST|DT': 0.01, 'FAST|DT_Out': 0.1}
    BaseDict = fastlib.paramsNoController(BaseDict)
    #BaseDict = fastlib.paramsStiff(BaseDict)
    #BaseDict = fastlib.paramsNoGen(BaseDict)
    PARAMS=[]
    for wsp,rpm in zip(WS,RPM): # NOTE: same length of WS and RPM otherwise do multiple for loops
        p=BaseDict.copy()
        p['EDFile|RotSpeed']       = rpm
        p['InflowFile|HWindSpeed'] = wsp
        p['InflowFile|WindType']   = 1 # Setting steady wind
        PARAMS.append(p)
    # --- Defining a function to name the files  based on the parameters
    def naming(p):
        return '{:03d}_ws{:04.1f}_om{:04.2f}'.format(p['__index__'],p['InflowFile|HWindSpeed'],p['EDFile|RotSpeed'])

    # --- Generating all files in a workdir
    fastfiles=fastlib.templateReplace(ref_dir,PARAMS,workdir=work_dir,name_function=naming,RemoveRefSubFiles=True,main_file=main_file)
    print(fastfiles)


    # --- Creating a batch script just in case
    with open(os.path.join(work_dir,'_RUN_ALL.bat'), 'w') as f:
        for l in [fastlib.FAST_EXE + ' '+ os.path.basename(f) for f in fastfiles]:
            f.write("%s\n" % l)

    # --- Creating a batch script just in case
    fastlib.writeBatch(os.path.join(work_dir,'_RUN_ALL.bat'),fastfiles,fastExe=FAST_EXE)
    # --- Running the simulations
    fastlib.run_fastfiles(fastfiles,fastExe=FAST_EXE,parallel=True,ShowOutputs=False,nCores=2)

    # --- Simple Postprocessing
    outFiles = [os.path.splitext(f)[0]+'.outb' for f in fastfiles]
    avg_results = fastlib.averagePostPro(outFiles,avgMethod='constantwindow',avgParam=10, ColMap = {'WS_[m/s]':'Wind1VelX_[m/s]'},ColSort='WS_[m/s]')
    print(avg_results)



if __name__=='__main__':
    ParametricExample()
