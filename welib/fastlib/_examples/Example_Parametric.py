import numpy as np
import os
try:
    import welib.fastlib.fastlib as fastlib
except:
    import fastlib


def ParametricExample():
    """ Example to run a set of OpenFAST simulations (parametric study)

    This script uses a reference directory (`ref_dir`) which contains a reference input file (.fst)
    1) The reference directory is copied to a working directory (`work_dir`).
    2) All the fast input files are generated in this directory based on a list of dictionaries (`PARAMS`).
    For each dictionary in this list:
       - The keys are "path" to a input parameter, e.g. `EDFile|RotSpeed`  or `FAST|TMax`.
         These should correspond to the variables used in the FAST inputs files.
       - The values are the values corresponding to this parameter
    For instance:
         PARAMS[0]['EDFile|RotSpeed']       = 5
         PARAMS[0]['InflowFile|HWindSpeed'] = 10

    3) The simulations are run, successively distributed on `nCores` CPUs.
    4) The output files are read, and averaged based on a method (e.g. average over a set of periods,
        see averagePostPro in fastlib for the different averaging methods).
       A pandas DataFrame is returned

    """
    # --- Parameters for this script
    ref_dir          = 'NREL5MW/'   # Folder where the fast input files are located (will be copied)
    work_dir         = 'NREL5MW_Parametric/'     # Output folder (will be created)
    main_file        = 'Main_Onshore_OF2.fst'  # Main file in ref_dir, used as a template
    FAST_EXE         = 'NREL5MW/OpenFAST2_x64s_ebra.exe' # Location of a FAST exe (and dll)

    # --- Defining the parametric study  (list of dictionnaries with keys as FAST parameters)
    WS = [3,5,6,7]
    RPM = [10,12,13,15]
    BaseDict = {'TMax': 10, 'DT': 0.01, 'DT_Out': 0.1}
    BaseDict = fastlib.paramsNoController(BaseDict)   # Remove the controller
    #BaseDict = fastlib.paramsControllerDLL(BaseDict) # Activate the controller
    #BaseDict = fastlib.paramsStiff(BaseDict)         # Make the turbine stiff (except generator)
    #BaseDict = fastlib.paramsNoGen(BaseDict)         # Remove the Generator DOF
    PARAMS=[]
    for i,(wsp,rpm) in enumerate(zip(WS,RPM)): # NOTE: same length of WS and RPM otherwise do multiple for loops
        p=BaseDict.copy()
        if i==0:
            p['AeroFile|TwrAero']       = True
        if i==2:
            p['EDFile|BldFile(1)|AdjBlMs'] =1.1
            p['EDFile|BldFile(2)|AdjBlMs'] =1.1
            p['EDFile|BldFile(3)|AdjBlMs'] =1.1

        p['EDFile|RotSpeed']       = rpm
        p['InflowFile|HWindSpeed'] = wsp
        p['InflowFile|WindType']   = 1 # Setting steady wind
        p['__name__']='{:03d}_ws{:04.1f}_om{:04.2f}'.format(i,p['InflowFile|HWindSpeed'],p['EDFile|RotSpeed'])
        PARAMS.append(p)
        i=i+1
    # --- Generating all files in a workdir
    fastfiles=fastlib.templateReplace(PARAMS,ref_dir,workdir=work_dir,RemoveRefSubFiles=True,main_file=main_file, oneSimPerDir=False)
    print(fastfiles)

    # --- Creating a batch script just in case
    fastlib.writeBatch(os.path.join(work_dir,'_RUN_ALL.bat'),fastfiles,fastExe=FAST_EXE)
    # --- Running the simulations
    fastlib.run_fastfiles(fastfiles,fastExe=FAST_EXE,parallel=True,ShowOutputs=False,nCores=4)

    # --- Simple Postprocessing
    outFiles = [os.path.splitext(f)[0]+'.outb' for f in fastfiles]
    avg_results = fastlib.averagePostPro(outFiles,avgMethod='periods',avgParam=1, ColMap = {'WS_[m/s]':'Wind1VelX_[m/s]'},ColSort='WS_[m/s]')
    print(avg_results)
    return avg_results



if __name__=='__main__':
    ParametricExample()
