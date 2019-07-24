import numpy as np
import os
try:
    import welib.fastlib.fastlib as fastlib
except:
    import fastlib


def PowerCurveParametricExample1():
    """ Example to run a set of FAST simulations to determine a power curve.
    In this example, the WS, RPM and Pitch are set within a for loop.
    If the controller and generator are active, these are just "initial conditions".
    Additional parameters may be set by adjusting the BaseDict.

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
    work_dir         = 'NREL5MW_ParametricPowerCurve1/'     # Output folder (will be created)
    main_file        = 'Main_Onshore_OF2.fst'  # Main file in ref_dir, used as a template
    FAST_EXE         = 'NREL5MW/OpenFAST2_x64s_ebra.exe' # Location of a FAST exe (and dll)

    # --- Defining the parametric study  (list of dictionnaries with keys as FAST parameters)
    WS    = [3,5,7,9 ,11,13,15]
    RPM   = [5,6,7,10,10,10,10] # initial conditions
    PITCH = [0,0,0,0 ,5 ,10,15] # initial conditions
    BaseDict = {'FAST|TMax': 100, 'FAST|DT': 0.01, 'FAST|DT_Out': 0.1}
    #BaseDict = fastlib.paramsNoController(BaseDict)
    #BaseDict = fastlib.paramsStiff(BaseDict)
    #BaseDict = fastlib.paramsNoGen(BaseDict)
    PARAMS=[]
    for wsp,rpm,pitch in zip(WS,RPM,PITCH): # NOTE: same length of WS and RPM otherwise do multiple for loops
        p=BaseDict.copy()
        p['EDFile|RotSpeed']       = rpm
        p['EDFile|BlPitch(1)']     = pitch
        p['EDFile|BlPitch(2)']     = pitch
        p['EDFile|BlPitch(3)']     = pitch
        p['InflowFile|HWindSpeed'] = wsp
        p['InflowFile|WindType']   = 1 # Setting steady wind
        PARAMS.append(p)
        # --- Defining a function to name the files  based on the parameters
        def naming(p):
            return 'ws{:04.1f}'.format(p['InflowFile|HWindSpeed'])

    # --- Generating all files in a workdir
    fastfiles=fastlib.templateReplace(ref_dir,PARAMS,workdir=work_dir,name_function=naming,RemoveRefSubFiles=True,main_file=main_file)
    print(fastfiles)

    # --- Creating a batch script just in case
    fastlib.writeBatch(os.path.join(work_dir,'_RUN_ALL.bat'), fastfiles,fastExe=FAST_EXE)
    # --- Running the simulations
    fastlib.run_fastfiles(fastfiles,fastExe=FAST_EXE,parallel=True,ShowOutputs=False,nCores=2)

    # --- Simple Postprocessing
    outFiles = [os.path.splitext(f)[0]+'.outb' for f in fastfiles]

    avg_results = fastlib.averagePostPro(outFiles,avgMethod='constantwindow',avgParam=10, ColMap = {'WS_[m/s]':'Wind1VelX_[m/s]'},ColSort='WS_[m/s]')
    print(avg_results)
    avg_results.to_csv('PowerCurve1.csv',sep='\t',index=False)


def PowerCurveParametricExample2():
    """ Example to run a set of FAST simulations to determine a power curve.
    In this example, the WS, RPM and Pitch are set within a for loop.
    If the controller and generator are active, these are just "initial conditions".
    Additional parameters may be set by adjusting the BaseDict.

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
    work_dir         = 'NREL5MW_ParametricPowerCurve2/'     # Output folder (will be created)
    main_file        = 'Main_Onshore_OF2.fst'  # Main file in ref_dir, used as a template
    FAST_EXE         = 'NREL5MW/OpenFAST2_x64s_ebra.exe' # Location of a FAST exe (and dll)
    out_Ext          = '.outb' # Output extension

    # --- Defining the parametric study  (list of dictionnaries with keys as FAST parameters)
    WS    = [3,5,7,9 ,11,13,15]
    RPM   = [5,6,7,10,10,10,10] # initial conditions
    PITCH = [0,0,0,0 ,5 ,10,15] # initial conditions
    BaseDict = {'FAST|TMax': 100, 'FAST|DT': 0.01, 'FAST|DT_Out': 0.1}
    PARAMS,_=fastlib.paramsWS_RPM_Pitch(WS,RPM,PITCH,BaseDict=BaseDict,FlatInputs=True)

    # --- Generating all files in a workdir
    fastfiles=fastlib.templateReplace(ref_dir,PARAMS,workdir=work_dir,RemoveRefSubFiles=True,RemoveAllowed=True,main_file=main_file)

    # --- Creating a batch script just in case
    fastlib.writeBatch(os.path.join(work_dir,'_RUN_ALL.bat'), fastfiles,fastExe=FAST_EXE)

    # --- Running the simulations
    fastlib.run_fastfiles(fastfiles,fastExe=FAST_EXE,parallel=True,ShowOutputs=False,nCores=2)

    # --- Simple Postprocessing
    outFiles = [os.path.splitext(f)[0]+out_Ext for f in fastfiles]
    avg_results = fastlib.averagePostPro(outFiles,avgMethod='constantwindow',avgParam=10, ColMap = {'WS_[m/s]':'Wind1VelX_[m/s]'},ColSort='WS_[m/s]')
    print(avg_results)
    avg_results.to_csv('PowerCurve2.csv',sep='\t',index=False)



if __name__=='__main__':
    PowerCurveParametricExample1()
    PowerCurveParametricExample2()
