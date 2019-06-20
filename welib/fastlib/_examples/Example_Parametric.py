import numpy as np
import os
import welib.fastlib as fastlib


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
        return '_ws{:02.0f}_rpm{:02.1f}'.format(p['InflowFile|HWindSpeed'],p['EDFile|RotSpeed'])

    # --- Generating all files in a workdir
    fastfiles=fastlib.templateReplace(ref_dir,PARAMS,workdir=work_dir,name_function=naming,RemoveRefSubFiles=True,main_file=main_file)
    print(fastfiles)


    # --- Creating a batch script just in case
    with open(os.path.join(work_dir,'_RUN_ALL.bat'), 'w') as f:
        for l in [fastlib.FAST_EXE + ' '+ os.path.basename(f) for f in fastfiles]:
            f.write("%s\n" % l)
    # --- Running the simulations
    fastlib.run_fastfiles(fastfiles,fastExe=FAST_EXE,parallel=True,ShowOutputs=False,nCores=2)

    # --- Simple Postprocessing
    outFiles = [os.path.splitext(f)[0]+'.outb' for f in fastfiles]
    avg_results = fastlib.averagePostPro(outFiles,TimeAvgWindow=10, ColMap = {'WS_[m/s]':'Wind1VelX_[m/s]'},ColSort='WS_[m/s]')
    print(avg_results)


def CPLambdaExample():
    """ Example to determine the CP-CT Lambda Pitch matrices of a turbine.
    This scrip uses the function CPCT_LambdaPitch which basically does the same as ParametricExample
    above.
    """
    ref_dir          = 'NREL5MW/'   # Folder where the fast input files are located (will be copied)
    main_file        = 'DLC120_ws07_ye000_s1_r1.fst'  # Main file in ref_dir, used as a template
    FAST_EXE         = 'NREL5MW/OpenFAST2_x64s_ebra.exe' # Location of a FAST exe (and dll)

    # --- Computing CP and CT matrices for range of lambda and pitches
    Lambda = np.linspace(0.1,10,3)
    Pitch  = np.linspace(-10,10,4)

    CP,CT,Lambda,Pitch,MaxVal,result = fastlib.CPCT_LambdaPitch(ref_dir,main_file,Lambda,Pitch,fastExe=FAST_EXE,ShowOutputs=False,nCores=4,TMax=10)

    print('CP max',MaxVal)

    # --- Plotting matrix of CP values
    from mpl_toolkits.mplot3d import Axes3D
    from matplotlib import cm
    import matplotlib.pyplot as plt
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    LAMBDA, PITCH = np.meshgrid(Lambda, Pitch)
    CP[CP<0]=0
    surf = ax.plot_surface(LAMBDA, PITCH, np.transpose(CP), cmap=cm.coolwarm, linewidth=0, antialiased=True,alpha=0.8)
    ax.scatter(MaxVal['lambda_opt'],MaxVal['pitch_opt'],MaxVal['CP_max'],c='k',marker='o',s=20)
    fig.colorbar(surf, shrink=0.5, aspect=5)
    plt.show()



if __name__=='__main__':
    ParametricExample()
    CPLambdaExample()
