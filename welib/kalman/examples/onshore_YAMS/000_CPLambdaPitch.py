import numpy as np
import os
import welib.weio as weio # https://github.com/ebranlard/weio
import wtDigiTwin.fast.fastlib as fastlib # latest fastlib is found at https://github.com/ebranlard/welib

def CPLambda():
    """ Determine the CP-CT Lambda Pitch matrices of a turbine.
    This scrip uses the function CPCT_LambdaPitch which basically does the same as ParametricExample
    above.
    """
    GE=True
    ReRun=False # we don't rerun simulations that were already run
    base             = '../../_data/NREL5MW'
    ref_dir          = '../../_data/NREL5MW/'  # Folder where the fast input files are located (will be copied)
    main_file        = '../../_data/NREL5MW/Main_Onshore.fst'  # Main file in ref_dir, used as a template
    FAST_EXE         = '../../_data/OpenFAST2_x64s_ebra.exe' # Location of a FAST exe (and dll)

    # --- Computing CP and CT matrices for range of lambda and pitches
    nLambda = 6
    nPitch  = 8
    Lambda = np.linspace(0.10 ,22,nLambda)
    Pitch  = np.linspace(-5,40,nPitch)

    CP,CT,Lambda,Pitch,MaxVal,result = fastlib.CPCT_LambdaPitch(ref_dir,main_file,Lambda,Pitch,fastExe=FAST_EXE,showOutputs=False,nCores=4,TMax=30,reRun=ReRun)

    print('CP max',MaxVal)

    np.savetxt(base+'_Lambda.csv',Lambda,delimiter = ',')
    np.savetxt(base+'_Pitch.csv' ,Pitch ,delimiter = ',')
    np.savetxt(base+'_CP.csv'    ,CP    ,delimiter = ',')
    np.savetxt(base+'_CT.csv'    ,CT    ,delimiter = ',')

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


# def focus_mesh(xmin,xmax,xfocus)

if __name__=='__main__':
    CPLambda()
