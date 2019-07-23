import numpy as np
import os
try:
    import welib.fastlib.fastlib as fastlib
except:
    import fastlib

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
    CPLambdaExample()
