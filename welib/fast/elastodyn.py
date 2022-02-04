import numpy as np
from welib.weio.fast_input_file import FASTInputFile
from welib.system.eva import eigMCK


def DTfreq(InFile):
    """ Returns ElastoDyn torsion drive train frequency
    INPUTS:
       - InFile: fst file
    OUTUPTS:
        - f0: natural frequency
        - fd: damped frequency
        - zeta: damping ratio
    """
    from welib.yams.windturbine import FASTWindTurbine
    WT     = FASTWindTurbine(InFile)
    nGear  = WT.ED['GBRatio']
    K_DT   = WT.ED['DTTorSpr']
    D_DT   = WT.ED['DTTorDmp']
    Jr_LSS = WT.rot.inertia[0,0]*1.000    # bld + hub, LSS
    Jg_LSS = WT.gen.inertia[0,0]          # gen, LSS

    M = np.array([[Jr_LSS+Jg_LSS,-Jg_LSS],[-Jg_LSS,Jg_LSS]])
    K = K_DT*np.array([[0, 0],[0, 1]])
    C = D_DT*np.array([[0, 0],[0, 1]])

    fd, zeta, Q, f0 = eigMCK(M, C, K)

    # f0 =  np.sqrt(K_DT/Jr_LSS + K_DT/Jg_LSS)/(2*np.pi))
    return f0, fd, zeta
