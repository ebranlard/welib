import numpy as np

def tuneTowerDamping(fstFile, zeta_target=None, modes=None, logdec_target=None, verbose=False):
    """ 
    Provide values of the tower damping ratios (zeta) necessary to achieve a target damping 
    for the global modes (full structure, onshore only for now)

    INPUTS:
     - fstFile: file name of a valid OpenFAST FST file, which reference valid ElastoDyn files
     - modes: list of mode names, in 'FA1', 'FA2', 'SS1', 'SS2'
            if None: modes = ['FA1','FA2','SS1','SS2']
     - zeta_target: damping ratios of the global mode names (not in %)
     OR
     - logdec_target: logarithmic decrements of the global modes

    OUTPUTS:
     - zeta_tower: dictionary with key the mode names and value the damping ratio of the tower
                   such that it should match the target global value. (not in %)

    EXAMPLES:
        zeta_tower = tuneTowerDamping('main.fst', zeta_target=[0.02, 0.03], modes=['FA1','SS1'])

        zeta_tower = tuneTowerDamping('main.fst', zeta_target=[0.02, 0.02, 0.03, 0.03])

    """
    from welib.yams.windturbine import FASTWindTurbine
    # --- Sanitation
    if logdec_target is None and zeta_target is None:
        raise Exception('Provide either zeta_target or logdec_target')
    if logdec_target is not None:
        logdec_target = np.asarray(logdec_target)
        zeta_target = logdec_target/(np.sqrt(4*np.pi**2 + logdec_target**2))
    zeta_target = np.asarray(zeta_target)
    if modes is None:
        modes = ['FA1','FA2','SS1','SS2']
    if len(zeta_target)!=len(modes):
        raise Exception('The number of modes and target damping should be the same')

    glM_Map    = {'fa1':0, 'fa2':1, 'ss1':2, 'ss2':3} # Maps of global Modes


    # --- Read wind turbine
    WT = FASTWindTurbine(fstFile, algo='OpenFAST') #, bldStartAtRotorCenter=False )

    # --- Compute target tower damping 
    zeta_tower = {}
    for zeta_target, mode in zip(zeta_target, modes):
        iDOF = glM_Map[mode.lower()]
        # --- Extracting relevant variables
        MT    = WT.RNA.mass     # Tower top mass
        zetat = WT.twr.damp_zeta[0]
        GMt   = WT.twr.MM [6+iDOF,6+iDOF] # Generalized mass matrix, isolated tower
        GMtt  = WT.twr.MM [6+iDOF,6+iDOF] + WT.RNA.mass
        GKt   = WT.twr.KK0[6+iDOF,6+iDOF] # Generalized stiffness matrix, isolated tower
        GKtt  = WT.twr.KK [6+iDOF,6+iDOF]  # Generalized stiffness matrix with stiffening
        om_t  = np.sqrt(GKt /GMt)    # without top mass and stiffening
        om_tt = np.sqrt(GKtt/(GMtt)) # With top mass and stiffening
        f_t   = om_t  /(2*np.pi) 
        f_tt  = om_tt /(2*np.pi) 
        #print('zetat ' ,zetat)
        #print('MT ' ,MT)
        #print('GMt ',GMt)
        #print('GMtt',GMtt)
        #print('GKt ',GKt)
        #print('GKtt',GKtt)
        #print('GDtt',GDtt)
        #print('ft  ={:8.5f}'.format(f_t  ))
        #print('ftt ={:8.5f}'.format(f_tt))
        zetatt = zetat*np.sqrt((GKt*GMt)/(GKtt*GMtt))
        zetat_target = zeta_target * np.sqrt((GKtt*GMtt)/(GKt*GMt))

        if verbose:
            print('Global {} mode: {:7.3f}Hz {:7.4f}%'.format(mode,f_tt, zetatt*100))
            print('             to get zeta_g={:7.4f}% , use zeta_t={:7.4f} [%]'.format(zeta_target*100, zetat_target*100))
        zeta_tower[mode] = zetat_target
    return zeta_tower
             




