import numpy as np

from welib.system.eva import eig


def CraigBampton(MM, KK, Ileader, nModesCB=None, Ifollow=None, F=None, DD=None): 
    """
    Performs the CraigBampton reduction of a system given some input master dofs index
    and a number of modes
        
    INPUTS
      Ileader : index of leader DOFs
      nModesCB: number of CB modes to keep
      MM, KK  : Maff and stiffness matrix
        
    INPUTS (Optional)
      nModesCB: number of CB modes to keep. Default: all
      Ifollow: indices of follower DOFs. Default: complementary set to Ileader
        
    OUTPUTS
      fc: critical frequency
      Mr,Kr,Fr,Dr: reduced mass, stiffness, force and damping  matrices
        
    AUTHOR: E. Branlard
    """
    
    # --- Input cleanup
    Ileader = np.asarray(Ileader).ravel()
    # --- Optional arguments
    if Ifollow is None:
        # Then we take the complementary to Ileader
        Iall    = np.arange(len(MM))
        Ifollow = [i for i in Iall if i not in Ileader]
    else:
        Ifollow = np.asarray(Ifollow).ravel()
    if nModesCB is None:
        nModesCB=len(Ifollow)

    # Partitioning - NOTE: leaders will be first in reduced matrix Mr and Kr
    Mll= MM[np.ix_(Ileader, Ileader)]
    Kll= KK[np.ix_(Ileader, Ileader)]
    Mff= MM[np.ix_(Ifollow, Ifollow)]
    Kff= KK[np.ix_(Ifollow, Ifollow)]
    Mlf= MM[np.ix_(Ileader, Ifollow)]
    Klf= KK[np.ix_(Ileader, Ifollow)]
    
    # --- Solve for Guyan modes
    Kff1Kfl = np.linalg.solve(Kff,(np.transpose(Klf))) # Kss1Ksm=Kss\(Kms');
    #Kff1Kfl = np.linalg.inv(Kff).dot(Klf.T)
    Kff1Kfl = np.linalg.lstsq(Kff,Klf.T, rcond=None)[0]
    Phi_G = - Kff1Kfl;

    # --- Solve EVP for constrained system
    Phi_CB, Lambda_CB = eig(Kff,Mff)
    Omega2 = np.diag(Lambda_CB).copy()
    Omega2[Omega2<0]=0.0
    f_CB  = np.sqrt(Omega2)/(2*np.pi)
    # --- Taking only thefirst few modes
    Phi_CB    = Phi_CB[:,:nModesCB]
    Lambda_CB = Lambda_CB[:,:nModesCB]
    f_CB      = f_CB[:nModesCB]
    # --- Using the T matrix:
    # # T=[eye(nm)  zeros(nm,nModesCB); -Kff1Kfl   Phi_CB];
    # # MM=[Mll Mlf; Mlf' Mff];
    # # KK=[Kll Klf; Klf' Kff];
    # # Mr=T' * MM * T;
    # # Kr=T' * KK * T;

    # --- Building reduced matrices
    #Mr11=Mmm-(Kss1Ksm')*Mms' - Mms*Kss1Ksm + (Kss1Ksm')*Mss*Kss1Ksm;
    #Kr11=Kmm-Kms*Kss1Ksm;
    #Mr12=(Mms-(Kss1Ksm')*Mss)*Psic;
    Mr11 = Mll - (np.transpose(Kff1Kfl)).dot(np.transpose(Mlf)) - Mlf.dot(Kff1Kfl) + (np.transpose(Kff1Kfl)).dot(Mff).dot(Kff1Kfl)
    Kr11 = Kll - Klf.dot(Kff1Kfl)
    Mr12 = (Mlf - (np.transpose(Kff1Kfl)).dot(Mff)).dot(Phi_CB)
    ZZ   = np.zeros((len(Ileader),nModesCB))

    # Building reduced matrix 
    Mr = np.block( [ [Mr11 , Mr12 ], [ Mr12.T, np.eye(nModesCB)       ] ])
    Kr = np.block( [ [Kr11  , ZZ  ], [ ZZ.T  ,  Lambda_CB[:nModesCB,:]] ])

    if DD is not None:
        raise Exception('Not done')
    if F is not None:
        raise Exception('Not done')

#     import pdb; pdb.set_trace()
    return Mr, Kr, Phi_G, Phi_CB, f_CB



if __name__=='__main__':
    np.set_printoptions(linewidth=500)
    L = 100
    EI = 1868211939147.334
    Maff = L * 8828.201296825122
    KK = EI / (L ** 3) * np.array([[12,6 * L,- 12,6 * L],[6 * L,4 * L ** 2,- 6 * L,2 * L ** 2],[- 12,- 6 * L,12,- 6 * L],[6 * L,2 * L ** 2,- 6 * L,4 * L ** 2]])
    MM = Maff / 420 * np.array([[156,22 * L,54,- 13 * L],[22 * L,4 * L ** 2,13 * L,- 3 * L ** 2],[54,13 * L,156,- 22 * L],[- 13 * L,- 3 * L ** 2,- 22 * L,4 * L ** 2]])
    print(MM)
    Mr,Kr,Phi_G,Phi_CB,f_CB = CraigBampton(MM,KK,[2], nModesCB=2)
    print(Mr)
    print(Kr)
    print(Phi_G)
    print(Phi_CB)
    print(f_CB)
    ## --- Solve EVA
    __,Lambda = eig(Kr,Mr)
    f= np.sqrt(np.sort(np.diag(Lambda)))/(2*np.pi)
    print(f)
#     f = np.sqrt(Omega2) / (2 * pi)
#     for i in np.arange(1,np.amin(8,Mr.shape[1-1])+1).reshape(-1):
#         print('f%d=%8.3f  Rayleigh Ratio=%.5f\n' % (i,f(i),(f(i) / fc) ** 2))


