import numpy as np

from welib.system.eva import eig


def CraigBampton(MM, KK, Ileader, nModesCB=None, Ifollow=None, F=None, DD=None, fullModesOut=False, discardIm=True): 
    """
    Performs the CraigBampton (CB) reduction of a system given some input master dofs index
    and a number of modes. Reduced matrices, and Guyan and Craig-Bampton modes are returned.
        
    INPUTS
      Ileader : index of leader DOFs
      nModesCB: number of CB modes to keep
      MM, KK  : Maff and stiffness matrix
        
    INPUTS (Optional)
      nModesCB: number of CB modes to keep. Default: all
      Ifollow: indices of follower DOFs. Default: complementary set to Ileader
      fullModesOut: if true, the Guyan and CB modes
      discardIm: if true, the imaginary part of the eigenvectors is discarded
        
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
    Phi_CB, Lambda_CB = eig(Kff,Mff, discardIm=discardIm)
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

    # --- Guyan frequencies
    Phi_G2, Lambda_G = eig(Kr11,Mr11, discardIm=discardIm)
    Omega2 = np.diag(Lambda_G).copy()
    Omega2[Omega2<0]=0.0
    f_G  = np.sqrt(Omega2)/(2*np.pi)

    # Building reduced matrix 
    Mr = np.block( [ [Mr11 , Mr12 ], [ Mr12.T, np.eye(nModesCB)       ] ])
    Kr = np.block( [ [Kr11  , ZZ  ], [ ZZ.T  ,  Lambda_CB[:nModesCB,:]] ])

    # --- Augmenting modes so that they have the same dimension as MM
    # Add "1" for Guyan modes, and "0" for CB modes
    if fullModesOut:
        Phi_G, Phi_CB = augmentModes(Ileader, Phi_G, Phi_CB, Ifollow=Ifollow)

    if DD is not None:
        raise NotImplementedError('Not done')
    if F is not None:
        raise NotImplementedError('Not done')

    I_G  = list(np.arange(len(Ileader)))
    I_CB = list(np.arange(len(Ifollow)) + len(I_G))

    return Mr, Kr, Phi_G, Phi_CB, f_G, f_CB, I_G, I_CB


def augmentModes(Ileader, Phi_G, Phi_CB, Ifollow=None):
    """ 
    Augment Guyan and Craig Bampton modes, so as to return full DOF vectors
    going back to the original size
    """
    # --- Augment modes so that they go back to same size after BC
    nl   = len(Ileader)
    nall = nl+Phi_G.shape[0]
    nf   = nall-nl
    if Ifollow is None:
        Iall    = np.arange(nall)
        Ifollow = list(np.setdiff1d(Iall, Ileader))
    # Guyan
    Phi_G_aug = np.zeros((nall, nl))
    Phi_G_aug[Ileader,:] = np.eye(nl)
    Phi_G_aug[Ifollow,:] = Phi_G
    # 
    Phi_CB_aug = np.zeros((nall, Phi_CB.shape[1]))
    Phi_CB_aug[Ileader,:] = 0
    Phi_CB_aug[Ifollow,:] = Phi_CB

    return Phi_G_aug, Phi_CB_aug



if __name__=='__main__':
    np.set_printoptions(linewidth=500)
    L = 100
    EI = 1868211939147.334
    Maff = L * 8828.201296825122
    KK = EI / (L ** 3) * np.array([[12,6 * L,- 12,6 * L],[6 * L,4 * L ** 2,- 6 * L,2 * L ** 2],[- 12,- 6 * L,12,- 6 * L],[6 * L,2 * L ** 2,- 6 * L,4 * L ** 2]])
    MM = Maff / 420 * np.array([[156,22 * L,54,- 13 * L],[22 * L,4 * L ** 2,13 * L,- 3 * L ** 2],[54,13 * L,156,- 22 * L],[- 13 * L,- 3 * L ** 2,- 22 * L,4 * L ** 2]])
    print(MM)
    Mr,Kr,Phi_G,Phi_CB,f_CB,f_G = CraigBampton(MM,KK,[2], nModesCB=2)
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


