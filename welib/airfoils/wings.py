"""
Aerodynamic coefficients for typical wings

References:
  [1] Tuck (1993):  Some accurate solutions of the lifting surface integral equation. J. Austral. Math. Soc. Ser. B 35127–144.
  [2] Katz, Plotkin (2001): A. Low-Speed Aerodynamics, 2nd Edition. Cambridge University Press
  [3] Branlard (2011): Wind turbine tip-loss corrections, review, implementation and investigation of new models. Appendix B. Figure B.1 
  [4] Branlard (2017): Wind turbine aerodynamics and vorticity based method, Springer. Chapter 45. Figure 45.3

"""
import numpy as np



def CLalpha_models(AR, wingtype, model):
    """ Lift slope for different wing types according to different model 
    INPUTS:
     - AR: aspect ratio of the wing (array like)
     - wingtype: string for the wing type, in ['rectangular_flatplate']
     - model: string for the model to use:
          if 'rectangular_flatplate', model in ['prandtl, 'helmbold', linear', 'tuck', 'jones', 'exp1']
    OUTPUTS:
     - AR (

    Example:
        AR    = np.linspace(0,15,100)
        AR, CLAp = CLalpha_models(AR  , 'rectangular_flatplate', 'Prandtl')

    """
    AR = np.asarray(AR)

    if wingtype.lower()=='rectangular_flatplate':
        if model.lower()=='prandtl':
            # Prandtl  [3]
            # TODO
            a0=5.989;
            #a0=2*np.pi
            xp=AR
            tau=0.0000055;
            Clalpha=a0/(1+a0/(np.pi*xp)*(1+tau));

        elif model.lower()=='helmbold':
            # Helmbold 1942 [3]
            a0=2*np.pi
            xh=AR
            Clalpha=a0/(np.sqrt(1+(a0/(np.pi*xh))**2)+ a0/(np.pi*xh));
        elif model.lower()=='linear':
            # Linear slope
            Clalpha=0.5*np.pi*AR

        elif model.lower()=='tuck':
            # Tuck 1993 [1]
            Clalpha=2*np.pi-np.pi*1/AR*(np.log(AR)+2.5620)+1.404*(1/AR)**2 *(np.log(AR)+3.645)
        elif model.lower()=='jones':
            # TODO I might have extracted that from a plot
            # [3]
            Jones =np.array([
                [1.1211,1.4700,1.8628,2.2775,2.6050,2.9435,3.3912,3.9700,4.6692,5.3248,5.7620,6.2867,6.7348,7.3252],
                [1.7463,2.0740,2.3799,2.6967,2.8935,3.0794,3.2981,3.5825,3.8344,4.0426,4.1414,4.2512,4.3392,4.4382]])
            AR      = Jones[0,:]+0.05
            Clalpha = Jones[1,:]
        elif model.lower()=='exp1':
            # [3]
            FlatPlateExp = np.array([
               [1.9937,1.9823,1.9941,2.6271,3.0199,3.0087,4.5055,5.0077,5.1065,4.9753,5.9918,6.0023,5.9692,0.5225,0.7067,0.9574,0.9577,0.9570],
               [2.5109,2.6526,2.3801,2.8063,3.1231,3.2103,3.7251,4.0203,3.9223,3.8785,4.1090,4.2508,4.3598,0.7968,1.3313,1.6261,1.5171,1.7461]])
            AR      = FlatPlateExp[0,:]+0.05
            Clalpha = FlatPlateExp[1,:]

        else:
            raise NotImplementedError('model', model)

    else:
        raise NotImplementedError('wingtype', wingtype)

    return AR,Clalpha



if __name__ == '__main__':
    import matplotlib.pyplot as plt
    AR = np.linspace(0,15,100)
    ARL = np.linspace(0,3,2)
    _, CLAp = CLalpha_models(AR  , 'rectangular_flatplate', 'Prandtl')
    _, CLAP = CLalpha_models(AR  , 'rectangular_flatplate', 'Prandtl')
    _, CLAH = CLalpha_models(AR  , 'rectangular_flatplate', 'Helmbold')
    _, CLAT = CLalpha_models(AR  , 'rectangular_flatplate', 'Tuck')
    _  , CLAL = CLalpha_models(ARL , 'rectangular_flatplate', 'linear')
    ARJ, CLAJ = CLalpha_models(AR  , 'rectangular_flatplate', 'Jones')
    ARE, CLAE = CLalpha_models(AR  , 'rectangular_flatplate', 'exp1')

    fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
    fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)

    ax.plot(ARL, CLAL , '-'   , label='Linear')
    ax.plot(AR , CLAP , '-.'   , label='Prandtl')
    ax.plot(AR , CLAH , '-'   , label='Helmbold')
    ax.plot(AR , CLAT , '--'   , label='Tuck')
    ax.plot(ARJ, CLAJ , '--'   , label='Jones')
    ax.plot(ARE, CLAE , '^'   , label='Exp.')
    ax.set_xlabel(r'AR [-]')
    ax.set_ylabel(r'$C_{l,\alpha}$ [-]')
    ax.set_ylim([0,6])
    ax.set_xlim([0,10])
    ax.legend()
    plt.show()
