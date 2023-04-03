""" 

COORDINATE SYSTEMS:
     NOTE:  xIEC, yIEC, zIEC =  xED, -zED, yED

     | ElastoDyn | IEC | 
     |   z*      |  i  | inertial
     |   a*      |  t  | tower
     |   t*      |  te | tower elements
     |   b*      |  p  | tower top
     |   d*      |  n  | nacelle (including yaw)
     |   c*      |  s  | non rotating shaft
     |   e*      |  a  | rotating shaft
     |   f*      |     | teetered
     |   g*      |  h  | hub (including delta 3 for 2-bladed
     |   g*'     |     | hub coordinate aligned with a given blade (rotated with azimuth between blade)
     |   i*      |  cK | coned for blade K    (rotatged with cone angle)
     |   j*      |  bK | pitched for blade K  (rotated with pitch angle)
     |   Lj*     |     |                      (rotated with structural twist)


BODIES:
     | ElastoDyn | IEC | 
     |    E      |  E  | : earth/inertial frame
     |    X      |  F  | : platform body
     |    N      |  N  | : nacelle body
     |    A      |     | : tail-furl body

POINTS:
     | ElastoDyn | IEC  | 
     |    Z      |  F   | : platform reference 
     |    Y      |  Gf  | : platform COG
     |    T0     |  T   | : tower-bottom
     |    O      |  N   | : tower-top / base-plate 
     |    U      |  Gn  | : Nacelle COG
     |    V      |      | : Rotor furl axis
     |    D      |      | : COG of structure that furls (excluding the rotor)
     |   IMU     |      | : Nacelle IMU
     |    P      |      | : Teeter pin
     |    Q      |      | : Apex of rotation, rotor center
     |    C      |  Gh  | : Hub COG
     |    W      |      | : specified point on the tail-furl axis
     |    I      |      | : tail boom COG
     |    J      |      | : tail fin COG


"""

import numpy as np
import os
import matplotlib.pyplot as plt
from welib.weio.fast_input_file import FASTInputFile
from welib.system.eva import eigMCK

# --------------------------------------------------------------------------------}
# --- GLOBAL CONSTANTS 
# --------------------------------------------------------------------------------{
#    INTEGER(IntKi), PARAMETER        :: MaxBl    =  3                                   ! Maximum number of blades allowed in simulation
#    INTEGER(IntKi), PARAMETER        :: NumBE    =  1                                   ! Number of blade-edge modes
#    INTEGER(IntKi), PARAMETER        :: NumBF    =  2                                   ! Number of blade-flap modes
DOF_Sg   =  0 # DOF index for platform surge
DOF_Sw   =  1 # DOF index for platform sway
DOF_Hv   =  2 # DOF index for platform heave
DOF_R    =  3 # DOF index for platform roll
DOF_P    =  4 # DOF index for platform pitch
DOF_Y    =  5 # DOF index for platform yaw
DOF_TFA1 =  6 # DOF index for 1st tower fore-aft mode
DOF_TSS1 =  7 # DOF index for 1st tower side-to-side mode
DOF_TFA2 =  8 # DOF index for 2nd tower fore-aft mode
DOF_TSS2 =  9 # DOF index for 2nd tower side-to-side mode
DOF_Yaw  = 10 # DOF index for nacelle-yaw
DOF_RFrl = 11 # DOF index for rotor-furl
DOF_GeAz = 12 # DOF index for the generator azimuth
DOF_DrTr = 13 # DOF index for drivetrain rotational-flexibility
DOF_TFrl = 14 # DOF index for tail-furl
#    INTEGER(IntKi), PARAMETER        :: DOF_BE (MaxBl,NumBE) = RESHAPE(  &              ! DOF indices for blade edge:
#                                                (/ 17, 20, 23 /),   (/MaxBl,NumBE/) )   !    1st blade edge mode for blades 1,2, and 3, respectively 17 + 3*(K-1)
#    INTEGER(IntKi), PARAMETER        :: DOF_BF (MaxBl,NumBF) = RESHAPE(  &              ! DOF indices for blade flap:
#                                                (/ 16, 19, 22,           &              !    1st blade flap mode for blades 1,2, and 3, respectively 16 + 3*(K-1)
#                                                   18, 21, 24 /),   (/MaxBl,NumBF/) )   !    2nd blade flap mode for blades 1,2, and 3, respectively 18 + 3*(K-1)
#    INTEGER(IntKi), PARAMETER        :: DOF_Teet = 22 !DOF_TFrl + 2*(NumBE+NumBF)+ 1    ! DOF index for rotor-teeter
ED_MaxDOFs  = 24
#    INTEGER(IntKi), PARAMETER        :: NPA      =  9                                   ! Number of DOFs that contribute to the angular velocity of the tail (body A) in the inertia frame.
#    INTEGER(IntKi), PARAMETER        :: NPB      =  7                                   ! Number of DOFs that contribute to the angular velocity of the tower top / baseplate (body B) in the inertia frame.
#    INTEGER(IntKi), PARAMETER        :: NPF      =  7                                   ! Number of DOFs that contribute to the angular velocity of the tower elements (body F) in the inertia frame                                           (body F) in the inertia frame.
#    INTEGER(IntKi), PARAMETER        :: NPG      = 10                                   ! Number of DOFs that contribute to the angular velocity of the generator (body G) in the inertia frame.
#    INTEGER(IntKi), PARAMETER        :: NPL      = 11                                   ! Number of DOFs that contribute to the angular velocity of the low-speed shaft (body L) in the inertia frame.
NPN      =  8                      # Number of DOFs that contribute to the angular velocity of the nacelle (body N) in the inertia frame.
#    INTEGER(IntKi), PARAMETER        :: NPR      =  9                                   ! Number of DOFs that contribute to the angular velocity of the structure that furls with the rotor (not including rotor) (body R) in the inertia frame.
NPX      =  3                      # Number of DOFs that contribute to the angular velocity of the platform (body X) in the inertia frame.
PX  = [ DOF_R, DOF_P, DOF_Y ] # Array of DOF indices (pointers) that contribute to the angular velocity of the platform                                                  (body X) in the inertia frame.
#    INTEGER(IntKi), PARAMETER        :: PF(NPF)  = (/ DOF_R, DOF_P, DOF_Y, DOF_TFA1, DOF_TSS1, DOF_TFA2, DOF_TSS2 /)                                                  ! Array of DOF indices (pointers) that contribute to the angular velocity of the tower elements                                            (body F) in the inertia frame.
#    INTEGER(IntKi), PARAMETER        :: PB(NPB)  = (/ DOF_R, DOF_P, DOF_Y, DOF_TFA1, DOF_TSS1, DOF_TFA2, DOF_TSS2 /)                                                  ! Array of DOF indices (pointers) that contribute to the angular velocity of the tower top / baseplate                                     (body B) in the inertia frame.
PN  = [ DOF_R, DOF_P, DOF_Y, DOF_TFA1, DOF_TSS1, DOF_TFA2, DOF_TSS2, DOF_Yaw ]  # Array of DOF indices (pointers) that contribute to the angular velocity of the nacelle                                                   (body N) in the inertia frame.
#    INTEGER(IntKi), PARAMETER        :: PR(NPR)  = (/ DOF_R, DOF_P, DOF_Y, DOF_TFA1, DOF_TSS1, DOF_TFA2, DOF_TSS2, DOF_Yaw, DOF_RFrl /)                               ! Array of DOF indices (pointers) that contribute to the angular velocity of the structure that furls with the rotor (not including rotor) (body R) in the inertia frame.
#    INTEGER(IntKi), PARAMETER        :: PL(NPL)  = (/ DOF_R, DOF_P, DOF_Y, DOF_TFA1, DOF_TSS1, DOF_TFA2, DOF_TSS2, DOF_Yaw, DOF_RFrl, DOF_GeAz, DOF_DrTr /)           ! Array of DOF indices (pointers) that contribute to the angular velocity of the low-speed shaft                                           (body L) in the inertia frame.
#    INTEGER(IntKi), PARAMETER        :: PG(NPG)  = (/ DOF_R, DOF_P, DOF_Y, DOF_TFA1, DOF_TSS1, DOF_TFA2, DOF_TSS2, DOF_Yaw, DOF_RFrl, DOF_GeAz /)                     ! Array of DOF indices (pointers) that contribute to the angular velocity of the generator                                                 (body G) in the inertia frame.
#    INTEGER(IntKi), PARAMETER        :: PA(NPA)  = (/ DOF_R, DOF_P, DOF_Y, DOF_TFA1, DOF_TSS1, DOF_TFA2, DOF_TSS2, DOF_Yaw, DOF_TFrl /)                               ! Array of DOF indices (pointers) that contribute to the angular velocity of the tail                                                      (body A) in the inertia frame.



def fitShapeFunction(x_bar, phi, exp=None, scale=True, plot=None):
    """ 
    Return polynomial fit for a given shapefunction
    The fit is such that phi_fit = a_i x_bar^e_i

    See also: from welib.yams.flexibility.polyshape
    INPUTS: 
      - x : dimensionless spanwise coordinate, from 0 to 1
            The points 0 and 1 need not be present.
      - exp: exponents of the polynomial. Should be length of coeff. 
            If None, exp = [2,3,4,5,6] as used in OpenFAST
      - scale: if True, scale the coefficients such that the sum is 1
    OUTPUTS:
      - pfit: fitted coefficients: [a_i] such that phi_fit = a_i x_bar^e_i
      - y_fit: fitted values phi_fit(x_bar)
      - fig: figure handle
      - fitter: object with dictionary fields 'coeffs', 'formula', 'fitted_function'
    """
    if exp is None:
        exp = np.arange(2,7)
    if  np.any(x_bar)<0 or np.any(x_bar)>1:
        raise Exception('`x_bar` should be between 0 and 1')

    from welib.tools.curve_fitting import model_fit
    phi_fit, pfit, fitter = model_fit('fitter: polynomial_discrete', x_bar, phi, exponents=exp)

    # --- Manipulation of the fitted model...
    pfit = np.around(pfit, 8) # sticking to 8 significant digits

    if scale:
        scale = np.sum(pfit)
        if np.abs(scale)<1e-8:
            print('[WARN] ElastoDyn: fitShapeFunction: Problem, sum of coefficient is close to 0')
        else:
            pfit = np.array(pfit)/scale
            pfit = np.around(pfit, 8) # sticking to 8 significant digits. NOTE: this might mess up the scale again..
    def fitted_function(xx):
        y=np.zeros(xx.shape)
        for i,(e,c) in enumerate(zip(exp, pfit)):
            y += c * xx**e
        return y
    phi_fit = fitted_function(x_bar)

    if plot:
        fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
        fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
        ax.plot(x_bar, phi          , label='Input Shape')
        ax.plot(x_bar, phi_fit, '--', label='Fitted Shape')
        ax.set_xlabel('x/L [-]')
        ax.set_ylabel('Shape function [-]')
        ax.legend(loc='upper left')
    else:
        fig=None


    return pfit, phi_fit, fig


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


def SHP(Fract, FlexL, ModShpAry, Deriv):
    """ SHP calculates the Derive-derivative of the shape function ModShpAry at Fract.
    NOTES: This function only works for Deriv = 0, 1, or 2. 
           Taken from ElastoDyn.f90
    """
    Swtch        = np.zeros((3, 1)); # Initialize Swtch(:) to 0
    Swtch[Deriv] = 1;
    shp          = 0.0;
    if Deriv==0:
        for i in np.arange(len(ModShpAry)):
            shp = shp + ModShpAry[i]*( Fract**(i+2) )
    else:
        for i in np.arange(len(ModShpAry)):
            I = i + 1
            J = I + 1;
            CoefTmp = Swtch[0] + Swtch[1]*J + Swtch[2]*I*J;
            if ( (J == 2) and (Deriv == 2) ):
                shp =       ModShpAry[i]*CoefTmp                         /( FlexL**Deriv );
            else:
                shp = shp + ModShpAry[i]*CoefTmp*( Fract**( J - Deriv ) )/( FlexL**Deriv );
    return shp


def getEDClass(class_or_filename):
    """
    Return ElastoDyn instance of FileCl
    INPUT: either
       - an instance of FileCl, as returned by reading the file, ED = weio.read(ED_filename)
       - a filepath to a ElastoDyn input file
       - a filepath to a main OpenFAST input file
    """
    if hasattr(class_or_filename,'startswith'): # if string
        ED = FASTInputFile(class_or_filename)
        if 'EDFile' in ED.keys(): # User provided a .fst file...
            parentDir=os.path.dirname(class_or_filename)
            EDfilename = os.path.join(parentDir, ED['EDFile'].replace('"',''))
            ED = FASTInputFile(EDfilename)
    else:
        ED = class_or_filename
    return ED



def ED_Parameters(fstFilename, split=False):

    fst = FASTInputFile(fstFilename)
    if 'EDFile' in fst.keys(): # User provided a .fst file...
        parentDir=os.path.dirname(fstFilename)
        EDfilename = os.path.join(parentDir, fst['EDFile'].replace('"',''))
        ED = getEDClass(EDfilename)
    else:
        raise NotImplementedError()

    prot, pbld, phub = rotorParameters(ED)
    ptwr = towerParameters(ED, gravity=fst['Gravity'], RotMass=prot['RotMass'])

    #pbld[0] = bladeDerivedParameters(pbld[0], inertiaAtBladeRoot=True)
    #ptwr = towerDerivedParameters(ptwr)

    # --- Platform
    pptfm={}
    pptfm['PtfmRefzt'] = ED['PtfmRefzt']
    pptfm['PtfmCMxt'] = ED['PtfmCMxt']
    pptfm['PtfmCMyt'] = ED['PtfmCMyt']


    pmisc = {}
    pmisc['rZT0zt']    = ED['TowerBsHt'] - ED['PtfmRefzt'] # zt-component of position vector rZT0.
    pmisc['RefTwrHt']  = ED['TowerHt']   - ED['PtfmRefzt'] # Vertical distance between ElastoDyn's undisplaced tower height (variable TowerHt) and ElastoDyn's inertia frame reference point (variable PtfmRef).
    pmisc['rZYzt']     = ED['PtfmCMzt'] - ED['PtfmRefzt']


    # --- Nacelle
    pnac = {}
    pnac['NacCMxn'] = ED['NacCMxn']
    pnac['NacCMyn'] = ED['NacCMyn']
    pnac['NacCMzn'] = ED['NacCMzn']

    # --- Generator
    pgen={}

    # --- Shaft
    psft={}

    # Calculate the turbine mass:
    #ptwr['TurbMass']  = ptwr['TwrTpMass'] + ptwr['TwrMass'];
    p=dict()
    if not split:
        p.update(pmisc)
        p.update(pptfm)
        p.update(ptwr)
        p.update(pnac)
        p.update(pgen)
        p.update(psft)
        p.update(phub)
        p.update(prot)
        p.update(pbld[0])
        return p
    else:
        p['misc'] = pmisc
        p['ptfm'] = pptfm
        p['twr']  = ptwr
        p['nac']  = pnac
        p['gen']  = pgen
        p['sft']  = psft
        p['hub']  = phub
        p['rot']  = prot
        p['bld']  = pbld
        return p




def rotorParameters(EDfilename, identicalBlades=True, pbld1=None):
    """ 
    Return rotor parameters computed like ElastoDyn
    p: rotor parametesr
    pbld: list of parameters for each blades.
    """
    ED = getEDClass(EDfilename)

    if pbld1 is None:
        if identicalBlades:
            pbld = [bladeParameters(EDfilename, 1)]*ED['NumBl']
        else:
            pbld=[bladeParameters(EDfilename, ibld+1) for ibld in  range(ED['NumBl'])]
    else:
        pbld = [pbld1]*ED['NumBl']

    p=dict()
    p['RotMass'] = sum([pbld[k]['BldMass'] for k in range(ED['NumBl'])])
    p['RotIner'] = sum([(pbld[k]['SecondMom'] + pbld[k]['BldMass']*ED['HubRad']*(2.0*pbld[k]['BldCG'] + ED['HubRad']))*(np.cos(pbld[k]['PreCone'])**2) for k in range(ED['NumBl'])])
    #if ( p.NumBl == 2 )   % 2-blader
    #    p.Hubg1Iner = ( InputFileData.HubIner - p.HubMass*( ( p.UndSling - p.HubCM )^2 ) )/( p.CosDel3^2 );
    #    p.Hubg2Iner = p.Hubg1Iner;
    #else                    % 3-blader
    #    p.Hubg1Iner = GetFASTPar(edDataOut, 'HubIner');
    #    p.Hubg2Iner = 0.0;
    #end
    phub=dict()
    phub['HubMass'] = ED['HubMass']
    phub['HubIner'] = ED['HubIner']

    p['RotMass'] += phub['HubMass']
    p['RotIner'] += phub['HubIner']

    p['TwrTpMass'] = p['RotMass'] + ED['NacMass'] + ED['YawBrMass'];

    return p, pbld, phub


def bladeParameters(EDfilename, ibld=1, RotSpeed=1, AdjBlMs=None):
    """
    Compute blade parameters in a way similar to OpenFAST
    See Routine Coeff from ElastoDyn.f90
    RotSpeed: used for rotational stiffening. Use 1 for unit contribution (proportioanl to omega**2) [rad/s]
    """
    from welib.yams.flexibility import polyshape
    from welib.yams.flexibility import GMBeam, GKBeam
    # --- Read inputs
    ED = getEDClass(EDfilename)
    try:
        EDbld   = os.path.join(os.path.dirname(ED.filename), ED['BldFile({})'.format(ibld)].replace('"',''))
    except:
        EDbld   = os.path.join(os.path.dirname(ED.filename), ED['BldFile{}'.format(ibld)].replace('"',''))
    bld     = FASTInputFile(EDbld)
    bldProp = bld.toDataFrame()

    # --- 
    p=dict()
    p['HubRad']   = ED['HubRad']
    p['BldNodes'] = ED['BldNodes']
    p['BldFlexL'] = ED['TipRad']- ED['HubRad'] # Length of the flexible portion of the blade.
    n=ED['BldNodes']
    p['DRNodes'] = np.ones(ED['BldNodes'])*p['BldFlexL']/ED['BldNodes']
    bld_fract    = np.arange(1./ED['BldNodes']/2., 1, 1./ED['BldNodes'])
    p['RNodes'] = bld_fract*p['BldFlexL']

    # Adjust mass
    if AdjBlMs is None:
        p['AdjBlMs'] = bld['AdjBlMs']
    else:
        p['AdjBlMs'] = AdjBlMs
    bldProp['BMassDen_[kg/m]'] *= p['AdjBlMs']


    # --- Interpolate the blade properties to this discretization:
    p['RNodesNorm'] = p['RNodes']/p['BldFlexL'];  # Normalized radius to analysis nodes relative to hub ( -1 < RNodesNorm(:) < 1 )
    p['Bl_s_span'] = np.concatenate(([0], p['RNodesNorm'], [1]))*p['BldFlexL'];  # Normalized radius to analysis nodes relative to hub ( -1 < RNodesNorm(:) < 1 )
    p['Bl_s_span_norm'] = p['Bl_s_span']/p['BldFlexL']
    p['BlFract']= bldProp['BlFract_[-]'].values
    StrcTwst    = bldProp['StrcTwst_[deg]'].values
    p['ThetaS']  = np.interp(p['RNodesNorm'], p['BlFract'], bldProp['StrcTwst_[deg]']) 
    p['MassB']   = np.interp(p['RNodesNorm'], p['BlFract'], bldProp['BMassDen_[kg/m]'])                     ;
    p['StiffBF'] = np.interp(p['RNodesNorm'], p['BlFract'], bldProp['FlpStff_[Nm^2]'])
    p['StiffBE'] = np.interp(p['RNodesNorm'], p['BlFract'], bldProp['EdgStff_[Nm^2]'])
    p['m_full']    = np.interp(p['Bl_s_span_norm'], p['BlFract'], bldProp['BMassDen_[kg/m]'])                     ;
    p['EI_F_full'] = np.interp(p['Bl_s_span_norm'], p['BlFract'], bldProp['FlpStff_[Nm^2]'])
    p['EI_E_full'] = np.interp(p['Bl_s_span_norm'], p['BlFract'], bldProp['EdgStff_[Nm^2]'])
    p['ThetaS']  = np.concatenate( ([StrcTwst[0]], p['ThetaS'] , [StrcTwst[-1]]) )
    # Set the blade damping and stiffness tuner
    try:
        p['BldFDamp'] = [bld['BldFlDmp(1)'], bld['BldFlDmp(2)'] ]
        p['BldEDamp'] = [bld['BldEdDmp(1)']]
        p['FStTunr']  = [bld['FlStTunr(1)'], bld['FlStTunr(2)'] ]
    except:
        p['BldFDamp'] = [bld['BldFlDmp1'], bld['BldFlDmp2'] ]
        p['BldEDamp'] = [bld['BldEdDmp1']]
        p['FStTunr']  = [bld['FlStTunr1'], bld['FlStTunr2'] ]
    # Set the mode shape coefficients 
    p['BldFl1Sh'] = [bld[c] for c in ['BldFl1Sh(2)', 'BldFl1Sh(3)', 'BldFl1Sh(4)', 'BldFl1Sh(5)', 'BldFl1Sh(6)']]
    p['BldFl2Sh'] = [bld[c] for c in ['BldFl2Sh(2)', 'BldFl2Sh(3)', 'BldFl2Sh(4)', 'BldFl2Sh(5)', 'BldFl2Sh(6)']]
    p['BldEdgSh'] = [bld[c] for c in ['BldEdgSh(2)', 'BldEdgSh(3)', 'BldEdgSh(4)', 'BldEdgSh(5)', 'BldEdgSh(6)']]
    p['CThetaS'] = np.cos(p['ThetaS']*np.pi/180);
    p['SThetaS'] = np.sin(p['ThetaS']*np.pi/180);

    p['PreCone']= ED['PreCone({:d})'.format(ibld)]*np.pi/180

    # --- Inertial properties
    # Initialize BldMass(), FirstMom(), and SecondMom() using TipMass() effects
    p['TipMass']   = ED['TipMass({:d})'.format(ibld)]                            
    p['BElmntMass'] = p['MassB']*p['DRNodes'] # Mass of blade element
    p['BldMass']   = sum(p['BElmntMass'])                + p['TipMass']
    p['FirstMom']  = sum(p['BElmntMass']*p['RNodes'])    + p['TipMass']*p['BldFlexL']               # wrt blade root    
    p['SecondMom'] = sum(p['BElmntMass']*p['RNodes']**2) + p['TipMass']*p['BldFlexL']*p['BldFlexL'] # wrt blade root
    p['FMomAbvNd']  = np.zeros(n)
    # Integrate to find FMomAbvNd:
    for J in  np.arange(ED['BldNodes']-1,-1,-1): # Loop through the blade nodes / elements in reverse
        p['FMomAbvNd'][J] = (0.5*p['BElmntMass'][J] )*(ED['HubRad'] + p['RNodes'][J] + 0.5*p['DRNodes'][J])
        if J == n-1: # Outermost blade element
           p['FMomAbvNd'][J] += p['TipMass'] * ED['TipRad']; # TipMass effects:
        else:  
           # Add to p['FMomAbvNd(K,J) the effects from the (not yet used) portion of element J+1
           p['FMomAbvNd'][J] += p['FMomAbvNd'][J+1] + (0.5*p['BElmntMass'][J+1])*( ED['HubRad'] + p['RNodes'][J+1] - 0.5*p['DRNodes'][J+1] );
    # Calculate BldCG() using FirstMom() and BldMass(); and calculate RotMass and RotIner:
    p['BldCG']= p['FirstMom']/p['BldMass'];
    p['MBF']       = np.zeros((2, 2))
    p['MBE']       = np.zeros((1, 1))
    p['KBFCent']   = np.zeros((2, 2))
    p['KBECent']   = np.zeros((1, 1))
    p['KBF']       = np.zeros((2, 2))
    p['KBE']       = np.zeros((1, 1))
    # Initialize the generalized blade masses using tip mass effects:
    p['MBF'][0,0] = p['TipMass'];
    p['MBF'][1,1] = p['TipMass'];
    p['MBE'][0,0] = p['TipMass'];

    # Shape functions and derivatives at all nodes and tip&root (2 additional points)
    exp = np.arange(2,7)
    p['ShapeF1_full'], p['dShapeF1_full'],p['ddShapeF1_full'] = polyshape(p['Bl_s_span'], coeff=p['BldFl1Sh'], exp=exp, x_max=p['BldFlexL'], doscale=False)
    p['ShapeF2_full'], p['dShapeF2_full'],p['ddShapeF2_full'] = polyshape(p['Bl_s_span'], coeff=p['BldFl2Sh'], exp=exp, x_max=p['BldFlexL'], doscale=False)
    p['ShapeE1_full'], p['dShapeE1_full'],p['ddShapeE1_full'] = polyshape(p['Bl_s_span'], coeff=p['BldEdgSh'], exp=exp, x_max=p['BldFlexL'], doscale=False)

    # Integrate to find the generalized mass of the blade (including tip mass effects).
    #   Ignore the cross-correlation terms of MBF (i.e. MBF(i,j) where i ~= j) since these terms will never be used.
    p['MBF'][0,0] = sum(p['BElmntMass']*p['ShapeF1_full'][1:-1]**2)
    p['MBF'][1,1] = sum(p['BElmntMass']*p['ShapeF2_full'][1:-1]**2)
    p['MBE'][0,0] = sum(p['BElmntMass']*p['ShapeE1_full'][1:-1]**2)

    ElmntStff      = p['StiffBF']*p['DRNodes']  # Flapwise stiffness of blade element J
    p['KBF'][0,0] = sum(ElmntStff*p['ddShapeF1_full'][1:-1]*p['ddShapeF1_full'][1:-1])
    p['KBF'][0,1] = sum(ElmntStff*p['ddShapeF1_full'][1:-1]*p['ddShapeF2_full'][1:-1])
    p['KBF'][1,0] = sum(ElmntStff*p['ddShapeF2_full'][1:-1]*p['ddShapeF1_full'][1:-1])
    p['KBF'][1,1] = sum(ElmntStff*p['ddShapeF2_full'][1:-1]*p['ddShapeF2_full'][1:-1])
    ElmntStff     = p['StiffBE']*p['DRNodes'] # Edgewise stiffness of blade element J
    p['KBE'][0,0] = sum(ElmntStff*p['ddShapeE1_full'][1:-1]*p['ddShapeE1_full'][1:-1])

    # Integrate to find the centrifugal-term of the generalized flapwise and edgewise
    #   stiffness of the blades.  Ignore the cross-correlation terms of KBFCent (i.e.
    #   KBFCent(i,j) where i ~= j) since these terms will never be used.
    ElmntStff      = p['FMomAbvNd']*p['DRNodes']*RotSpeed**2 # Centrifugal stiffness of blade element J
    p['KBFCent'][0,0] = sum(ElmntStff*p['dShapeF1_full'][1:-1]**2)
    p['KBFCent'][1,1] = sum(ElmntStff*p['dShapeF2_full'][1:-1]**2)
    p['KBECent'][0,0] = sum(ElmntStff*p['dShapeE1_full'][1:-1]**2)

    # Calculate the 2nd derivatives of the twisted shape functions (include root and tip):
    p['TwistedSF'] = np.zeros((2, 3, ED['BldNodes']+2, 3)); # x/y, BF1/BF2/BE, node, deriv
    p['TwistedSF'][0,0,:,2] =  p['ddShapeF1_full'][:]*p['CThetaS'][:] # 2nd deriv. of Phi1(J) for blade K
    p['TwistedSF'][1,0,:,2] = -p['ddShapeF1_full'][:]*p['SThetaS'][:] # 2nd deriv. of Psi1(J) for blade K
    p['TwistedSF'][0,1,:,2] =  p['ddShapeF2_full'][:]*p['CThetaS'][:] # 2nd deriv. of Phi2(J) for blade K
    p['TwistedSF'][1,1,:,2] = -p['ddShapeF2_full'][:]*p['SThetaS'][:] # 2nd deriv. of Psi2(J) for blade K
    p['TwistedSF'][0,2,:,2] =  p['ddShapeE1_full'][:]*p['SThetaS'][:] # 2nd deriv. of Phi3(J) for blade K
    p['TwistedSF'][1,2,:,2] =  p['ddShapeE1_full'][:]*p['CThetaS'][:] # 2nd deriv. of Psi3(J) for blade K
    # Integrate to find the 1st derivatives of the twisted shape functions:
    for J  in np.arange(n):    # Loop through the blade nodes / elements
        TwstdSF= np.zeros((2, 3, 2))
        order=1;
        for I in [0,1]:   # Loop through Phi and Psi
            for L in [0,1,2]:  # Loop through all blade DOFs
                TwstdSF[I,L,order] = p['TwistedSF'][I,L,J+1,order+1]*0.5*p['DRNodes'][J];
                p['TwistedSF'][I,L,J+1,order] = TwstdSF[I,L,order];
        if  J != 0:    # All but the innermost blade element
            # Add the effects from the (not yet used) portion of element J-1
            for I in [0,1]:  # Loop through Phi and Psi
                for L in [0,1,2]: # Loop through all blade DOFs
                    p['TwistedSF'][I,L,J+1,order] += p['TwistedSF'][I,L,J,order] + TwstdSFOld[I,L,order];
        # Store the TwstdSF and AxRdBld terms of the current element (these will be used for the next element)
#         TwstdSFOld = TwstdSF;
#     for J  in np.arange(n):    # Loop through the blade nodes / elements
#         TwstdSF= np.zeros((2, 3, 2))
        # Integrate to find the twisted shape functions themselves (i.e., their zeroeth derivative):
        order = 0
        for I in [0,1]:   # Loop through Phi and Psi
            for L in [0,1,2]:  # Loop through all blade DOFs
                TwstdSF[I,L, order] = p['TwistedSF'][I,L,J+1, order+1]*0.5*p['DRNodes'][J];
                p['TwistedSF'][I,L,J+1,order] = TwstdSF[ I,L, order ];
        if  J != 0:    # All but the innermost blade element
            # Add the effects from the (not yet used) portion of element J-1
            for I in [0,1]:   # Loop through Phi and Psi
                for L in [0,1,2]:  # Loop through all blade DOFs
                    p['TwistedSF'][I,L,J+1,order] += p['TwistedSF'][I,L,J,order]  + TwstdSFOld[I,L, order];
        TwstdSFOld = TwstdSF;
    # Integrate to find the 1st and zeroeth derivatives of the twisted shape functions at the tip:
    for I in [0,1]:   # Loop through Phi and Psi
        for L in [0,1,2]:  # Loop through all blade DOFs
            p['TwistedSF'][I,L,-1,1] = p['TwistedSF'][I,L,-2,1] + TwstdSFOld[I,L,1]
            p['TwistedSF'][I,L,-1,0] = p['TwistedSF'][I,L,-2,0] + TwstdSFOld[I,L,0]
    # Blade root
    p['TwistedSF'][:,:,0,1] = 0.0;
    p['TwistedSF'][:,:,0,0] = 0.0;

    # Integrate to find the blade axial reduction shape functions:
    p['AxRedBld']  = np.zeros((3, 3, ED['BldNodes']+2)); #
    for J  in np.arange(n):    # Loop through the blade nodes / elements
        AxRdBld= np.zeros((3, 3))
        for I in [0,1,2]:     # Loop through all blade DOFs
            for L in [0,1,2]:  # Loop through all blade DOFs
                AxRdBld[I,L] = 0.5*p['DRNodes'][J]*( p['TwistedSF'][0,I,J+1,1]*p['TwistedSF'][0,L,J+1,1] + p['TwistedSF'][1,I,J+1,1]*p['TwistedSF'][1,L,J+1,1] );
                p['AxRedBld'][I,L,J+1] = AxRdBld[I,L]
        if  J != 0:    # All but the innermost blade element
            # Add the effects from the (not yet used) portion of element J-1
            for I in [0,1,2]:     # Loop through all blade DOFs
                for L in [0,1,2]:  # Loop through all blade DOFs
                    p['AxRedBld'][I,L,J+1] += p['AxRedBld'][I,L,J]  + AxRdBldOld[I,L]
        AxRdBldOld = AxRdBld;
    # Integrate to find the blade axial reduction shape functions at the tip:
    for I in [0,1,2]:     # Loop through all blade DOFs
        for L in [0,1,2]:  # Loop through all blade DOFs
            p['AxRedBld'][I,L,-1] = p['AxRedBld'][I,L,-2] + AxRdBldOld[I,L]
    # Blade root
    p['AxRedBld'] [:,:,0  ] = 0.0;

    # Apply the flapwise modal stiffness tuners of the blades to KBF():
    for I in [0,1]:     # Loop through flap DOFs
        for L in [0,1]:  # Loop through flap DOFs
            p['KBF'][I,L] = np.sqrt( p['FStTunr'][I]*p['FStTunr'][L] )*p['KBF'][I,L];
    # Calculate the blade natural frequencies:
    p['FreqBF'] = np.zeros((2,3))
    p['FreqBE'] = np.zeros((1,3))
    for I in [0,1]:     # Loop through flap DOFs
        p['FreqBF'][I,0] = (1/2/np.pi)*np.sqrt(   p['KBF'][I,I]                     /( p['MBF'][I,I] - p['TipMass']) )# Natural blade I-flap frequency w/o centrifugal stiffening nor     tip mass effects
        p['FreqBF'][I,1] = (1/2/np.pi)*np.sqrt(   p['KBF'][I,I]                     /  p['MBF'][I,I]                 )# Natural blade I-flap frequency w/o centrifugal stiffening, but w/ tip mass effects
        p['FreqBF'][I,2] = (1/2/np.pi)*np.sqrt( ( p['KBF'][I,I] + p['KBFCent'][I,I])/  p['MBF'][I,I]                 )# Natural blade I-flap frequency w/  centrifugal stiffening and     tip mass effects
    I=0
    p['FreqBE'][I,0] =     (1/2/np.pi)*np.sqrt(   p['KBE'][I,I]                     /( p['MBE'][I,I] - p['TipMass']) )# Natural blade 1-edge frequency w/o centrifugal stiffening nor      tip mass effects
    p['FreqBE'][I,1] =     (1/2/np.pi)*np.sqrt(   p['KBE'][I,I]                     /  p['MBE'][I,I]                 )# Natural Blade 1-edge frequency w/o  centrifugal stiffening, but w/ tip mass effects
    p['FreqBE'][I,2] =     (1/2/np.pi)*np.sqrt( ( p['KBE'][I,I] + p['KBECent'][I,I])/  p['MBE'][I,I]                 )# Natural Blade 1-edge frequency w/  centrifugal stiffening and      tip mass effects
    # Calculate the generalized damping of the blades:
    p['CBF'] = np.zeros((2,2))
    p['CBE'] = np.zeros((1,1))
    for I in [0,1]:     # Loop through flap DOFs
        for L in [0,1]:     # Loop through flap DOFs
            p['CBF'][I,L] = ( 0.01*p['BldFDamp'][L] )*p['KBF'][I,L]/( np.pi*p['FreqBF'][L,0] );
    L=0; I=0;
    p['CBE'][I,L] = ( 0.01*p['BldEDamp'][L] )*p['KBE'][I,L]/( np.pi*p['FreqBE'][L,0] );

    nq = 3

    # --- Twisted and untwisted shape functions
    nNodes=n+2
    p['Ut'] = np.zeros((nq, 3, nNodes))
    p['Vt'] = np.zeros((nq, 3, nNodes))
    p['Kt'] = np.zeros((nq, 3, nNodes))
    p['U']  = np.zeros((nq, 3, nNodes))
    p['V']  = np.zeros((nq, 3, nNodes))
    p['K']  = np.zeros((nq, 3, nNodes))
    for j,idir,name in zip(range(0,nq), (0,0,1), ('F1','F2','E1')): # direction is x, x, y
        p['Ut'][j][0,:] = p['TwistedSF'][0, j, :, 0]  # x
        p['Ut'][j][1,:] = p['TwistedSF'][1, j, :, 0]  # y
        p['Vt'][j][0,:] = p['TwistedSF'][0, j, :, 1]  # x
        p['Vt'][j][1,:] = p['TwistedSF'][1, j, :, 1]  # y
        p['Kt'][j][0,:] = p['TwistedSF'][0, j, :, 2]  # x
        p['Kt'][j][1,:] = p['TwistedSF'][1, j, :, 2]  # y
        p['U'][j][idir,:]  = p['Shape'+name+'_full']  
        p['V'][j][idir,:]  = p['dShape'+name+'_full'] 
        p['K'][j][idir,:]  = p['ddShape'+name+'_full']

    # --- Parameters consistent with "YAMS" flexibility module
    p['s_span'] = p['Bl_s_span']
    p['m'] = p['m_full']
    p['EI'] = np.zeros((3,nNodes))
    p['EI'][0,:] = p['EI_F_full']
    p['EI'][1,:] = p['EI_E_full']

    p['s_G0'] = np.zeros((3, len(p['Bl_s_span'])))
    p['s_G0'][2,:] = p['s_span']  # TODO add hub radius

    #KK0 = GKBeam(s_span, EI, PhiK, bOrth=False)
    #if bStiffening:
    #    KKg     = GKBeamStiffnening(s_span, PhiV, gravity, m, Mtop, Omega, main_axis=main_axis)
    #    KKg_self= GKBeamStiffnening(s_span, PhiV, gravity, m, Mtop, Omega, main_axis=main_axis, bSelfWeight=True , bMtop=False, bRot=False)
    #    KKg_Mtop= GKBeamStiffnening(s_span, PhiV, gravity, m, Mtop, Omega, main_axis=main_axis, bSelfWeight=False, bMtop=True,  bRot=False)
    #    KKg_rot = GKBeamStiffnening(s_span, PhiV, gravity, m, Mtop, Omega, main_axis=main_axis, bSelfWeight=False, bMtop=False, bRot=True)
#     MM, IT = GMBeam(s_G0, p['s_span'], p['m_full'], p['Ut'], rot_terms=True, method='OpenFAST', main_axis='z', U_untwisted=p['U']) 
    return p


def bladeDerivedParameters(p, inertiaAtBladeRoot=True):
    """ Compute blade derived parameters, suitable for SID:
      - Inertial matrices
      - Stiffness matrices
    The parameters are computed "by hand" (as opposed to using GMBeam)
    NOTE: this function is mostly for debugging purposes. 
          Use GMBeam with method='OpenFAST' instead
    """

    nq = 3
    if inertiaAtBladeRoot:
        rh = 0
    else:
        rh = p['HubRad'] # Hub Radius # TODO make this an option if from blade root or not

    # --- Rigid mass matrix terms
    # Mxt term 
    p['mdCM'] = np.zeros((3))
    p['mdCM'][2]=  sum(p['BElmntMass'][:]*(p['RNodes']+rh));
    # 
    p['mdCM_M1'] = np.zeros((3,nq))
    for j in np.arange(nq):
        p['mdCM_M1'][0,j]= sum(p['TwistedSF'][0, j, 1:-1, 0]*p['BElmntMass'])
        p['mdCM_M1'][1,j]= sum(p['TwistedSF'][1, j, 1:-1, 0]*p['BElmntMass'])

    p['J']    = np.zeros((3,3))
    p['J'][0,0] = sum(p['BElmntMass'][:]*(p['RNodes']+rh)**2)
    p['J'][1,1] = sum(p['BElmntMass'][:]*(p['RNodes']+rh)**2)
    # --- Elastic matrices
    p['Ke'] = np.zeros((nq,nq))
    p['Ke0'] = np.zeros((nq,nq)) # Without any stiffening
    p['De'] = np.zeros((nq,nq))
    p['Me'] = np.zeros((nq,nq))
    # Me
    p['Me'][0,0]= p['MBF'][0, 0]
    p['Me'][0,1]= p['MBF'][0, 1]
    p['Me'][1,0]= p['MBF'][1, 0]
    p['Me'][1,1]= p['MBF'][1, 1]
    p['Me'][2,2]= p['MBE'][0, 0]
    # Ke
    p['Ke0'][0,0]= p['KBF'][0, 0]
    p['Ke0'][0,1]= p['KBF'][0, 1]
    p['Ke0'][1,0]= p['KBF'][1, 0]
    p['Ke0'][1,1]= p['KBF'][1, 1]
    p['Ke0'][2,2]= p['KBE'][0, 0]
    p['Ke'] = p['Ke0'] # OpenFAST does not put the stiffening in this
    # De
    p['De'][0,0]= p['CBF'][0, 0]
    p['De'][0,1]= p['CBF'][0, 1]
    p['De'][1,0]= p['CBF'][1, 0]
    p['De'][1,1]= p['CBF'][1, 1]
    p['De'][2,2]= p['CBE'][0, 0]

    # KgOm
    p['Kg_Om'] = np.zeros((nq,nq))
    p['Kg_Om'][0,0] = p['KBFCent'][0,0]
    p['Kg_Om'][1,1] = p['KBFCent'][1,1]
    p['Kg_Om'][2,2] = p['KBECent'][0,0]

    # --- Elastic mass matrix terms
    # Ct
    p['Ct'] = np.zeros((nq,3))
    p['Ct'][0,0] = sum(np.squeeze(p['TwistedSF'][0, 0, 1:-1, 0])*p['BElmntMass'][:])# 1st mode acceleration in x
    p['Ct'][0,1] = sum(np.squeeze(p['TwistedSF'][1, 0, 1:-1, 0])*p['BElmntMass'][:])# 1st mode acceleration in y
    p['Ct'][1,0] = sum(np.squeeze(p['TwistedSF'][0, 1, 1:-1, 0])*p['BElmntMass'][:])# 2nd mode acceleration in x
    p['Ct'][1,1] = sum(np.squeeze(p['TwistedSF'][1, 1, 1:-1, 0])*p['BElmntMass'][:])# 2nd mode acceleration in y
    p['Ct'][2,0] = sum(np.squeeze(p['TwistedSF'][0, 2, 1:-1, 0])*p['BElmntMass'][:])# 3nd mode acceleration in x
    p['Ct'][2,1] = sum(np.squeeze(p['TwistedSF'][1, 2, 1:-1, 0])*p['BElmntMass'][:])# 3nd mode acceleration in y

    # Cr
    p['Cr'] = np.zeros((nq,3))
    p['Cr'][0,0]= -sum((p['RNodes']+rh)*np.squeeze(p['TwistedSF'][1, 0, 1:-1, 0])*p['BElmntMass'][:])# 1st mode acceleration about x axis (1) -> movement in negative y
    p['Cr'][0,1]=  sum((p['RNodes']+rh)*np.squeeze(p['TwistedSF'][0, 0, 1:-1, 0])*p['BElmntMass'][:])# 1st mode acceleration about y axis (2) -> movement in x
    p['Cr'][1,0]= -sum((p['RNodes']+rh)*np.squeeze(p['TwistedSF'][1, 1, 1:-1, 0])*p['BElmntMass'][:])# 2nd mode acceleration about x axis (1) -> movement in negative y
    p['Cr'][1,1]=  sum((p['RNodes']+rh)*np.squeeze(p['TwistedSF'][0, 1, 1:-1, 0])*p['BElmntMass'][:])# 2nd mode acceleration about y axis (2) -> movement in x
    p['Cr'][2,0]= -sum((p['RNodes']+rh)*np.squeeze(p['TwistedSF'][1, 2, 1:-1, 0])*p['BElmntMass'][:])# 3nd mode acceleration about x axis (1) -> movement in negative y
    p['Cr'][2,1]=  sum((p['RNodes']+rh)*np.squeeze(p['TwistedSF'][0, 2, 1:-1, 0])*p['BElmntMass'][:])# 3nd mode acceleration about y axis (2) -> movement in x

    # --- Oe,  Oe_j = \int [~Phi_j] [~s] = { \int [~s] [~Phi_j] }^t = -1/2 *(Gr_j)^t
    # M0: nq x 6
    # Straight blade: Oe6_j= [0, 0, 0, 0, szy, szx]
    p['Oe6'] = np.zeros((nq,6))
    for j in np.arange(nq):
        szx = sum((p['RNodes']+rh) * p['TwistedSF'][0, j, 1:-1, 0]*p['BElmntMass'])
        szy = sum((p['RNodes']+rh) * p['TwistedSF'][1, j, 1:-1, 0]*p['BElmntMass'])
        p['Oe6'][j,4] = szy
        p['Oe6'][j,5] = szx

    # M1: nq x 6 x nq
    #Oe = np.zeros((nf,3,3))
    #Oe6= np.zeros((nf,6))
    o=0 # derivative order 0=shape
    Oe_M1 = np.zeros((nq,6,nq))
    SS=np.zeros((nq,nq,3,3))
    for j in np.arange(nq):
        for k in np.arange(nq):
            s=np.zeros((3,3))
            for i1 in [0,1]: # NOTE: z,2 is 0
                for i2 in [0,1]: # NOTE: z,2 is 0
                    s[i1,i2] = sum(np.squeeze(p['TwistedSF'][i1, k, 1:-1, o])*np.squeeze(p['TwistedSF'][i2, j, 1:-1, o])*p['BElmntMass'])
            SS[j,k,:,:]=s
            #Oe6_j=  [       -(syy+szz),        -(sxx+szz),       -(sxx+syy),       sxy+syx,        syz+szy,    sxz+szx] 
            Oe_M1[j,:,k] = [-(s[1,1]+s[2,2]), -(s[0,0]+s[2,2]), -(s[0,0]+s[1,1]), s[0,1]+s[1,0], s[1,2]+s[2,1], s[0,2]+s[2,0] ]
    p['SS_M1']=SS
    # Centrifugal stiffening
    # TODO TODO
    # TODO TODO TODO
    # Below is for flap1 + edge1 not flap1 flap2 edge
    Oe_M1g = np.zeros((nq,6,nq))
    Oe_M1g[0,0,0]= p['KBFCent'][0,0];
    Oe_M1g[0,1,0]= p['KBFCent'][0,0];
    Oe_M1g[2,0,2]= p['KBECent'][0,0];
    Oe_M1g[2,1,2]= p['KBECent'][0,0];

    p['Oe_M1']  = Oe_M1
    p['Oe_M1g'] = Oe_M1g


    return p


def towerParameters(EDfilename, gravity, RotMass=None, noInertialCouplings=True):
    """
    Compute tower parameters exactly like OpenFAST
    See Routine Coeff from ElastoDyn.f90
    """
    from welib.yams.flexibility import polyshape
    # --- Read inputs
    ED = getEDClass(EDfilename)
    EDtwr   = os.path.join(os.path.dirname(ED.filename), ED['TwrFile'].replace('"',''))
    twr     = FASTInputFile(EDtwr)
    twrProp = twr.toDataFrame()
    n = ED['TwrNodes']

    # --- 
    p=dict()

    # --- Main dimensions
    p['TwrFlexL'] = ED['TowerHt'] - ED['TowerBsHt'] # Height / length of the flexible portion of the tower.
    p['TwrNodes'] = ED['TwrNodes']
    n             = ED['TwrNodes']
    nModesPerDir  = 2
    nDeriv        = 3
    # --- Spanwise nodes
    twr_fract     = np.arange(1./n/2., 1, 1./n)
    p['DHNodes'] = np.ones(n)*p['TwrFlexL']/n
    p['HNodes']     = twr_fract * p['TwrFlexL']  # NOTE: we don't add the TowerBsHt!
    p['HNodesNorm'] = p['HNodes']/p['TwrFlexL']
    p['Twr_s_span'] = np.concatenate(([0], p['HNodesNorm'], [1]))*p['TwrFlexL']; # Midpoints + 0 and L
    p['Twr_s_span_norm'] = p['Twr_s_span']/p['TwrFlexL']
    # --- Interpolate properties to new nodal positions
    HtFract= twrProp['HtFract_[-]'].values
    p['MassT']     = np.interp(p['HNodesNorm'], HtFract, twrProp['TMassDen_[kg/m]'])     ;
    p['StiffTFA']  = np.interp(p['HNodesNorm'], HtFract, twrProp['TwFAStif_[Nm^2]'])     ;
    p['StiffTSS']  = np.interp(p['HNodesNorm'], HtFract, twrProp['TwSSStif_[Nm^2]'])     ;
    p['m_full']    = np.interp(p['Twr_s_span_norm'],HtFract, twrProp['TMassDen_[kg/m]'])     ;
    p['EI_FA_full'] = np.interp(p['Twr_s_span_norm'], HtFract, twrProp['TwFAStif_[Nm^2]'])
    p['EI_SS_full'] = np.interp(p['Twr_s_span_norm'], HtFract, twrProp['TwSSStif_[Nm^2]'])
    # Shape coefficients
    p['TwFAM1Sh'] = [twr[c] for c in ['TwFAM1Sh(2)', 'TwFAM1Sh(3)', 'TwFAM1Sh(4)', 'TwFAM1Sh(5)', 'TwFAM1Sh(6)']]
    p['TwFAM2Sh'] = [twr[c] for c in ['TwFAM2Sh(2)', 'TwFAM2Sh(3)', 'TwFAM2Sh(4)', 'TwFAM2Sh(5)', 'TwFAM2Sh(6)']]
    p['TwSSM1Sh'] = [twr[c] for c in ['TwSSM1Sh(2)', 'TwSSM1Sh(3)', 'TwSSM1Sh(4)', 'TwSSM1Sh(5)', 'TwSSM1Sh(6)']]
    p['TwSSM2Sh'] = [twr[c] for c in ['TwSSM2Sh(2)', 'TwSSM2Sh(3)', 'TwSSM2Sh(4)', 'TwSSM2Sh(5)', 'TwSSM2Sh(6)']]
    p['FAStTunr'] = [twr['FAStTunr(1)'], twr['FAStTunr(2)']]
    p['SSStTunr'] = [twr['SSStTunr(1)'], twr['SSStTunr(2)']]
    p['TwrFADmp'] = [twr['TwrFADmp(1)'], twr['TwrFADmp(2)']]
    p['TwrSSDmp'] = [twr['TwrSSDmp(1)'], twr['TwrSSDmp(2)']]
    # Calculate the tower-top mass:
    p['TwrTpMass'] = RotMass + ED['NacMass'] + ED['YawBrMass'];
    p['TElmntMass'] = p['MassT']*p['DHNodes']  # Mass of tower element J
    p['TwrMass'] = sum(p['TElmntMass'])
    # Mass above node
    p['TMssAbvNd'] = np.zeros(n)
    for J in np.arange(n-1,-1,-1): 
        p['TMssAbvNd'][J] = 0.5*p['TElmntMass'][J];
        if J == n-1:   # Uppermost tower element
            # Add the TwrTpMass effects:
            p['TMssAbvNd'][J] = p['TMssAbvNd'][J] + p['TwrTpMass']
        else:         # All other tower elements
            # Add to TMssAbvNd'][J] the effects from the (not yet used) portion of element J+1
            p['TMssAbvNd'][J] = 0.5*p['TElmntMass'][J+1] + p['TMssAbvNd'][J] + p['TMssAbvNd'][J+1];

    # --- Tower shape functions (all derivatives) for mode 1&2 FA and SS
    p['TwrFASF'] = np.zeros((nModesPerDir, n+2, nDeriv)) # NOTE: full (+2)
    p['TwrSSSF'] = np.zeros((nModesPerDir, n+2, nDeriv)) # NOTE: full (+2)
    p['TwrFASF'][0,:,0],p['TwrFASF'][0,:,1],p['TwrFASF'][0,:,2] = polyshape(p['Twr_s_span'], coeff=p['TwFAM1Sh'], x_max=p['TwrFlexL'], doscale=False)
    p['TwrFASF'][1,:,0],p['TwrFASF'][1,:,1],p['TwrFASF'][1,:,2] = polyshape(p['Twr_s_span'], coeff=p['TwFAM2Sh'], x_max=p['TwrFlexL'], doscale=False)
    p['TwrSSSF'][0,:,0],p['TwrSSSF'][0,:,1],p['TwrSSSF'][0,:,2] = polyshape(p['Twr_s_span'], coeff=p['TwSSM1Sh'], x_max=p['TwrFlexL'], doscale=False)
    p['TwrSSSF'][1,:,0],p['TwrSSSF'][1,:,1],p['TwrSSSF'][1,:,2] = polyshape(p['Twr_s_span'], coeff=p['TwSSM2Sh'], x_max=p['TwrFlexL'], doscale=False)
    
    # --- Generalized mass
    p['MTFA'] = np.zeros((2, 2))
    p['MTSS'] = np.zeros((2, 2))
    for I in [0,1]:  # Loop through all tower modes in a single direction
        p['MTFA'][I,I] = p['TwrTpMass'];
        p['MTSS'][I,I] = p['TwrTpMass'];
    for I in [0,1]:    # Loop through all tower DOFs in one direction
        for L in [0,1]:  # Loop through all tower DOFs in one direction
            p['MTFA'][I,L] += sum(p['TElmntMass']*p['TwrFASF'][I,1:-1,0]*p['TwrFASF'][L,1:-1,0])
            p['MTSS'][I,L] += sum(p['TElmntMass']*p['TwrSSSF'][I,1:-1,0]*p['TwrSSSF'][L,1:-1,0])
    if noInertialCouplings:
        p['MTFA'][0,1]=0
        p['MTFA'][1,0]=0
        p['MTSS'][0,1]=0
        p['MTSS'][1,0]=0

    # --- Generalized stiffness
    p['KTFA'] = np.zeros((nModesPerDir, nModesPerDir))
    p['KTSS'] = np.zeros((nModesPerDir, nModesPerDir))
    for I in [0,1]:    # Loop through all tower DOFs in one direction
        for L in [0,1]:  # Loop through all tower DOFs in one direction
            p['KTFA'][I,L] = sum(p['StiffTFA']*p['DHNodes']*p['TwrFASF'][I,1:-1,2]*p['TwrFASF'][L,1:-1,2])
            p['KTSS'][I,L] = sum(p['StiffTSS']*p['DHNodes']*p['TwrSSSF'][I,1:-1,2]*p['TwrSSSF'][L,1:-1,2])
    # --- Self-weight geometrical stiffness
    p['KTFAGrav'] = np.zeros((nModesPerDir, nModesPerDir))
    p['KTSSGrav'] = np.zeros((nModesPerDir, nModesPerDir))
    for I in [0,1]:    # Loop through all tower DOFs in one direction
        p['KTFAGrav'][I,I] = sum( - p['TMssAbvNd']*p['DHNodes']*p['TwrFASF'][I,1:-1,1]**2)*gravity
        p['KTSSGrav'][I,I] = sum( - p['TMssAbvNd']*p['DHNodes']*p['TwrSSSF'][I,1:-1,1]**2)*gravity
    # --- Tower top geometric stiffness
    p['KTFAGravTT'] = np.zeros((nModesPerDir, nModesPerDir))
    p['KTSSGravTT'] = np.zeros((nModesPerDir, nModesPerDir))
    for I in [0,1]:    # Loop through all tower DOFs in one direction
        # TODO CHECK SIGN
        p['KTFAGravTT'][I,I] = sum(- p['TwrTpMass']* p['DHNodes']*p['TwrFASF'][I,1:-1,1]**2)*gravity
        p['KTSSGravTT'][I,I] = sum(- p['TwrTpMass']* p['DHNodes']*p['TwrSSSF'][I,1:-1,1]**2)*gravity
    # --- Integrate to find the tower axial reduction shape functions:
    p['AxRedTFA']= np.zeros((nModesPerDir, nModesPerDir, n+2)) # NOTE: full (+2)
    p['AxRedTSS']= np.zeros((nModesPerDir, nModesPerDir, n+2)) # NOTE: full (+2)
    for J in np.arange(n):   # Loop through the tower nodes / elements
        AxRdTFA= np.zeros((2, 2))
        AxRdTSS= np.zeros((2, 2))
        for I in [0,1]:     # Loop through all tower DOFs in one direction
            for L in [0,1]:  # Loop through all tower DOFs in one direction
                AxRdTFA [I,L] = 0.5*p['DHNodes'][J]*p['TwrFASF'][I,J+1,1]*p['TwrFASF'][L,J,1];
                AxRdTSS [I,L] = 0.5*p['DHNodes'][J]*p['TwrSSSF'][I,J+1,1]*p['TwrSSSF'][L,J,1];
                p['AxRedTFA'][I,L,J+1] = AxRdTFA[I,L]
                p['AxRedTSS'][I,L,J+1] = AxRdTSS[I,L]
        if J != 0:    # All but the lowermost tower element
            # Add the effects from the (not yet used) portion of element J-1
            for I in [0,1]:     # Loop through all tower DOFs in one direction
                for L in [0,1]:  # Loop through all tower DOFs in one direction
                    p['AxRedTFA'][I,L,J+1] += p['AxRedTFA'][I,L,J]+ AxRdTFAOld[I,L]
                    p['AxRedTSS'][I,L,J+1] += p['AxRedTSS'][I,L,J]+ AxRdTSSOld[I,L]
        # Store the AxRdTFA and AxRdTSS terms of the current element (these will be used for the next element)
        AxRdTFAOld = AxRdTFA;
        AxRdTSSOld = AxRdTSS;
    # Integrate to find the tower axial reduction shape functions at the tower-top:
    for I in [0,1]:     # Loop through all tower DOFs in one direction
        for L in [0,1]:  # Loop through all tower DOFs in one direction
            p['AxRedTFA'][I,L,-1] = p['AxRedTFA'][I,L,-2] + AxRdTFAOld[I,L]
            p['AxRedTSS'][I,L,-1] = p['AxRedTSS'][I,L,-2] + AxRdTSSOld[I,L]
    # Apply the modal stiffness tuners of the tower to KTFA() and KTSS():
    for I in [0,1]:     # Loop through all tower DOFs in one direction
        for L in [0,1]:  # Loop through all tower DOFs in one direction
            p['KTFA'][I,L] = np.sqrt( p['FAStTunr'][I]*p['FAStTunr'][L] )*p['KTFA'][I,L]
            p['KTSS'][I,L] = np.sqrt( p['SSStTunr'][I]*p['SSStTunr'][L] )*p['KTSS'][I,L]
    # Calculate the tower natural frequencies:
    p['FreqTFA'] = np.zeros((nModesPerDir, 3)) # NOTE: third frequency not computed by OpenFAST
    p['FreqTSS'] = np.zeros((nModesPerDir, 3))
    p['CTFA']    = np.zeros((nModesPerDir, 2))
    p['CTSS']    = np.zeros((nModesPerDir, 2))
    for I in [0,1]:     # Loop through all tower DOFs in one direction
        allKTFAGrav= p['KTFAGrav'][I,I] + p['KTFAGravTT'][I,I]
        allKTSSGrav= p['KTSSGrav'][I,I] + p['KTSSGravTT'][I,I]
        p['FreqTFA'][I,0] = (1/2/np.pi)*np.sqrt(   p['KTFA'][I,I]                       /( p['MTFA'][I,I] - p['TwrTpMass'] ) ) # Natural tower I-fore-aft frequency w/o gravitational destiffening nor tower-top mass effects
        p['FreqTFA'][I,1] = (1/2/np.pi)*np.sqrt( ( p['KTFA'][I,I] + p['KTFAGrav'][I,I] )/  p['MTFA'][I,I]                    ) # Natural tower I-fore-aft frequency w/  gravitational destiffening and tower-top mass effects
        p['FreqTFA'][I,2] = (1/2/np.pi)*np.sqrt( ( p['KTFA'][I,I] + allKTFAGrav         )/ p['MTFA'][I,I]                    ) # Natural tower I-fore-aft frequency w/  gravitational destiffening and tower-top mass effects
        p['FreqTSS'][I,0] = (1/2/np.pi)*np.sqrt(   p['KTSS'][I,I]                       /( p['MTSS'][I,I] - p['TwrTpMass'] ) ) # Natural tower I-side-to-side frequency w/o gravitational destiffening nor tower-top mass effects
        p['FreqTSS'][I,1] = (1/2/np.pi)*np.sqrt( ( p['KTSS'][I,I] + p['KTSSGrav'][I,I]  )/ p['MTSS'][I,I]                    ) # Natural tower I-side-to-side frequency w/  gravitational destiffening and tower-top mass effects
        p['FreqTSS'][I,2] = (1/2/np.pi)*np.sqrt( ( p['KTSS'][I,I] + allKTSSGrav         )/ p['MTSS'][I,I]                    ) # Natural tower I-side-to-side frequency w/  gravitational destiffening and tower-top mass effects
    # Calculate the generalized damping of the tower:
    for I in [0,1]:     # Loop through all tower DOFs in one direction
        for L in [0,1]:  # Loop through all tower DOFs in one direction
            p['CTFA'][I,L] = ( 0.01*p['TwrFADmp'][L] )*p['KTFA'][I,L]/( np.pi*p['FreqTFA'][L,0] );
            p['CTSS'][I,L] = ( 0.01*p['TwrSSDmp'][L] )*p['KTSS'][I,L]/( np.pi*p['FreqTSS'][L,0] );


    # ---  Shape functions (FA1 FA2 SS1 SS2)
    nq=4
    nNodes=n+2
    p['U']  = np.zeros((nq, 3, nNodes))
    p['V']  = np.zeros((nq, 3, nNodes))
    p['K']  = np.zeros((nq, 3, nNodes))
    j=0
    for idir, name in zip((0,1),('FA','SS')):
        for jm in [0,1]:
            p['U'][j][idir,:] = p['Twr{}SF'.format(name)][jm, :, 0]
            p['V'][j][idir,:] = p['Twr{}SF'.format(name)][jm, :, 1]
            p['K'][j][idir,:] = p['Twr{}SF'.format(name)][jm, :, 2]
            j+=1
    # --- Parameters consistent with "YAMS" flexibility module
    p['s_span'] = p['Twr_s_span']
    p['m'] = p['m_full']
    p['EI'] = np.zeros((3,n+2))
    p['EI'][0,:] = p['EI_FA_full']
    p['EI'][1,:] = p['EI_SS_full']

    p['s_G0'] = np.zeros((3, len(p['Twr_s_span'])))
    p['s_G0'][2,:] = p['Twr_s_span']  # TODO add hub radius
#     for k,v in p.items():
#         if hasattr(v, '__len__'):
#             v = np.asarray(v)
#             if len(v.shape)>=3:
#                 print('{:15s}:({})'.format(k,v.shape))
#             elif len(v.shape)==2:
#                 print('{:15s}:\n {} ({})'.format(k,v, v.shape))
#             else:
#                 n=len(v)
#                 print('{:15s}:{} ({})'.format(k,v,n))
#         else:
#             print('{:15s}:{}'.format(k,v))


    # --- "purely" elastic mass
    p['MTFA_e'] = p['MTFA']
    p['MTSS_e'] = p['MTSS']
    for I in [0,1]: #Loop through all tower modes in a single direction
        p['MTFA_e'][I,I] = p['MTFA'][I,I] - p['TwrTpMass']
        p['MTSS_e'][I,I] = p['MTSS'][I,I] - p['TwrTpMass']

    # TODO TODO TODO Check signs
    p['KTFAGrav_nd']   = -p['KTFAGrav']/gravity                  # reverting to a dimensionless factor
    p['KTSSGrav_nd']   = -p['KTSSGrav']/gravity                  # reverting to a dimensionless factor
    p['KTFAGravTT_nd'] = -p['KTFAGravTT']/gravity/p['TwrTpMass'] # reverting to a dimensionless factor
    p['KTSSGravTT_nd'] = -p['KTSSGravTT']/gravity/p['TwrTpMass'] # reverting to a dimensionless factor


    return p

def towerDerivedParameters(p):
    """ Compute blade derived parameters, suitable for SID:
      - Inertial matrices
      - Stiffness matrices
    The parameters are computed "by hand" (as opposed to using GMBeam)
    NOTE: this function is mostly for debugging purposes. 
          Use GMBeam with method='OpenFAST' instead

    Order of shape functions: FA1 FA2 SS1 SS2
    """

    # param.tower_Ct1_1_1_3= p.KTFAGrav(1, 1);
    # param.tower_Ct1_2_2_3= p.KTSSGrav(1, 1);

    # param.tower_frame_11_origin1_1_1_1= 1;
    # param.tower_frame_11_origin1_2_2_1= 1;
    # param.tower_frame_11_phi1_1_3_1 = p.KTFAGravTT(1, 1)   ;
    # param.tower_frame_11_phi1_2_3_2 = p.KTSSGravTT(1, 1)   ;
    # param.tower_frame_11_psi0_1_2   = -p.TwrFASF(1, end, 2);
    # param.tower_frame_11_psi0_2_1   = p.TwrSSSF(1, end, 2) ;

    nq = 4
    # --- Rigid mass matrix terms
    p['mdCM'] = np.zeros((3))
    p['mdCM'][2]=  sum(p['TElmntMass'][:]*(p['HNodes']));
    p['J']    = np.zeros((3,3))
    p['J'][0,0] = sum(p['TElmntMass'][:]*(p['HNodes'])**2)
    p['J'][1,1] = sum(p['TElmntMass'][:]*(p['HNodes'])**2)
    # --- Elastic matrices
    p['Ke0'] = np.zeros((nq,nq)) # Without stiffening
    p['De'] = np.zeros((nq,nq))
    p['Me'] = np.zeros((nq,nq))
    for I in [0,1]:     # Loop through all tower DOFs in one direction
        for L in [0,1]:  # Loop through all tower DOFs in one direction
            p['Me'][I  ,L]   = p['MTFA_e'][I, L]
            p['Me'][I+2,L+2] = p['MTSS_e'][I, L]
            p['Ke0'][I  ,L]   = p['KTFA'][I, L]
            p['Ke0'][I+2,L+2] = p['KTSS'][I, L]
            p['De'][I  ,L]   = p['CTFA'][I, L]
            p['De'][I+2,L+2] = p['CTSS'][I, L]

    # --- Self-weight, and top mass geometrical stiffness
    p['Kg_SW'] = np.zeros((nq,nq))
    p['Kg_TM'] = np.zeros((nq,nq))
    for I in [0,1]: 
        for L in [0,1]:  # Loop through all tower DOFs in one direction
            p['Kg_SW'][I,   L]   = p['KTFAGrav'][I,L]
            p['Kg_SW'][I+2, L+2] = p['KTSSGrav'][I,L]
            p['Kg_TM'][I,   L]   = p['KTFAGravTT'][I,L]
            p['Kg_TM'][I+2, L+2] = p['KTSSGravTT'][I,L]

    p['Ke'] = p['Ke0'] + p['Kg_SW'] # + p['Kg_TM'] # NOTE: TM  not included, needs further validation
    # --- Elastic mass matrix terms
    # Ct
    p['Ct'] = np.zeros((nq,3))
    p['Ct'][0,0] = sum(np.squeeze(p['TwrFASF'][0, 1:-1, 0])*p['TElmntMass'][:])
    p['Ct'][1,0] = sum(np.squeeze(p['TwrFASF'][1, 1:-1, 0])*p['TElmntMass'][:])
    p['Ct'][2,1] = sum(np.squeeze(p['TwrSSSF'][0, 1:-1, 0])*p['TElmntMass'][:])
    p['Ct'][3,1] = sum(np.squeeze(p['TwrSSSF'][1, 1:-1, 0])*p['TElmntMass'][:])
    # Cr
    p['Cr'] = np.zeros((nq,3))
    p['Cr'][0,1]=  sum((p['HNodes'])*np.squeeze(p['TwrFASF'][0, 1:-1, 0])*p['TElmntMass'][:])
    p['Cr'][1,1]=  sum((p['HNodes'])*np.squeeze(p['TwrFASF'][1, 1:-1, 0])*p['TElmntMass'][:])
    p['Cr'][2,0]= -sum((p['HNodes'])*np.squeeze(p['TwrSSSF'][0, 1:-1, 0])*p['TElmntMass'][:])
    p['Cr'][3,0]= -sum((p['HNodes'])*np.squeeze(p['TwrSSSF'][1, 1:-1, 0])*p['TElmntMass'][:])
    return p



def ED_qDict2q(qDict):
    """ 
    From a dictionary of values, return the flat array of degrees of freedom ordered liked ElastoDyn.
    The dictionary is a convenient way to just specify few DOFs
    Example:
       qDict={'Hv':0}
    """
    q = np.zeros(ED_MaxDOFs)

    qMAP = {'Sg':DOF_Sg,'Sw':DOF_Sw,'Hv':DOF_Hv, 'R':DOF_R, 'P':DOF_P, 'Y':DOF_Y}
    qMAP.update({'TFA1':DOF_TFA1, 'TSS1':DOF_TSS1, 'TFA2':DOF_TFA2, 'TSS2':DOF_TSS2})
    qMAP.update({'Yaw':DOF_Yaw, 'RFrl':DOF_RFrl, 'GeAz':DOF_GeAz, 'DrTr':DOF_DrTr})

    for k,v in qDict.items():
        if k not in qMAP.keys():
            raise Exception('Key {} not supported by qMAP'.format(k))
        if k=='TSS1' or k=='TSS2':
            q[qMAP[k]] = -v # <<<< NOTE: DOF has a negative convention
        else:
            q[qMAP[k]] = v
    return q


def EDVec2IEC(v, offset=None):
    """ Convert a vector/array of coordinates from ElastoDyn coordinate system to IEC"""
    if offset is None:
        offset = np.array([0,0,0])
    v = np.asarray(v)
    if len(v.shape)==1:
        vIEC = np.array([ v[0]+offset[0], -v[2]+offset[1],  v[1]+offset[2]])
    elif len(v.shape)==2 and v.shape[1]==3:
        vIEC = np.column_stack( [ v[:,0]+offset[0], -v[:,2]+offset[1],  v[:,1]+offset[2]] )
    else:
        raise NotImplementedError()
    return vIEC

def EDSysVectoIECDCM(a1, a2, a3):
    """ Return DCM from ElastoDyn coordsys vectors (in a different coordinate system)"""
    R_g2t = np.zeros((3,3))
    R_g2t[:,0] = [ a1[0], -a3[0],  a2[0]]
    R_g2t[:,1] = [-a1[2],  a3[2], -a2[2]]
    R_g2t[:,2] = [ a1[1], -a3[1],  a2[1]]
    return R_g2t

def EDSysVecRot(R,z1,z2,z3):
    """ 
    Apply transformation matrix to vectors to obtain other vectors
    """
    a1 = R[0,0]*z1 + R[0,1]*z2 + R[0,2]*z3
    a2 = R[1,0]*z1 + R[1,1]*z2 + R[1,2]*z3
    a3 = R[2,0]*z1 + R[2,1]*z2 + R[2,2]*z3
    return a1, a2 ,a3




def ED_CoordSys(q=None, qDict=None, TwrFASF=None, TwrSSSF=None):
    """ 
    Return ElastoDyn coordinate systems in ElastoDyn and IEC convention:

    See ElastoDyn line 5978 
    See ElastoDyn line 1678

    INPUTS:
     - q 
       OR
     - qDict
    """
    # 
    if qDict is not None:
        q = ED_qDict2q(qDict)


    from welib.yams.rotations import smallRot_OF
    # --- Inertial frame coordinate system:
    z1 = np.array([ 1, 0, 0]) # Vector / direction z1 (=  xi from the IEC coord. system).
    z2 = np.array([ 0, 1, 0]) # Vector / direction z2 (=  zi from the IEC coord. system).
    z3 = np.array([ 0, 0, 1]) # Vector / direction z3 (= -yi from the IEC coord. system).
    R_g2g = np.eye(3)
    # --- Tower base / platform coordinate system:
    # Vector / direction a1 (=  xt from the IEC coord. system)
    # Vector / direction a2 (=  zt from the IEC coord. system)
    # Vector / direction a3 (= -yt from the IEC coord. system)
    R = smallRot_OF(q[DOF_R], q[DOF_Y], -q[DOF_P])
    a1, a2, a3 = EDSysVecRot(R, z1, z2, z3)
    # IEC Coordinate system
    R_g2t = EDSysVectoIECDCM(a1, a2, a3)

    # --- Tower element-fixed coordinate system:
    # Vector / direction t1 for tower node J (=  Lxt from the IEC coord. system).
    # Vector / direction t2 for tower node J (=  Lzt from the IEC coord. system).
    # Vector / direction t3 for tower node J (= -Lyt from the IEC coord. system).
    if TwrFASF is not None:
        TwrNodes = TwrFASF.shape[1] - 2
        t1 = np.zeros((TwrNodes,3))
        t2 = np.zeros((TwrNodes,3))
        t3 = np.zeros((TwrNodes,3))
        R_g2Ts = np.zeros((TwrNodes,3,3))# List of transformations from global to tower elements
        for j in range(TwrNodes):
            # Slope V:    Mode 1   V             Mode 2   V 
            ThetaFA = -TwrFASF[0,j,1]*q[DOF_TFA1] - TwrFASF[1,j,1]*q[DOF_TFA2]
            ThetaSS =  TwrSSSF[0,j,1]*q[DOF_TSS1] + TwrSSSF[1,j,1]*q[DOF_TSS2]
            R = smallRot_OF(ThetaSS, 0, ThetaFA)
            t1[j,:], t2[j,:], t3[j,:] = EDSysVecRot(R, a1, a2, a3)
            R_g2Ts[j,:,:] =EDSysVectoIECDCM(t1[j,:], t2[j,:], t3[j,:])

        # --- Tower-top / base plate coordinate system:
        # Slope V:    Mode 1   V             Mode 2   V 
        #print('>>> Coupling coeffs FA', -TwrFASF[0,-1,1], -TwrFASF[0,-1,1])
        #print('>>> Coupling coeffs SS',  TwrSSSF[0,-1,1],  TwrSSSF[0,-1,1])
        ThetaFA = -TwrFASF[0,-1,1]*q[DOF_TFA1] - TwrFASF[1,-1,1]*q[DOF_TFA2]
        ThetaSS =  TwrSSSF[0,-1,1]*q[DOF_TSS1] + TwrSSSF[1,-1,1]*q[DOF_TSS2]
        R = smallRot_OF(ThetaSS, 0, ThetaFA)
        b1, b2, b3 = EDSysVecRot(R, a1, a2, a3)
        R_g2b = EDSysVectoIECDCM(b1, b2, b3)
    else:
        t1, t2, t3 = None, None, None
        R_g2Ts = None

        b1, b2, b3  = a1, a2, a3
        R_g2b = R_g2t

    # --- Nacelle / yaw coordinate system:
    CNacYaw  = np.cos(q[DOF_Yaw])
    SNacYaw  = np.sin(q[DOF_Yaw])
    d1 = CNacYaw*b1 - SNacYaw*b3     # Vector / direction d1 (=  xn from the IEC coord. system).
    d2 = b2                          # Vector / direction d2 (=  zn from the IEC coord. system).
    d3 = SNacYaw*b1 + CNacYaw*b3     # Vector / direction d3 (= -yn from the IEC coord. system).
    R_g2n = EDSysVectoIECDCM(d1, d2, d3)

    # --- To be continued
    # ...
    #     ! Rotor-furl coordinate system:
    # 
    #    CRotFurl = COS( x%QT(DOF_RFrl) )
    #    SRotFurl = SIN( x%QT(DOF_RFrl) )
    # 
    #    CoordSys%rf1 = ( (   1.0 - p%CRFrlSkw2*p%CRFrlTlt2 )*CRotFurl   + p%CRFrlSkw2*p%CRFrlTlt2          )*CoordSys%d1 &
    #                 + ( p%CRFrlSkew*p%CSRFrlTlt*( 1.0 -     CRotFurl ) - p%SRFrlSkew*p%CRFrlTilt*SRotFurl )*CoordSys%d2 &
    #                 + ( p%CSRFrlSkw*p%CRFrlTlt2*( CRotFurl - 1.0     ) -             p%SRFrlTilt*SRotFurl )*CoordSys%d3
    #    CoordSys%rf2 = ( p%CRFrlSkew*p%CSRFrlTlt*( 1.0 -     CRotFurl ) + p%SRFrlSkew*p%CRFrlTilt*SRotFurl )*CoordSys%d1 &
    #                 + (             p%CRFrlTlt2*            CRotFurl   +             p%SRFrlTlt2          )*CoordSys%d2 &
    #                 + ( p%SRFrlSkew*p%CSRFrlTlt*( CRotFurl - 1.0     ) + p%CRFrlSkew*p%CRFrlTilt*SRotFurl )*CoordSys%d3
    #    CoordSys%rf3 = ( p%CSRFrlSkw*p%CRFrlTlt2*( CRotFurl - 1.0     ) +             p%SRFrlTilt*SRotFurl )*CoordSys%d1 &
    #                 + ( p%SRFrlSkew*p%CSRFrlTlt*( CRotFurl - 1.0     ) - p%CRFrlSkew*p%CRFrlTilt*SRotFurl )*CoordSys%d2 &
    #                 + ( (   1.0 - p%SRFrlSkw2*p%CRFrlTlt2 )*CRotFurl   + p%SRFrlSkw2*p%CRFrlTlt2          )*CoordSys%d3
    #    CoordSys%rfa = p%CRFrlSkew*p%CRFrlTilt*CoordSys%d1 + p%SRFrlTilt*CoordSys%d2 - p%SRFrlSkew*p%CRFrlTilt*CoordSys%d3
    # 
    # 
    #       ! Shaft coordinate system:
    #    CoordSys%c1 =  p%CShftSkew*p%CShftTilt*CoordSys%rf1 + p%SShftTilt*CoordSys%rf2 - p%SShftSkew*p%CShftTilt*CoordSys%rf3  ! Vector / direction c1 (=  xs from the IEC coord. system).
    #    CoordSys%c2 = -p%CShftSkew*p%SShftTilt*CoordSys%rf1 + p%CShftTilt*CoordSys%rf2 + p%SShftSkew*p%SShftTilt*CoordSys%rf3  ! Vector / direction c2 (=  zs from the IEC coord. system).
    #    CoordSys%c3 =  p%SShftSkew*            CoordSys%rf1                            + p%CShftSkew*            CoordSys%rf3  ! Vector / direction c3 (= -ys from the IEC coord. system).
    #       ! Azimuth coordinate system:
    #    CAzimuth = COS( x%QT(DOF_DrTr) + x%QT(DOF_GeAz) )
    #    SAzimuth = SIN( x%QT(DOF_DrTr) + x%QT(DOF_GeAz) )
    #    CoordSys%e1 =  CoordSys%c1                                  ! Vector / direction e1 (=  xa from the IEC coord. system).
    #    CoordSys%e2 =  CAzimuth*CoordSys%c2 + SAzimuth*CoordSys%c3  ! Vector / direction e2 (=  ya from the IEC coord. system).
    #    CoordSys%e3 = -SAzimuth*CoordSys%c2 + CAzimuth*CoordSys%c3  ! Vector / direction e3 (=  za from the IEC coord. system).
    #       ! Teeter coordinate system:
    #       ! Lets define TeetAng, which is the current teeter angle (= QT(DOF_Teet) for
    #       !   2-blader or 0 for 3-blader) and is used in place of QT(DOF_Teet)
    #       !   throughout SUBROUTINE RtHS().  Doing it this way, we can run the same
    #       !   equations of motion for both the 2 and 3-blader configurations even
    #       !   though a 3-blader does not have a teetering DOF.
    #    IF ( p%NumBl == 2 )  THEN ! 2-blader
    #       dat%TeetAng    = x%QT (DOF_Teet)
    #       dat%TeetAngVel = x%QDT(DOF_Teet)
    #    ELSE                    ! 3-blader
    #       dat%TeetAng    = 0.0  ! Teeter is not an available DOF for a 3-blader
    #       dat%TeetAngVel = 0.0  ! Teeter is not an available DOF for a 3-blader
    #    ENDIF
    #    CTeetAng = COS( dat%TeetAng )
    #    STeetAng = SIN( dat%TeetAng )
    #    CoordSys%f1 = CTeetAng*CoordSys%e1 - STeetAng*CoordSys%e3       ! Vector / direction f1.
    #    CoordSys%f2 = CoordSys%e2                                       ! Vector / direction f2.
    #    CoordSys%f3 = STeetAng*CoordSys%e1 + CTeetAng*CoordSys%e3       ! Vector / direction f3.
    #       ! Hub / delta-3 coordinate system:
    #    CoordSys%g1 =  CoordSys%f1                                      ! Vector / direction g1 (=  xh from the IEC coord. system).
    #    CoordSys%g2 =  p%CosDel3*CoordSys%f2 + p%SinDel3*CoordSys%f3    ! Vector / direction g2 (=  yh from the IEC coord. system).
    #    CoordSys%g3 = -p%SinDel3*CoordSys%f2 + p%CosDel3*CoordSys%f3    ! Vector / direction g3 (=  zh from the IEC coord. system).
    #    DO K = 1,p%NumBl ! Loop through all blades
    #       ! Hub (Prime) coordinate system rotated to match blade K.
    #        gRotAng = p%TwoPiNB*(K-1)
    #       CgRotAng = COS( gRotAng )
    #       SgRotAng = SIN( gRotAng )
    #       g1Prime =  CoordSys%g1
    #       g2Prime =  CgRotAng*CoordSys%g2 + SgRotAng*CoordSys%g3
    #       g3Prime = -SgRotAng*CoordSys%g2 + CgRotAng*CoordSys%g3
    #       ! Coned coordinate system:
    #       CoordSys%i1(K,:) = p%CosPreC(K)*g1Prime - p%SinPreC(K)*g3Prime  ! i1(K,:) = vector / direction i1 for blade K (=  xcK from the IEC coord. system).
    #       CoordSys%i2(K,:) = g2Prime                                      ! i2(K,:) = vector / direction i2 for blade K (=  ycK from the IEC coord. system).
    #       CoordSys%i3(K,:) = p%SinPreC(K)*g1Prime + p%CosPreC(K)*g3Prime  ! i3(K,:) = vector / direction i3 for blade K (=  zcK from the IEC coord. system).
    #       ! Blade / pitched coordinate system:
    #       CosPitch = COS( REAL(BlPitch(K),R8Ki) )
    #       SinPitch = SIN( REAL(BlPitch(K),R8Ki) )
    #       CoordSys%j1(K,:) = CosPitch*CoordSys%i1(K,:) - SinPitch*CoordSys%i2(K,:)      ! j1(K,:) = vector / direction j1 for blade K (=  xbK from the IEC coord. system).
    #       CoordSys%j2(K,:) = SinPitch*CoordSys%i1(K,:) + CosPitch*CoordSys%i2(K,:)      ! j2(K,:) = vector / direction j2 for blade K (=  ybK from the IEC coord. system).
    #       CoordSys%j3(K,:) = CoordSys%i3(K,:)                                           ! j3(K,:) = vector / direction j3 for blade K (=  zbK from the IEC coord. system).
    #       DO J = 0,p%TipNode ! Loop through the blade nodes / elements
    #       ! Blade coordinate system aligned with local structural axes (not element fixed):
    #          Lj1 = p%CThetaS(K,J)*CoordSys%j1(K,:) - p%SThetaS(K,J)*CoordSys%j2(K,:)  ! vector / direction Lj1 at node J for blade K
    #          Lj2 = p%SThetaS(K,J)*CoordSys%j1(K,:) + p%CThetaS(K,J)*CoordSys%j2(K,:)  ! vector / direction Lj2 at node J for blade K
    #          Lj3 = CoordSys%j3(K,:)                                                   ! vector / direction Lj3 at node J for blade K
    #       ! Blade element-fixed coordinate system aligned with local structural axes:
    #          ThetaOoP =   p%TwistedSF(K,1,1,J,1)*x%QT( DOF_BF(K,1) ) &
    #                     + p%TwistedSF(K,1,2,J,1)*x%QT( DOF_BF(K,2) ) &
    #                     + p%TwistedSF(K,1,3,J,1)*x%QT( DOF_BE(K,1) )
    #          ThetaIP  = - p%TwistedSF(K,2,1,J,1)*x%QT( DOF_BF(K,1) ) &
    #                     - p%TwistedSF(K,2,2,J,1)*x%QT( DOF_BF(K,2) ) &
    #                     - p%TwistedSF(K,2,3,J,1)*x%QT( DOF_BE(K,1) )
    #          ThetaLxb = p%CThetaS(K,J)*ThetaIP - p%SThetaS(K,J)*ThetaOoP
    #          ThetaLyb = p%SThetaS(K,J)*ThetaIP + p%CThetaS(K,J)*ThetaOoP
    #          CALL SmllRotTrans( 'blade deflection (ElastoDyn SetCoordSy)', ThetaLxb, ThetaLyb, 0.0_R8Ki, TransMat, TRIM(Num2LStr(t))//' s', ErrStat2, ErrMsg2 ) ! Get the transformation matrix, TransMat, from blade coordinate system aligned with local structural axes (not element fixed) to blade element-fixed coordinate system aligned with local structural axes.
    #          CoordSys%n1(K,J,:) = TransMat(1,1)*Lj1 + TransMat(1,2)*Lj2 + TransMat(1,3)*Lj3   ! Vector / direction n1 for node J of blade K (= LxbK from the IEC coord. system).
    #          CoordSys%n2(K,J,:) = TransMat(2,1)*Lj1 + TransMat(2,2)*Lj2 + TransMat(2,3)*Lj3   ! Vector / direction n2 for node J of blade K (= LybK from the IEC coord. system).
    #          CoordSys%n3(K,J,:) = TransMat(3,1)*Lj1 + TransMat(3,2)*Lj2 + TransMat(3,3)*Lj3   ! Vector / direction n3 for node J of blade K (= LzbK from the IEC coord. system).
    #       ! skip these next CoordSys variables at the root and the tip; they are required only for AD14:
    #          if (j == 0 .or. j==p%TipNode) cycle  
    #       ! Blade element-fixed coordinate system used for calculating and returning
    #       !    aerodynamics loads:
    #       ! This coordinate system is rotated about positive n3 by the angle
    #       !    BlPitch(K) + ThetaS(K,J) and is coincident with the i-vector triad
    #       !    when the blade is undeflected.
    #          CPitPTwstS = CosPitch*p%CThetaS(K,J) - SinPitch*p%SThetaS(K,J)  ! = COS( BlPitch(K) + ThetaS(K,J) ) found using the sum of angles formulae of cosine.
    #          SPitPTwstS = CosPitch*p%SThetaS(K,J) + SinPitch*p%CThetaS(K,J)  ! = SIN( BlPitch(K) + ThetaS(K,J) ) found using the sum of angles formulae of   sine.
    #          CoordSys%m1(K,J,:)  =  CPitPTwstS*CoordSys%n1(K,J,:) + SPitPTwstS*CoordSys%n2(K,J,:)   ! m1(K,J,:) = vector / direction m1 for node J of blade K (used to calc. and return aerodynamic loads from AeroDyn).
    #          CoordSys%m2(K,J,:)  = -SPitPTwstS*CoordSys%n1(K,J,:) + CPitPTwstS*CoordSys%n2(K,J,:)   ! m2(K,J,:) = vector / direction m2 for node J of blade K (used to calc. and return aerodynamic loads from AeroDyn).
    #          CoordSys%m3(K,J,:)  =  CoordSys%n3(K,J,:)                                              ! m3(K,J,:) = vector / direction m3 for node J of blade K (used to calc. and return aerodynamic loads from AeroDyn).
    #       ! Calculate the trailing edge coordinate system used in noise calculations.
    #       ! This coordinate system is blade element-fixed and oriented with the local
    #       !   aerodynamic axes (te2 points toward trailing edge, te1 points toward
    #       !   suction surface):
    #          CPitPTwstA = CosPitch*p%CAeroTwst(J) - SinPitch*p%SAeroTwst(J)  ! = COS( BlPitch(K) + AeroTwst(J) ) found using the sum of angles formulae of cosine.
    #          SPitPTwstA = CosPitch*p%SAeroTwst(J) + SinPitch*p%CAeroTwst(J)  ! = SIN( BlPitch(K) + AeroTwst(J) ) found using the sum of angles formulae of   sine.
    #          CoordSys%te1(K,J,:) =  CPitPTwstA*CoordSys%m1(K,J,:) - SPitPTwstA*CoordSys%m2(K,J,:)   ! te1(K,J,:) = vector / direction te1 for node J of blade K (used to calc. noise and to calc. and return aerodynamic loads from AeroDyn).
    #          CoordSys%te2(K,J,:) =  SPitPTwstA*CoordSys%m1(K,J,:) + CPitPTwstA*CoordSys%m2(K,J,:)   ! te2(K,J,:) = vector / direction te2 for node J of blade K (used to calc. noise and to calc. and return aerodynamic loads from AeroDyn).
    #          CoordSys%te3(K,J,:) =  CoordSys%m3(K,J,:)                                              ! te3(K,J,:) = vector / direction te3 for node J of blade K (used to calc. noise and to calc. and return aerodynamic loads from AeroDyn).
    #       ENDDO ! J - Blade nodes / elements
    #    ENDDO ! K - Blades
    #       ! Tail-furl coordinate system:
    #    CTailFurl = COS( x%QT(DOF_TFrl) )
    #    STailFurl = SIN( x%QT(DOF_TFrl) )
    #    CoordSys%tf1 = ( ( 1.0 - p%CTFrlSkw2*p%CTFrlTlt2 )*CTailFurl  + p%CTFrlSkw2*p%CTFrlTlt2           )*CoordSys%d1 &
    #                 + ( p%CTFrlSkew*p%CSTFrlTlt*(  1.0 - CTailFurl ) - p%STFrlSkew*p%CTFrlTilt*STailFurl )*CoordSys%d2 &
    #                 + ( p%CSTFrlSkw*p%CTFrlTlt2*( CTailFurl - 1.0  ) -             p%STFrlTilt*STailFurl )*CoordSys%d3
    #    CoordSys%tf2 = ( p%CTFrlSkew*p%CSTFrlTlt*(  1.0 - CTailFurl ) + p%STFrlSkew*p%CTFrlTilt*STailFurl )*CoordSys%d1 &
    #                 + (             p%CTFrlTlt2*         CTailFurl +               p%STFrlTlt2           )*CoordSys%d2 &
    #                 + ( p%STFrlSkew*p%CSTFrlTlt*( CTailFurl - 1.0  ) + p%CTFrlSkew*p%CTFrlTilt*STailFurl )*CoordSys%d3
    #    CoordSys%tf3 = ( p%CSTFrlSkw*p%CTFrlTlt2*( CTailFurl - 1.0  ) +             p%STFrlTilt*STailFurl )*CoordSys%d1 &
    #                 + ( p%STFrlSkew*p%CSTFrlTlt*( CTailFurl - 1.0  ) - p%CTFrlSkew*p%CTFrlTilt*STailFurl )*CoordSys%d2 &
    #                 + ( ( 1.0 - p%STFrlSkw2*p%CTFrlTlt2 )*CTailFurl  + p%STFrlSkw2*p%CTFrlTlt2           )*CoordSys%d3
    #    CoordSys%tfa = p%CTFrlSkew*p%CTFrlTilt*CoordSys%d1 + p%STFrlTilt*CoordSys%d2 - p%STFrlSkew*p%CTFrlTilt*CoordSys%d3
    #     R_

    CoordSys = dict()
    # ED
    CoordSys['z1'] = z1
    CoordSys['z2'] = z2
    CoordSys['z3'] = z3
    CoordSys['a1'] = a1
    CoordSys['a2'] = a2
    CoordSys['a3'] = a3
    CoordSys['t1'] = t1
    CoordSys['t2'] = t2
    CoordSys['t3'] = t3
    CoordSys['d1'] = d1
    CoordSys['d2'] = d2
    CoordSys['d3'] = d3
    # IEC
    CoordSys['R_g2f']  = R_g2t
    CoordSys['R_g2t']  = R_g2t
    CoordSys['R_g2Ts'] = R_g2Ts # To tower elements
    CoordSys['R_g2n']  = R_g2n  # To nacelle (including nacelle yaw)
    return CoordSys


def ED_Positions(q=None, qDict=None, CoordSys=None, p=None, dat=None, IEC=None):
    """ 
    See ElastoDyn.f90 CalculatePositions
    # !> This routine is used to calculate the positions stored in other states that are used in both the
    # !! CalcOutput and CalcContStateDeriv routines.
    # SUBROUTINE CalculatePositions( p, x, CoordSys, dat )
    """
    if qDict is not None:
        QT = ED_qDict2q(qDict)
    else:
        QT = q

    if dat is None:
        dat = dict()
    TTopNode=-1
    # Define the position vectors between the various points on the wind turbine
    #   that are not dependent on the distributed tower or blade parameters:
    dat['rZ']    = QT[DOF_Sg]* CoordSys['z1'] + QT[DOF_Hv] * CoordSys['z2'] - QT[DOF_Sw]* CoordSys['z3']    # Position vector from inertia frame origin to platform reference (point Z).
    dat['rZY']   = p['rZYzt']*  CoordSys['a2'] + p['PtfmCMxt']*CoordSys['a1'] - p['PtfmCMyt']*CoordSys['a3'] # Position vector from platform reference (point Z) to platform mass center (point Y).      
    dat['rZT0']  = p['rZT0zt']* CoordSys['a2']                                                               # Position vector from platform reference (point Z) to tower base (point T(0))
    dat['rZO']   = ( QT[DOF_TFA1] + QT[DOF_TFA2] )*CoordSys['a1'] # Position vector from platform reference (point Z) to tower-top / base plate (point O).
    dat['rZO']   +=  ( p['RefTwrHt'] - 0.5*(      p['AxRedTFA'][0,0,TTopNode]*QT[DOF_TFA1]*QT[DOF_TFA1] \
                                            +     p['AxRedTFA'][1,1,TTopNode]*QT[DOF_TFA2]*QT[DOF_TFA2] \
                                            + 2.0*p['AxRedTFA'][0,1,TTopNode]*QT[DOF_TFA1]*QT[DOF_TFA2] \
                                            +     p['AxRedTSS'][0,0,TTopNode]*QT[DOF_TSS1]*QT[DOF_TSS1] \
                                            +     p['AxRedTSS'][1,1,TTopNode]*QT[DOF_TSS2]*QT[DOF_TSS2] \
                                            + 2.0*p['AxRedTSS'][0,1,TTopNode]*QT[DOF_TSS1]*QT[DOF_TSS2] ))*CoordSys['a2'] 
    dat['rZO']   +=  ( QT[DOF_TSS1] + QT[DOF_TSS2])*CoordSys['a3']
    dat['rOU']   =   p['NacCMxn']*CoordSys['d1']  +  p['NacCMzn']  *CoordSys['d2']  -  p['NacCMyn']  *CoordSys['d3']  # Position vector from tower-top / base plate (point O) to nacelle center of mass (point U).
#     dat['rOV']   = p['RFrlPnt_n(1)*CoordSys['d1  +  p['RFrlPnt_n(3)*CoordSys['d2  -  p['RFrlPnt_n(2)*CoordSys['d3                            ! Position vector from tower-top / base plate (point O) to specified point on rotor-furl axis (point V).
#     dat['rVIMU'] =   p['rVIMUxn*CoordSys['rf1 +  p['rVIMUzn  *CoordSys['rf2 -   p['rVIMUyn *CoordSys['rf3                           ! Position vector from specified point on rotor-furl axis (point V) to nacelle IMU (point IMU).
#     dat['rVD']   =     p['rVDxn*CoordSys['rf1 +    p['rVDzn  *CoordSys['rf2 -     p['rVDyn *CoordSys['rf3                           ! Position vector from specified point on rotor-furl axis (point V) to center of mass of structure that furls with the rotor (not including rotor) (point D).
#     dat['rVP']   =     p['rVPxn*CoordSys['rf1 +    p['rVPzn  *CoordSys['rf2 -     p['rVPyn *CoordSys['rf3 + p['OverHang*CoordSys['c1  ! Position vector from specified point on rotor-furl axis (point V) to teeter pin (point P).
#     dat['rPQ']   = -p['UndSling*CoordSys['g1                                                                                    ! Position vector from teeter pin (point P) to apex of rotation (point Q).
#     dat['rQC']   =     p['HubCM*CoordSys['g1                                                                                    ! Position vector from apex of rotation (point Q) to hub center of mass (point C).
#     dat['rOW']   = p['TFrlPnt_n(1)*CoordSys['d1  + p['TFrlPnt_n(3) *CoordSys['d2 -  p['TFrlPnt_n(2)*CoordSys['d3                             ! Position vector from tower-top / base plate (point O) to specified point on  tail-furl axis (point W).
#     dat['rWI']   =     p['rWIxn*CoordSys['tf1 +      p['rWIzn*CoordSys['tf2 -     p['rWIyn*CoordSys['tf3                            ! Position vector from specified point on  tail-furl axis (point W) to tail boom center of mass     (point I).
#     dat['rWJ']   =     p['rWJxn*CoordSys['tf1 +      p['rWJzn*CoordSys['tf2 -     p['rWJyn*CoordSys['tf3                            ! Position vector from specified point on  tail-furl axis (point W) to tail fin  center of mass     (point J).
#     dat['rPC']   = dat['rPQ'] + dat['rQC']                                                                                # Position vector from teeter pin (point P) to hub center of mass (point C).
    dat['rT0O']  = dat['rZO'] - dat['rZT0']   # Position vector from the tower base (point T(0)) to tower-top / base plate (point O).
    dat['rO']    = dat['rZ']  + dat['rZO']    # Position vector from inertial frame origin to tower-top / base plate (point O).
#     dat['rV']    = dat['rO']  + dat['rOV']                                                                                # Position vector from inertial frame origin to specified point on rotor-furl axis (point V)
#     !dat['rP']   = dat['rO']  + dat['rOV'] + dat['rVP']                                                                   # Position vector from inertial frame origin to teeter pin (point P).
#     dat['rP']    = dat['rV']  + dat['rVP']                                                                                # Position vector from inertial frame origin to teeter pin (point P).
#     dat['rQ']    = dat['rP']  + dat['rPQ']                                                                                # Position vector from inertial frame origin to apex of rotation (point Q).
#     dat['rJ']    = dat['rO']  + dat['rOW']+ dat['rWJ']                                                                    # Position vector from inertial frame origin to tail fin center of mass (point J).'
# 
# 
#     DO K = 1,p['NumBl ! Loop through all blades
#        ! Calculate the position vector of the tip:
#        dat['rS0S(:,K,p['TipNode) = ( p['TwistedSF(K,1,1,p['TipNode,0)*QT( DOF_BF(K,1) ) &                                       ! Position vector from the blade root (point S(0)) to the blade tip (point S(p['BldFlexL)).
#                                      + p['TwistedSF(K,1,2,p['TipNode,0)*QT( DOF_BF(K,2) ) &
#                                      + p['TwistedSF(K,1,3,p['TipNode,0)*QT( DOF_BE(K,1) )                     )*CoordSys['j1(K,:) &
#                                    + ( p['TwistedSF(K,2,1,p['TipNode,0)*QT( DOF_BF(K,1) ) &
#                                      + p['TwistedSF(K,2,2,p['TipNode,0)*QT( DOF_BF(K,2) ) &
#                                      + p['TwistedSF(K,2,3,p['TipNode,0)*QT( DOF_BE(K,1) )                     )*CoordSys['j2(K,:) &
#                                    + ( p['BldFlexL - 0.5* &
#                                    (      p['AxRedBld(K,1,1,p['TipNode)*QT( DOF_BF(K,1) )*QT( DOF_BF(K,1) ) &
#                                      +    p['AxRedBld(K,2,2,p['TipNode)*QT( DOF_BF(K,2) )*QT( DOF_BF(K,2) ) &
#                                      +    p['AxRedBld(K,3,3,p['TipNode)*QT( DOF_BE(K,1) )*QT( DOF_BE(K,1) ) &
#                                      + 2.*p['AxRedBld(K,1,2,p['TipNode)*QT( DOF_BF(K,1) )*QT( DOF_BF(K,2) ) &
#                                      + 2.*p['AxRedBld(K,2,3,p['TipNode)*QT( DOF_BF(K,2) )*QT( DOF_BE(K,1) ) &
#                                      + 2.*p['AxRedBld(K,1,3,p['TipNode)*QT( DOF_BF(K,1) )*QT( DOF_BE(K,1) ) ) )*CoordSys['j3(K,:)
#        dat['rQS (:,K,p['TipNode) = dat['rS0S(:,K,p['TipNode) + p['HubRad*CoordSys['j3(K,:)                                      ! Position vector from apex of rotation (point Q) to the blade tip (point S(p['BldFlexL)).
#        dat['rS  (:,K,p['TipNode) = dat['rQS (:,K,p['TipNode) + dat['rQ                                                     ! Position vector from inertial frame origin      to the blade tip (point S(p['BldFlexL)).
#        
#        ! position vectors for blade root node:
#        dat['rQS (:,K,0) = p['HubRad*CoordSys['j3(K,:)    
#        dat['rS  (:,K,0) = p['HubRad*CoordSys['j3(K,:) + dat['rQ
#           ! Calculate the position vector from the teeter pin to the blade root:
#        dat['rPS0(:,K) = dat['rPQ + p['HubRad*CoordSys['j3(K,:)   ! Position vector from teeter pin (point P) to blade root (point S(0)).
#        DO J = 1,p['BldNodes ! Loop through the blade nodes / elements
#        ! Calculate the position vector of the current node:
#           dat['rS0S(:,K,J) = (  p['TwistedSF(K,1,1,J,0)*QT( DOF_BF(K,1) ) &                                                   ! Position vector from the blade root (point S(0)) to the current node (point S(RNodes(J)).
#                                  + p['TwistedSF(K,1,2,J,0)*QT( DOF_BF(K,2) ) &
#                                  + p['TwistedSF(K,1,3,J,0)*QT( DOF_BE(K,1) )                          )*CoordSys['j1(K,:) &
#                              + (   p['TwistedSF(K,2,1,J,0)*QT( DOF_BF(K,1) ) &
#                                  + p['TwistedSF(K,2,2,J,0)*QT( DOF_BF(K,2) ) &
#                                  + p['TwistedSF(K,2,3,J,0)*QT( DOF_BE(K,1) )                          )*CoordSys['j2(K,:) &
#                              + (  p['RNodes(J) - 0.5* &
#                                (      p['AxRedBld(K,1,1,J)*QT( DOF_BF(K,1) )*QT( DOF_BF(K,1) ) &
#                                 +     p['AxRedBld(K,2,2,J)*QT( DOF_BF(K,2) )*QT( DOF_BF(K,2) ) &
#                                 +     p['AxRedBld(K,3,3,J)*QT( DOF_BE(K,1) )*QT( DOF_BE(K,1) ) &
#                                 + 2.0*p['AxRedBld(K,1,2,J)*QT( DOF_BF(K,1) )*QT( DOF_BF(K,2) ) &
#                                 + 2.0*p['AxRedBld(K,2,3,J)*QT( DOF_BF(K,2) )*QT( DOF_BE(K,1) ) &
#                                 + 2.0*p['AxRedBld(K,1,3,J)*QT( DOF_BF(K,1) )*QT( DOF_BE(K,1) )    ) )*CoordSys['j3(K,:)
#           dat['rQS (:,K,J) = dat['rS0S(:,K,J) + p['HubRad*CoordSys['j3(K,:)                                                ! Position vector from apex of rotation (point Q) to the current node (point S(RNodes(J)).
#           dat['rS  (:,K,J) = dat['rQS (:,K,J) + dat['rQ                                                               ! Position vector from inertial frame origin      to the current node (point S(RNodes(J)).
#        END DO !J = 1,p['BldNodes ! Loop through the blade nodes / elements
#     END DO !K = 1,p['NumBl
    # --- Tower element positions
    TwrNodes = p['TwrFASF'].shape[1] - 2
    dat['rZT'] = np.zeros((TwrNodes+1,3))
    dat['rT0T'] = np.zeros((TwrNodes,3))
    dat['rT'] = np.zeros((TwrNodes,3))
    dat['rZT'][0,:] = dat['rZT0']
    # TODO vectorize
    # TODO TODO TODO check that this is correct due to messy ED indexing (some starting at 0 some at 1)
    for j in range(TwrNodes):
        jj=j+1
        # Calculate the position vector of the current node:
        dat['rT0T'][j,:] = ( p['TwrFASF'][0,jj,0]*QT[DOF_TFA1] + p['TwrFASF'][1,jj,0]*QT[DOF_TFA2] )*CoordSys['a1']# Position vector from base of flexible portion of tower (point T(0)) to current node (point T(j)).
        dat['rT0T'][j,:] +=                     + ( p['HNodes'][j] - 0.5*(     p['AxRedTFA'][0,0,j]*QT[DOF_TFA1]*QT[DOF_TFA1] \
                                               +     p['AxRedTFA'][1,1,jj]*QT[DOF_TFA2]*QT[DOF_TFA2] \
                                               + 2.0*p['AxRedTFA'][0,1,jj]*QT[DOF_TFA1]*QT[DOF_TFA2] \
                                               +     p['AxRedTSS'][0,0,jj]*QT[DOF_TSS1]*QT[DOF_TSS1] \
                                               +     p['AxRedTSS'][1,1,jj]*QT[DOF_TSS2]*QT[DOF_TSS2] \
                                               + 2.0*p['AxRedTSS'][0,1,jj]*QT[DOF_TSS1]*QT[DOF_TSS2]  ) )*CoordSys['a2']
        dat['rT0T'][j,:] += ( p['TwrSSSF'][0,jj,0]*QT[DOF_TSS1] + p['TwrSSSF'][1,jj,0]*QT[DOF_TSS2] )*CoordSys['a3']
        dat['rZT'][jj,:] = dat['rZT0'] + dat['rT0T'][j,:]  # Position vector from platform reference (point Z) to the current node (point T(HNodes(j)).
        dat['rT'][j,:]  = dat['rZ']   + dat['rZT'][jj,:]   # Position vector from inertial frame origin        to the current node (point T(HNodes(j)).


    # --- IEC
    if IEC is None:
        IEC = dict()
    IEC['r_F0'] = np.array([0,0,p['PtfmRefzt']])
    IEC['r_F']  = EDVec2IEC(dat['rZ']) + IEC['r_F0']    # ED ReftPtfm point
    IEC['r_FT'] = EDVec2IEC(dat['rZT0'])  # Tower base from Floater
    IEC['r_T']  = IEC['r_F'] + IEC['r_FT']
    IEC['r_TN'] = EDVec2IEC(dat['rT0O'])  # Tower top / base plate/ nacelle origin from tower base
    IEC['r_N']  = EDVec2IEC(dat['rO'])+ IEC['r_F0']    # Tower top / base plate/ nacelle origin
    IEC['r_NGn'] = EDVec2IEC(dat['rOU'])   # From Tower top to nacelle COG
    IEC['r_Gn'] = IEC['r_N'] + IEC['r_NGn'] # Nacelle COG
    IEC['r_Ts'] = EDVec2IEC(dat['rT'], IEC['r_F0']  )

    return dat, IEC


def ED_AngPosVelPAcc(q=None, qd=None, qDict=None, qdDict=None, CoordSys=None, p=None, dat=None, IEC=None):
    """ 
    See ElastoDyn.f90 CalculateAngularPosVelPAcc
    This routine is used to calculate the angular positions, velocities, and partial accelerations stored in other states that are used in
    both the CalcOutput and CalcContStateDeriv routines.
    SUBROUTINE CalculateAngularPosVelPAcc( p, x, CoordSys, dat )

    """
    if qDict is not None:
        QT = ED_qDict2q(qDict)
    else:
        QT = q
    if qdDict is not None:
        QDT = ED_qDict2q(qdDict)
    else:
        QDT = qd

    if dat is None:
        dat = dict()
#    REAL(ReKi)                   :: AngVelHM  (3)                                   ! Angular velocity of eleMent J of blade K (body M) in the hub (body H).
#    REAL(ReKi)                   :: AngAccELt (3)                                   ! Portion of the angular acceleration of the low-speed shaft (body L) in the inertia frame (body E for earth) associated with everything but the QD2T()'s.
    #NDOF=ED_MaxDOFs
    NDOF=11 # TODO
    TwrNodes = p['TwrFASF'].shape[1] - 2
    NumBl=3
    # These angular velocities are allocated to start numbering a dimension with 0 instead of 1:
    dat['PAngVelEB']=np.zeros((NDOF,2,3))
    dat['PAngVelER']=np.zeros((NDOF,2,3))
    dat['PAngVelEX']=np.zeros((NDOF,2,3))
    dat['PAngVelEA']=np.zeros((NDOF,2,3))
    dat['PAngVelEN']=np.zeros((NDOF,2,3))
    #dat['PAngVelEG']=np.zeros((NDOF,2,3))
    #dat['PAngVelEH']=np.zeros((NDOF,2,3))
    #dat['PAngVelEL']=np.zeros((NDOF,2,3))
    dat['PAngVelEF']=np.zeros((TwrNodes+1, NDOF,2,3))
    #dat['PAngVelEM']=np.zeros((NumBl,0:p%TipNode,NDOF,2,3))

    dat['AngPosXF']=np.zeros((TwrNodes+1,3))
    dat['AngPosEF']=np.zeros((TwrNodes+1,3))
    dat['AngVelEF']=np.zeros((TwrNodes+1,3))



    # --------------------------------------------------------------------------------}
    # --- Angular and partial angular velocities
    # --------------------------------------------------------------------------------{
    #   Define the angular and partial angular velocities of all of the rigid bodies in the inertia frame:
    #   NOTE: PAngVelEN(I,D,:) = the Dth-derivative of the partial angular velocity of DOF I for body N in body E.
    # --- Platform "X"
    dat['PAngVelEX'][       :,0,:] = 0.0
    dat['PAngVelEX'][DOF_R   ,0,:] =  CoordSys['z1']
    dat['PAngVelEX'][DOF_P   ,0,:] = -CoordSys['z3']
    dat['PAngVelEX'][DOF_Y   ,0,:] =  CoordSys['z2']
    dat['AngVelEX'] = QDT[DOF_R   ]*dat['PAngVelEX'][DOF_R,0,:] + QDT[DOF_P]*dat['PAngVelEX'][DOF_P,0,:]  + QDT[DOF_Y]*dat['PAngVelEX'][DOF_Y,0,:]
    dat['AngPosEX'] = QT [DOF_R   ]*dat['PAngVelEX'][DOF_R,0,:] + QT [DOF_P]*dat['PAngVelEX'][DOF_P,0,:]  + QT [DOF_Y]*dat['PAngVelEX'][DOF_Y,0,:]
    # --- Tower top "B"
    dat['PAngVelEB'][       :,0,:] = dat['PAngVelEX'][:,0,:]
    dat['PAngVelEB'][DOF_TFA1,0,:] = -p['TwrFASF'][0,-1,1]*CoordSys['a3']
    dat['PAngVelEB'][DOF_TSS1,0,:] =  p['TwrSSSF'][0,-1,1]*CoordSys['a1']
    dat['PAngVelEB'][DOF_TFA2,0,:] = -p['TwrFASF'][1,-1,1]*CoordSys['a3']
    dat['PAngVelEB'][DOF_TSS2,0,:] =  p['TwrSSSF'][1,-1,1]*CoordSys['a1']
    dat['AngVelEB']                =  dat['AngVelEX'].copy()
    dat['AngVelEB']                += QDT[DOF_TFA1]*dat['PAngVelEB'][DOF_TFA1,0,:]
    dat['AngVelEB']                += QDT[DOF_TSS1]*dat['PAngVelEB'][DOF_TSS1,0,:]
    dat['AngVelEB']                += QDT[DOF_TFA2]*dat['PAngVelEB'][DOF_TFA2,0,:]
    dat['AngVelEB']                += QDT[DOF_TSS2]*dat['PAngVelEB'][DOF_TSS2,0,:]
    dat['AngPosXB']                 = QT [DOF_TFA1]*dat['PAngVelEB'][DOF_TFA1,0,:]
    dat['AngPosXB']                += QT [DOF_TSS1]*dat['PAngVelEB'][DOF_TSS1,0,:]
    dat['AngPosXB']                += QT [DOF_TFA2]*dat['PAngVelEB'][DOF_TFA2,0,:]
    dat['AngPosXB']                += QT [DOF_TSS2]*dat['PAngVelEB'][DOF_TSS2,0,:]
    # --- Nacelle "N" (yawed tower top)
    dat['PAngVelEN'][       :,0,:]= dat['PAngVelEB'][:,0,:]
    dat['PAngVelEN'][DOF_Yaw ,0,:]= CoordSys['d2']
    dat['AngVelEN']               = dat['AngVelEB'] + QDT[DOF_Yaw ]*dat['PAngVelEN'][DOF_Yaw ,0,:]
# 
#     dat['PAngVelER(       :,0,:)= dat['PAngVelEN(:,0,:)
#     dat['PAngVelER(DOF_RFrl,0,:)= CoordSys['rfa
#     dat['AngVelER               = dat['AngVelEN + QDT(DOF_RFrl)*dat['PAngVelER(DOF_RFrl,0,:)
# 
#     dat['PAngVelEL(       :,0,:)= dat['PAngVelER(:,0,:)
#     dat['PAngVelEL(DOF_GeAz,0,:)= CoordSys['c1
#     dat['PAngVelEL(DOF_DrTr,0,:)= CoordSys['c1
#     dat['AngVelEL               = dat['AngVelER + QDT(DOF_GeAz)*dat['PAngVelEL(DOF_GeAz,0,:) &
#                                                             + QDT(DOF_DrTr)*dat['PAngVelEL(DOF_DrTr,0,:)
# 
#     dat['PAngVelEH(       :,0,:)= dat['PAngVelEL(:,0,:)
#     dat['AngVelEH               = dat['AngVelEL
# IF ( p['NumBl == 2 )  THEN ! 2-blader
#    dat['PAngVelEH(DOF_Teet,0,:)= CoordSys['f2
#    dat['AngVelEH               = dat['AngVelEH + QDT(DOF_Teet)*dat['PAngVelEH(DOF_Teet,0,:)
# ENDIF
#  
#     dat['PAngVelEG(       :,0,:) = dat['PAngVelER(:,0,:)
#     dat['PAngVelEG(DOF_GeAz,0,:) = p['GBRatio*CoordSys['c1
#     dat['AngVelEG                = dat['AngVelER + QDT(DOF_GeAz)*dat['PAngVelEG(DOF_GeAz,0,:)
#  
#     dat['PAngVelEA(       :,0,:) = dat['PAngVelEN(:,0,:)
#     dat['PAngVelEA(DOF_TFrl,0,:) = CoordSys['tfa
#     dat['AngVelEA                = dat['AngVelEN + QDT(DOF_TFrl)*dat['PAngVelEA(DOF_TFrl,0,:)
#  
#  
#  
#     ! Define the 1st derivatives of the partial angular velocities of all
#     !   of the rigid bodies in the inertia frame and the portion of the angular
#     !   acceleration of the rigid bodies in the inertia frame associated with
#     !   everything but the QD2T()'s:
#     dat['PAngVelEX(       :,1,:) = 0.0
#     dat['AngAccEXt               = 0.0
#     dat['PAngVelEB(       :,1,:) =                  dat['PAngVelEX(:,1,:)
#     dat['PAngVelEB(DOF_TFA1,1,:) = CROSS_PRODUCT(   dat['AngVelEX,                   dat['PAngVelEB(DOF_TFA1,0,:) )
#     dat['PAngVelEB(DOF_TSS1,1,:) = CROSS_PRODUCT(   dat['AngVelEX,                   dat['PAngVelEB(DOF_TSS1,0,:) )
#     dat['PAngVelEB(DOF_TFA2,1,:) = CROSS_PRODUCT(   dat['AngVelEX,                   dat['PAngVelEB(DOF_TFA2,0,:) )
#     dat['PAngVelEB(DOF_TSS2,1,:) = CROSS_PRODUCT(   dat['AngVelEX,                   dat['PAngVelEB(DOF_TSS2,0,:) )
#     dat['AngAccEBt               =                  dat['AngAccEXt + QDT(DOF_TFA1)*dat['PAngVelEB(DOF_TFA1,1,:) &
#                                                                          + QDT(DOF_TSS1)*dat['PAngVelEB(DOF_TSS1,1,:) &
#                                                                          + QDT(DOF_TFA2)*dat['PAngVelEB(DOF_TFA2,1,:) &
#                                                                          + QDT(DOF_TSS2)*dat['PAngVelEB(DOF_TSS2,1,:)
#     dat['PAngVelEN(       :,1,:) =                 dat['PAngVelEB(:,1,:)
#     dat['PAngVelEN(DOF_Yaw ,1,:) = CROSS_PRODUCT(  dat['AngVelEB,                    dat['PAngVelEN(DOF_Yaw ,0,:) )
#     dat['AngAccENt               =                 dat['AngAccEBt  + QDT(DOF_Yaw )*dat['PAngVelEN(DOF_Yaw ,1,:)
#  
#     dat['PAngVelER(       :,1,:) =                 dat['PAngVelEN(:,1,:)
#     dat['PAngVelER(DOF_RFrl,1,:) = CROSS_PRODUCT(  dat['AngVelEN,                    dat['PAngVelER(DOF_RFrl,0,:) )
#     dat['AngAccERt               =                 dat['AngAccENt  + QDT(DOF_RFrl)*dat['PAngVelER(DOF_RFrl,1,:)
#  
#     dat['PAngVelEL(       :,1,:) =                 dat['PAngVelER(:,1,:)
#     dat['PAngVelEL(DOF_GeAz,1,:) = CROSS_PRODUCT(  dat['AngVelER,                    dat['PAngVelEL(DOF_GeAz,0,:) )
#     dat['PAngVelEL(DOF_DrTr,1,:) = CROSS_PRODUCT(  dat['AngVelER,                    dat['PAngVelEL(DOF_DrTr,0,:) )
#             AngAccELt               =                 dat['AngAccERt  + QDT(DOF_GeAz)*dat['PAngVelEL(DOF_GeAz,1,:) &
#                                                                          + QDT(DOF_DrTr)*dat['PAngVelEL(DOF_DrTr,1,:)
#  
#     dat['PAngVelEH(       :,1,:) = dat['PAngVelEL(:,1,:)
#     dat['AngAccEHt               =                  AngAccELt
#  IF ( p['NumBl == 2 )  THEN ! 2-blader
#     dat['PAngVelEH(DOF_Teet,1,:) = CROSS_PRODUCT(  dat['AngVelEH,                    dat['PAngVelEH(DOF_Teet,0,:) )
#     dat['AngAccEHt               =                 dat['AngAccEHt   + QDT(DOF_Teet)*dat['PAngVelEH(DOF_Teet,1,:)
#  ENDIF
#  
#     dat['PAngVelEG(       :,1,:) = dat['PAngVelER(:,1,:)
#     dat['PAngVelEG(DOF_GeAz,1,:) = CROSS_PRODUCT(  dat['AngVelER,                    dat['PAngVelEG(DOF_GeAz,0,:) )
#     dat['AngAccEGt               =                 dat['AngAccERt  + QDT(DOF_GeAz)*dat['PAngVelEG(DOF_GeAz,1,:)
#  
#     dat['PAngVelEA(       :,1,:) = dat['PAngVelEN(:,1,:)
#     dat['PAngVelEA(DOF_TFrl,1,:) = CROSS_PRODUCT(  dat['AngVelEN,                    dat['PAngVelEA(DOF_TFrl,0,:) )
#     dat['AngAccEAt               =                 dat['AngAccENt  + QDT(DOF_TFrl)*dat['PAngVelEA(DOF_TFrl,1,:)

#     DO K = 1,p['NumBl ! Loop through all blades
#        DO J = 0,p['TipNode ! Loop through the blade nodes / elements
#        ! Define the partial angular velocities of the current node (body M(RNodes(J))) in the inertia frame:
#        ! NOTE: PAngVelEM(K,J,I,D,:) = the Dth-derivative of the partial angular velocity
#        !   of DOF I for body M of blade K, element J in body E.
#  
#           dat['PAngVelEM(K,J,          :,0,:) = dat['PAngVelEH(:,0,:)
#           dat['PAngVelEM(K,J,DOF_BF(K,1),0,:) = - p['TwistedSF(K,2,1,J,1)*CoordSys['j1(K,:) &
#                                                    + p['TwistedSF(K,1,1,J,1)*CoordSys['j2(K,:)
#           dat['PAngVelEM(K,J,DOF_BF(K,2),0,:) = - p['TwistedSF(K,2,2,J,1)*CoordSys['j1(K,:) &
#                                                    + p['TwistedSF(K,1,2,J,1)*CoordSys['j2(K,:)
#           dat['PAngVelEM(K,J,DOF_BE(K,1),0,:) = - p['TwistedSF(K,2,3,J,1)*CoordSys['j1(K,:) &
#                                                    + p['TwistedSF(K,1,3,J,1)*CoordSys['j2(K,:)
#                                        AngVelHM  =     QDT(DOF_BF(K,1))*dat['PAngVelEM(K,J,DOF_BF(K,1),0,:) &
#                                                      + QDT(DOF_BF(K,2))*dat['PAngVelEM(K,J,DOF_BF(K,2),0,:) &
#                                                      + QDT(DOF_BE(K,1))*dat['PAngVelEM(K,J,DOF_BE(K,1),0,:)
#            dat['AngVelEM(:,J,K              ) =  dat['AngVelEH + AngVelHM
#            dat['AngPosHM(:,K,J              ) =     QT (DOF_BF(K,1))*dat['PAngVelEM(K,J,DOF_BF(K,1),0,:) &
#                                                      + QT (DOF_BF(K,2))*dat['PAngVelEM(K,J,DOF_BF(K,2),0,:) &
#                                                      + QT (DOF_BE(K,1))*dat['PAngVelEM(K,J,DOF_BE(K,1),0,:)
#           dat['AngAccEKt(:,J              ,K) =  dat['AngAccEHt + QDT(DOF_BF(K,1))*dat['PAngVelEM(K,J,DOF_BF(K,1),1,:) & 
#                                                                       + QDT(DOF_BF(K,2))*dat['PAngVelEM(K,J,DOF_BF(K,2),1,:) & 
#                                                                       + QDT(DOF_BE(K,1))*dat['PAngVelEM(K,J,DOF_BE(K,1),1,:)   
#        ! Define the 1st derivatives of the partial angular velocities of the current node (body M(RNodes(J))) in the inertia frame:
#     ! NOTE: These are currently unused by the code, therefore, they need not
#     !       be calculated.  Thus, they are currently commented out.  If it
#     !       turns out that they are ever needed (i.e., if inertias of the
#     !       blade elements are ever added, etc...) simply uncomment out these computations:
#     !      dat['PAngVelEM(K,J,          :,1,:) = dat['PAngVelEH(:,1,:)
#     !      dat['PAngVelEM(K,J,DOF_BF(K,1),1,:) = CROSS_PRODUCT(   dat['AngVelEH, PAngVelEM(K,J,DOF_BF(K,1),0,:) )
#     !      dat['PAngVelEM(K,J,DOF_BF(K,2),1,:) = CROSS_PRODUCT(   dat['AngVelEH, PAngVelEM(K,J,DOF_BF(K,2),0,:) )
#     !      dat['PAngVelEM(K,J,DOF_BE(K,1),1,:) = CROSS_PRODUCT(   dat['AngVelEH, PAngVelEM(K,J,DOF_BE(K,1),0,:) )
#        END DO !J = 1,p['BldNodes ! Loop through the blade nodes / elements
#     END DO !K = 1,p['NumBl
    # --- Tower values:
    for J in range(TwrNodes+1):
#        ! Define the partial angular velocities (and their 1st derivatives) of the
#        !   current node (body F(HNodes(J))  in the inertia frame.
#        ! Also define the overall angular velocity of the current node in the inertia frame.
#        !   Also, define the portion of the angular acceleration of the current node
#        !   in the inertia frame associated with everything but the QD2T()'s:
#        ! NOTE: PAngVelEF(J,I,D,:) = the Dth-derivative of the partial angular velocity
#        !   of DOF I for body F of element J in body E.
       dat['PAngVelEF'] [J,       :,0,:] = dat['PAngVelEX'][:,0,:]
       dat['PAngVelEF'] [J,DOF_TFA1,0,:] = -p['TwrFASF'][0,J,1]*CoordSys['a3'] # Local slope
       dat['PAngVelEF'] [J,DOF_TSS1,0,:] =  p['TwrSSSF'][0,J,1]*CoordSys['a1']
       dat['PAngVelEF'] [J,DOF_TFA2,0,:] = -p['TwrFASF'][1,J,1]*CoordSys['a3']
       dat['PAngVelEF'] [J,DOF_TSS2,0,:] =  p['TwrSSSF'][1,J,1]*CoordSys['a1']
       dat['PAngVelEF'] [J,       :,1,:] = dat['PAngVelEX'][:,1,:]
       dat['PAngVelEF'] [J,DOF_TFA1,1,:] = np.cross(  dat['AngVelEX']  ,  dat['PAngVelEF'][J,DOF_TFA1,0,:] )
       dat['PAngVelEF'] [J,DOF_TSS1,1,:] = np.cross(  dat['AngVelEX']  ,  dat['PAngVelEF'][J,DOF_TSS1,0,:] )
       dat['PAngVelEF'] [J,DOF_TFA2,1,:] = np.cross(  dat['AngVelEX']  ,  dat['PAngVelEF'][J,DOF_TFA2,0,:] )
       dat['PAngVelEF'] [J,DOF_TSS2,1,:] = np.cross(  dat['AngVelEX']  ,  dat['PAngVelEF'][J,DOF_TSS2,0,:] )
       dat['AngVelEF']  [J,:]            =  dat['AngVelEX'].copy()
       dat['AngVelEF']  [J,:]           += QDT[DOF_TFA1]*dat['PAngVelEF'][J,DOF_TFA1,0,:]
       dat['AngVelEF']  [J,:]           += QDT[DOF_TSS1]*dat['PAngVelEF'][J,DOF_TSS1,0,:]
       dat['AngVelEF']  [J,:]           += QDT[DOF_TFA2]*dat['PAngVelEF'][J,DOF_TFA2,0,:]
       dat['AngVelEF']  [J,:]           += QDT[DOF_TSS2]*dat['PAngVelEF'][J,DOF_TSS2,0,:]
       dat['AngPosXF']  [J,:]           =  QT [DOF_TFA1]*dat['PAngVelEF'][J,DOF_TFA1,0,:]
       dat['AngPosXF']  [J,:]           += QT [DOF_TSS1]*dat['PAngVelEF'][J,DOF_TSS1,0,:]
       dat['AngPosXF']  [J,:]           += QT [DOF_TFA2]*dat['PAngVelEF'][J,DOF_TFA2,0,:]
       dat['AngPosXF']  [J,:]           += QT [DOF_TSS2]*dat['PAngVelEF'][J,DOF_TSS2,0,:]
       dat['AngPosEF']  [J,:]           =  dat['AngPosEX']  + dat['AngPosXF'][J,:]
#        dat['AngAccEFt'][J,:]          =  dat['AngAccEXt']
#        dat['AngAccEFt'][J,:]         += QDT[DOF_TFA1]*dat['PAngVelEF'][J,DOF_TFA1,1,:]
#        dat['AngAccEFt'][J,:]         += QDT[DOF_TSS1]*dat['PAngVelEF'][J,DOF_TSS1,1,:]
#        dat['AngAccEFt'][J,:]         += QDT[DOF_TFA2]*dat['PAngVelEF'][J,DOF_TFA2,1,:]
#        dat['AngAccEFt'][J,:]         += QDT[DOF_TSS2]*dat['PAngVelEF'][J,DOF_TSS2,1,:]
    # ---
    if IEC is None:
        IEC = dict()
    IEC['theta_f'] = EDVec2IEC(dat['AngPosEX'])
    IEC['omega_f'] = EDVec2IEC(dat['AngVelEX'])
    IEC['omega_t'] = IEC['omega_f'].copy()
    IEC['theta_fn'] = EDVec2IEC(dat['AngPosXB']) # TODO TODO ADD YAW
    IEC['theta_n'] = IEC['theta_f'] + IEC['theta_fn']
    IEC['omega_n'] = EDVec2IEC(dat['AngVelEN'])

    IEC['omega_Ts'] = EDVec2IEC(dat['AngVelEF']) # Tower nodes ang vel
    IEC['theta_fTs'] = EDVec2IEC(dat['AngPosXF']) # Tower nodes ang pos from platform
    IEC['theta_Ts'] = EDVec2IEC(dat['AngPosEF']) # Tower nodes ang pos from platform

    return dat, IEC


def ED_LinVelPAcc(q=None, qd=None, qDict=None, qdDict=None, CoordSys=None, p=None, dat=None, IEC=None):
    """ 
    See ElastoDyn.f90 CalculateLinearVelPAcc
    !> This routine is used to calculate the linear velocities and accelerations stored in other states that are used in
    !! both the CalcOutput and CalcContStateDeriv routines.
    SUBROUTINE CalculateLinearVelPAcc( p, x, CoordSys, dat )

    BODIES:
        E: the earth/inertial frame
        X: the platform body
        N: the nacelle body
        A: the tail-furl body

    """
    if qDict is not None:
        QT = ED_qDict2q(qDict)
    else:
        QT = q
    if qdDict is not None:
        QDT = ED_qDict2q(qdDict)
    else:
        QDT = qd

    if dat is None:
        dat = dict()
#    ALLOCATE( RtHS%LinVelES( Dims, 0:p%TipNode, p%NumBl ), &
#              RtHS%AngVelEM( Dims, 0:p%TipNode, p%NumBl ), STAT=ErrStat )
#    ! These linear velocities are allocated to start numbering a dimension with 0 instead of 1:
#    ALLOCATE ( RtHS%PLinVelEIMU(p%NDOF,0:1,Dims) , STAT=ErrStat )
#    ALLOCATE ( RtHS%PLinVelEO(p%NDOF,0:1,Dims) , STAT=ErrStat )
#    ALLOCATE ( RtHS%PLinVelES(p%NumBl,0:p%TipNode,p%NDOF,0:1,Dims) , STAT=ErrStat )
#    ALLOCATE ( RtHS%PLinVelET(0:p%TwrNodes,p%NDOF,0:1,Dims) , STAT=ErrStat )
#    ALLOCATE ( RtHS%PLinVelEC(p%NDOF,0:1,3) , STAT=ErrStat )
#    ALLOCATE ( RtHS%PLinVelED(p%NDOF,0:1,3) , STAT=ErrStat )
#    ALLOCATE ( RtHS%PLinVelEI(p%NDOF,0:1,3) , STAT=ErrStat )
#    ALLOCATE ( RtHS%PLinVelEJ(p%NDOF,0:1,3) , STAT=ErrStat )
#    ALLOCATE ( RtHS%PLinVelEP(p%NDOF,0:1,3) , STAT=ErrStat )
#    ALLOCATE ( RtHS%PLinVelEQ(p%NDOF,0:1,3) , STAT=ErrStat )
#               RtHS%PLinVelEW(p%NDOF,0:1,3) , &
#               RtHS%PLinVelEY(p%NDOF,0:1,3) , STAT=ErrStat )
    NDOF=11 # TODO
    TwrNodes = p['TwrFASF'].shape[1] - 2
    NumBl=3

    dat['PLinVelEZ'] = np.zeros((NDOF,2, 3))
    dat['PLinVelEO'] = np.zeros((NDOF,2, 3))
    dat['PLinVelEU'] = np.zeros((NDOF,2, 3))
    dat['PLinVelET'] = np.zeros((TwrNodes+1,NDOF,2, 3))
#    !-------------------------------------------------------------------------------------------------
#    ! Partial linear velocities and accelerations
#    !-------------------------------------------------------------------------------------------------
#       ! Define the partial linear velocities (and their 1st derivatives) of all of
#       !   the points on the wind turbine in the inertia frame that are not
#       !   dependent on the distributed tower or blade parameters.  Also, define
#       !   the portion of the linear acceleration of the points in the inertia
#       !   frame associated with everything but the QD2T()'s:
#       ! NOTE: PLinVelEX(I,D,:) = the Dth-derivative of the partial linear velocity
#       !   of DOF I for point X in body E.
#     EwXXrZY   = CROSS_PRODUCT( dat['AngVelEX, dat['rZY   ) !
#     EwRXrVD   = np.cross( dat['AngVelER'], dat['rVD']   ) # Cross products
#     EwRXrVP   = np.cross( dat['AngVelER'], dat['rVP']   ) # in the following
#     EwHXrPQ   = np.cross( dat['AngVelEH'], dat['rPQ']   ) # DO...LOOPs
#     EwHXrQC   = np.cross( dat['AngVelEH'], dat['rQC']   ) #
#     EwNXrOW   = np.cross( dat['AngVelEN'], dat['rOW']   ) #
#     EwAXrWI   = np.cross( dat['AngVelEA'], dat['rWI']   ) #
#     EwAXrWJ   = np.cross( dat['AngVelEA'], dat['rWJ']   ) #

    # --- Platform reference point "Z"
    dat['PLinVelEZ'][       :,:,:] = 0.0
    dat['PLinVelEZ'][DOF_Sg  ,0,:] =  CoordSys['z1']
    dat['PLinVelEZ'][DOF_Sw  ,0,:] = -CoordSys['z3']
    dat['PLinVelEZ'][DOF_Hv  ,0,:] =  CoordSys['z2']
    dat['LinVelEZ']                =   QDT[DOF_Sg]*dat['PLinVelEZ'][DOF_Sg,0,:]+ QDT[DOF_Sw]*dat['PLinVelEZ'][DOF_Sw,0,:]+ QDT[DOF_Hv]*dat['PLinVelEZ'][DOF_Hv,0,:]

    # --- Platform COG "Y"
#    dat['PLinVelEY(       :,:,:) = dat['PLinVelEZ(:,:,:)
#    DO I = 1,NPX   ! Loop through all DOFs associated with the angular motion of the platform (body X)
#       TmpVec0 = CROSS_PRODUCT( dat['PAngVelEX(PX(I)   ,0,:), dat['rZY  )
#       TmpVec1 = CROSS_PRODUCT( dat['PAngVelEX(PX(I)   ,0,:),     EwXXrZY  )
#       dat['PLinVelEY(PX(I),0,:) = TmpVec0   +                       dat['PLinVelEY(PX(I)   ,0,:)
#       dat['PLinVelEY(PX(I),1,:) = TmpVec1   +                       dat['PLinVelEY(PX(I)   ,1,:)
#       dat['LinAccEYt           = dat['LinAccEYt + x%QDT(PX(I) )*dat['PLinVelEY(PX(I)   ,1,:)

    # --- Tower Top "O"
    dat['PLinVelEO']               = dat['PLinVelEZ'][:,:,:].copy()
    dat['PLinVelEO'][DOF_TFA1,0,:] = CoordSys['a1'] - (p['AxRedTFA'][0,0,-1]* QT[DOF_TFA1] + p['AxRedTFA'][0,1,-1]* QT[DOF_TFA2])*CoordSys['a2']
    dat['PLinVelEO'][DOF_TSS1,0,:] = CoordSys['a3'] - (p['AxRedTSS'][0,0,-1]* QT[DOF_TSS1] + p['AxRedTSS'][0,1,-1]* QT[DOF_TSS2])*CoordSys['a2']
    dat['PLinVelEO'][DOF_TFA2,0,:] = CoordSys['a1'] - (p['AxRedTFA'][1,1,-1]* QT[DOF_TFA2] + p['AxRedTFA'][0,1,-1]* QT[DOF_TFA1])*CoordSys['a2']
    dat['PLinVelEO'][DOF_TSS2,0,:] = CoordSys['a3'] - (p['AxRedTSS'][1,1,-1]* QT[DOF_TSS2] + p['AxRedTSS'][0,1,-1]* QT[DOF_TSS1])*CoordSys['a2']
    TmpVec1 = np.cross(   dat['AngVelEX']   , dat['PLinVelEO'][DOF_TFA1,0,:] )
    TmpVec2 = np.cross(   dat['AngVelEX']   , dat['PLinVelEO'][DOF_TSS1,0,:] )
    TmpVec3 = np.cross(   dat['AngVelEX']   , dat['PLinVelEO'][DOF_TFA2,0,:] )
    TmpVec4 = np.cross(   dat['AngVelEX']   , dat['PLinVelEO'][DOF_TSS2,0,:] )
    dat['PLinVelEO'][DOF_TFA1,1,:] = TmpVec1 - ( p['AxRedTFA'][0,0,-1]*QDT[DOF_TFA1] + p['AxRedTFA'][0,1,-1]*QDT[DOF_TFA2] )*CoordSys['a2']
    dat['PLinVelEO'][DOF_TSS1,1,:] = TmpVec2 - ( p['AxRedTSS'][0,0,-1]*QDT[DOF_TSS1] + p['AxRedTSS'][0,1,-1]*QDT[DOF_TSS2] )*CoordSys['a2']
    dat['PLinVelEO'][DOF_TFA2,1,:] = TmpVec3 - ( p['AxRedTFA'][1,1,-1]*QDT[DOF_TFA2] + p['AxRedTFA'][0,1,-1]*QDT[DOF_TFA1] )*CoordSys['a2']
    dat['PLinVelEO'][DOF_TSS2,1,:] = TmpVec4 - ( p['AxRedTSS'][1,1,-1]*QDT[DOF_TSS2] + p['AxRedTSS'][0,1,-1]*QDT[DOF_TSS1] )*CoordSys['a2']
    LinVelXO               =              QDT[DOF_TFA1]*dat['PLinVelEO'][DOF_TFA1,0,:] \
                                        + QDT[DOF_TSS1]*dat['PLinVelEO'][DOF_TSS1,0,:] \
                                        + QDT[DOF_TFA2]*dat['PLinVelEO'][DOF_TFA2,0,:] \
                                        + QDT[DOF_TSS2]*dat['PLinVelEO'][DOF_TSS2,0,:]
    dat['LinAccEOt']       =              QDT[DOF_TFA1]*dat['PLinVelEO'][DOF_TFA1,1,:] \
                                        + QDT[DOF_TSS1]*dat['PLinVelEO'][DOF_TSS1,1,:] \
                                        + QDT[DOF_TFA2]*dat['PLinVelEO'][DOF_TFA2,1,:] \
                                        + QDT[DOF_TSS2]*dat['PLinVelEO'][DOF_TSS2,1,:]

    dat['LinVelEO'] = LinVelXO + dat['LinVelEZ'].copy()
    EwXXrZO   = np.cross( dat['AngVelEX'], dat['rZO']   ) #
    for I in range(NPX): # Loop through all DOFs associated with the angular motion of the platform (body X)
        TmpVec0 = np.cross( dat['PAngVelEX'][PX[I]   ,0,:], dat['rZO']              )
        TmpVec1 = np.cross( dat['PAngVelEX'][PX[I]   ,0,:], EwXXrZO + LinVelXO      )
        dat['PLinVelEO'][PX[I],0,:] = TmpVec0    +                       dat['PLinVelEO'][PX[I]   ,0,:]
        dat['PLinVelEO'][PX[I],1,:] = TmpVec1    +                       dat['PLinVelEO'][PX[I]   ,1,:]
        dat['LinVelEO']           +=  QDT[PX[I] ]*dat['PLinVelEO'][PX[I],0,:]
        dat['LinAccEOt']          +=  QDT[PX[I] ]*dat['PLinVelEO'][PX[I],1,:]

    # --- Nacelle COG "U"
    dat['PLinVelEU'] = dat['PLinVelEO'].copy()
    dat['LinVelEU']  = dat['LinVelEZ'].copy()        # NOTE: MANU
    dat['LinAccEUt']  = 0
    EwNXrOU   = np.cross( dat['AngVelEN'], dat['rOU']   ) #
    for I in range(NPN): # Loop through all DOFs associated with the angular motion of the nacelle (body N)
       TmpVec0 = np.cross( dat['PAngVelEN'][PN[I]   ,0,:], dat['rOU']               )
       TmpVec1 = np.cross( dat['PAngVelEN'][PN[I]   ,0,:],  EwNXrOU                 )
       TmpVec2 = np.cross( dat['PAngVelEN'][PN[I]   ,1,:], dat['rOU']               )
       dat['PLinVelEU'][PN[I],0,:] = TmpVec0    +               dat['PLinVelEU'][PN[I]   ,0,:]
       dat['PLinVelEU'][PN[I],1,:] = TmpVec1    + TmpVec2 +     dat['PLinVelEU'][PN[I]   ,1,:]
       dat['LinVelEU']             +=  QDT[PN[I] ]*dat['PLinVelEU'][PN[I],0,:]  # NOTE: TODO TODO THIS IS MANU TRIAL
       dat['LinAccEUt']            +=  QDT[PN[I] ]*dat['PLinVelEU'][PN[I],1,:]

    # --- Rotor Furl "V"
#     dat['PLinVelEV'] = dat['PLinVelEO'].copy()
#     EwNXrOV   = np.cross( dat['AngVelEN'], dat['rOV']   ) #
#     for I in range(NPN): # Loop through all DOFs associated with the angular motion of the nacelle (body N)
#        TmpVec0 = np.cross( dat['PAngVelEN(PN[I]   ,0,:), dat['rOV                 )
#        TmpVec1 = np.cross( dat['PAngVelEN(PN[I]   ,0,:),     EwNXrOV                 )
#        TmpVec2 = np.cross( dat['PAngVelEN(PN[I]   ,1,:), dat['rOV                 )
#        dat['PLinVelEV'][PN[I],0,:] = TmpVec0    +               dat['PLinVelEV[PN[I],0,:]
#        dat['PLinVelEV'][PN[I],1,:] = TmpVec1    + TmpVec2 +     dat['PLinVelEV[PN[I],1,:]
#        LinAccEVt                   += QDT(PN[I] )*dat['PLinVelEV'][PN[I]   ,1,:]
#  
#     dat['PLinVelED(       :,:,:) = dat['PLinVelEV'].copy()
#     for I in range(NPR): # Loop through all DOFs associated with the angular motion of the structure that furls with the rotor (not including rotor) (body R)
#        TmpVec0 = CROSS_PRODUCT( dat['PAngVelER(PR[I]   ,0,:), dat['rVD                 )
#        TmpVec1 = CROSS_PRODUCT( dat['PAngVelER(PR[I]   ,0,:),     EwRXrVD                 )
#        TmpVec2 = CROSS_PRODUCT( dat['PAngVelER(PR[I]   ,1,:), dat['rVD                 )
#        dat['PLinVelED(PR[I],0,:) = TmpVec0    +                       dat['PLinVelED'][PR[I]   ,0,:]
#        dat['PLinVelED(PR[I],1,:) = TmpVec1    + TmpVec2 +             dat['PLinVelED'][PR[I]   ,1,:]
#        dat['LinAccEDt            += QDT(PR[I] )*dat['PLinVelED'][PR[I]   ,1,:]
#  
    # --- Nacelle IMU
#     dat['PLinVelEIMU'] = dat['PLinVelEV'].copy()
#     dat['LinVelEIMU']  =  dat['LinVelEZ'].copy()
#     EwRXrVIMU = np.cross( dat['AngVelER'], dat['rVIMU'] ) # that are used
#     for I in range(NPR): # Loop through all DOFs associated with the angular motion of the structure that furls with the rotor (not including rotor) (body R)
#        TmpVec0 = CROSS_PRODUCT( dat['PAngVelER(PR[I]   ,0,:), dat['rVIMU               )
#        TmpVec1 = CROSS_PRODUCT( dat['PAngVelER(PR[I]   ,0,:),     EwRXrVIMU               )
#        TmpVec2 = CROSS_PRODUCT( dat['PAngVelER(PR[I]   ,1,:), dat['rVIMU               )
#        dat['PLinVelEIMU(PR[I],0,:) = TmpVec0    +                         dat['PLinVelEIMU(PR[I] ,0,:)
#        dat['PLinVelEIMU(PR[I],1,:) = TmpVec1    + TmpVec2 +               dat['PLinVelEIMU(PR[I] ,1,:)
#        dat['LinVelEIMU']           +=  QDT(PR[I] )*dat['PLinVelEIMU(PR[I] ,0,:)
#        dat['LinAccEIMUt']          +=  QDT(PR[I] )*dat['PLinVelEIMU(PR[I] ,1,:)
#  
#     dat['PLinVelEP(       :,:,:) = dat['PLinVelEV(:,:,:)
#     for I in range(NPR): # Loop through all DOFs associated with the angular motion of the structure that furls with the rotor (not including rotor) (body R)
#        TmpVec0 = CROSS_PRODUCT(             dat['PAngVelER(PR[I]   ,0,:),     dat['rVP                 )
#        TmpVec1 = CROSS_PRODUCT(             dat['PAngVelER(PR[I]   ,0,:), EwRXrVP                 )
#        TmpVec2 = CROSS_PRODUCT(             dat['PAngVelER(PR[I]   ,1,:),     dat['rVP                 )
#        dat['PLinVelEP(PR[I],0,:) = TmpVec0    +               dat['PLinVelEP(PR[I]   ,0,:)
#        dat['PLinVelEP(PR[I],1,:) = TmpVec1    + TmpVec2 +     dat['PLinVelEP(PR[I]   ,1,:)
#         LinAccEPt           =  LinAccEPt + QDT(PR[I] )*dat['PLinVelEP(PR[I]   ,1,:)
#  
    # --- Rotor center "Q"
#     dat['PLinVelEQ'] = dat['PLinVelEP'].copy()
#     dat['LinVelEQ']  =  dat['LinVelEZ'].copy()
#     DO I = 1,p['NPH   ! Loop through all DOFs associated with the angular motion of the hub (body H)
#        TmpVec0 = CROSS_PRODUCT( dat['PAngVelEH(p['PH(I)   ,0,:),   dat['rPQ  )
#        TmpVec1 = CROSS_PRODUCT( dat['PAngVelEH(p['PH(I)   ,0,:),       EwHXrPQ  )
#        TmpVec2 = CROSS_PRODUCT( dat['PAngVelEH(p['PH(I)   ,1,:),   dat['rPQ  )
#        dat['PLinVelEQ(p['PH(I),0,:) = TmpVec0    +                 dat['PLinVelEQ(p['PH(I)   ,0,:)
#        dat['PLinVelEQ(p['PH(I),1,:) = TmpVec1    + TmpVec2 +       dat['PLinVelEQ(p['PH(I)   ,1,:)
#        dat['LinVelEQ               +=  QDT(p['PH(I) )*dat['PLinVelEQ(p['PH(I)   ,0,:)
#        LinAccEQt                   +=  QDT(p['PH(I) )*dat['PLinVelEQ(p['PH(I)   ,1,:)
#  
    # --- Hub COG "C"
#     dat['PLinVelEC(       :,:,:) = dat['PLinVelEQ(:,:,:)
#     DO I = 1,p['NPH   ! Loop through all DOFs associated with the angular motion of the hub (body H)
#        TmpVec0 = CROSS_PRODUCT( dat['PAngVelEH(p['PH(I)   ,0,:), dat['rQC )
#        TmpVec1 = CROSS_PRODUCT( dat['PAngVelEH(p['PH(I)   ,0,:),     EwHXrQC )
#        TmpVec2 = CROSS_PRODUCT( dat['PAngVelEH(p['PH(I)   ,1,:), dat['rQC )
#        dat['PLinVelEC(p['PH(I),0,:) = TmpVec0    +                         dat['PLinVelEC(p['PH(I)   ,0,:)
#        dat['PLinVelEC(p['PH(I),1,:) = TmpVec1    + TmpVec2 +               dat['PLinVelEC(p['PH(I)   ,1,:)
#        dat['LinAccECt              =  dat['LinAccECt + QDT(p['PH(I) )*dat['PLinVelEC(p['PH(I)   ,1,:)
#  
#  
#     DO K = 1,p['NumBl ! Loop through all blades
#        DO J = 0,p['TipNode ! Loop through the blade nodes / elements
#        ! Define the partial linear velocities (and their 1st derivatives) of the
#        !   current node (point S(RNodes(J))) in the inertia frame.  Also define
#        !   the overall linear velocity of the current node in the inertia frame.
#        !   Also, define the portion of the linear acceleration of the current node
#        !   in the inertia frame associated with everything but the QD2T()'s:
#  
#           EwHXrQS = CROSS_PRODUCT(  dat['AngVelEH, dat['rQS(:,K,J) )
#  
#           dat['PLinVelES(K,J,          :,:,:) = dat['PLinVelEQ(:,:,:)
#           dat['PLinVelES(K,J,DOF_BF(K,1),0,:) = p['TwistedSF(K,1,1,J,0)                          *CoordSys['j1(K,:) &  !bjj: this line can be optimized
#                                                  + p['TwistedSF(K,2,1,J,0)                          *CoordSys['j2(K,:) &
#                                                  - (   p['AxRedBld(K,1,1,J)*QT ( DOF_BF(K,1) ) &
#                                                      + p['AxRedBld(K,1,2,J)*QT ( DOF_BF(K,2) ) &
#                                                      + p['AxRedBld(K,1,3,J)*QT ( DOF_BE(K,1) )   )*CoordSys['j3(K,:)
#           dat['PLinVelES(K,J,DOF_BE(K,1),0,:) = p['TwistedSF(K,1,3,J,0)                          *CoordSys['j1(K,:) &
#                                                  + p['TwistedSF(K,2,3,J,0)                          *CoordSys['j2(K,:) &
#                                                  - (   p['AxRedBld(K,3,3,J)*QT ( DOF_BE(K,1) ) &
#                                                      + p['AxRedBld(K,2,3,J)*QT ( DOF_BF(K,2) ) &
#                                                      + p['AxRedBld(K,1,3,J)*QT ( DOF_BF(K,1) )   )*CoordSys['j3(K,:)
#           dat['PLinVelES(K,J,DOF_BF(K,2),0,:) = p['TwistedSF(K,1,2,J,0)                          *CoordSys['j1(K,:) &
#                                                  + p['TwistedSF(K,2,2,J,0)                          *CoordSys['j2(K,:) &
#                                                  - (   p['AxRedBld(K,2,2,J)*QT ( DOF_BF(K,2) ) &
#                                                      + p['AxRedBld(K,1,2,J)*QT ( DOF_BF(K,1) ) &
#                                                      + p['AxRedBld(K,2,3,J)*QT ( DOF_BE(K,1) )   )*CoordSys['j3(K,:)
#  
#           TmpVec1 = CROSS_PRODUCT( dat['AngVelEH, dat['PLinVelES(K,J,DOF_BF(K,1),0,:) )
#           TmpVec2 = CROSS_PRODUCT( dat['AngVelEH, dat['PLinVelES(K,J,DOF_BE(K,1),0,:) )
#           TmpVec3 = CROSS_PRODUCT( dat['AngVelEH, dat['PLinVelES(K,J,DOF_BF(K,2),0,:) )
#  
#           dat['PLinVelES(K,J,DOF_BF(K,1),1,:) = TmpVec1 &
#                                                  - (   p['AxRedBld(K,1,1,J)*QDT( DOF_BF(K,1) ) &
#                                                      + p['AxRedBld(K,1,2,J)*QDT( DOF_BF(K,2) ) &
#                                                      + p['AxRedBld(K,1,3,J)*QDT( DOF_BE(K,1) )   )*CoordSys['j3(K,:)
#           dat['PLinVelES(K,J,DOF_BE(K,1),1,:) = TmpVec2 &
#                                                  - (   p['AxRedBld(K,3,3,J)*QDT( DOF_BE(K,1) ) &
#                                                      + p['AxRedBld(K,2,3,J)*QDT( DOF_BF(K,2) ) &
#                                                      + p['AxRedBld(K,1,3,J)*QDT( DOF_BF(K,1) )   )*CoordSys['j3(K,:)
#           dat['PLinVelES(K,J,DOF_BF(K,2),1,:) = TmpVec3 &
#                                                  - (   p['AxRedBld(K,2,2,J)*QDT( DOF_BF(K,2) ) &
#                                                      + p['AxRedBld(K,1,2,J)*QDT( DOF_BF(K,1) ) &
#                                                      + p['AxRedBld(K,2,3,J)*QDT( DOF_BE(K,1) )   )*CoordSys['j3(K,:)
#  
#           LinVelHS                 = QDT( DOF_BF(K,1) )*dat['PLinVelES(K,J,DOF_BF(K,1),0,:) &
#                                    + QDT( DOF_BE(K,1) )*dat['PLinVelES(K,J,DOF_BE(K,1),0,:) &
#                                    + QDT( DOF_BF(K,2) )*dat['PLinVelES(K,J,DOF_BF(K,2),0,:)
#           dat['LinAccESt(:,K,J) = QDT( DOF_BF(K,1) )*dat['PLinVelES(K,J,DOF_BF(K,1),1,:) &
#                                    + QDT( DOF_BE(K,1) )*dat['PLinVelES(K,J,DOF_BE(K,1),1,:) &
#                                    + QDT( DOF_BF(K,2) )*dat['PLinVelES(K,J,DOF_BF(K,2),1,:)
#  
#           dat['LinVelES(:,J,K)  = LinVelHS + dat['LinVelEZ
#           DO I = 1,p['NPH   ! Loop through all DOFs associated with the angular motion of the hub (body H)
#  
#              TmpVec0 = CROSS_PRODUCT(   dat['PAngVelEH(p['PH(I),0,:), dat['rQS(:,K,J)            )  !bjj: this line can be optimized
#              TmpVec1 = CROSS_PRODUCT(   dat['PAngVelEH(p['PH(I),0,:),     EwHXrQS        + LinVelHS )  !bjj: this line can be optimized
#              TmpVec2 = CROSS_PRODUCT(   dat['PAngVelEH(p['PH(I),1,:), dat['rQS(:,K,J)            )  !bjj: this line can be optimized
#  
#              dat['PLinVelES(K,J,p['PH(I),0,:) = dat['PLinVelES(K,J,p['PH(I),0,:) + TmpVec0            !bjj: this line can be optimized
#              dat['PLinVelES(K,J,p['PH(I),1,:) = dat['PLinVelES(K,J,p['PH(I),1,:) + TmpVec1 + TmpVec2  !bjj: this line can be optimized
#  
#              dat['LinVelES(:,J,K)          = dat['LinVelES(:,J,K)   + QDT(p['PH(I))*dat['PLinVelES(K,J,p['PH(I),0,:)  !bjj: this line can be optimized
#              dat['LinAccESt(:,K,J)         = dat['LinAccESt(:,K,J)  + QDT(p['PH(I))*dat['PLinVelES(K,J,p['PH(I),1,:)  !bjj: this line can be optimized
#  
#           END DO ! I - all DOFs associated with the angular motion of the hub (body H)
#  
#        END DO !J = 0,p['TipNodes ! Loop through the blade nodes / elements
#        
#        
#     !JASON: USE TipNode HERE INSTEAD OF BldNodes IF YOU ALLOCATE AND DEFINE n1, n2, n3, m1, m2, AND m3 TO USE TipNode.  THIS WILL REQUIRE THAT THE AERODYNAMIC AND STRUCTURAL TWISTS, AeroTwst() AND ThetaS(), BE KNOWN AT THE TIP!!!
#        !IF (.NOT. p['BD4Blades) THEN
#        !   dat['LinVelESm2(K) = DOT_PRODUCT( dat['LinVelES(:,p['TipNode,K), CoordSys['m2(K,p['BldNodes,:) )
#        !END IF
#              
#     END DO !K = 1,p['NumBl
#  
#  
#     dat['PLinVelEW(       :,:,:) = dat['PLinVelEO(:,:,:)
#     DO I = 1,NPN   ! Loop through all DOFs associated with the angular motion of the nacelle (body N)
#  
#        TmpVec0 = CROSS_PRODUCT( dat['PAngVelEN(PN(I)   ,0,:), dat['rOW                 )
#        TmpVec1 = CROSS_PRODUCT( dat['PAngVelEN(PN(I)   ,0,:),     EwNXrOW                 )
#        TmpVec2 = CROSS_PRODUCT( dat['PAngVelEN(PN(I)   ,1,:), dat['rOW                 )
#  
#        dat['PLinVelEW(PN(I),0,:) = TmpVec0    +               dat['PLinVelEW(PN(I)   ,0,:)
#        dat['PLinVelEW(PN(I),1,:) = TmpVec1    + TmpVec2 +     dat['PLinVelEW(PN(I)   ,1,:)
#  
#         LinAccEWt                   =  LinAccEWt + QDT(PN(I) )*dat['PLinVelEW(PN(I)   ,1,:)
#  
#     ENDDO          ! I - all DOFs associated with the angular motion of the nacelle (body N)
#  
#  
#     ! Velocities of point I (tail boom center of mass)
#     dat['PLinVelEI(       :,:,:) = dat['PLinVelEW(:,:,:)
#     DO I = 1,NPA   ! Loop through all DOFs associated with the angular motion of the tail (body A)
#  
#        TmpVec0 = CROSS_PRODUCT( dat['PAngVelEA(PA(I)   ,0,:), dat['rWI                 )
#        TmpVec1 = CROSS_PRODUCT( dat['PAngVelEA(PA(I)   ,0,:),     EwAXrWI                 )
#        TmpVec2 = CROSS_PRODUCT( dat['PAngVelEA(PA(I)   ,1,:), dat['rWI                 )
#  
#        dat['PLinVelEI(PA(I),0,:) = TmpVec0    +                       dat['PLinVelEI(PA(I)   ,0,:)
#        dat['PLinVelEI(PA(I),1,:) = TmpVec1    + TmpVec2 +             dat['PLinVelEI(PA(I)   ,1,:)
#  
#        dat['LinAccEIt            =  dat['LinAccEIt + QDT(PA(I) )*dat['PLinVelEI(PA(I)   ,1,:)
#  
#     ENDDO          ! I - all DOFs associated with the angular motion of the tail (body A)
#  
#  
#     ! Velocities of point J (tail fin center of mass)
#     dat['PLinVelEJ(       :,:,:) = dat['PLinVelEW(:,:,:)
#     dat['LinVelEJ                = dat['LinVelEZ
#     DO I = 1,NPA   ! Loop through all DOFs associated with the angular motion of the tail (body A)
#  
#        TmpVec0 = CROSS_PRODUCT( dat['PAngVelEA(PA(I)   ,0,:), dat['rWJ                 )
#        TmpVec1 = CROSS_PRODUCT( dat['PAngVelEA(PA(I)   ,0,:),     EwAXrWJ                 )
#        TmpVec2 = CROSS_PRODUCT( dat['PAngVelEA(PA(I)   ,1,:), dat['rWJ                 )
#  
#        dat['PLinVelEJ(PA(I),0,:) = TmpVec0    +               dat['PLinVelEJ(PA(I)   ,0,:)
#        dat['PLinVelEJ(PA(I),1,:) = TmpVec1    + TmpVec2 +     dat['PLinVelEJ(PA(I)   ,1,:)
#  
#         dat['LinVelEJ            =  dat['LinVelEJ  + QDT(PA(I) )*dat['PLinVelEJ(PA(I)   ,0,:)
#         dat['LinAccEJt           =  dat['LinAccEJt + QDT(PA(I) )*dat['PLinVelEJ(PA(I)   ,1,:)
#  
#     ENDDO          ! I - all DOFs associated with the angular motion of the tail (body A)
#  
#  
    # --- Tower nodes
    dat['LinAccETt'] = np.zeros((TwrNodes+1,3))
    dat['LinVelET'] = np.zeros((TwrNodes+1,3))
    for J in range(TwrNodes+1): #Loop through the tower nodes / elements
       # Define the partial linear velocities (and their 1st derivatives) of the current node (point T(HNodes(J))) in the inertia frame.
       #  Also define the overall linear velocity of the current node in the inertia frame.
       #  Also, define the portion of the linear acceleration of the current node in the inertia frame associated with
       #    everything but the QD2T()'s:
       EwXXrZT                   = np.cross(  dat['AngVelEX'], dat['rZT'][J,:] )
       dat['PLinVelET'][J,       :,:,:] = dat['PLinVelEZ'].copy()
       dat['PLinVelET'][J,DOF_TFA1,0,:] = p['TwrFASF'][0,J,0]*CoordSys['a1'] - (   p['AxRedTFA'][0,0,J]* QT[DOF_TFA1]  + p['AxRedTFA'][0,1,J]* QT[DOF_TFA2]   )*CoordSys['a2']
       dat['PLinVelET'][J,DOF_TSS1,0,:] = p['TwrSSSF'][0,J,0]*CoordSys['a3'] - (   p['AxRedTSS'][0,0,J]* QT[DOF_TSS1]  + p['AxRedTSS'][0,1,J]* QT[DOF_TSS2]   )*CoordSys['a2']
       dat['PLinVelET'][J,DOF_TFA2,0,:] = p['TwrFASF'][1,J,0]*CoordSys['a1'] - (   p['AxRedTFA'][1,1,J]* QT[DOF_TFA2]  + p['AxRedTFA'][0,1,J]* QT[DOF_TFA1]   )*CoordSys['a2']
       dat['PLinVelET'][J,DOF_TSS2,0,:] = p['TwrSSSF'][1,J,0]*CoordSys['a3'] - (   p['AxRedTSS'][1,1,J]* QT[DOF_TSS2]  + p['AxRedTSS'][0,1,J]* QT[DOF_TSS1]   )*CoordSys['a2']
       TmpVec1 = np.cross( dat['AngVelEX'], dat['PLinVelET'][J,DOF_TFA1,0,:] )
       TmpVec2 = np.cross( dat['AngVelEX'], dat['PLinVelET'][J,DOF_TSS1,0,:] )
       TmpVec3 = np.cross( dat['AngVelEX'], dat['PLinVelET'][J,DOF_TFA2,0,:] )
       TmpVec4 = np.cross( dat['AngVelEX'], dat['PLinVelET'][J,DOF_TSS2,0,:] )
       dat['PLinVelET'][J,DOF_TFA1,1,:] = TmpVec1 - (   p['AxRedTFA'][0,0,J]*QDT[DOF_TFA1]  + p['AxRedTFA'][0,1,J]*QDT[DOF_TFA2]   )*CoordSys['a2']
       dat['PLinVelET'][J,DOF_TSS1,1,:] = TmpVec2 - (   p['AxRedTSS'][0,0,J]*QDT[DOF_TSS1]  + p['AxRedTSS'][0,1,J]*QDT[DOF_TSS2]   )*CoordSys['a2']
       dat['PLinVelET'][J,DOF_TFA2,1,:] = TmpVec3 - (   p['AxRedTFA'][1,1,J]*QDT[DOF_TFA2]  + p['AxRedTFA'][0,1,J]*QDT[DOF_TFA1]   )*CoordSys['a2']
       dat['PLinVelET'][J,DOF_TSS2,1,:] = TmpVec4 - (   p['AxRedTSS'][1,1,J]*QDT[DOF_TSS2]  + p['AxRedTSS'][0,1,J]*QDT[DOF_TSS1]   )*CoordSys['a2']
       LinVelXT               = QDT[DOF_TFA1]*dat['PLinVelET'][J,DOF_TFA1,0,:] \
                              + QDT[DOF_TSS1]*dat['PLinVelET'][J,DOF_TSS1,0,:] \
                              + QDT[DOF_TFA2]*dat['PLinVelET'][J,DOF_TFA2,0,:] \
                              + QDT[DOF_TSS2]*dat['PLinVelET'][J,DOF_TSS2,0,:]
       dat['LinAccETt'][J,:] =  QDT[DOF_TFA1]*dat['PLinVelET'][J,DOF_TFA1,1,:] \
                              + QDT[DOF_TSS1]*dat['PLinVelET'][J,DOF_TSS1,1,:] \
                              + QDT[DOF_TFA2]*dat['PLinVelET'][J,DOF_TFA2,1,:] \
                              + QDT[DOF_TSS2]*dat['PLinVelET'][J,DOF_TSS2,1,:]
       dat['LinVelET'][J,:] = LinVelXT + dat['LinVelEZ']
       for I in range(NPX): # Loop through all DOFs associated with the angular motion of the platform (body X)
          TmpVec0   = np.cross( dat['PAngVelEX'][PX[I],0,:], dat['rZT'][J,:]          )
          TmpVec1   = np.cross( dat['PAngVelEX'][PX[I],0,:], EwXXrZT      + LinVelXT )
          dat['PLinVelET'][J,PX[I],0,:] = dat['PLinVelET'][J,PX[I],0,:] + TmpVec0
          dat['PLinVelET'][J,PX[I],1,:] = dat['PLinVelET'][J,PX[I],1,:] + TmpVec1
          dat['LinVelET'][ J,        :] += QDT[PX[I]]*dat['PLinVelET'][J,PX[I],0,:]
          dat['LinAccETt'][J,        :] += QDT[PX[I]]*dat['PLinVelET'][J,PX[I],1,:]
    # ---
    if IEC is None:
        IEC = dict()
    IEC['v_F'] = EDVec2IEC(dat['LinVelEZ'])
    IEC['v_N'] = EDVec2IEC(dat['LinVelEO'])
    IEC['ud_N'] = EDVec2IEC(LinVelXO) # Elastic velocity of tower top
    IEC['v_Gn'] = EDVec2IEC(dat['LinVelEU']) # velocity of nacelle COG
    IEC['v_Ts'] = EDVec2IEC(dat['LinVelET']) # velocity of tower nodes

    return dat, IEC


def ED_CalcOutputs(x, p, noAxRed=False):
    """ 
    INPUTS: 
      - x: states (for now structural only)
      - p: parameters, e.g. as returned by p = ED_Parameters(fstSim)
    OUTPUTS:
      - CS: dictionary of Coordinate system info
      - dat: dictionary of structural variables in "ElastoDyn internal" coordinate system
      - IEC: dictionary of structural variables in OpenFAST/IEC coordinate system

    Example:
        p = ED_Parameters(fstSim)
        x['qDict']   = {'Sg': 10.0, 'Sw':20.0, 'Hv': 5.0, 'R':0.0, 'P':0.3, 'Y':0, 'TFA1':1.0, 'TSS1':10.0, 'Yaw':np.pi/8}
        x['qdDict']  = {'Sg':  1.0, 'Sw': 2.0, 'Hv': 3.0, 'R':0.1, 'P':0.3, 'Y':0, 'TFA1':2.0, 'TSS1':4.0,  'Yaw':0.0}
        
    """
    if noAxRed:
        p['AxRedTFA']*=0 #
        p['AxRedTSS']*=0 #
    qDict  = x['qDict']
    qdDict = x['qdDict']
    CS = ED_CoordSys(qDict=qDict, TwrFASF=p['TwrFASF'], TwrSSSF=p['TwrSSSF'])
    dat, IEC = ED_Positions(qDict=qDict, CoordSys=CS, p=p)
    dat, IEC = ED_AngPosVelPAcc(qDict=qDict, qdDict=qdDict, CoordSys=CS, p=p, dat=dat, IEC=IEC)
    dat, IEC = ED_LinVelPAcc   (qDict=qDict, qdDict=qdDict, CoordSys=CS, p=p, dat=dat, IEC=IEC)
    return CS, dat, IEC

if __name__ == '__main__':
#     EDfilename='../yams/_Jens/FEMBeam_NewFASTCoeffs/data/NREL5MW_ED_Onshore.dat'
#     EDfilename='../../data/NREL5MW/onshore/NREL5MW_ED_Onshore.dat'
#     bladeParameters(EDfilename)

    EDfilename  = '../../data/NREL5MW/onshore/NREL5MW_ED_Onshore.dat'
    FSTfilename = '../../data/NREL5MW/Main_Onshore.fst'

    fst = FASTInputFile(FSTfilename)

    p,pbld = rotorParameters(EDfilename)

    ptwr=towerParameters(EDfilename, gravity=fst['Gravity'], RotMass=p['RotMass'])

    # Calculate the turbine mass:
    ptwr['TurbMass']  = ptwr['TwrTpMass'] + ptwr['TwrMass'];


