import numpy as np
import os
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


def SHP(Fract, FlexL, ModShpAry, Deriv):
    """ SHP calculates the Derive-derivative of the shape function ModShpAry at Fract.
    #  NOTE: This function only works for Deriv = 0, 1, or 2. """
    if ( Deriv < 0 or Deriv > 2 ):
        raise Exception('Function SHP input Deriv={} is invalid. Deriv must be 0, 1, or 2.'.format(Deriv))
    elif ( Fract < 0.0 or Fract > 1.0 ):
        raise Exception('Function SHP input Fract={} does not meet the condition 0<=Fract<=1.'.format(Fract))

    Swtch        = np.zeros((3, 1)); # Initialize Swtch(:) to 0
    Swtch[Deriv] = 1;
    shp          = 0.0;

    for i in np.arange(len(ModShpAry)):
        I = i + 1
        J = I + 1;
        CoefTmp = Swtch[0] + Swtch[1]*J + Swtch[2]*I*J;

        if ( (J == 2) and (Deriv == 2) ):
            shp =       ModShpAry[i]*CoefTmp                         /( FlexL**Deriv );
        else:
            shp = shp + ModShpAry[i]*CoefTmp*( Fract**( J - Deriv ) )/( FlexL**Deriv );
    return shp


def rotorParameters(EDfilename):
    ED = FASTInputFile(EDfilename)
    pbld=[bladeParameters(EDfilename, ibld+1) for ibld in  range(ED['NumBl'])]

    p['RotMass'] = sum([pbld[k]['BldMass'] for k in range(ED['NumBl'])])
    p['RotIner'] = sum([(pbld[k]['SecondMom'] + pbld[k]['BldMass']*ED['HubRad']*(2.0*pbld[k]['BldCG'] + ED['HubRad']))*(pbld[k]['CosPreC'](K)**2) for k in range(ED['NumBl'])])

def bladeParameters(EDfilename, ibld=1, RotSpeed=1):
    """
    Compute blade parameters in a way similar to OpenFAST
    See Routine Coeff from ElastoDyn.f90
    RotSpeed: used for rotational stiffening. Use 1 for unit contribution (proportioanl to omega**2) [rad/s]
    """
    # --- Read inputs
    ED      = FASTInputFile(EDfilename)
    EDbld   = os.path.join(os.path.dirname(EDfilename), ED['BldFile({})'.format(ibld)].replace('"',''))
    bld     = FASTInputFile(EDbld)
    bldProp = bld.toDataFrame()

    # --- 
    p=dict()
    BD4Blades = False
    p['BldNodes'] = ED['BldNodes']
    p['BldFlexL'] = ED['TipRad']- ED['HubRad'] # Length of the flexible portion of the blade.
    p['BldMass'] =0 
    n=ED['BldNodes']
    if not BD4Blades:
        p['DRNodes'] = np.ones(ED['BldNodes'])*p['BldFlexL']/ED['BldNodes']
        bld_fract    = np.arange(1./ED['BldNodes']/2., 1, 1./ED['BldNodes'])
        p['RNodes'] = bld_fract*p['BldFlexL']
    # --- Interpolate the blade properties to this discretization:
    p['RNodesNorm'] = p['RNodes']/p['BldFlexL'];  # Normalized radius to analysis nodes relative to hub ( -1 < RNodesNorm(:) < 1 )
    p['BlFract']= bldProp['BlFract_[-]'].values
    StrcTwst    = bldProp['StrcTwst_[deg]'].values
    p['ThetaS']  = np.interp(p['RNodesNorm'], p['BlFract'], bldProp['StrcTwst_[deg]']) 
    p['MassB']   = np.interp(p['RNodesNorm'], p['BlFract'], bldProp['BMassDen_[kg/m]'])                     ;
    p['StiffBF'] = np.interp(p['RNodesNorm'], p['BlFract'], bldProp['FlpStff_[Nm^2]'])
    p['StiffBE'] = np.interp(p['RNodesNorm'], p['BlFract'], bldProp['EdgStff_[Nm^2]'])
    p['ThetaS']  = np.concatenate( ([StrcTwst[0]], p['ThetaS'] , [StrcTwst[-1]]) )
    # Set the blade damping and stiffness tuner
    p['BldFDamp'] = [bld['BldFlDmp(1)'], bld['BldFlDmp(2)'] ]
    p['BldEDamp'] = [bld['BldEdDmp(1)']]
    p['FStTunr']  = [bld['FlStTunr(1)'], bld['FlStTunr(2)'] ]
    # Set the mode shape coefficients 
    p['BldFl1Sh'] = [bld[c] for c in ['BldFl1Sh(2)', 'BldFl1Sh(3)', 'BldFl1Sh(4)', 'BldFl1Sh(5)', 'BldFl1Sh(6)']]
    p['BldFl2Sh'] = [bld[c] for c in ['BldFl2Sh(2)', 'BldFl2Sh(3)', 'BldFl2Sh(4)', 'BldFl2Sh(5)', 'BldFl2Sh(6)']]
    p['BldEdgSh'] = [bld[c] for c in ['BldEdgSh(2)', 'BldEdgSh(3)', 'BldEdgSh(4)', 'BldEdgSh(5)', 'BldEdgSh(6)']]
    p['CThetaS'] = np.cos(p['ThetaS']*np.pi/180);
    p['SThetaS'] = np.sin(p['ThetaS']*np.pi/180);

    # --- Inertial properties
    # Initialize BldMass(), FirstMom(), and SecondMom() using TipMass() effects
    p['TipMass']   = ED['TipMass({:d})'.format(ibld)]                            
    p['BldMass']   = p['TipMass']                            
    p['FirstMom']  = p['TipMass']*p['BldFlexL']             
    p['SecondMom'] = p['TipMass']*p['BldFlexL']*p['BldFlexL']
    p['BElmntMass'] = np.zeros(n) # Mass of blade element FMomAbvNd J
    p['FMomAbvNd']  = np.zeros(n)
    for J in  np.arange(ED['BldNodes']-1,-1,-1): # Loop through the blade nodes / elements in reverse
        # Calculate the mass of the current element
        p['BElmntMass'][J] = p['MassB'][J]*p['DRNodes'][J] # Mass of blade elementp['FMomAbvNd J
        # Integrate to find some blade properties which will be output in .fsm
        p['BldMass']   += p['BElmntMass'][J];
        p['FirstMom']  += p['BElmntMass'][J]*p['RNodes'][J];
        p['SecondMom'] += p['BElmntMass'][J]*p['RNodes'][J]*p['RNodes'][J];
        # Integrate to find FMomAbvNd:
        p['FMomAbvNd'][J] = (0.5*p['BElmntMass'][J] )*(ED['HubRad'] + p['RNodes'][J] + 0.5*p['DRNodes'][J])
        if J == n-1: # Outermost blade element
           # Add the TipMass() effects:
           p['FMomAbvNd'][J] += p['TipMass'] * ED['TipRad'];
        else:  
           # Add to p['FMomAbvNd(K,J) the effects from the (not yet used) portion of element J+1
           p['FMomAbvNd'][J] += p['FMomAbvNd'][J+1] + (0.5*p['BElmntMass'][J+1])*( ED['HubRad'] + p['RNodes'][J+1] - 0.5*p['DRNodes'][J+1] );
    if not BD4Blades:
        # Calculate BldCG() using FirstMom() and BldMass(); and calculate RotMass and RotIner:
        p['BldCG']= p['FirstMom']/p['BldMass'];

    p['MBF']       = np.zeros((2, 2))                 ;
    p['MBE']       = np.zeros((1, 1))                 ;
    p['KBFCent']   = np.zeros((2, 2))                 ;
    p['KBECent']   = np.zeros((1, 1))                ;
    p['KBF']       = np.zeros((2, 2))                ;
    p['KBE']       = np.zeros((1, 1))                ;
    p['TwistedSF'] = np.zeros((2, 3, ED['BldNodes']+2, 3)); # x/y, BF1/BF2/BE, node, deriv
    p['AxRedBld']  = np.zeros((3, 3, ED['BldNodes']+2)); #
    p['ShapeF1']   = np.zeros(n);
    p['ShapeF2']   = np.zeros(n);
    p['ShapeE1']   = np.zeros(n);
    p['dShapeF1']   = np.zeros(n);
    p['dShapeF2']   = np.zeros(n);
    p['dShapeE1']   = np.zeros(n);
    p['ddShapeF1']   = np.zeros(n);
    p['ddShapeF2']   = np.zeros(n);
    p['ddShapeE1']   = np.zeros(n);
    # Initialize the generalized blade masses using tip mass effects:
    p['MBF'][0,0] = p['TipMass'];
    p['MBF'][1,1] = p['TipMass'];
    p['MBE'][0,0] = p['TipMass'];
    for J  in np.arange(n):    # Loop through the blade nodes / elements
        # Integrate to find the generalized mass of the blade (including tip mass effects).
        #   Ignore the cross-correlation terms of MBF (i.e. MBF(i,j) where i ~= j) since
        #   these terms will never be used.
        p['ShapeF1'][J] = SHP( p['RNodesNorm'][J], p['BldFlexL'], p['BldFl1Sh'], 0);
        p['ShapeF2'][J] = SHP( p['RNodesNorm'][J], p['BldFlexL'], p['BldFl2Sh'], 0);
        p['ShapeE1'][J] = SHP( p['RNodesNorm'][J], p['BldFlexL'], p['BldEdgSh'], 0);
        p['MBF'][0,0]  += p['BElmntMass'][J]*p['ShapeF1'][J]*p['ShapeF1'][J];
        p['MBF'][1,1]  += p['BElmntMass'][J]*p['ShapeF2'][J]*p['ShapeF2'][J];
        p['MBE'][0,0]  += p['BElmntMass'][J]*p['ShapeE1'][J]*p['ShapeE1'][J];
        # Integrate to find the generalized stiffness of the blade (not including centrifugal effects).
        ElmntStff      = p['StiffBF'][J]*p['DRNodes'][J]  # Flapwise stiffness of blade element J
        p['ddShapeF1'][J] = SHP( p['RNodesNorm'][J], p['BldFlexL'], p['BldFl1Sh'], 2);
        p['ddShapeF2'][J] = SHP( p['RNodesNorm'][J], p['BldFlexL'], p['BldFl2Sh'], 2);
        p['KBF'][0,0] += ElmntStff*p['ddShapeF1'][J]*p['ddShapeF1'][J]
        p['KBF'][0,1] += ElmntStff*p['ddShapeF1'][J]*p['ddShapeF2'][J]
        p['KBF'][1,0] += ElmntStff*p['ddShapeF2'][J]*p['ddShapeF1'][J]
        p['KBF'][1,1] += ElmntStff*p['ddShapeF2'][J]*p['ddShapeF2'][J]
        ElmntStff      = p['StiffBE'][J]*p['DRNodes'][J] # Edgewise stiffness of blade element J
        p['ddShapeE1'][J] = SHP( p['RNodesNorm'][J], p['BldFlexL'], p['BldEdgSh'], 2);
        p['KBE'][0,0] += ElmntStff*p['ddShapeE1'][J]*p['ddShapeE1'][J];
        # Integrate to find the centrifugal-term of the generalized flapwise and edgewise
        #   stiffness of the blades.  Ignore the cross-correlation terms of KBFCent (i.e.
        #   KBFCent(i,j) where i ~= j) since these terms will never be used.
        ElmntStff      = p['FMomAbvNd'][J]*p['DRNodes'][J]*RotSpeed**2 # Centrifugal stiffness of blade element J
        p['dShapeF1'][J]  = SHP( p['RNodesNorm'][J], p['BldFlexL'], p['BldFl1Sh'], 1)
        p['dShapeF2'][J]  = SHP( p['RNodesNorm'][J], p['BldFlexL'], p['BldFl2Sh'], 1)
        p['dShapeE1'][J]  = SHP( p['RNodesNorm'][J], p['BldFlexL'], p['BldEdgSh'], 1)
        p['KBFCent'][0,0] += ElmntStff*p['dShapeF1'][J]*p['dShapeF1'][J]
        p['KBFCent'][1,1] += ElmntStff*p['dShapeF2'][J]*p['dShapeF2'][J]
        p['KBECent'][0,0] += ElmntStff*p['dShapeE1'][J]*p['dShapeE1'][J]
        # Calculate the 2nd derivatives of the twisted shape functions:
        Shape             = SHP( p['RNodesNorm'][J], p['BldFlexL'], p['BldFl1Sh'], 2);
        p['TwistedSF'][0,0,J+1,2] =  Shape*p['CThetaS'][J+1] # 2nd deriv. of Phi1(J) for blade K
        p['TwistedSF'][1,0,J+1,2] = -Shape*p['SThetaS'][J+1] # 2nd deriv. of Psi1(J) for blade K
        Shape  = SHP( p['RNodesNorm'][J], p['BldFlexL'], p['BldFl2Sh'], 2);
        p['TwistedSF'][0,1,J+1,2] =  Shape*p['CThetaS'][J+1] # 2nd deriv. of Phi2(J) for blade K
        p['TwistedSF'][1,1,J+1,2] = -Shape*p['SThetaS'][J+1] # 2nd deriv. of Psi2(J) for blade K
        Shape  = SHP( p['RNodesNorm'][J], p['BldFlexL'], p['BldEdgSh'], 2);
        p['TwistedSF'][0,2,J+1,2] =  Shape*p['SThetaS'][J+1] # 2nd deriv. of Phi3(J) for blade K
        p['TwistedSF'][1,2,J+1,2] =  Shape*p['CThetaS'][J+1] # 2nd deriv. of Psi3(J) for blade K
        # Integrate to find the 1st derivatives of the twisted shape functions:
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
        # Integrate to find the blade axial reduction shape functions:
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
        # Store the TwstdSF and AxRdBld terms of the current element (these will be used for the next element)
        TwstdSFOld = TwstdSF;
        AxRdBldOld = AxRdBld;
    # End loop on J

    if BD4Blades:
        # the 1st and zeroeth derivatives of the twisted shape functions at the blade root:
        p['TwistedSF'][:,:,:,1] = 0.0;
        p['TwistedSF'][:,:,:,0] = 0.0;
        p['AxRedBld'] [:,:,:  ] = 0.0;
    else:
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
        # Calculate the 2nd derivatives of the twisted shape functions at the blade root:
        Shape  = SHP( 0.0, p['BldFlexL'], p['BldFl1Sh'], 2)
        p['TwistedSF'][0,0,0,2] =  Shape*p['CThetaS'][0] # 2nd deriv. of Phi1(0) for blade K
        p['TwistedSF'][1,0,0,2] = -Shape*p['SThetaS'][0] # 2nd deriv. of Psi1(0) for blade K
        Shape  = SHP( 0.0, p['BldFlexL'], p['BldFl2Sh'], 2)
        p['TwistedSF'][0,1,0,2] =  Shape*p['CThetaS'][0] # 2nd deriv. of Phi2(0) for blade K
        p['TwistedSF'][1,1,0,2] = -Shape*p['SThetaS'][0] # 2nd deriv. of Psi2(0) for blade K
        Shape  = SHP( 0.0, p['BldFlexL'], p['BldEdgSh'], 2)
        p['TwistedSF'][0,2,0,2] =  Shape*p['SThetaS'][0] # 2nd deriv. of Phi3(0) for blade K
        p['TwistedSF'][1,2,0,2] =  Shape*p['CThetaS'][0] # 2nd deriv. of Psi3(0) for blade K
        # Calculate the 2nd derivatives of the twisted shape functions at the tip:
        Shape  = SHP( 1.0, p['BldFlexL'], p['BldFl1Sh'], 2)
        p['TwistedSF'][0,0,-1,2] =  Shape*p['CThetaS'][-1] # 2nd deriv. of Phi1(p['TipNode) for blade K
        p['TwistedSF'][1,0,-1,2] = -Shape*p['SThetaS'][-1] # 2nd deriv. of Psi1(p['TipNode) for blade K
        Shape  = SHP( 1.0, p['BldFlexL'], p['BldFl2Sh'], 2)
        p['TwistedSF'][0,1,-1,2] =  Shape*p['CThetaS'][-1] # 2nd deriv. of Phi2(p['TipNode) for blade K
        p['TwistedSF'][1,1,-1,2] = -Shape*p['SThetaS'][-1] # 2nd deriv. of Psi2(p['TipNode) for blade K
        Shape  = SHP( 1.0, p['BldFlexL'], p['BldEdgSh'], 2)
        p['TwistedSF'][0,2,-1,2] =  Shape*p['SThetaS'][-1] # 2nd deriv. of Phi3(p['TipNode) for blade K
        p['TwistedSF'][1,2,-1,2] =  Shape*p['CThetaS'][-1] # 2nd deriv. of Psi3(p['TipNode) for blade K
        # Integrate to find the 1st and zeroeth derivatives of the twisted shape functions at the tip:
        for I in [0,1]:   # Loop through Phi and Psi
            for L in [0,1,2]:  # Loop through all blade DOFs
                p['TwistedSF'][I,L,-1,1] = p['TwistedSF'][I,L,-2,1] + TwstdSFOld[I,L,1]
                p['TwistedSF'][I,L,-1,0] = p['TwistedSF'][I,L,-2,0] + TwstdSFOld[I,L,0]
        # the 1st and zeroeth derivatives of the twisted shape functions at the blade root:
        p['TwistedSF'][:,:,0,1] = 0.0;
        p['TwistedSF'][:,:,0,0] = 0.0;
        p['AxRedBld'] [:,:,0  ] = 0.0;
        # Integrate to find the blade axial reduction shape functions at the tip:
        for I in [0,1,2]:     # Loop through all blade DOFs
            for L in [0,1,2]:  # Loop through all blade DOFs
                p['AxRedBld'][I,L,-1] = p['AxRedBld'][I,L,-2] + AxRdBldOld[I,L]

    return p

def towerParameters():
    pass
    #p.TwrFlexL  = p.TowerHt   - p.TowerBsHt;                                         % Height / length of the flexible portion of the tower.
 # 
# # --------------------------------------------------------------------------------}
# ## --- Tower 
# # --------------------------------------------------------------------------------{
# # setdiff(fieldnames(p), fields)
# # fields= fieldnames(p);
# # Calculate the tower-top mass:
# 
# # p['TwrTpMass = p['RotMass + p['RFrlMass + p['BoomMass + p['TFinMass + p['NacMass + p['YawBrMass;
# p['TwrTpMass = p['RotMass + p['NacMass + p['YawBrMass;
# 
# 
# for J = p['TwrNodes:-1:1 # Loop through the tower nodes / elements in reverse
# 
# 
#     # Calculate the mass of the current element
# 
#     p['TElmntMass(J)    = p['MassT(J)*p['DHNodes(J);     # Mass of tower element J
# 
# 
#     # Integrate to find the tower mass which will be output in .fsm
# 
#     p['TwrMass      = p['TwrMass + p['TElmntMass(J);
# 
# 
#     # Integrate to find TMssAbvNd:
# 
#     TMssAbvNd   (J) = 0.5*p['TElmntMass(J);
# 
#     if ( J == p['TwrNodes )   # Uppermost tower element
#         # Add the TwrTpMass effects:
# 
# #         TMssAbvNd(J) = TMssAbvNd(J) + p['TwrTpMass;
#         TMssAbvNd(J) = TMssAbvNd(J);
#     else                       # All other tower elements
#         # Add to TMssAbvNd(J) the effects from the (not yet used) portion of element J+1
# 
#         TMssAbvNd(J) = 0.5*p['TElmntMass(J+1) + TMssAbvNd(J) + TMssAbvNd(J+1);
#     end
# 
# 
# end # J - Tower nodes / elements in reverse
# 
# 
# 
# # Initialize the generalized tower masses using tower-top mass effects:
# p['MTFA= zeros(2, 2);
# p['MTSS= zeros(2, 2);
# for I = 1:2  # Loop through all tower modes in a single direction
#     p['MTFA(I,I) = p['TwrTpMass;
#     p['MTSS(I,I) = p['TwrTpMass;
# end       # I - All tower modes in a single direction
# 
# # set values for tower base (note that we haven't corrctly defined the values for (:,0,2) in the arrays below):
# p['TwrFASF(   1:2,1,1:2) = 0.0;
# p['TwrSSSF(   1:2,1,1:2) = 0.0;
# p['AxRedTFA(1:2,1:2,1)   = 0.0;
# p['AxRedTSS(1:2,1:2,1)   = 0.0;
# 
# p['KTFA= zeros(2, 2);
# p['KTSS= zeros(2, 2);
# p['KTFAGrav= zeros(2, 2);
# p['KTSSGrav= zeros(2, 2);
# 
# for J = 1:p['TwrNodes    # Loop through the tower nodes / elements
# 
# 
#     # Calculate the tower shape functions (all derivatives):
# 
#     p['TwrFASF(1,J+1,3) = SHP( p['HNodesNorm(J), p['TwrFlexL, p['TwFAM1Sh(:), 2);
#     p['TwrFASF(2,J+1,3) = SHP( p['HNodesNorm(J), p['TwrFlexL, p['TwFAM2Sh(:), 2);
#     p['TwrFASF(1,J+1,2) = SHP( p['HNodesNorm(J), p['TwrFlexL, p['TwFAM1Sh(:), 1);
#     p['TwrFASF(2,J+1,2) = SHP( p['HNodesNorm(J), p['TwrFlexL, p['TwFAM2Sh(:), 1);
#     p['TwrFASF(1,J+1,1) = SHP( p['HNodesNorm(J), p['TwrFlexL, p['TwFAM1Sh(:), 0);
#     p['TwrFASF(2,J+1,1) = SHP( p['HNodesNorm(J), p['TwrFlexL, p['TwFAM2Sh(:), 0);
# 
#     p['TwrSSSF(1,J+1,3) = SHP( p['HNodesNorm(J), p['TwrFlexL, p['TwSSM1Sh(:), 2);
#     p['TwrSSSF(2,J+1,3) = SHP( p['HNodesNorm(J), p['TwrFlexL, p['TwSSM2Sh(:), 2);
#     p['TwrSSSF(1,J+1,2) = SHP( p['HNodesNorm(J), p['TwrFlexL, p['TwSSM1Sh(:), 1);
#     p['TwrSSSF(2,J+1,2) = SHP( p['HNodesNorm(J), p['TwrFlexL, p['TwSSM2Sh(:), 1);
#     p['TwrSSSF(1,J+1,1) = SHP( p['HNodesNorm(J), p['TwrFlexL, p['TwSSM1Sh(:), 0);
#     p['TwrSSSF(2,J+1,1) = SHP( p['HNodesNorm(J), p['TwrFlexL, p['TwSSM2Sh(:), 0);
# 
# 
#     # Integrate to find the generalized mass of the tower (including tower-top mass effects).
#     #   Ignore the cross-correlation terms of p['MTFA (i.e. p['MTFA(i,j) where i ~= j) and p['MTSS
#     #   since these terms will never be used.
# 
# 
#     for I = 1:2     # Loop through all tower DOFs in one direction
#         p['MTFA  (I,I) = p['MTFA  (I,I) + p['TElmntMass(J)*p['TwrFASF(I,J+1,1)^2;
#         p['MTSS  (I,I) = p['MTSS  (I,I) + p['TElmntMass(J)*p['TwrSSSF(I,J+1,1)^2;
#     end          # I - through all tower DOFs in one direction
# 
# 
#     # Integrate to find the generalized stiffness of the tower (not including gravitational
#     #    effects).
# 
#     ElStffFA       = p['StiffTFA(J)*p['DHNodes(J);                        # Fore-aft stiffness of tower element J
#     ElStffSS       = p['StiffTSS(J)*p['DHNodes(J);                        # Side-to-side stiffness of tower element J
# 
#     for I = 1:2     # Loop through all tower DOFs in one direction
#         for L = 1:2  # Loop through all tower DOFs in one direction
#             p['KTFA (I,L) = p['KTFA    (I,L) + ElStffFA *p['TwrFASF(I,J+1,3)*p['TwrFASF(L,J,3);
#             p['KTSS (I,L) = p['KTSS    (I,L) + ElStffSS *p['TwrSSSF(I,J+1,3)*p['TwrSSSF(L,J,3);
#         end       # L - All tower DOFs in one direction
#     end          # I - through all tower DOFs in one direction
# 
# 
#     # Integrate to find the gravitational-term of the generalized stiffness of the tower.
#     #   Ignore the cross-correlation terms of p['KTFAGrav (i.e. p['KTFAGrav(i,j) where i ~= j)
#     #   and p['KTSSGrav since these terms will never be used.
# 
# #     ElmntStff      = -TMssAbvNd(J)*p['DHNodes(J)*p['Gravity;              # Gravitational stiffness of tower element J
#     ElmntStff      = -TMssAbvNd(J)*p['DHNodes(J)*1;              # stiffness of tower element J due to unity acceleration
# 
#     for I = 1:2     # Loop through all tower DOFs in one direction
#         p['KTFAGrav(I,I) = p['KTFAGrav(I,I) + ElmntStff*p['TwrFASF(I,J+1,2)^2;
#         p['KTSSGrav(I,I) = p['KTSSGrav(I,I) + ElmntStff*p['TwrSSSF(I,J+1,2)^2;
#     end
# 
# 
#     # Integrate to find the tower axial reduction shape functions:
#     AxRdTFA= zeros(2, 2);
#     AxRdTSS= zeros(2, 2);
#     for I = 1:2     # Loop through all tower DOFs in one direction
#         for L = 1:2  # Loop through all tower DOFs in one direction
#             AxRdTFA (I,L) = 0.5*p['DHNodes(J)*p['TwrFASF(I,J+1,2)*p['TwrFASF(L,J,2);
#             AxRdTSS (I,L) = 0.5*p['DHNodes(J)*p['TwrSSSF(I,J+1,2)*p['TwrSSSF(L,J,2);
# 
#             p['AxRedTFA(I,L,J+1) = AxRdTFA(I,L);
#             p['AxRedTSS(I,L,J+1) = AxRdTSS(I,L);
#         end       # L - All tower DOFs in one direction
#     end
# 
#     if ( J ~= 1 )    # All but the lowermost tower element
#         # Add the effects from the (not yet used) portion of element J-1
# 
#         for I = 1:2     # Loop through all tower DOFs in one direction
#             for L = 1:2  # Loop through all tower DOFs in one direction
#                 p['AxRedTFA(I,L,J+1) = p['AxRedTFA(I,L,J+1) + p['AxRedTFA(I,L,J)+ AxRdTFAOld(I,L);
#                 p['AxRedTSS(I,L,J+1) = p['AxRedTSS(I,L,J+1) + p['AxRedTSS(I,L,J)+ AxRdTSSOld(I,L);
#             end       # L - All tower DOFs in one direction
#         end
#     end
# 
# 
#     # Store the AxRdTFA and AxRdTSS terms of the current element (these will be used for the next element)
# 
#     AxRdTFAOld = AxRdTFA;
#     AxRdTSSOld = AxRdTSS;
# 
# 
# end # J - Tower nodes / elements
# 
# ElmntStff      = -1*p['TwrFlexL*1;              # due to unitiy mass under unity acceleration stiffness of tower element J due to unity acceleration
# for I= 1:2
#     p['KTFAGravTT(I,I) = ElmntStff*p['TwrFASF(I,p['TwrNodes+1,2)^2;
#     p['KTSSGravTT(I,I) = ElmntStff*p['TwrSSSF(I,p['TwrNodes+1,2)^2;
# end
# 
# # Apply the modal stiffness tuners of the tower to KTFA() and KTSS():
# 
# for I = 1:2     # Loop through all tower DOFs in one direction
#     for L = 1:2  # Loop through all tower DOFs in one direction
#         p['KTFA(I,L) = sqrt( p['FAStTunr(I)*p['FAStTunr(L) )*p['KTFA(I,L);
# 
#         p['KTSS(I,L) = sqrt( p['SSStTunr(I)*p['SSStTunr(L) )*p['KTSS(I,L);
#     end       # L - All tower DOFs in one direction
# end          # I - through all tower DOFs in one direction
# 
# 
# # Calculate the tower natural frequencies:
# 
# for I = 1:2     # Loop through all tower DOFs in one direction
#     allKTFAGrav= (p['KTFAGrav(I,I) + p['TwrTpMass*p['KTFAGravTT(I,I))*p['Gravity;
#     allKTSSGrav= (p['KTSSGrav(I,I) + p['TwrTpMass*p['KTSSGravTT(I,I))*p['Gravity;
#     p['FreqTFA(I,1) = (1/2/pi)*sqrt(   p['KTFA(I,I)                  /( p['MTFA(I,I) - p['TwrTpMass ) );  # Natural tower I-fore-aft frequency w/o gravitational destiffening nor tower-top mass effects
#     p['FreqTFA(I,2) = (1/2/pi)*sqrt( ( p['KTFA(I,I) + allKTFAGrav )/  p['MTFA(I,I)               );  # Natural tower I-fore-aft frequency w/  gravitational destiffening and tower-top mass effects
#     p['FreqTSS(I,1) = (1/2/pi)*sqrt(   p['KTSS(I,I)                  /( p['MTSS(I,I) - p['TwrTpMass ) );  # Natural tower I-side-to-side frequency w/o gravitational destiffening nor tower-top mass effects
#     p['FreqTSS(I,2) = (1/2/pi)*sqrt( ( p['KTSS(I,I) + allKTSSGrav )/  p['MTSS(I,I)               );  # Natural tower I-side-to-side frequency w/  gravitational destiffening and tower-top mass effects
# end          # I - All tower DOFs in one direction
# 
# 
# # Calculate the generalized damping of the tower:
# 
# for I = 1:2     # Loop through all tower DOFs in one direction
#     for L = 1:2  # Loop through all tower DOFs in one direction
#         p['CTFA(I,L) = ( 0.01*p['TwrFADmp(L) )*p['KTFA(I,L)/( pi*p['FreqTFA(L,1) );
# 
#         p['CTSS(I,L) = ( 0.01*p['TwrSSDmp(L) )*p['KTSS(I,L)/( pi*p['FreqTSS(L,1) );
#     end       # L - All tower DOFs in one direction
# end          # I - All tower DOFs in one direction
# 
# 
# # Calculate the tower shape functions (all derivatives) at the tower-top:
# 
# p['TwrFASF(1,p['TTopNode+1,3) = SHP( 1.0, p['TwrFlexL, p['TwFAM1Sh(:), 2);
# p['TwrFASF(2,p['TTopNode+1,3) = SHP( 1.0, p['TwrFlexL, p['TwFAM2Sh(:), 2);
# p['TwrFASF(1,p['TTopNode+1,2) = SHP( 1.0, p['TwrFlexL, p['TwFAM1Sh(:), 1);
# p['TwrFASF(2,p['TTopNode+1,2) = SHP( 1.0, p['TwrFlexL, p['TwFAM2Sh(:), 1);
# p['TwrFASF(1,p['TTopNode+1,1) = SHP( 1.0, p['TwrFlexL, p['TwFAM1Sh(:), 0);
# p['TwrFASF(2,p['TTopNode+1,1) = SHP( 1.0, p['TwrFlexL, p['TwFAM2Sh(:), 0);
# 
# p['TwrSSSF(1,p['TTopNode+1,3) = SHP( 1.0, p['TwrFlexL, p['TwSSM1Sh(:), 2);
# p['TwrSSSF(2,p['TTopNode+1,3) = SHP( 1.0, p['TwrFlexL, p['TwSSM2Sh(:), 2);
# p['TwrSSSF(1,p['TTopNode+1,2) = SHP( 1.0, p['TwrFlexL, p['TwSSM1Sh(:), 1);
# p['TwrSSSF(2,p['TTopNode+1,2) = SHP( 1.0, p['TwrFlexL, p['TwSSM2Sh(:), 1);
# p['TwrSSSF(1,p['TTopNode+1,1) = SHP( 1.0, p['TwrFlexL, p['TwSSM1Sh(:), 0);
# p['TwrSSSF(2,p['TTopNode+1,1) = SHP( 1.0, p['TwrFlexL, p['TwSSM2Sh(:), 0);
# 
# 
# # Integrate to find the tower axial reduction shape functions at the tower-top:
# 
# for I = 1:2     # Loop through all tower DOFs in one direction
#     for L = 1:2  # Loop through all tower DOFs in one direction
#         p['AxRedTFA(I,L,p['TTopNode+1) = p['AxRedTFA(I,L,p['TwrNodes+1)+ AxRdTFAOld(I,L);
#         p['AxRedTSS(I,L,p['TTopNode+1) = p['AxRedTSS(I,L,p['TwrNodes+1)+ AxRdTSSOld(I,L);
#     end       # L - All tower DOFs in one direction
# end
# 
# 
# # Calculate the turbine mass:
# 
# p['TurbMass  = p['TwrTpMass + p['TwrMass;
# 
# # setdiff(fieldnames(p), fields)
# 
#   pass


    for k,v in p.items():
        if hasattr(v, '__len__'):
            v = np.asarray(v)
            if len(v.shape)>=3:
                print('{:15s}:({})'.format(k,v.shape))
            elif len(v.shape)==2:
                print('{:15s}:\n {} ({})'.format(k,v, v.shape))
            else:
                n=len(v)
                print('{:15s}:{} ({})'.format(k,v,n))
        else:
            print('{:15s}:{}'.format(k,v))


if __name__ == '__main__':
    #EDfilename='../yams/_Jens/FEMBeam_NewFASTCoeffs/data/NREL5MW_ED_Onshore.dat'
    EDfilename='../../data/NREL5MW/Main_Onshore.fst'
    EDfilename='../../data/NREL5MW/onshore/NREL5MW_ED_Onshore.dat'
    bladeParameters(EDfilename)

