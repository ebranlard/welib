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


def rotorParameters(EDfilename, identicalBlades=True, pbld1=None):
    ED = FASTInputFile(EDfilename)

    if pbld1 is None:
        if identicalBlades:
            pbld = [bladeParameters(EDfilename, ibld+1)]*ED['NumBl']
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
    p['RotMass'] += ED['HubMass']
    p['RotIner'] += ED['HubIner']

    return p, pbld


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
    p['HubRad'] = ED['HubRad']
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
    p['s_span'] = np.concatenate(([0], p['RNodesNorm'], [1]))*p['BldFlexL'];  # Normalized radius to analysis nodes relative to hub ( -1 < RNodesNorm(:) < 1 )
    p['BlFract']= bldProp['BlFract_[-]'].values
    StrcTwst    = bldProp['StrcTwst_[deg]'].values
    p['ThetaS']  = np.interp(p['RNodesNorm'], p['BlFract'], bldProp['StrcTwst_[deg]']) 
    p['MassB']   = np.interp(p['RNodesNorm'], p['BlFract'], bldProp['BMassDen_[kg/m]'])                     ;
    p['StiffBF'] = np.interp(p['RNodesNorm'], p['BlFract'], bldProp['FlpStff_[Nm^2]'])
    p['StiffBE'] = np.interp(p['RNodesNorm'], p['BlFract'], bldProp['EdgStff_[Nm^2]'])
    p['m_full']    = np.interp(p['s_span'], p['BlFract'], bldProp['BMassDen_[kg/m]'])                     ;
    p['EI_F_full'] = np.interp(p['s_span'], p['BlFract'], bldProp['FlpStff_[Nm^2]'])
    p['EI_E_full'] = np.interp(p['s_span'], p['BlFract'], bldProp['EdgStff_[Nm^2]'])
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

    p['PreCone']= ED['PreCone({:d})'.format(ibld)]*np.pi/180

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


    nq = 3
    rh = p['HubRad'] # Hub Radius # TODO make this an option if from blade root or not

    # --- Shape functions
    nNodes=n+2
    p['U']=np.zeros((nq, 3, nNodes))
    p['V']=np.zeros((nq, 3, nNodes))
    p['K']=np.zeros((nq, 3, nNodes))
    for j in range(0,nq):
        p['U'][j][0,:] = p['TwistedSF'][0, j, :, 0]  # x
        p['U'][j][1,:] = p['TwistedSF'][1, j, :, 0]  # y
        p['V'][j][0,:] = p['TwistedSF'][0, j, :, 1]  # x
        p['V'][j][1,:] = p['TwistedSF'][1, j, :, 1]  # y
        p['K'][j][0,:] = p['TwistedSF'][0, j, :, 2]  # x
        p['K'][j][1,:] = p['TwistedSF'][1, j, :, 2]  # y

    from welib.yams.flexibility import GMBeam, GKBeam
    # TODO TODO TODO
    #MM, Gr, Ge, Oe, Oe6 = GMBeam(s_G0, s_span, m, PhiU, jxxG=jxxG, bUseIW=True, main_axis=main_axis, bAxialCorr=bAxialCorr, bOrth=False, rot_terms=True)
    #KK0 = GKBeam(s_span, EI, PhiK, bOrth=False)
    #if bStiffening:
    #    KKg     = GKBeamStiffnening(s_span, PhiV, gravity, m, Mtop, Omega, main_axis=main_axis)
    #    KKg_self= GKBeamStiffnening(s_span, PhiV, gravity, m, Mtop, Omega, main_axis=main_axis, bSelfWeight=True , bMtop=False, bRot=False)
    #    KKg_Mtop= GKBeamStiffnening(s_span, PhiV, gravity, m, Mtop, Omega, main_axis=main_axis, bSelfWeight=False, bMtop=True,  bRot=False)
    #    KKg_rot = GKBeamStiffnening(s_span, PhiV, gravity, m, Mtop, Omega, main_axis=main_axis, bSelfWeight=False, bMtop=False, bRot=True)
    s_G0 = np.zeros((3, len(p['s_span'])))
    s_G0[2,:] = p['s_span'] + rh # TODO hub radius option
    MM, Gr, Ge, Oe, Oe6 = GMBeam(s_G0, p['s_span'], p['m_full'], p['U'], rot_terms=True) 
    #, jxxG=jxxG, bUseIW=True, main_axis=main_axis, bAxialCorr=bAxialCorr, bOrth=False, rot_terms=True)

    # --- Rigid mass matrix terms
    p['J']    = np.zeros((3,3))
    p['mdCM'] = np.zeros((3,3))
    p['mdCM'][2,0]= sum(p['BElmntMass'][:]*(p['RNodes']+rh));
    #sid.md1_1_1_1= sum(squeeze(p.TwistedSF(1, 1, 1, 2:end-1, 1)).*p.BElmntMass(:, 1));
    #sid.md1_1_2_1= sum(squeeze(p.TwistedSF(1, 2, 1, 2:end-1, 1)).*p.BElmntMass(:, 1));
    #sid.md1_2_1_1= sum(squeeze(p.TwistedSF(1, 1, 3, 2:end-1, 1)).*p.BElmntMass(:, 1));
    #sid.md1_2_2_1= sum(squeeze(p.TwistedSF(1, 2, 3, 2:end-1, 1)).*p.BElmntMass(:, 1));

    p['J'][0,0] = sum(p['BElmntMass'][:]*(p['RNodes']+rh)**2) # TODO better integration
    p['J'][1,1] = sum(p['BElmntMass'][:]*(p['RNodes']+rh)**2) # TODO


    # --- Elastic matrices
    p['Ke'] = np.zeros((nq,nq))
    p['De'] = np.zeros((nq,nq))
    p['Me'] = np.zeros((nq,nq))
    # Me
    p['Me'][0,0]= p['MBF'][0, 0]
    p['Me'][0,1]= p['MBF'][0, 1]
    p['Me'][1,0]= p['MBF'][1, 0]
    p['Me'][1,1]= p['MBF'][1, 1]
    p['Me'][2,2]= p['MBE'][0, 0]
    # Ke
    p['Ke'][0,0]= p['KBF'][0, 0]
    p['Ke'][0,1]= p['KBF'][0, 1]
    p['Ke'][1,0]= p['KBF'][1, 0]
    p['Ke'][1,1]= p['KBF'][1, 1]
    p['Ke'][2,2]= p['KBE'][0, 0]
    # De
    p['De'][0,0]= p['CBF'][0, 0]
    p['De'][0,1]= p['CBF'][0, 1]
    p['De'][1,0]= p['CBF'][1, 0]
    p['De'][1,1]= p['CBF'][1, 1]
    p['De'][2,2]= p['CBE'][0, 0]

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
    # M1: nq x 6 x nq
    #Oe = np.zeros((nf,3,3))
    #Oe6= np.zeros((nf,6))
    # TODO TODO TODO
    # Below is for flap1 + edge1 not flap1 flap2 edge
    # and ordering is nq x nq x 6
    o=0 # derivative order 0=shape
    Oe_M1_1_1_1= - sum(np.squeeze(p['TwistedSF'][1, 0, 1:-1, o])*np.squeeze(p['TwistedSF'][1, 0, 1:-1, o])*p['BElmntMass']) + p['KBFCent'][0, 0];
    Oe_M1_1_1_2= - sum(np.squeeze(p['TwistedSF'][0, 0, 1:-1, o])*np.squeeze(p['TwistedSF'][0, 0, 1:-1, o])*p['BElmntMass']) + p['KBFCent'][0, 0];
    Oe_M1_1_1_4= 2*sum(np.squeeze(p['TwistedSF'][0, 0, 1:-1, o])*np.squeeze(p['TwistedSF'][1, 0, 1:-1, o])*p['BElmntMass'])
    Oe_M1_1_2_1= - sum(np.squeeze(p['TwistedSF'][1, 0, 1:-1, o])*np.squeeze(p['TwistedSF'][1, 2, 1:-1, o])*p['BElmntMass'])
    Oe_M1_1_2_2= - sum(np.squeeze(p['TwistedSF'][0, 0, 1:-1, o])*np.squeeze(p['TwistedSF'][0, 2, 1:-1, o])*p['BElmntMass'])
    Oe_M1_1_2_4=   sum(np.squeeze(p['TwistedSF'][0, 0, 1:-1, o])*np.squeeze(p['TwistedSF'][1, 2, 1:-1, o])*p['BElmntMass']) + sum(np.squeeze(p['TwistedSF'][1, 0, 1:-1, o])*np.squeeze(p['TwistedSF'][0, 2, 1:-1, o])*p['BElmntMass'])
    Oe_M1_2_1_1= - sum(np.squeeze(p['TwistedSF'][1, 0, 1:-1, o])*np.squeeze(p['TwistedSF'][1, 2, 1:-1, o])*p['BElmntMass'])
    Oe_M1_2_1_2= - sum(np.squeeze(p['TwistedSF'][0, 0, 1:-1, o])*np.squeeze(p['TwistedSF'][0, 2, 1:-1, o])*p['BElmntMass'])
    Oe_M1_2_1_4=   sum(np.squeeze(p['TwistedSF'][0, 0, 1:-1, o])*np.squeeze(p['TwistedSF'][1, 2, 1:-1, o])*p['BElmntMass']) + sum(np.squeeze(p['TwistedSF'][1, 0, 1:-1, o])*np.squeeze(p['TwistedSF'][0, 2, 1:-1, o])*p['BElmntMass'])
    Oe_M1_2_2_1= - sum(np.squeeze(p['TwistedSF'][1, 2, 1:-1, o])*np.squeeze(p['TwistedSF'][1, 2, 1:-1, o])*p['BElmntMass']) + p['KBECent'][0,0];
    Oe_M1_2_2_2= - sum(np.squeeze(p['TwistedSF'][0, 2, 1:-1, o])*np.squeeze(p['TwistedSF'][0, 2, 1:-1, o])*p['BElmntMass']) + p['KBECent'][0,0];
    Oe_M1_2_2_4= 2*sum(np.squeeze(p['TwistedSF'][0, 2, 1:-1, o])*np.squeeze(p['TwistedSF'][1, 2, 1:-1, o])*p['BElmntMass']);

    #for j in range(nf):
    #    sxx = trapzs(s_G[0,:]*U[j][0,:]*m)
    #    sxy = trapzs(s_G[0,:]*U[j][1,:]*m)
    #    sxz = trapzs(s_G[0,:]*U[j][2,:]*m)
    #    syx = trapzs(s_G[1,:]*U[j][0,:]*m)
    #    syy = trapzs(s_G[1,:]*U[j][1,:]*m)
    #    syz = trapzs(s_G[1,:]*U[j][2,:]*m)
    #    szx = trapzs(s_G[2,:]*U[j][0,:]*m)
    #    szy = trapzs(s_G[2,:]*U[j][1,:]*m)
    #    szz = trapzs(s_G[2,:]*U[j][2,:]*m)
    #    Gr[j][0,:] = 2*np.array([ syy+szz, -syx  , -szx     ])
    #    Gr[j][1,:] = 2*np.array([ -sxy   ,sxx+szz, -szy     ])
    #    Gr[j][2,:] = 2*np.array([ -sxz   , -syz  , sxx+syy  ])

    #    Oe[j] = -0.5*Gr[j].T
    #    Oe6[j][0] = Oe[j][0,0]
    #    Oe6[j][1] = Oe[j][1,1]
    #    Oe6[j][2] = Oe[j][2,2]
    #    Oe6[j][3] = Oe[j][0,1] + Oe[j][1,0] 
    #    Oe6[j][4] = Oe[j][1,2] + Oe[j][2,1]
    #    Oe6[j][5] = Oe[j][0,2] + Oe[j][2,0]
    #import pdb; pdb.set_trace()
    return p

def towerParameters(EDfilename, gravity, RotMass=None):
    """
    Compute tower parameters exactly like OpenFAST
    See Routine Coeff from ElastoDyn.f90
    """
    # --- Read inputs
    ED      = FASTInputFile(EDfilename)
    EDtwr   = os.path.join(os.path.dirname(EDfilename), ED['TwrFile'].replace('"',''))
    twr     = FASTInputFile(EDtwr)
    twrProp = twr.toDataFrame()
    n = ED['TwrNodes']

    # --- 
    p=dict()
    p['TwrFlexL']  = ED['TowerHt'] - ED['TowerBsHt'] # Height / length of the flexible portion of the tower.

    # New nodal positions
    n = ED['TwrNodes']
    twr_fract    = np.arange(1./n/2., 1, 1./n)
    p['DHNodes'] = np.ones(n)*p['TwrFlexL']/n
    p['HNodes']  = twr_fract*(ED['TowerHt']-ED['TowerBsHt']) + ED['TowerBsHt']
    p['HNodesNorm'] = p['HNodes']/p['TwrFlexL']
    # Interpolate properties to new nodal positions
    HtFract= twrProp['HtFract_[-]'].values
    p['MassT']    = np.interp(p['HNodesNorm'], HtFract, twrProp['TMassDen_[kg/m]']);
    p['StiffTFA'] = np.interp(p['HNodesNorm'], HtFract, twrProp['TwFAStif_[Nm^2]']);
    p['StiffTSS'] = np.interp(p['HNodesNorm'], HtFract, twrProp['TwSSStif_[Nm^2]']);
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
    p['TElmntMass'] = np.zeros(n)
    p['TMssAbvNd'] = np.zeros(n)
    p['TwrMass']=0
    for J in np.arange(n-1,-1,-1): # Loop through the tower nodes / elements in reverse
        # Calculate the mass of the current element
        p['TElmntMass'][J]    = p['MassT'][J]*p['DHNodes'][J];     # Mass of tower element J
        # Integrate to find the tower mass which will be output in .fsm
        p['TwrMass'] += p['TElmntMass'][J];
        # Integrate to find TMssAbvNd:
        p['TMssAbvNd'][J] = 0.5*p['TElmntMass'][J];
        if J == n-1:   # Uppermost tower element
            # Add the TwrTpMass effects:
            p['TMssAbvNd'][J] = p['TMssAbvNd'][J];
        else:         # All other tower elements
            # Add to TMssAbvNd'][J] the effects from the (not yet used) portion of element J+1
            p['TMssAbvNd'][J] = 0.5*p['TElmntMass'][J+1] + p['TMssAbvNd'][J] + p['TMssAbvNd'][J+1];
# end loop on J - Tower nodes / elements in reverse
# Initialize the generalized tower masses using tower-top mass effects:
    p['MTFA']= np.zeros((2, 2))
    p['MTSS']= np.zeros((2, 2))
    for I in [0,1]:  # Loop through all tower modes in a single direction
        p['MTFA'][I,I] = p['TwrTpMass'];
        p['MTSS'][I,I] = p['TwrTpMass'];
    nModesPerDir=2
    nDeriv =3
    p['TwrFASF'] = np.zeros((nModesPerDir, n+2, nDeriv)) # NOTE: full (+2)
    p['TwrSSSF'] = np.zeros((nModesPerDir, n+2, nDeriv)) # NOTE: full (+2)
    p['AxRedTFA']= np.zeros((nModesPerDir, nModesPerDir, n+2)) # NOTE: full (+2)
    p['AxRedTSS']= np.zeros((nModesPerDir, nModesPerDir, n+2)) # NOTE: full (+2)
    p['KTFA']       = np.zeros((nModesPerDir, nModesPerDir))
    p['KTSS']       = np.zeros((nModesPerDir, nModesPerDir))
    p['KTFAGrav']   = np.zeros((nModesPerDir, nModesPerDir))
    p['KTSSGrav']   = np.zeros((nModesPerDir, nModesPerDir))
    p['KTFAGravTT'] = np.zeros((nModesPerDir, nModesPerDir))
    p['KTSSGravTT'] = np.zeros((nModesPerDir, nModesPerDir))
    # set values for tower base (note that we haven't corrctly defined the values for (:,0,2) in the arrays below):
    #p['TwrFASF'] [:,0  ,:2] = 0.0;
    #p['TwrSSSF'] [:,0  ,:2] = 0.0;
    #p['AxRedTFA'][:,:2 ,0]  = 0.0;
    #p['AxRedTSS'][:,:2 ,0]  = 0.0;
    
    for J in np.arange(n):   # Loop through the tower nodes / elements
        # Calculate the tower shape functions (all derivatives):
        p['TwrFASF'][0,J+1,2] = SHP( p['HNodesNorm'][J], p['TwrFlexL'], p['TwFAM1Sh'], 2);
        p['TwrFASF'][1,J+1,2] = SHP( p['HNodesNorm'][J], p['TwrFlexL'], p['TwFAM2Sh'], 2);
        p['TwrFASF'][0,J+1,1] = SHP( p['HNodesNorm'][J], p['TwrFlexL'], p['TwFAM1Sh'], 1);
        p['TwrFASF'][1,J+1,1] = SHP( p['HNodesNorm'][J], p['TwrFlexL'], p['TwFAM2Sh'], 1);
        p['TwrFASF'][0,J+1,0] = SHP( p['HNodesNorm'][J], p['TwrFlexL'], p['TwFAM1Sh'], 0);
        p['TwrFASF'][1,J+1,0] = SHP( p['HNodesNorm'][J], p['TwrFlexL'], p['TwFAM2Sh'], 0);
        p['TwrSSSF'][0,J+1,2] = SHP( p['HNodesNorm'][J], p['TwrFlexL'], p['TwSSM1Sh'], 2);
        p['TwrSSSF'][1,J+1,2] = SHP( p['HNodesNorm'][J], p['TwrFlexL'], p['TwSSM2Sh'], 2);
        p['TwrSSSF'][0,J+1,1] = SHP( p['HNodesNorm'][J], p['TwrFlexL'], p['TwSSM1Sh'], 1);
        p['TwrSSSF'][1,J+1,1] = SHP( p['HNodesNorm'][J], p['TwrFlexL'], p['TwSSM2Sh'], 1);
        p['TwrSSSF'][0,J+1,0] = SHP( p['HNodesNorm'][J], p['TwrFlexL'], p['TwSSM1Sh'], 0);
        p['TwrSSSF'][1,J+1,0] = SHP( p['HNodesNorm'][J], p['TwrFlexL'], p['TwSSM2Sh'], 0);
        # Integrate to find the generalized mass of the tower (including tower-top mass effects).
        #   Ignore the cross-correlation terms of p['MTFA (i.e. p['MTFA(i,j) where i ~= j) and p['MTSS
        #   since these terms will never be used.
        for I in [0,1]:     # Loop through all tower DOFs in one direction
            p['MTFA'][I,I] += p['TElmntMass'][J]*p['TwrFASF'][I,J+1,0]**2;
            p['MTSS'][I,I] += p['TElmntMass'][J]*p['TwrSSSF'][I,J+1,0]**2;
        # Integrate to find the generalized stiffness of the tower (not including gravitational effects).
        ElStffFA = p['StiffTFA'][J]*p['DHNodes'][J] # Fore-aft stiffness of tower element J
        ElStffSS = p['StiffTSS'][J]*p['DHNodes'][J] # Side-to-side stiffness of tower element J
        for I in [0,1]:    # Loop through all tower DOFs in one direction
            for L in [0,1]:  # Loop through all tower DOFs in one direction
                p['KTFA'][I,L] += ElStffFA *p['TwrFASF'][I,J+1,2]*p['TwrFASF'][L,J,2];
                p['KTSS'][I,L] += ElStffSS *p['TwrSSSF'][I,J+1,2]*p['TwrSSSF'][L,J,2];
        # Integrate to find the gravitational-term of the generalized stiffness of the tower.
        #   Ignore the cross-correlation terms of p['KTFAGrav (i.e. p['KTFAGrav(i,j) where i ~= j)
        #   and p['KTSSGrav since these terms will never be used.
        ElmntStff = -p['TMssAbvNd'][J]*p['DHNodes'][J]*1  # stiffness of tower element J due to unity acceleration
        for I in [0,1]:    # Loop through all tower DOFs in one direction
            p['KTFAGrav'][I,I] += ElmntStff*p['TwrFASF'][I,J+1,1]**2
            p['KTSSGrav'][I,I] += ElmntStff*p['TwrSSSF'][I,J+1,1]**2
        # Integrate to find the tower axial reduction shape functions:
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
    #end # J - Tower nodes / elements
    ElmntStff      = -1*p['TwrFlexL']*1; # due to unitiy mass under unity acceleration stiffness of tower element J due to unity acceleration
    for I in [0,1]:
        p['KTFAGravTT'][I,I] = ElmntStff*p['TwrFASF'][I,-2,1]**2;
        p['KTSSGravTT'][I,I] = ElmntStff*p['TwrSSSF'][I,-2,1]**2;
    # Apply the modal stiffness tuners of the tower to KTFA() and KTSS():
    for I in [0,1]:     # Loop through all tower DOFs in one direction
        for L in [0,1]:  # Loop through all tower DOFs in one direction
            p['KTFA'][I,L] = np.sqrt( p['FAStTunr'][I]*p['FAStTunr'][L] )*p['KTFA'][I,L]
            p['KTSS'][I,L] = np.sqrt( p['SSStTunr'][I]*p['SSStTunr'][L] )*p['KTSS'][I,L]
    # Calculate the tower natural frequencies:
    p['FreqTFA'] = np.zeros((nModesPerDir, 2))
    p['FreqTSS'] = np.zeros((nModesPerDir, 2))
    p['CTFA']    = np.zeros((nModesPerDir, 2))
    p['CTSS']    = np.zeros((nModesPerDir, 2))
    for I in [0,1]:     # Loop through all tower DOFs in one direction
        allKTFAGrav= (p['KTFAGrav'][I,I] + p['TwrTpMass']*p['KTFAGravTT'][I,I])*gravity
        allKTSSGrav= (p['KTSSGrav'][I,I] + p['TwrTpMass']*p['KTSSGravTT'][I,I])*gravity
        p['FreqTFA'][I,0] = (1/2/np.pi)*np.sqrt(   p['KTFA'][I,I]                /( p['MTFA'][I,I] - p['TwrTpMass'] ) ) # Natural tower I-fore-aft frequency w/o gravitational destiffening nor tower-top mass effects
        p['FreqTFA'][I,1] = (1/2/np.pi)*np.sqrt( ( p['KTFA'][I,I] + allKTFAGrav )/  p['MTFA'][I,I]                    ) # Natural tower I-fore-aft frequency w/  gravitational destiffening and tower-top mass effects
        p['FreqTSS'][I,0] = (1/2/np.pi)*np.sqrt(   p['KTSS'][I,I]                /( p['MTSS'][I,I] - p['TwrTpMass'] ) ) # Natural tower I-side-to-side frequency w/o gravitational destiffening nor tower-top mass effects
        p['FreqTSS'][I,1] = (1/2/np.pi)*np.sqrt( ( p['KTSS'][I,I] + allKTSSGrav )/  p['MTSS'][I,I]                    ) # Natural tower I-side-to-side frequency w/  gravitational destiffening and tower-top mass effects
    # Calculate the generalized damping of the tower:
    for I in [0,1]:     # Loop through all tower DOFs in one direction
        for L in [0,1]:  # Loop through all tower DOFs in one direction
            p['CTFA'][I,L] = ( 0.01*p['TwrFADmp'][L] )*p['KTFA'][I,L]/( np.pi*p['FreqTFA'][L,0] );
            p['CTSS'][I,L] = ( 0.01*p['TwrSSDmp'][L] )*p['KTSS'][I,L]/( np.pi*p['FreqTSS'][L,0] );
    # Calculate the tower shape functions (all derivatives) at the tower-top:
    p['TwrFASF'][0,-1,2] = SHP( 1.0, p['TwrFlexL'], p['TwFAM1Sh'], 2);
    p['TwrFASF'][1,-1,2] = SHP( 1.0, p['TwrFlexL'], p['TwFAM2Sh'], 2);
    p['TwrFASF'][0,-1,1] = SHP( 1.0, p['TwrFlexL'], p['TwFAM1Sh'], 1);
    p['TwrFASF'][1,-1,1] = SHP( 1.0, p['TwrFlexL'], p['TwFAM2Sh'], 1);
    p['TwrFASF'][0,-1,0] = SHP( 1.0, p['TwrFlexL'], p['TwFAM1Sh'], 0);
    p['TwrFASF'][1,-1,0] = SHP( 1.0, p['TwrFlexL'], p['TwFAM2Sh'], 0);
    p['TwrSSSF'][0,-1,2] = SHP( 1.0, p['TwrFlexL'], p['TwSSM1Sh'], 2);
    p['TwrSSSF'][1,-1,2] = SHP( 1.0, p['TwrFlexL'], p['TwSSM2Sh'], 2);
    p['TwrSSSF'][0,-1,1] = SHP( 1.0, p['TwrFlexL'], p['TwSSM1Sh'], 1);
    p['TwrSSSF'][1,-1,1] = SHP( 1.0, p['TwrFlexL'], p['TwSSM2Sh'], 1);
    p['TwrSSSF'][0,-1,0] = SHP( 1.0, p['TwrFlexL'], p['TwSSM1Sh'], 0);
    p['TwrSSSF'][1,-1,0] = SHP( 1.0, p['TwrFlexL'], p['TwSSM2Sh'], 0);
    # Integrate to find the tower axial reduction shape functions at the tower-top:
    for I in [0,1]:     # Loop through all tower DOFs in one direction
        for L in [0,1]:  # Loop through all tower DOFs in one direction
            p['AxRedTFA'][I,L,-1] = p['AxRedTFA'][I,L,-2] + AxRdTFAOld[I,L]
            p['AxRedTSS'][I,L,-1] = p['AxRedTSS'][I,L,-2] + AxRdTSSOld[I,L]
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
    return p


if __name__ == '__main__':
#     EDfilename='../yams/_Jens/FEMBeam_NewFASTCoeffs/data/NREL5MW_ED_Onshore.dat'
#     EDfilename='../../data/NREL5MW/onshore/NREL5MW_ED_Onshore.dat'
#     bladeParameters(EDfilename)

    EDfilename  = '../../data/NREL5MW/onshore/NREL5MW_ED_Onshore.dat'
    FSTfilename = '../../data/NREL5MW/Main_Onshore.fst'

    fst = FASTInputFile(FSTfilename)

    p,pbld = rotorParameters(EDfilename)
    print(p)
    print(fst['Gravity'])

    ptwr=towerParameters(EDfilename, gravity=fst['Gravity'], RotMass=p['RotMass'])

    # Calculate the turbine mass:
    ptwr['TurbMass']  = ptwr['TwrTpMass'] + ptwr['TwrMass'];


