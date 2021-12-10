import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# Local 
from welib.weio.fast_input_file import FASTInputFile
from welib.tools.clean_exceptions import *
from welib.tools.tictoc import Timer


class SubDyn(object):
    """
    Tools for SubDyn

    - Setup a FEM model, compute Guyan and CB modes
    - Get a dataframe with properties
    - More todo

    """

    def __init__(self, sdFilename=None, sdData=None):
        super(SubDyn, self).__init__()

        # Read SubDyn file
        if sdFilename is not None:
            self.File = FASTInputFile(sdFilename)
        if sdData is not None:
            self.File = sdData

        self.M_tip=None
        # internal
        self._graph=None

    def __repr__(self):
        s='<{} object>:\n'.format(type(self).__name__)
        s+='|properties:\n'
        s+='|- File: (input file data)\n'
        s+='|methods:\n'
        s+='|- setTopMass\n'
        s+='|- beamDataFrame, beamFEM, beamModes\n'
        s+='|- toYAMS\n'
        return s


    def setTopMass(self):
        # TODO
        # Add an optional top mass and ineria
        if TopMass:
            # NOTE: you can use welib.yams.windturbine to compute RNA mass and inertia
            Mtop = 50000  # Top mass [kg]
            M_tip= rigidBodyMassMatrixAtP(m=Mtop, J_G=None, Ref2COG=None)
        else:
            M_tip=None

    @property
    def graph(self):
        import copy
        if self._graph is None:
            with Timer('Setting up graph'):
                print('>>> NDIV is ', self.File['NDiv'])
                self._graph = self.File.toGraph().divideElements(self.File['NDiv'],
                        excludeDataKey='Type', excludeDataList=['Cable','Rigid'], method='insert',
                        keysNotToCopy=['IBC','RBC','addedMassMatrix'] # Very important
                        )
            # Satinization
            #idBC_Fixed    =  0 # Fixed BC
            #idBC_Internal = 10 # Free/Internal BC
            #idBC_Leader   = 20 # Leader DOF
            MapIBC={0:0, 1:20} # 1 in the input file is leader
            MapRBC={0:10, 1:0} # 1 in the input file is fixed
            for n in self._graph.Nodes:
                if 'IBC' in n.data.keys():
                    IBC = n.data['IBC']
                    n.data['IBC'] = [MapIBC[i] for i in IBC[:6]]
                if 'RBC' in n.data.keys():
                    RBC = n.data['RBC']
                    n.data['RBC'] = [MapRBC[i] for i in RBC[:6]]
                    if any(RBC)==0:
                        raise NotImplementedError('SSI')

        return copy.deepcopy(self._graph)

    # --------------------------------------------------------------------------------}
    # --- Functions for beam-like structure (Spar, Monopile)
    # --------------------------------------------------------------------------------{
    def beamDataFrame(self, equispacing=False):
        """ """
        # --- Parameters
        UseSubDynModel=True
        TopMass = False

        # Convert to "welib.fem.Graph" class to easily handle the model (overkill for a monopile)
        df = self.graph.sortNodesBy('z').nodalDataFrame()

        if equispacing:
            from welib.tools.pandalib import pd_interp1
            # Interpolate dataframe to equispaced values
            xOld  = df['z']    # NOTE: FEM uses "x" as main axis
            nSpan = len(xOld)
            x = np.linspace(np.min(xOld),np.max(xOld), nSpan)
            df = pd_interp1(x, 'z', df)

        x   = df['z'] # NOTE: FEM uses "x" as main axis
        D   = df['D'] # Diameter [m]
        t   = df['t'] # thickness [m]
        # Derive section properties for a hollow cylinder based on diameter and thickness
        A        = np.pi*( (D/2)**2 - (D/2-t)**2) # Area for annulus [m^2]
        I        = np.pi/64*(D**4-(D-2*t)**4)     # Second moment of area for annulus (m^4)
        Kt       = I                              # Torsion constant, same as I for annulus [m^4]
        Ip       = 2*I                            # Polar second moment of area [m^4]
        df['A']  = A
        df['I']  = I
        df['Kt'] = Kt
        df['Ip'] = Ip
        df['m']  = df['rho'].values*A

        return df

    def beamFEM(self, df=None):
        """ return FEM model for beam-like structures, like Spar/Monopile"""
        import welib.FEM.fem_beam as femb

        BC       = 'clamped-free' # TODO Boundary condition: free-free or clamped-free
        element  = 'frame3d'      # Type of element used in FEM

        if df is None:
            df = self.beamDataFrame()
        x   = df['z']              # NOTE: FEM uses "x" as main axis
        E   = df['E']              # Young modules [N/m^2]
        G   = df['G']              # Shear modules [N/m^2]
        rho = df['rho']            # material density [kg/m^3]
        Ip  = df['Ip']
        I   = df['I']
        A   = df['A']
        Kt  = df['Kt']

        # --- Compute FEM model and mode shapes
        with Timer('Setting up FEM model'):
            FEM=femb.cbeam(x,m=rho*A,EIx=E*Ip,EIy=E*I,EIz=E*I,EA=E*A,A=A,E=E,G=G,Kt=Kt,
                        element=element, BC=BC, M_tip=self.M_tip)
        return FEM

    def beamModes(self, nCB=8, FEM = None):
        """ Returns mode shapes for beam-like structures, like Spar/Monopile """
        import welib.FEM.fem_beam as femb
        element  = 'frame3d'      # Type of element used in FEM
        if FEM is None:
            FEM = self.beamFEM()
        # --- Perform Craig-Bampton reduction, fixing the top node of the beam
        with Timer('FEM eigenvalue analysis'):
            Q_G,_Q_CB, df_G, df_CB, Modes_G, Modes_CB, CB = femb.CB_topNode(FEM, nCB=nCB, element=element, main_axis='x')
        # df_CB.to_csv('_CB.csv',index=False)
        # df_G.to_csv('_Guyan.csv',index=False)
        return  Q_G,_Q_CB, df_G, df_CB, Modes_G, Modes_CB, CB 

    def beamModesPlot(self):
        """ """
        # TODO
        nModesPlot=8
        # --- Show frequencies to screen
        print('Mode   Frequency  Label ')
        for i in np.arange(8):
            print('{:4d} {:10.3f}   {:s}'.format(i+1,FEM['freq'][i],FEM['modeNames'][i]))

        # --- Plot mode components for first few modes
        print(x.shape)
        #Q=FEM['Q'] ; modeNames = FEM['modeNames']
        #Q=Q_CB ;modeNames = names_CB
        Modes=Modes_CB
        nModesPlot=min(len(Modes),nModesPlot)

        fig,axes = plt.subplots(1, nModesPlot, sharey=False, figsize=(12.4,2.5))
        fig.subplots_adjust(left=0.04, right=0.98, top=0.91, bottom=0.11, hspace=0.40, wspace=0.30)
        for i in np.arange(nModesPlot):
            key= list(Modes.keys())[i]

            axes[i].plot(x, Modes[key]['comp'][:,0]  ,'-'  , label='ux')
            axes[i].plot(x, Modes[key]['comp'][:,1]  ,'-'  , label='uy')
            axes[i].plot(x, Modes[key]['comp'][:,2]  ,'-'  , label='uz')
            axes[i].plot(x, Modes[key]['comp'][:,3]  ,':'  , label='vx')
            axes[i].plot(x, Modes[key]['comp'][:,4]  ,':'  , label='vy')
            axes[i].plot(x, Modes[key]['comp'][:,5]  ,':'  , label='vz')
            axes[i].set_xlabel('')
            axes[i].set_ylabel('')
            axes[i].set_title(Modes[key]['label'])
            if i==0:
                axes[i].legend()

    # --------------------------------------------------------------------------------}
    # --- Functions for general FEM model (jacket, flexible floaters)
    # --------------------------------------------------------------------------------{
    def FEM(self, TP=(0,0,0), gravity = 9.81):
        """ return FEM model for beam-like structures, like Spar/Monopile

        TP: position of transition point
        """
        import welib.FEM.fem_beam as femb
        import welib.FEM.fem_model as femm
        BC       = 'clamped-free' # TODO Boundary condition: free-free or clamped-free
        element  = 'frame3d'      # Type of element used in FEM

        FEMMod = self.File['FEMMod']
        if FEMMod==1:
            element='frame3d'
        elif FEMMod==2:
            element='frame3dlin'
        elif FEMMod==3:
            element='timoshenko'
        else:
            raise NotImplementedError()
        print('>>> Main Element', FEMMod)
        # Get graph
        graph = self.graph
        print('>>> graph\n',graph)
        #graph.toJSON('_GRAPH.json')
        # Convert to FEM model
        with Timer('From graph'):
            model = femm.FEMModel.from_graph(self.graph, mainElementType=element, TP=TP, gravity=gravity)
        #model.toJSON('_MODEL.json')
        with Timer('Assembly'):
            model.assembly()
        with Timer('Internal constraints'):
            model.applyInternalConstraints()
            model.partition()
        with Timer('BC'):
            model.applyFixedBC()
        with Timer('EIG'):
            model.eig(normQ='byMax')
        with Timer('CB'):
            model.CraigBampton(nModesCB = self.File['Nmodes'], zeta=self.File['JDampings']/100 )

        model.rigidBodyEquivalent()

        #print(model)
        #print(model.DOF2Nodes)
        #print(graph)
#         n2e = graph.nodes2Elements
#         conn = graph.connectivity
#         print(n2e)
#         print(conn)

        ## --- Compute FEM model and mode shapes
        #with Timer('Setting up FEM model'):
        #    FEM=femb.cbeam(x,m=rho*A,EIx=E*Ip,EIy=E*I,EIz=E*I,EA=E*A,A=A,E=E,G=G,Kt=Kt,
        #                element=element, BC=BC, M_tip=self.M_tip)


        return model




    def toYAMSData(self, shapes=[0,4], main_axis='z'):
        """ 
        Convert to Data needed to setup a Beam Model in YAMS (see bodies.py in yams)
        """
        from welib.mesh.gradient import gradient_regular

        # --- Perform Craig-Bampton reduction, fixing the top node of the beam
        # Get beam data frame
        df = self.beamDataFrame(equispacing=True)
        if np.any(df['y']!=0): 
            raise NotImplementedError('FASTBeamBody for substructure only support monopile, structure not fully vertical in file: {}'.format(self.File.filename))
        if np.any(df['x']!=0): 
            raise NotImplementedError('FASTBeamBody for substructure only support monopile, structure not fully vertical in file: {}'.format(self.File.filename))

        FEM = self.beamFEM(df)
        Q_G,_Q_CB, df_G, df_CB, Modes_G, Modes_CB, CB = self.beamModes(nCB=0, FEM=FEM)

        x     = df['z'].values
        nSpan = len(x)

        # TODO TODO finda way to use these matrices instead of the ones computed with flexibility
        #print('CB MM\n',CB['MM'])
        #print('CB KK\n',CB['KK'])

        # --- Setup shape functions
        if main_axis=='x':
            raise NotImplementedError('')
        else:
            pass
            # we need to swap the CB modes
        nShapes=len(shapes)
        PhiU = np.zeros((nShapes,3,nSpan)) # Shape
        PhiV = np.zeros((nShapes,3,nSpan)) # Shape
        PhiK = np.zeros((nShapes,3,nSpan)) # Shape
        dx=np.unique(np.around(np.diff(x),4))
        if len(dx)>1:
            print(x)
            print(dx)
            raise NotImplementedError()
        for iShape, idShape in enumerate(shapes):
            if idShape==0:
                # shape 0 "ux"  (uz in FEM)
                PhiU[iShape][0,:] = df_G['G3_uz'].values
                PhiV[iShape][0,:] =-df_G['G3_ty'].values
                PhiK[iShape][0,:] = gradient_regular(PhiV[iShape][0,:],dx=dx[0],order=4)
            elif idShape==1:
                # shape 1,  "uy"
                PhiU[iShape][1,:] = df_G['G2_uy'].values
                PhiV[iShape][1,:] = df_G['G2_tz'].values
                PhiK[iShape][1,:] = gradient_regular(PhiV[iShape][1,:],dx=dx[0],order=4)
            elif idShape==4:
                # shape 4,  "vy"  (vz in FEM)
                PhiU[iShape][0,:] = df_G['G6_uy'].values
                PhiV[iShape][0,:] = df_G['G6_tz'].values
                PhiK[iShape][0,:] = gradient_regular(PhiV[iShape][0,:],dx=dx[0],order=4)
            else:
                raise NotImplementedError()

        # --- Dictionary structure for YAMS
        p=dict()
        p['s_span']=x-np.min(x)
        p['s_P0']=np.zeros((3,nSpan))
        if main_axis=='z':
            p['s_P0'][2,:]=x-np.min(x)
            p['r_O']   = (df['x'].values[0], df['y'].values[0], df['z'].values[0])
            p['R_b2g'] = np.eye(3)
        p['m']  = df['m'].values
        p['EI'] = np.zeros((3,nSpan))
        if main_axis=='z':
            p['EI'][0,:]=df['E'].values*df['I'].values
            p['EI'][1,:]=df['E'].values*df['I'].values
        p['jxxG']  = df['rho']*df['Ip']          # TODO verify
        p['s_min'] = p['s_span'][0]
        p['s_max'] = p['s_span'][-1]
        p['PhiU']  = PhiU
        p['PhiV']  = PhiV
        p['PhiK']  = PhiK

        # --- Damping
        damp_zeta     = None
        RayleighCoeff = None
        DampMat       = None
        if self.File['GuyanDampMod']==1:
            # Rayleigh Damping
            RayleighCoeff=self.File['RayleighDamp']
            #if RayleighCoeff[0]==0:
            #    damp_zeta=omega*RayleighCoeff[1]/2. 
        elif self.File['GuyanDampMod']==2:
            # Full matrix
            DampMat = self.File['GuyanDampMatrix']
            DampMat=DampMat[np.ix_(shapes,shapes)]

        return p, damp_zeta, RayleighCoeff, DampMat



