"""
A FEM Model is a special kind of Graph
   Graphs have: Nodes, Elements, Node properties, and a "Model" storing it all

"""
import numpy as np
import pandas as pd
import copy
from scipy.optimize import OptimizeResult as OdeResultsClass 

from welib.FEM.utils import DCM, rigidTransformationMatrix
from welib.FEM.fem_elements import *   # Elements used
from welib.FEM.fem_core import insertFixedBCinModes
from welib.FEM.graph import NodeProperty
from welib.FEM.graph import GraphModel 

# class MaterialProperty(NodeProperty):
#     def __init__(self):
#         Property.__init__(self)
#         pass

# --------------------------------------------------------------------------------}
# --- Main FEM model class  
# --------------------------------------------------------------------------------{
class FEMModel(GraphModel):
    @classmethod
    def from_graph(cls, graph, ndiv=None, mainElementType='frame3d', refPoint=None, gravity=9.81):
        """ 
        Returns FEM model from graph
        DOFs number are attributed to graph
        """
        # Copying graph to avoid any access of memory by reference
        g = copy.deepcopy(graph)
        # Creating "FEM" Nodes without DOFs
        for i,n in enumerate(g.Nodes):
            g.Nodes[i] = FEMNode(n.ID, n.x, n.y, n.z, DOFs=[], **n.data)
        # Creating "FEM" Elements without DOFs
        for i,e in enumerate(g.Elements):
            if e.data['Type']=='Beam':
                if mainElementType=='frame3d':
                    g.Elements[i] = SubDynFrame3dElement(e.ID, e.nodeIDs, None, e.propset, e.propIDs, None, **e.data)
                elif mainElementType=='timoshenko':
                    g.Elements[i] = SubDynTimoshenko3dElement(e.ID, e.nodeIDs, None, e.propset, e.propIDs, None, **e.data)
                else:
                    raise NotImplementedError()
            elif e.data['Type']=='Cable':
                if mainElementType in ['frame3d', 'timoshenko']:
                    g.Elements[i] = SubDynCable3dElement(e.ID, e.nodeIDs, None, e.propset, e.propIDs, None, **e.data)
                else:
                    raise NotImplementedError()
            elif e.data['Type']=='Rigid':
                if mainElementType in ['frame3d', 'timoshenko']:
                    g.Elements[i] = SubDynRigid3dElement(e.ID, e.nodeIDs, None, e.propset, e.propIDs, None, **e.data)
                else:
                    raise NotImplementedError()
            else:
                raise NotImplementedError()
        # Update connectivity
        g.connecticityHasChanged()
        g.updateConnectivity()
        # Distribute DOF to nodes and elements
        g = distributeDOF(g, mainElementType=mainElementType)
        # Create instance of FEM Model
        self = cls(Elements=g.Elements, Nodes=g.Nodes, NodePropertySets=g.NodePropertySets, 
                ElemPropertySets=g.ElemPropertySets, MiscPropertySets=g.MiscPropertySets, refPoint=refPoint, gravity=gravity)
        #print(self)
        return self

    @classmethod
    def from_cbeam(cls, FEM, CB=None):
        """ 
        Convert a simple fem_beam representation as returned by cbeam
        into a more advanced FEM Model

        - FEM: dictionary, output from cbeam
        - CB: dictionary, output from CB_topNode
        """
        #         FEM={'xNodes':xNodes, 'MM':MM, 'KK':KK, 'MM_full':MM_,'KK_full':KK_,'Tr':Tr,
        #             'IFull2BC':IFull2BC, 'IBC2Full':IBC2Full,
        #             'Elem2Nodes':Elem2Nodes, 'Nodes2DOF':Nodes2DOF, 'Elem2DOF':Elem2DOF,
        #             'Q':Q,'freq':freq, 'modeNames':modeNames}
        Nodes2DOF = FEM['Nodes2DOF']
        Elem2Nodes = FEM['Elem2Nodes']
        # --- Create nodes and elements
        Nodes=[]
        for i,n in enumerate(FEM['xNodes'].T):
            Nodes.append(FEMNode(i, n[0], n[1], n[2], DOFs=Nodes2DOF[i]))
        nNodes = len(Nodes)
        Elements = []
        for ie in range(nNodes-1):
            nodes = [Nodes[Elem2Nodes[ie,0]] , Nodes[Elem2Nodes[ie,1]]]
            Elements.append(FEMElement(ie, nodeIDs=Elem2Nodes[ie], nodes=nodes))
        # --- Creating FEM Model
        self = cls(Elements=Elements, Nodes=Nodes) #refPoint=refPoint, gravity=gravity)
        # Providing data directly
        self.MM       = FEM['MM_full'] # Mass matrix   without boundary conditions
        self.KK       = FEM['KK_full'] # Stiff matrix  without boundary conditions
        self.CC       = None # Damping matrix ..
        self.FF       = None # 
        self.freq     = FEM['freq'] # Frequency of the Full FEM system after Boundary conditions and internal constraints
        self.Q        = FEM['Q']   # Modes of the full FEM system after Boundary conditions and internal constraints
        #self.T_c = FEM['Tr']
        self.T_c = np.eye(self.MM.shape[0]) # NOTE: Q is full already
        # Set CB data directly
        if CB is not None:
            self.MM_CB  = CB['MM']
            self.KK_CB  = CB['KK']
            self.Phi_G  = CB['Phi_G']
            self.Phi_CB = CB['Phi_CB']
            self.f_G    = CB['f_G']
            self.f_CB   = CB['f_CB']
            self.Q_G    = CB['Q_G']
            self.Q_CB   = CB['Q_CB']
        self.setModes()
        return self



    def __init__(self, Elements=None, Nodes=None, NodePropertySets=None, ElemPropertySets=None, MiscPropertySets=None,
            mainElementType='frame3d', gravity=9.81, main_axis='z', refPoint=None): # FEM specific
        """ 

         - refPoint: reference point from which a rigid body connection to the interface is computed
        """

        GraphModel.__init__(self, Elements=Elements, Nodes=Nodes, NodePropertySets=NodePropertySets,
                ElemPropertySets=ElemPropertySets, MiscPropertySets=MiscPropertySets)

        self.refPoint        = refPoint
        self.gravity         = gravity
        self.main_axis       = main_axis
        self.mainElementType = mainElementType
        # Main data generated
        self.MM_init  = None # Initial mass matrix, before internal constaints
        self.KK_init  = None # Initial stiffness matrix, before internal constaints
        self.FF_init  = None # Initial load vector     , before internal constaints
        self.MM       = None # Mass matrix   without boundary conditions
        self.KK       = None # Stiff matrix  without boundary conditions
        self.CC       = None # Damping matrix ..
        self.FF       = None # 
        self.FExt     = None # 
        self.freq     = None # Frequency of the Full FEM system after Boundary conditions and internal constraints
        self.Q        = None # Modes of the full FEM system after Boundary conditions and internal constraints
                             # NOTE: augmented to include the boundary condition in them
         # CB reduction data
        self.MM_CB  = None
        self.KK_CB  = None
        self.Phi_G  = None
        self.Phi_CB = None
        self.f_G    = None
        self.f_CB   = None
        self.Q_G    = None
        self.Q_CB   = None

        # internal data
        self._nDOF    = None

        # 
        if Nodes is not None:
            #print('>>>', self.Nodes[0].data)
            if 'DOFs' not in self.Nodes[0].data.keys():
                print('>>> Attributing DOFs to Nodes and Elements')
                g = distributeDOF(self, mainElementType=mainElementType)
            elif len(self.Nodes[0].data['DOFs'])==0:
                print('>>> Attributing DOFs to Nodes and Elements')
                g = distributeDOF(self, mainElementType=mainElementType)

    # --------------------------------------------------------------------------------}
    # --- Handling of DOFs nodes elements 
    # --------------------------------------------------------------------------------{
    @property
    def nDOF(self):
        if self._nDOF is None:
            self._nDOF = np.max(self.Nodes[-1].data['DOFs'])+1
        return self._nDOF
    @property
    def nDOFc(self):
        DOFs = []
        for n in self.Nodes:
            DOFs += n.data['DOFs_c']
        return np.max(DOFs)+1


    @property
    def DOF2Nodes(self):
        """ 
        Return mapping from DOF to nodes.
        Using "DOF" for system *before* internal constraints have been eliminated

        Return array nDOF x 4, with columns: 
          0: DOF index (not ID),   (starts at 0)
          1: Node index (not ID),  (starts at 0)
          2: Number of DOF on node,
          3: DOF number out of all the DOF on node (starts at 1)
        """
        DOF2Nodes=np.zeros((self.nDOF,4),int)
        for iN,node in enumerate(self.Nodes):
            for iiDOF,iDOF in enumerate(node.data['DOFs']):
                DOF2Nodes[iDOF,0] = iDOF
                DOF2Nodes[iDOF,1] = iN
                DOF2Nodes[iDOF,2] = len(node.data['DOFs'])
                DOF2Nodes[iDOF,3] = iiDOF+1
        return DOF2Nodes

    @property
    def DOFc2Nodes(self):
        """ 
        Return mapping from DOF to nodes.
        Using "DOF" for system where internal constraints have been eliminated

        Return array nDOF x 4, with columns: 
          0: DOF index (not ID),   (starts at 0)
          1: Node index (not ID),  (starts at 0)
          2: Number of DOF on node,
          3: DOF number out of all the DOF on node (starts at 1)
        """
        DOF2Nodes=np.zeros((self.nDOFc,4),int)
        for iN,node in enumerate(self.Nodes):
            DOFs = node.data['DOFs_c']
            #if len(DOFs)==0:
            #    if 'ID_link' in node.data.keys():
            #        print('>>>>>>  node' , node.ID, 'linked to', node.data['ID_link'] )
            #        DOFs = self.getNode(node.data['ID_link']).data['DOFs_c']
            for iiDOF,iDOF in enumerate(DOFs):
                DOF2Nodes[iDOF,0] = iDOF
                DOF2Nodes[iDOF,1] = iN
                DOF2Nodes[iDOF,2] = len(DOFs)
                DOF2Nodes[iDOF,3] = iiDOF+1
        return DOF2Nodes


    # Interface, reaction, internal nodes TODO, make this cleaner
    @property
    def interfaceNodes(self): return [n for n in self.Nodes if 'IBC' in n.data]
    @property
    def reactionNodes(self): return [n for n in self.Nodes if 'RBC' in n.data]
    @property
    def internalNodes(self): return [n for n in self.Nodes if 'RBC' not in n.data and 'IBC' not in n.data]
    
    # --------------------------------------------------------------------------------}
    # --- FEM
    # --------------------------------------------------------------------------------{
    def assembly(self,  gravity=None, Elements=None,):
        if Elements is None:
            Elements = self.Elements
        if gravity is None:
            gravity=self.gravity

        nDOF=self.nDOF # np.max(self.Nodes[-1].data['DOFs'])+1
        MM = np.zeros((nDOF,nDOF))
        KK = np.zeros((nDOF,nDOF))
        FF = np.zeros(nDOF)
        
        # loop over all elements, compute element matrices and assemble into global matrices
        for e in Elements:
            # --- Element mass, stiffness, gravity force and other force
            Ke   = e.Ke()
            Me   = e.Me()
            Fe_g = e.Fe_g(gravity)
            Fe_o = e.Fe_o()
            IDOF = e.data['DOFs']

            # --- Assembly in global unconstrained system
            IDOF = e.data['DOFs']
            KK[np.ix_(IDOF,IDOF)] += Ke
            MM[np.ix_(IDOF,IDOF)] += Me
            FF[IDOF] += Fe_g + Fe_o

        # Add concentrated masses to Mass matrix and gravity vector
        for n in self.Nodes:
            if 'addedMassMatrix' in n.data.keys():
                IDOF = n.data['DOFs']
                MM[np.ix_(IDOF,IDOF)] += n.data['addedMassMatrix']
                FF[IDOF[2]] -= n.data['addedMassMatrix'][0,0]*gravity # gravity along z DOF index "2"

        self.KK_init= KK
        self.MM_init= MM 
        self.FF_init= FF
        self.resetExternalLoads()
        return self

    def rayleighDamping(self, alpha, beta):
        # TODO DD_BC?
        if self.CC is None:
            self.CC = alpha * self.MM + beta * self.KK
        if self.MM_BC is not None:
            self.CC_BC = alpha * self.MM_BC + beta * self.KK_BC

    # --------------------------------------------------------------------------------}
    # --- Direct elimination of Internal Constraints (Rigid links and rotational joints)
    # --------------------------------------------------------------------------------{
    def applyInternalConstraints(self):
        """ 
        Apply internal constraints such as rigid links and rotational joints
        using a direct elminination technique

        - Tc: reduction matrix such that x_init= Tc.x     dimension: n x n_c
              where x is the reduced vector of DOF and x_init is the intial vector (more DOF)

        """
        rotJoints  = [n.ID for n in self.Nodes    if n.data['Type']!=idJointCantilever]
        rigidLinks = [e.ID for e in self.Elements if e.data['TypeID']==idMemberRigid    ]
        if len(rotJoints)>0 or len(rigidLinks)>0:
            print('Number of Rotational joints:',len(rotJoints))
            print('Number of Rigid Links      :',len(rigidLinks))
            from .direct_elimination import nDOF_c, buildTMatrix, rigidLinkAssemblies
            RA = rigidLinkAssemblies(self)
            self.T_c = buildTMatrix(self, RA)

        else:
            self.T_c = np.eye(self.MM_init.shape[0])
            # Store new DOF indices
            for n in self.Nodes:
                n.data['DOFs_c'] = list(n.data['DOFs'])
                
        # Distribute Nodal DOFs to Element DOFs
        for e in self.Elements:
            e.data['DOFs_c'] =[]
            for n in e.nodes:
                e.data['DOFs_c'] +=n.data['DOFs_c']

        self.MM = (self.T_c.T).dot(self.MM_init).dot(self.T_c)
        self.KK = (self.T_c.T).dot(self.KK_init).dot(self.T_c)
        self.FF = (self.T_c.T).dot(self.FF_init)  # TODO verify order

        # --- Creating a convenient Map from DOF to Nodes
        #p%DOFred2Nodes=-999
        #do iNode=1,p%nNodes
        #   nDOFPerNode = len(p%NodesDOFred(iNode))
        #   do iiDOF = 1, nDOFPerNode
        #      iDOF = p%NodesDOFred(iNode)%List(iiDOF)
        #      p%DOFred2Nodes(iDOF,1) = iNode       ! First column is Node index
        #      p%DOFred2Nodes(iDOF,2) = nDOFPerNode ! Second column is number of DOF per node
        #      p%DOFred2Nodes(iDOF,3) = iiDOF       ! Third column is number of DOF per node


    def partition(self):
        """
        Partition DOFs into leader/follower and fixed, following the order convention of SubDyn.
        
        Intermediate variables:
            Partition DOFs and Nodes into sets: 
            Nodes are partitioned into the I,C,L (and R) sets, Nodes_I, Nodes_C, Nodes_L, with:
                    I="Interface" nodes
                    C="Reaction" nodes
                    L=Interior nodes
                    R=I+C
            DOFs indices are partitioned into B, F, L
                    B=Leader DOFs (Rbar in SubDyn documentation)
                    F=Fixed DOFS
                    L=Interior DOFs
            Subpartitions of both categories use the convention: "NodePartition_DOFPartition"
               e.g. C_F : "reaction" nodes DOFs that are fixed
                    C_L : "reaction" nodes DOFs that will be counted as internal
                    I_B : "interface" nodes DOFs that are leader DOFs
        """
        reactionNodes  = self.reactionNodes
        interfaceNodes = self.interfaceNodes
        # --- Count DOFs - NOTE: we count node by node
        nDOF___  = sum([len(n.data['DOFs_c'])                                 for n in self.Nodes])
        # Interface DOFs
        nDOFI__  = sum([len(n.data['DOFs_c'])                       for n in interfaceNodes])
        nDOFI_B = sum([sum(np.array(n.data['IBC'])==idDOF_Leader)   for n in interfaceNodes])
        # DOFs of reaction nodes
        nDOFC__ = sum([len(n.data['DOFs_c'])                        for n in reactionNodes]) 
        nDOFC_B = sum([sum(np.array(n.data['RBC'])==idDOF_Leader)   for n in reactionNodes])
        # Total number of DOFs in each category:
        self.nDOFR__ = nDOFI__ + nDOFC__ # Total number, used to be called "nDOFR"
        self.nDOF__B = nDOFC_B + nDOFI_B

        # --- Distibutes the I, L, C nodal DOFs into  B, F, L sub-categories 
        # NOTE: order is important for compatibility with SubDyn
        # TODO: list comprehension
        IDI__ = []
        IDI_B = []
        IDI_F = []
        for n in interfaceNodes:
            IDI__ += n.data['DOFs_c'] # NOTE: respects order
            IDI_B += [dof for i,dof in enumerate(n.data['DOFs_c']) if n.data['IBC'][i]==idDOF_Leader]
            IDI_F += [dof for i,dof in enumerate(n.data['DOFs_c']) if n.data['IBC'][i]==idDOF_Fixed ]
        IDI__ = IDI_B+IDI_F
        IDC__ = []
        IDC_B = []
        IDC_L = []
        IDC_F = []
        for n in reactionNodes:
            IDC__ += n.data['DOFs_c'] # NOTE: respects order
            IDC_B += [dof for i,dof in enumerate(n.data['DOFs_c']) if n.data['RBC'][i]==idDOF_Leader  ]
            IDC_L += [dof for i,dof in enumerate(n.data['DOFs_c']) if n.data['RBC'][i]==idDOF_Internal]
            IDC_F += [dof for i,dof in enumerate(n.data['DOFs_c']) if n.data['RBC'][i]==idDOF_Fixed   ]
        IDR=IDC__+IDI__
        IDL_L = []
        for n in self.internalNodes:
            IDL_L += n.data['DOFs_c']

        # --- Total indices per partition B, F, L
        self.DOFc_Leader     =         IDC_B + IDI_B  # boundary/retained/leader DOFs
        self.DOFc_Fixed      =         IDC_F + IDC_B  # Fixed DOFs
        self.DOFc_Follower   = IDL_L + IDC_L          # internal DOFs
        self.DOFc_Boundary   = IDR    # I+C Boundary nodes for rigid body equivalent
        self.DOFc_Internal   = IDL_L  # L   Internal nodes for rigid body equivalent
        self.DOFc_Interface  = IDI__  # I   Interface
        ID_ALL = self.DOFc_Leader + self.DOFc_Fixed + self.DOFc_Follower
        for i in np.arange(nDOF___):
            if i not in ID_ALL:
                raise Exception('DOF {} not found in DOF list after partition'.format(i))

    # --------------------------------------------------------------------------------}
    # ---  
    # --------------------------------------------------------------------------------{
    def applyFixedBC(self, IFixed=None):
        """ 
        Apply boundary conditions. (Fixed boundary conditions only)
        IFixed: Index of fixed DOF
        """
        # NOTE: we use the matrices where internal constraints have been eliminated 
        MM = self.MM
        KK = self.KK
        CC = self.CC
        nDOF_tot = MM.shape[0]
        IDOF_All = np.arange(0,nDOF_tot)

        # Tip and root degrees of freedom
        #IDOF_root = Nodes2DOF[Elem2Nodes[0,:][0] ,:]
        #IDOF_tip  = Nodes2DOF[Elem2Nodes[-1,:][1],:]

        # --- Boundary condition transformation matrix (removes row/columns)
        Tr=np.eye(nDOF_tot)
        
        # Root and Tip BC
        if IFixed is None:
            IFixed=[]
            for n in self.reactionNodes:
                I = n.data['RBC'][:6]
                IFixed += [n.data['DOFs_c'][ii] for ii,i in enumerate(I) if int(i)==idDOF_Fixed]

        Tr = np.delete(Tr, IFixed, axis=1) # removing columns

        Mr = (Tr.T).dot(MM).dot(Tr)
        Kr = (Tr.T).dot(KK).dot(Tr)
        if CC is not None:
            Cr = (Tr.T).dot(CC).dot(Tr)
        else:
            Cr = None

        # --- Create mapping from M to Mr
        nDOF_r = Mr.shape[0]
        IDOF_BC = list(np.setdiff1d(IDOF_All, IFixed))
        IFull2BC = np.zeros(nDOF_tot,dtype=int)
        IBC2Full = np.zeros(nDOF_r,dtype=int)
        k=0
        for i in IDOF_All:
            if i in IFixed:
                IFull2BC[i]=-1
            else:
                IFull2BC[i]=k
                IBC2Full[k]=i
                k+=1

        self.MM_BC = Mr
        self.KK_BC = Kr
        self.CC_BC = Cr
        self.FF_BC = Tr.T.dot(self.FF)
        self.T_BC  = Tr
        self.T_Full2BC  = Tr.T.dot(self.T_c.T)
        #
        self.IFull2BC = IFull2BC
        self.IBC2Full = IBC2Full
        self.DOFr_Leader   = [IFull2BC[i] for i in self.DOFc_Leader]
        self.DOFr_Follower = [IFull2BC[i] for i in self.DOFc_Follower]
        self.DOFr_Fixed    = [IFull2BC[i] for i in self.DOFc_Fixed]

        return Mr, Kr, Tr, IFull2BC, IBC2Full


    def eig(self, normQ='byMax'):
        from welib.system.eva import eig
        KK = self.KK_BC
        MM = self.MM_BC

        # --- Compute modes and frequencies
        [Q, freq]= eig(KK, MM, freq_out=True, normQ=normQ, discardIm=True)

        Q   = insertFixedBCinModes(Q, self.T_BC)
        self.freq = freq
        self.Q    = Q
        return Q, freq


    # --------------------------------------------------------------------------------}
    # --- Reference point 
    # --------------------------------------------------------------------------------{
    @property
    def T_refPoint(self):
        """ Rigid body transformation matrix from interface DOFs to refpoint"""
        if self.refPoint is None: 
            raise Exception('Cannot compute T_refPoint, refPoint is None')
        return rigidTransformationMatrix(self.DOFc_Interface, self.refPoint, self.DOFc2Nodes, self.points)
    # --------------------------------------------------------------------------------}
    # --- General FEM Utils
    # --------------------------------------------------------------------------------{
    def setFullMatrices(self,MM,KK,CC=None):
        self.MM=MM
        self.KK=KK
        if CC is not None:
            self.CC=CC

    def CraigBampton(self, nModesCB=None, BC_before_CB=True):
        """
        Perform Craig-Bampton reduction

        nModesCB: number of CB modes to retain
        zeta :  damping ratios for CB modes 
        BC_before_CB: if true, using the matrices where the fixed BC have been applied 
        """
        from welib.FEM.reduction import CraigBampton
        if BC_before_CB:
            M, K = self.MM_BC, self.KK_BC
            Ileader, Ifollow = self.DOFr_Leader, self.DOFr_Follower
            if nModesCB is None:
                nModesCB=M.shape[0] - len(Ileader)
            # NOTE: we return all CB modes at first
            Mr, Kr, Phi_G, Phi_CB, f_G, f_CB, I1, I2 = CraigBampton(M, K, Ileader=Ileader, Ifollow=Ifollow, nModesCB=None, discardIm=True)
            # Small cleanup
            Phi_G [np.abs(Phi_G )<1e-11] = 0
            Phi_CB[np.abs(Phi_CB)<1e-11] = 0
            Mr    [np.abs(Mr    )<1e-11] = 0
            Kr    [np.abs(Kr    )<1e-11] = 0

            self.DOF_Leader_CB   = I1
            self.DOF_Follower_CB = I2
            self.MM_CB  = Mr
            self.KK_CB  = Kr
            self.Phi_G  = Phi_G
            self.Phi_CB = Phi_CB
            self.f_G    = np.real(f_G)
            self.f_CB   = f_CB[:nModesCB]
            omega_CB = 2*np.pi*f_CB[:nModesCB]

            self.nModesCB= nModesCB
        else:
            raise NotImplementedError()
        #if Ifixed is not None:
        #    M,K = self.applyFixBC()
        #else:
        #    M,K = self.MM, self.KK

    def rigidBodyEquivalent(self):
        """ 
        Compute Rigid body equivalent mass matrix at origin, 
        NOTE: Without SSI mass
        """
        # 
        from welib.yams.utils import identifyRigidBodyMM, rigidBodyMassMatrixAtP
        # Transformation matrix from leader DOFs to Origin
        TIR= rigidTransformationMatrix(self.DOFc_Boundary, (0,0,0), self.DOFc2Nodes, self.points)
        # Compute Rigid body mass matrix (without Soil, and using both Interface and Reactions nodes as leader DOF)
        if self.nDOFR__!=self.nDOF__B: # Most likely the case
            MBB = self.rigidBody() 
        else:
            MBB = model.MM_CB[np.ix_(model.DOF_Leader_CB, model.DOF_Leader_CB)]
        M_O = TIR.T.dot(MBB).dot(TIR)
        # Clean up for values that ought to be 0
        M_O[0,1:3] = 0; 
        M_O[1,0  ] = 0; M_O[1,2] = 0; M_O[1,4] = 0;
        M_O[2,0:1] = 0; M_O[2,5] = 0;
        M_O[3,0  ] = 0; M_O[4,1] = 0; M_O[5,2] = 0;
        M_O[np.abs(M_O)<1e-6]=0
        self.M_O            = M_O
        mass, J_G, Ref2COG  = identifyRigidBodyMM(M_O) # NOTE ref is (0,0,0)
        self.mass           = mass
        self.center_of_mass = Ref2COG
        self.M_G    = rigidBodyMassMatrixAtP(mass, J_G, (0,0,0))
        if self.refPoint is not None:
            self.M_ref = rigidBodyMassMatrixAtP(mass, J_G, -np.array(self.refPoint)+np.array(Ref2COG))

    def rigidBody(self):
        """ Extract rigid body mass without SSI
        Both "interface" and "reaction" nodes are fixed
        NOTE: effectively performs a Guyan reduction """
        from welib.FEM.reduction import CraigBampton
        # --- Remove SSI from Mass and stiffness matrix (NOTE: use NodesDOFred, reduced matrix)
        #CALL InsertSoilMatrices(Init%M, Init%K, p%NodesDOFred, Init, p, ErrStat2, ErrMsg2, Substract=.True.);
        # --- Perform Guyan reduction to get MBB
        Ileader = self.DOFc_Boundary
        Ifollow = self.DOFc_Internal
        print('leader',Ileader)
        print('follow',Ifollow)
        Mr, Kr, Phi_G, Phi_CB, f_G, f_CB, I1, I2 = CraigBampton(self.MM, self.KK, Ileader=Ileader, Ifollow=Ifollow, nModesCB=0)
        #! --- Insert SSI from Mass and stiffness matrix again
        #CALL InsertSoilMatrices(Init%M, Init%K, p%NodesDOFred, Init, p, ErrStat2, ErrMsg2, Substract=.False.); if(Failed()) return
        return Mr[np.ix_(I1,I1)] # return Guyan only

    
    # --------------------------------------------------------------------------------}
    # --- IO 
    # --------------------------------------------------------------------------------{
    def nodesDisp(self, UDOF_c, IDOF=None, scale=True, maxAmplitude=None, sortDim=None,):
        """ 
        Returns nNodes x 3 x nShapes array of nodal displacements 

        INPUTS:
          - UDOF: nDOF_c x nModes: array of DOF "displacements" for each mode
                  in the system where internal constraints have been eliminated
          - IDOF: Optional array of subset/reordered DOF. 1:nDOF_c if not provided
          - scale: if True, modes are shapes according based on `maxAmplitude`
          - maxAmplitude: if provided, scale used for the mode scaling. If not provided,
                     maxAmplitude is set to 10% of the maximum dimension of the structure
          - sortDim: 0,1,2: sort by x,y,z
        """
        if IDOF is None:
            IDOF = list(np.arange(self.nDOF))
        if maxAmplitude is None:
            maxAmplitude = self.maxDimension * 0.1 # 10% of max dimension of the structure
        if True:
            DOF2Nodes = self.DOF2Nodes
            UDOF      = self.T_c.dot(UDOF_c)
        else:
            DOF2Nodes = self.DOF2Nodes
            UDOF      = UDOF_c
        # dimension: n x n_c
        #self.extent
        #self.points
        # --- 
        INodes = list(np.sort(np.unique(DOF2Nodes[IDOF,1]))) # Sort nodes
        nShapes = UDOF.shape[1]
        disp = np.empty((len(INodes),3,nShapes)); disp.fill(np.nan)
        pos  = np.empty((len(INodes),3))         ; pos.fill(np.nan)

        # --- METHOD 1 - Loop through DOFs KEEP ME
        #print('>>> fem_model: TODO problem with pos, order appears wrong')
        #for i,iDOF in enumerate(IDOF):
        #    iNode       = DOF2Nodes[iDOF,1]
        #    nDOFPerNode = DOF2Nodes[iDOF,2]
        #    nodeDOF     = DOF2Nodes[iDOF,3]
        #    iiNode      = INodes.index(iNode)
        #    node = self.Nodes[iNode-1]
        #    if nodeDOF<=3:
        #        pos[iiNode, 0]= node.x
        #        pos[iiNode, 1]= node.y
        #        pos[iiNode, 2]= node.z
        #        for iShape in np.arange(nShapes):
        #            disp[iiNode, nodeDOF-1, iShape] = UDOF[i, iShape]
        # --- METHOD 2 - Loop through Nodes
        Ix=[]; Iy=[]; Iz=[]
        for i,n in enumerate(self.Nodes):
            pos[i, 0]= n.x
            pos[i, 1]= n.y
            pos[i, 2]= n.z
            Ix.append(n.data['DOFs'][0])
            Iy.append(n.data['DOFs'][1])
            Iz.append(n.data['DOFs'][2])
        for iShape in np.arange(nShapes):
            disp[:, 0, iShape] = UDOF[Ix, iShape]
            disp[:, 1, iShape] = UDOF[Iy, iShape]
            disp[:, 2, iShape] = UDOF[Iz, iShape]

        # Scaling 
        if scale:
            for iShape in np.arange(nShapes):
                maxDisp=np.nanmax(np.abs(disp[:, :, iShape]))
                if maxDisp>1e-5:
                    disp[:, :, iShape] *= maxAmplitude/maxDisp
        # Sorting according to a dimension
        if sortDim is not None: 
            I=np.argsort(pos[:,sortDim])
            INodes = np.array(INodes)[I]
            disp   = disp[I,:,:]
            pos    = pos[I,:]
        return disp, pos, INodes


    def getModes(self, scale=True, maxAmplitude=None, sortDim=None):
        """ return Guyan and CB modes

          - maxAmplitude: if provided, scale used for the mode scaling. If not provided,
                     maxAmplitude is set to 10% of the maximum dimension of the structure
          - sortDim: 0,1,2: sort by x,y,z

        """
        if maxAmplitude is None:
            maxAmplitude = self.maxDimension * 0.1 # 10% of max dimension of the structure

        # CB modes
        PhiM     = self.Phi_CB
        PhiM_aug = np.zeros((self.nDOFc, PhiM.shape[1]))
        PhiM_aug[self.DOFc_Follower, : ] = PhiM
        dispCB, posCB, INodesCB = self.nodesDisp(PhiM_aug, scale=scale, maxAmplitude=maxAmplitude, sortDim=sortDim)

        # Guyan modes
        PhiR     = self.Phi_G
        PhiR_aug = np.zeros((self.nDOFc, PhiR.shape[1]))
        for i in np.arange(len(self.DOFc_Leader)):
            PhiR_aug[self.DOFc_Leader[i] , i] = 1
        PhiR_aug[self.DOFc_Follower, : ] = PhiR
        PhiR_Intf = PhiR_aug.dot(self.T_refPoint) # nDOF x 6 (since TI is nGY x 6)
        dispGy, posGy, INodesGy = self.nodesDisp(PhiR_Intf, scale=scale, maxAmplitude=maxAmplitude, sortDim=sortDim)

        return dispGy, posGy, INodesGy, dispCB, posCB, INodesCB


    def setModes(self, nModesFEM=30, nModesCB=None):

        # FEM Modes
        if self.Q is not None:
            dispFEM, posFEM, INodesFEM = self.nodesDisp(self.Q)
            for iMode in range(min(dispFEM.shape[2], nModesFEM)):
                self.addMode(displ=dispFEM[:,:,iMode], name='FEM{:d}'.format(iMode+1), freq=self.freq[iMode], group='FEM')

        # GY CB Modes
        if self.f_CB is not None:
            if nModesCB is None:
                nModesCB = len(self.f_CB)

            if self.Q_G is not None and self.Q_CB is not None:
                dispGy, posGy, INodesGy = self.nodesDisp(self.Q_G)
                dispCB, posCB, INodesCB = self.nodesDisp(self.Q_CB)
            else:
                dispGy, posGy, InodesGy, dispCB, posCB, InodesCB = self.getModes(sortDim=None) 
            for iMode in range(dispGy.shape[2]):
                self.addMode(displ=dispGy[:,:,iMode], name='GY{:d}'.format(iMode+1), freq=self.f_G[iMode], group='GY')

            for iMode in range(min(len(self.f_CB), nModesCB)):
                self.addMode(displ=dispCB[:,:,iMode], name='CB{:d}'.format(iMode+1), freq=self.f_CB[iMode], group='CB') 


    def toYAML(model, filename=None, more=False, ioff=1):
        """ 
        Write a YAML summary file, similar to SubDyn
        see welib.fast.subdyn. subdyntoYAMLSum

        """
        def yaml_array(var, M, Fmt='{:15.6e}', comment=''):
            M = np.atleast_2d(M)
            if len(comment)>0:
                s='{}: # {} x {} {}\n'.format(var, M.shape[0], M.shape[1], comment)
            else:
                s='{}: # {} x {}\n'.format(var, M.shape[0], M.shape[1])

            if M.shape[0]==1:
                if M.shape[1]==0:
                    s+= '  - [ ]\n'
                else:
                    for l in M:
                        s+= '  - [' + ','.join([Fmt.format(le) for le in l]) + ',]\n'
            else:
                for l in M:
                    s+= '  - [' + ','.join([Fmt.format(le) for le in l]) + ']\n'
            s = s.replace('e+','E+').replace('e-','E-')
            return s

        # --- Helper functions
        def nodeID(nodeID):
            if hasattr(nodeID,'__len__'):
                return [model.Nodes.index(model.getNode(n))+ioff for n in nodeID]
            else:
                return model.Nodes.index(model.getNode(nodeID))+ioff

        def elemID(elemID):
            #e=model.getElement(elemID)
            for ie,e in enumerate(model.Elements):
                if e.ID==elemID:
                    return ie+ioff
        def elemType(elemType):
            from welib.FEM.fem_elements import idMemberBeam, idMemberCable, idMemberRigid
            return {'SubDynBeam3d':idMemberBeam, 'SubDynFrame3d':idMemberBeam, 'Beam':idMemberBeam, 'Frame3d':idMemberBeam,
                    'SubDynTimoshenko3d':idMemberBeam,
                    'SubDynCable3d':idMemberCable, 'Cable':idMemberCable,
                    'Rigid':idMemberRigid,
                    'SubDynRigid3d':idMemberRigid}[elemType]

        def propID(propID, propset):
            prop = model.NodePropertySets[propset]
            for ip, p in enumerate(prop):
                if p.ID == propID:
                    return ip+ioff
# 
#     SD_Vars = subdynPartitionVars(model)

        # --- Helper functions
        s=''
        s += '#____________________________________________________________________________________________________\n'
        s += '# RIGID BODY EQUIVALENT DATA\n'
        s += '#____________________________________________________________________________________________________\n'
        if not hasattr(model, 'M_O'):
            model.rigidBodyEquivalent()
        s0 = 'Mass: {:15.6e} # Total Mass\n'.format(model.M_O[0,0])
        s += s0.replace('e+','E+').replace('e-','E-')
        s0 = 'CM_point: [{:15.6e},{:15.6e},{:15.6e},] # Center of mass coordinates (Xcm,Ycm,Zcm)\n'.format(model.center_of_mass[0],model.center_of_mass[1],model.center_of_mass[2])
        s += s0.replace('e+','E+').replace('e-','E-')
        if model.refPoint is not None:
            s0 = 'TP_point: [{:15.6e},{:15.6e},{:15.6e},] # Transition piece reference point\n'.format(model.refPoint[0],model.refPoint[1],model.refPoint[2])
            s += s0.replace('e+','E+').replace('e-','E-')
        s += yaml_array('MRB',  model.M_O,  comment = 'Rigid Body Equivalent Mass Matrix w.r.t. (0,0,0).')
        if model.refPoint is not None:
            s += yaml_array('M_P' , model.M_ref,comment = 'Rigid Body Equivalent Mass Matrix w.r.t. TP Ref point')
        s += yaml_array('M_G' , model.M_G,  comment = 'Rigid Body Equivalent Mass Matrix w.r.t. CM (Xcm,Ycm,Zcm).')
        s += '#____________________________________________________________________________________________________\n'
        s += '# GUYAN MATRICES at the TP reference point\n'
        s += '#____________________________________________________________________________________________________\n'
        #s += yaml_array('KBBt' , model.KBBt,  comment = '')
        #s += yaml_array('MBBt' , model.MBBt,  comment = '')
        #s += yaml_array('CBBt' , model.CBBt,  comment = '(user Guyan Damping + potential joint damping from CB-reduction)')
        s += '#____________________________________________________________________________________________________\n'
        s += '# SYSTEM FREQUENCIES\n'
        s += '#____________________________________________________________________________________________________\n'
        s += '#Eigenfrequencies [Hz] for full system, with reaction constraints (+ Soil K/M + SoilDyn K0) \n'
        s += yaml_array('Full_frequencies', model.freq)
#             s += '#Frequencies of Guyan modes [Hz]\n'
#             s += yaml_array('GY_frequencies', model.f_G)
#             s += '#Frequencies of Craig-Bampton modes [Hz]\n'
#             s += yaml_array('CB_frequencies', model.f_CB)
        s += '#____________________________________________________________________________________________________\n'
        s += '# Internal FEM representation\n'
        s += '#____________________________________________________________________________________________________\n'
        s += 'nNodes_I: {:7d} # Number of Nodes: "interface" (I)\n'.format(len(model.interfaceNodes))
        s += 'nNodes_C: {:7d} # Number of Nodes: "reactions" (C)\n'.format(len(model.reactionNodes))
        s += 'nNodes_L: {:7d} # Number of Nodes: "internal"  (L)\n'.format(len(model.internalNodes))
        s += 'nNodes  : {:7d} # Number of Nodes: total   (I+C+L)\n'.format(len(model.Nodes))
#             if more:
#                 s += 'nDOFI__ : {:7d} # Number of DOFs: "interface"          (I__)\n'.format(len(SD_Vars['IDI__']))
#                 s += 'nDOFI_B : {:7d} # Number of DOFs: "interface" retained (I_B)\n'.format(len(SD_Vars['IDI_B']))
#                 s += 'nDOFI_F : {:7d} # Number of DOFs: "interface" fixed    (I_F)\n'.format(len(SD_Vars['IDI_F']))
#                 s += 'nDOFC__ : {:7d} # Number of DOFs: "reactions"          (C__)\n'.format(len(SD_Vars['IDC__']))
#                 s += 'nDOFC_B : {:7d} # Number of DOFs: "reactions" retained (C_B)\n'.format(len(SD_Vars['IDC_B']))
#                 s += 'nDOFC_L : {:7d} # Number of DOFs: "reactions" internal (C_L)\n'.format(len(SD_Vars['IDC_L']))
#                 s += 'nDOFC_F : {:7d} # Number of DOFs: "reactions" fixed    (C_F)\n'.format(len(SD_Vars['IDC_F']))
#                 s += 'nDOFR__ : {:7d} # Number of DOFs: "intf+react"         (__R)\n'.format(len(SD_Vars['IDR__']))
#                 s += 'nDOFL_L : {:7d} # Number of DOFs: "internal"  internal (L_L)\n'.format(len(SD_Vars['IDL_L']))
#             s += 'nDOF__B : {:7d} # Number of DOFs:             retained (__B)\n'.format(SD_Vars['nDOF__B'])
#             s += 'nDOF__L : {:7d} # Number of DOFs:             internal (__L)\n'.format(SD_Vars['nDOF__L'])
#             s += 'nDOF__F : {:7d} # Number of DOFs:             fixed    (__F)\n'.format(SD_Vars['nDOF__F'])
#             s += 'nDOF_red: {:7d} # Number of DOFs: total\n'                     .format(SD_Vars['nDOF___'])
        s += yaml_array('Nodes_I', nodeID([n.ID for n in model.interfaceNodes]), Fmt='{:7d}', comment='"interface" nodes"');
        s += yaml_array('Nodes_C', nodeID([n.ID for n in model.reactionNodes ]), Fmt='{:7d}', comment='"reaction" nodes"');
        s += yaml_array('Nodes_L', nodeID([n.ID for n in model.internalNodes ]), Fmt='{:7d}', comment='"internal" nodes"');
#             if more:
#                 s += yaml_array('DOF_I__', np.array(SD_Vars['IDI__'])+1,   Fmt='{:7d}', comment = '"interface"           DOFs"')
#                 s += yaml_array('DOF_I_B', np.array(SD_Vars['IDI_B'])+1,   Fmt='{:7d}', comment = '"interface" retained  DOFs')
#                 s += yaml_array('DOF_I_F', np.array(SD_Vars['IDI_F'])+1,   Fmt='{:7d}', comment = '"interface" fixed     DOFs')
#                 s += yaml_array('DOF_C__', np.array(SD_Vars['IDC__'])+1,   Fmt='{:7d}', comment = '"reaction"            DOFs"')
#                 s += yaml_array('DOF_C_B', np.array(SD_Vars['IDC_B'])+1,   Fmt='{:7d}', comment = '"reaction"  retained  DOFs')
#                 s += yaml_array('DOF_C_L', np.array(SD_Vars['IDC_L'])+1,   Fmt='{:7d}', comment = '"reaction"  internal  DOFs')
#                 s += yaml_array('DOF_C_F', np.array(SD_Vars['IDC_F'])+1,   Fmt='{:7d}', comment = '"reaction"  fixed     DOFs')
#                 s += yaml_array('DOF_L_L', np.array(SD_Vars['IDL_L'])+1,   Fmt='{:7d}', comment = '"internal"  internal  DOFs')
#                 s += yaml_array('DOF_R_' , np.array(SD_Vars['IDR__'])+1,   Fmt='{:7d}', comment = '"interface&reaction"  DOFs')
#             s += yaml_array('DOF___B', np.array(model.DOFc_Leader  )+1, Fmt='{:7d}',  comment='all         retained  DOFs');
#             s += yaml_array('DOF___F', np.array(model.DOFc_Fixed   )+1, Fmt='{:7d}',  comment='all         fixed     DOFs');
#             s += yaml_array('DOF___L', np.array(model.DOFc_Follower)+1, Fmt='{:7d}',  comment='all         internal  DOFs');
#             s += '\n'
        s += '#Index map from DOF to nodes\n'
        s += '#     Node No.,  DOF/Node,   NodalDOF\n'
        s += 'DOF2Nodes: # {} x 3 (nDOFRed x 3, for each constrained DOF, col1: node index, col2: number of DOF, col3: DOF starting from 1)\n'.format(model.nDOFc)
        DOFc2Nodes = model.DOFc2Nodes
        for l in DOFc2Nodes:
            s +='  - [{:7d},{:7d},{:7d}] # {}\n'.format(l[1]+1, l[2], l[3], l[0]+1 )
        s += '#     Node_[#]          X_[m]           Y_[m]           Z_[m]       JType_[-]       JDirX_[-]       JDirY_[-]       JDirZ_[-]  JStff_[Nm/rad]\n'
        s += 'Nodes: # {} x 9\n'.format(len(model.Nodes))
        for n in model.Nodes:
            s += '  - [{:7d}.,{:15.3f},{:15.3f},{:15.3f},{:14d}.,   0.000000E+00,   0.000000E+00,   0.000000E+00,   0.000000E+00]\n'.format(nodeID(n.ID), n.x, n.y, n.z, int(n.data['Type']) )
        # NOTE: propIDs might not exist
#         s += '#    Elem_[#]    Node_1   Node_2   Prop_1   Prop_2     Type     Length_[m]      Area_[m^2]  Dens._[kg/m^3]        E_[N/m2]        G_[N/m2]       shear_[-]       Ixx_[m^4]       Iyy_[m^4]       Jzz_[m^4]          T0_[N]\n'
#         s += 'Elements: # {} x 16\n'.format(len(model.Elements))
#         for e in model.Elements:
#             I = e.inertias
#             s0='  - [{:7d}.,{:7d}.,{:7d}.,{:7d}.,{:7d}.,{:7d}.,{:15.3f},{:15.3f},{:15.3f},{:15.6e},{:15.6e},{:15.6e},{:15.6e},{:15.6e},{:15.6e},{:15.6e}]\n'.format(
#                 elemID(e.ID), nodeID(e.nodeIDs[0]), nodeID(e.nodeIDs[1]), propID(e.propIDs[0], e.propset), propID(e.propIDs[1], e.propset), elemType(e.data['Type']), 
#                 e.length, e.area, e.rho, e.E, e.G, e.kappa, I[0], I[1], I[2], e.T0)
#             s += s0.replace('e+','E+').replace('e-','E-')

#             s += '#____________________________________________________________________________________________________\n'
#             s += '#User inputs\n'
#             s += '\n'
#             s += '#Number of properties (NProps):{:6d}\n'.format(len(model.NodePropertySets['Beam']))
#             s += '#Prop No         YoungE         ShearG        MatDens          XsecD          XsecT\n'
#             for ip,p in enumerate(model.NodePropertySets['Beam']):
#                 s0='#{:8d}{:15.6e}{:15.6e}{:15.6e}{:15.6e}{:15.6e}\n'.format(p.ID, p['E'],p['G'],p['rho'],p['D'],p['t'])
#                 s +=  s0.replace('e+','E+').replace('e-','E-')
#             s +='\n'
#             s += '#No. of Reaction DOFs:{:6d}\n'.format(len(SD_Vars['IDC__']) )
#             s += '#React. DOF_ID    BC\n'
#             s += '\n'.join(['#{:10d}{:10s}'.format(idof+1,'     Fixed' ) for idof in SD_Vars['IDC_F']])
#             s += '\n'.join(['#{:10d}{:10s}'.format(idof+1,'     Free'  ) for idof in SD_Vars['IDC_L']])
#             s += '\n'.join(['#{:10d}{:10s}'.format(idof+1,'     Leader') for idof in SD_Vars['IDC_B']])
#             s += '\n\n'
#             s += '#No. of Interface DOFs:{:6d}\n'.format(len(SD_Vars['IDI__']))
#             s += '#Interf. DOF_ID    BC\n'
#             s += '\n'.join(['#{:10d}{:10s}'.format(idof+1,'    Fixed' ) for idof in SD_Vars['IDI_F']])
#             s += '\n'.join(['#{:10d}{:10s}'.format(idof+1,'    Leader') for idof in SD_Vars['IDI_B']])
#             s += '\n\n'
#             CM = []
#             from welib.yams.utils import identifyRigidBodyMM
#             for n in model.Nodes:
#                 if 'addedMassMatrix' in n.data:
#                     mass, J_G, ref2COG = identifyRigidBodyMM(n.data['addedMassMatrix'])
#                     CM.append( (n.ID, mass, J_G, ref2COG) )
#             s += '#Number of concentrated masses (NCMass):{:6d}\n'.format(len(CM))
#             s += '#JointCMas           Mass            JXX            JYY            JZZ            JXY            JXZ            JYZ           MCGX           MCGY           MCGZ\n'
#             for cm in CM:
#                 s0 = '# {:9.0f}.{:15.6e}{:15.6e}{:15.6e}{:15.6e}{:15.6e}{:15.6e}{:15.6e}{:15.6e}{:15.6e}{:15.6e}\n'.format( nodeID(cm[0]),  cm[1], cm[2][0,0], cm[2][1,1], cm[2][2,2], cm[2][0,1], cm[2][0,2], cm[2][1,2],cm[3][0],cm[3][1],cm[3][2] )
#                 s += s0.replace('e+','E+').replace('e-','E-')
#             s += '\n'
#             #s += '#Number of members    18\n'
#             #s += '#Number of nodes per member:     2\n'
#             #s += '#Member I Joint1_ID Joint2_ID    Prop_I    Prop_J           Mass         Length     Node IDs...\n'
#             #s += '#       77        61        60        11        11   1.045888E+04   2.700000E+00       19    18\n'
#             #s += '#____________________________________________________________________________________________________\n'
#             #s += '#Direction Cosine Matrices for all Members: GLOBAL-2-LOCAL. No. of 3x3 matrices=    18\n'
#             #s += '#Member I        DC(1,1)        DC(1,2)        DC(1,3)        DC(2,1)        DC(2,2)        DC(2,3)        DC(3,1)        DC(3,2)        DC(3,3)\n'
#             #s += '#       77  1.000E+00  0.000E+00  0.000E+00  0.000E+00 -1.000E+00  0.000E+00  0.000E+00  0.000E+00 -1.000E+00\n'
        s += '#____________________________________________________________________________________________________\n'
        s += '#FEM Eigenvectors ({} x {}) [m or rad], full system with reaction constraints (+ Soil K/M + SoilDyn K0)\n'.format(*model.Q.shape)
        s += yaml_array('Full_Modes', model.Q)
#             s += '#____________________________________________________________________________________________________\n'
#             s += '#CB Matrices (PhiM,PhiR) (reaction constraints applied)\n'
#             s += yaml_array('PhiM', model.Phi_CB[:,:model.nModesCB] ,comment='(CB modes)')
#             s += yaml_array('PhiR', model.Phi_G,  comment='(Guyan modes)')
#             s += '\n'
        if more:
            s += '#____________________________________________________________________________________________________\n'
            s += '# ADDITIONAL DEBUGGING INFORMATION\n'
            s += '#____________________________________________________________________________________________________\n'
            s +=  ''
            e = model.Elements[0]
#                 rho=e.rho
#                 A = e.area
#                 L = e.length
#                 t= rho*A*L
#                 s0 = '{:15.6e}{:15.6e}{:15.6e}{:15.6e}{:15.6e}{:15.6e}{:15.6e}{:15.6e}{:15.6e}{:15.6e}{:15.6e}\n'.format(model.gravity,e.area, e.length, e.inertias[0], e.inertias[1], e.inertias[2], e.kappa, e.E, e.G, e.rho, t)
#                 s0 = s0.replace('e+','E+').replace('e-','E-')
#                 s += s0
            s += yaml_array('KeLocal' +str(), model.Elements[0].Ke(local=True))
# 
#                 for ie,e in enumerate(model.Elements):
#                     s += yaml_array('DC' +str(ie+1), e.DCM.transpose())
#                     s += yaml_array('Ke' +str(ie+1), e.Ke())
#                     s += yaml_array('Me' +str(ie+1), e.Me())
#                     s += yaml_array('FGe'+str(ie+1), e.Fe_g(model.gravity))
#                     s += yaml_array('FCe'+str(ie+1), e.Fe_o())
# 
#                     s += yaml_array('KeLocal' +str(ie+1), e.Ke(local=True))
#                     s += yaml_array('MeLocal' +str(ie+1), e.Me(local=True))
#                     s += yaml_array('FGeLocal'+str(ie+1), e.Fe_g(model.gravity, local=True))
#                     s += yaml_array('FCeLocal'+str(ie+1), e.Fe_o(local=True))
# 
        s += '#____________________________________________________________________________________________________\n'
        e = model.Elements[0]
        s += yaml_array('Ke', e.Ke(local=True), comment='First element stiffness matrix'); # TODO not in local
        s += yaml_array('Me', e.Me(local=True), comment='First element mass matrix');
        s += yaml_array('FGe', e.Fe_g(model.gravity,local=True), comment='First element gravity vector');
        s += yaml_array('FCe', e.Fe_o(local=True), comment='First element cable pretension');
        s += '#____________________________________________________________________________________________________\n'
        s += '#FULL FEM K and M matrices. TOTAL FEM TDOFs:    {}\n'.format(model.nDOF); # NOTE: wrong in SubDyn, should be nDOFc
        s += yaml_array('K', model.KK, comment='Stiffness matrix');
        s += yaml_array('M', model.MM, comment='Mass matrix');
        s += '#____________________________________________________________________________________________________\n'
        s += '#Gravity and cable loads applied at each node of the system (before DOF elimination with T matrix)\n'
        s += yaml_array('FG', model.FF_init, comment=' ');
#                 s += '#____________________________________________________________________________________________________\n'
#                 s += '#Additional CB Matrices (MBB,MBM,KBB) (constraint applied)\n'
#                 s += yaml_array('MBB'    , model.MBB, comment='');
#                 s += yaml_array('MBM'    , model.MBM[:,:model.nModesCB], comment='');
#                 s += yaml_array('CMMdiag', model.CMM, comment='(2 Zeta OmegaM)');
#                 s += yaml_array('KBB'    , model.KBB, comment='');
#                 s += yaml_array('KMM'    , np.diag(model.KMM), comment='(diagonal components, OmegaL^2)');
#                 s += yaml_array('KMMdiag', np.diag(model.KMM)[:model.nModesCB], comment='(diagonal components, OmegaL^2)');
#                 s += yaml_array('PhiL'   , model.Phi_CB, comment='');
#                 s += 'PhiLOm2-1: # 18 x 18 \n'
#                 s += 'KLL^-1: # 18 x 18 \n'
        s += '#____________________________________________________________________________________________________\n'
        s += yaml_array('T_red', model.T_c, Fmt = '{:9.2e}', comment='(Constraint elimination matrix)');
#             s += 'AA: # 16 x 16 (State matrix dXdx)\n'
#             s += 'BB: # 16 x 48 (State matrix dXdu)\n'
#             s += 'CC: # 6 x 16 (State matrix dYdx)\n'
#             s += 'DD: # 6 x 48 (State matrix dYdu)\n'
        if model.refPoint is not None:
            s += '#____________________________________________________________________________________________________\n'
            s += yaml_array('TI', model.T_refPoint,  Fmt = '{:9.2e}',comment='(TP refpoint Transformation Matrix TI)');
        if filename is not None:
            with open(filename, 'w') as f:
                f.write(s);
        else:
            return s

    # --------------------------------------------------------------------------------}
    # --- LOADS
    # --------------------------------------------------------------------------------{
    def resetExternalLoads(self):
        self.FExt = np.zeros(len(self.FF_init))

    def setExternalLoadOnNode(self, node, Load):
        DOFs = node.data['DOFs'] # Before internal constraints
        self.FExt[DOFs] = Load # NOTE: Load vector should be in harmony with number of DOF

    def loadVector_BC(self, Fext, internal):
        """ Returns load vector suitable for system with BC included"""
        #FExt_c = (self.T_c.T).dot(Fext)
        #FExt_BC = (self.T_BC.T).dot(FExt_c)
        FExt_BC = (self.T_Full2BC).dot(Fext)
        if internal:
            F_BC = FExt_BC + self.FF_BC
        else:
            F_BC = FExt_BC
        return F_BC
    # --------------------------------------------------------------------------------}
    # --- STATIC ANALYSES
    # --------------------------------------------------------------------------------{
    def staticDisplacements(self, F=None, internal=True):
        if F is None:
            F =  self.FExt
        F_BC = self.loadVector_BC(F, internal=internal)

        # Solve for displacements
        U_BC = np.linalg.solve(self.KK_BC, F_BC)

        U_c = self.T_BC.dot(U_BC)
        U   = self.T_c.dot(U_c)
        return U, U_c, U_BC


    # --------------------------------------------------------------------------------}
    # --- DYNAMIC ANALYSES
    # --------------------------------------------------------------------------------{
    def accelerations(self, gz, gzp, F):
        """ 
        gz:  DOFs
        gzp: velocities
        """
        gz  = gz.flatten()
        gzp = gzp.flatten()
        F   = F.flatten()
        RHS  = -self.KK_BC.dot(gz) -self.CC_BC.dot(gzp) + F
        #     gzpp = M_U\(M_L\RHS)  ;
        #gzpp = np.linalg.solve(self.MM_BC,RHS) # TODO LU Decomposition
        gzpp = self.MM_BC_inv.dot(RHS) # TODO LU Decomposition
        return gzpp

    def dqdt_firstOrder(self, t, q):
        q = q.flatten()
        nq  = len(q)
        nx  = int(nq/2)
        if np.mod(int(t*1000000),10000)==0:
            print('>>> dqdt {:15.6f} {:15.3f} {:15.3f}'.format(t, min(np.max(q),1e8), min(np.max(q[nx-6]),1e8)))
        dqdt_ = np.zeros(q.shape)
        x  = q[0:nx]
        xd = q[nx:]

        # Loads
        self.resetExternalLoads()
        if self._force_fn is not None:
            Fext = self._force_fn(t, x, xd, self)
        else:
            Fext = self.Fext
        F_BC = self.loadVector_BC(Fext, internal=True)
        xdd = self.accelerations(x, xd, F_BC)
        dqdt_[0:nx] =  xd
        dqdt_[nx:]  =  xdd
        return dqdt_




    def integrate(self, t_eval, gz0=None, gzp0=None, fext=None, algorithm ='solve_ivp', method='RK45', **options):
        print('----------------- TIME INTEGRATION --------------------------')
        self.MM_BC_inv=np.linalg.inv(self.MM_BC)

        if gz0 is None:
            gz0 = np.zeros(self.KK_BC.shape[0])
        if gzp0 is None:
            gzp0 = np.zeros(self.KK_BC.shape[0])
        q0 = np.concatenate((gz0.flatten(), gzp0.flatten()))

        self._force_fn = fext

        if algorithm=='solve_ivp':
            from scipy.integrate import  solve_ivp #odeint
            print(q0)
            print(q0.shape)
            odefun = self.dqdt_firstOrder
            res = solve_ivp(fun=odefun, t_span=[t_eval[0], t_eval[-1]], y0=q0, t_eval=t_eval, method=method, vectorized=True, max_step=0.01, atol=1, rtol=1)
        elif algorithm=='welib':
            from welib.ode import integrate
            res = integrate(self.dqdt_firstOrder, t_eval, q0, method=method)
        else:
            raise NotImplementedError()
        return res
#     % [M_L,M_U] = lu(Mr);
#     % [T,Y] = ode23(@(t,y)fFirstOrder(t,y,@fAccTowerFEM,[],Mr,Kr,Dr,M_L,M_U),vt,zeros(2*nDOF_tot,1));
#     % [T,Y] = ode23(@(t,y)fYDotFEM(t,y,Mr,Kr,Dr,M_U,M_L),vt,zeros(2*nDOF_tot,1));
#     [T,Y] = fodeNewmark(@(t,x,xp)fMDKR_FEM(t,x,xp,Mr,Kr,Dr),vt,zeros(2*nDOF_tot,1),options);
#     %% Computing tip-displacement

# % --- Sub Function for ode solver
# function [M,D,K,F] = fMDKR_FEM(t,gz,gzp,M,K,D)
#     global FORCE_FUNCTION
#     global FORCE_POINT
#     nDOF_tot=size(M,1);
#     % --- Prescribed force
#     [Fz] = FORCE_FUNCTION(t);
#     F=zeros(nDOF_tot,1);
#     if ~isequal(FORCE_POINT,'top'); error('Only top force for now'); end
#     % Last Translational DOF gets force
#     F(end-1)=Fz;


    def __repr__(self):
        s=super(FEMModel, self).__repr__()
        s+='\nFEM Model specific\n'
        s+='- gravity  : {}\n'.format(self.gravity)
        s+='- main_axis: {}\n'.format(self.main_axis)
        def printshape(M,sM,s1='-'):
            s=''
            if M is not None:
                M = np.asarray(M)
                if len(M.shape)==2:
                    s=s1+' {:5s}: shape {} x {}\n'.format(sM, M.shape[0], M.shape[1])
                else:
                    s=s1+' {:5s}: shape {}\n'.format(sM, M.shape[0])
            return s
        s+=printshape(self.MM, 'MM')
        s+=printshape(self.KK, 'KK')
        s+=printshape(self.CC, 'CC')
        s+=printshape(self.FF, 'FF')
        s+=printshape(self.MM_BC, 'MM_BC')
        s+=printshape(self.KK_BC, 'KK_BC')
        s+=printshape(self.FF_BC, 'FF_BC')
        s+=printshape(self.freq, 'freq')
        s+=printshape(self.Q, 'Q')
        s+='* nDOF: {}\n'.format(self.nDOF)
        s+='* nDOFc: {}\n'.format(self.nDOFc)
        s+=printshape(self.DOF2Nodes, 'DOF2Nodes','*')
        s+=printshape(self.DOFc2Nodes, 'DOFc2Nodes','*')
        s+='* interfaceNodes: {}\n'.format(self.interfaceNodes)
        s+='* reactionNodes: {}\n'.format(self.reactionNodes)
        s+=printshape(self.internalNodes, 'internalNodes','*')
        return s

# --------------------------------------------------------------------------------}
# --- Helper functions 
# --------------------------------------------------------------------------------{
def distributeDOF(g, mainElementType='frame3d'):
    """ 
    Given a list of Nodes and Elements, distribute degrees of freedom (DOFs) 
    to each node and elements.
     For Cantilever Joint -> Condensation into 3 translational and 3 rotational DOFs
     For other joint type -> Condensation of the 3 translational DOF
                          -> Keeping 3 rotational DOF for each element connected to the joint

    INPUTS:
       g: graph
    """

    if mainElementType in ['frame3d', 'timoshenko','timoshenko3d','euler-bernoulli3d']:
        nDOFPerElem  = 12
        nRotPerNode  = 3
        nDispPerNode = 3
        nDOFPerNode  = 6
    else:
        raise NotImplementedError()

    # --- Attributing DOF to each elements
    NodesDOF={}
    ElemsDOF={}
    for e in g.Elements:
        ElemsDOF[e.ID] = [-1]*nDOFPerElem

    iPrev=0
    for iNode, n in enumerate(g.Nodes):
        # --- Distribute to joints iPrev + 1:6, or, iPrev + 1:(3+3m)
        elems = g.nodeIDs2Elements[n.ID] # elements connected to current node

        if int(n.data['Type']) == idJointCantilever:
            nRot=nRotPerNode
        else:
            nRot= nRotPerNode*len(elems) 
        NodesDOF[n.ID] = list(np.arange(0,3+nRot) + iPrev)
        # --- Distribute to elements
        for e in elems:
            nodeIDlist = g.elementIDs2NodeIDs[e.ID]
            if nodeIDlist.index(n.ID) ==0: # then current node is elem node 1
                iOff = 0
            else:                          # current node is elem node 2
                iOff = 6
            # Displacements
            ElemsDOF[e.ID][iOff:iOff+nDispPerNode] =  NodesDOF[n.ID][0:nDispPerNode]
            # Rotations
            if int(n.data['Type']) == idJointCantilever:
                ElemsDOF[e.ID][iOff+nDispPerNode:iOff+nDOFPerNode] =  NodesDOF[n.ID][nDispPerNode:nDOFPerNode]
            else:
                ElemsDOF[e.ID][iOff+nDispPerNode:iOff+nDOFPerNode] = NodesDOF[n.ID][nRotPerNode*iElem:nRotPerNode*(iElem+1)] # TODO verify

        iPrev = iPrev + len(NodesDOF[n.ID])

    # --- Safety check
    for e in g.Elements:
        if any(ElemsDOF[e.ID])<0:
            Exception('Implementation error in Distribute DOF, some element DOF were not allocated')

    # --- Store DOFs
    for n in g.Nodes:
        n.data['DOFs'] = NodesDOF[n.ID]
    for e in g.Elements:
        e.data['DOFs'] = ElemsDOF[e.ID]
    return g


if __name__=='__main__':
    np.set_printoptions(linewidth=500)
    mdl=SubDynModel()
    mdl.fromSummaryFile('../../data/Monopile/Pendulum.SD.sum.yaml')

