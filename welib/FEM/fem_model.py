"""
A FEM Model is a special kind of Graph
   Graphs have: Nodes, Elements, Node properties, and a "Model" storing it all

"""
import numpy as np
import pandas as pd


from welib.FEM.utils import DCM, rigidTransformationMatrix
from welib.FEM.fem_elements import *   # Elements used
from welib.FEM.fem_core import insertFixedBCinModes

from welib.tools.clean_exceptions import *
# from welib.FEM.graph import Node as GraphNode
# from welib.FEM.graph import Element as GraphElement
# from welib.FEM.graph import NodeProperty
# from welib.FEM.graph import GraphModel
from welib.weio.tools.graph import NodeProperty
from welib.weio.tools.graph import GraphModel 



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
        import copy
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
            else:
                raise NotImplementedError()
        #print('----------------------- GRAPH')
        #print(graph)
        #print('----------------------- G')
        #print(g)
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

    def __init__(self, Elements=[], Nodes=[], NodePropertySets=dict(), ElemPropertySets=dict(), MiscPropertySets=dict(),
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
        self.freq     = None # Frequency of the Full FEM system after Boundary conditions and internal constraints
        self.Q        = None # Modes of the full FEM system after Boundary conditions and internal constraints
                             # NOTE: augmented to include the boundary condition in them

        # internal data
        self._nDOF    = None

        # 
        if 'DOFs' not in self.Nodes[0].data.keys():
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
    def nDOF_c(self):
        self._nDOF_c = np.max(self.Nodes[-1].data['DOFs_c'])+1
        return self._nDOF


    @property
    def DOF2Nodes(self):
        DOF2Nodes=np.zeros((self.nDOF,4),int)
        for iN,node in enumerate(self.Nodes):
            for iiDOF,iDOF in enumerate(node.data['DOFs']):
                DOF2Nodes[iDOF,0] = iDOF
                DOF2Nodes[iDOF,1] = iN
                DOF2Nodes[iDOF,2] = len(node.data['DOFs'])
                DOF2Nodes[iDOF,3] = iiDOF+1
        return DOF2Nodes

    @property
    def DOF_c2Nodes(self):
        DOF2Nodes=np.zeros((self.nDOF_c,4),int)
        for iN,node in enumerate(self.Nodes):
            for iiDOF,iDOF in enumerate(node.data['DOFs_c']):
                DOF2Nodes[iDOF,0] = iDOF
                DOF2Nodes[iDOF,1] = iN
                DOF2Nodes[iDOF,2] = len(node.data['DOFs_c'])
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
        return self

    # --------------------------------------------------------------------------------}
    # --- Direct elimination of Internal Constraints (Rigid links and rotational joints)
    # --------------------------------------------------------------------------------{
    def applyInternalConstraints(self):
        """ 
        Apply internal constraints such as rigid links and rotational joints
        using a direct elminination technique

        - Tc: reduction matrix such that x_init= Tc.x
              where x is the reduced vector of DOF and x_init is the intial vector (more DOF)

        """
        rotJoints  = [n.ID for n in self.Nodes    if n.data['Type']!=idJointCantilever]
        rigidLinks = [e.ID for e in self.Elements if e.data['Type']==idMemberRigid    ]
        if len(rotJoints)>0 or len(rigidLinks)>0:
            print('Number of Rotational joints:',len(rotJoints))
            print('Number of Rigid Links      :',len(rigidLinks))
            from .direct_elimination import nDOF_c, buildTMatrix, rigidLinkAssemblies
            raise NotImplementedError('Direct elimination')

        else:
            self.T_c = np.eye(self.MM_init.shape[0])
            # Store new DOF indices
            for n in self.Nodes:
                n.data['DOFs_c'] = list(n.data['DOFs'])
            for e in self.Elements:
                e.data['DOFs_c'] =[]
                for n in e.nodes:
                    e.data['DOFs_c'] +=n.data['DOFs_c']

        self.MM = (self.T_c.T).dot(self.MM_init).dot(self.T_c)
        self.KK = (self.T_c.T).dot(self.KK_init).dot(self.T_c)

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
        Partition DOFs into leader/follower and fixed, following the order convetion of SubDyn.
        
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
        self.DOF_Leader     =         IDC_B + IDI_B  # boundary/retained/leader DOFs
        self.DOF_Fixed      =         IDC_F + IDC_B  # Fixed DOFs
        self.DOF_Follower   = IDL_L + IDC_L          # internal DOFs
        self.DOF_Boundary   = IDR    # I+C Boundary nodes for rigid body equivalent
        self.DOF_Internal   = IDL_L  # L   Internal nodes for rigid body equivalent
        self.DOF_Interface  = IDI__  # I   Interface
        ID_ALL = self.DOF_Leader + self.DOF_Fixed + self.DOF_Follower
        for i in np.arange(nDOF___):
            if i not in ID_ALL:
                raise Exception('DOF {} not found in DOF list after partition'.format(i))

    # --------------------------------------------------------------------------------}
    # ---  
    # --------------------------------------------------------------------------------{
    def applyFixedBC(self, IFixed=None):
        """ 
        Apply boundary conditions. (Fixed boundary conditions only)
        """
        # NOTE: we use the matrices where internal constraints have been eliminated 
        MM = self.MM
        KK = self.KK
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
                IFixed += [n.data['DOFs'][ii] for ii,i in enumerate(I) if int(i)==idDOF_Fixed]

        Tr = np.delete(Tr, IFixed, axis=1) # removing columns

        Mr = (Tr.T).dot(MM).dot(Tr)
        Kr = (Tr.T).dot(KK).dot(Tr)

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
        self.T_BC  = Tr
        #
        self.DOF_Leader_r   = [IFull2BC[i] for i in self.DOF_Leader]
        self.DOF_Follower_r = [IFull2BC[i] for i in self.DOF_Follower]
        self.DOF_Fixed_r    = [IFull2BC[i] for i in self.DOF_Fixed]

        return Mr, Kr, Tr, IFull2BC, IBC2Full


    def eig(self, normQ='byMax'):
        from welib.system.eva import eig
        KK = self.KK_BC
        MM = self.MM_BC

        # --- Compute modes and frequencies
        [Q, freq]= eig(KK, MM, freq_out=True, normQ=normQ)

        self.freq = freq
        self.Q   = insertFixedBCinModes(Q, self.T_BC)


    # --------------------------------------------------------------------------------}
    # --- Reference point 
    # --------------------------------------------------------------------------------{
    @property
    def T_refPoint(self):
        """ Rigid body transformation matrix from interface DOFs to refpoint"""
        if self.refPoint is None: 
            raise Exception('Cannot compute T_refPoint, refPoint is None')
        return rigidTransformationMatrix(self.DOF_Interface, self.refPoint, self.DOF_c2Nodes, self.points)
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
            Ileader, Ifollow = self.DOF_Leader_r, self.DOF_Follower_r
            if nModesCB is None:
                nModesCB=M.shape[0] - len(Ileader)
            # NOTE: we return all CB modes at first
            Mr, Kr, Phi_G, Phi_CB, f_G, f_CB, I1, I2 = CraigBampton(M, K, Ileader=Ileader, Ifollow=Ifollow, nModesCB=None)
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
        TIR= rigidTransformationMatrix(self.DOF_Boundary, (0,0,0), self.DOF_c2Nodes, self.points)
        # Compute Rigid body mass matrix (without Soil, and using both Interface and Reactions nodes as leader DOF)
        if self.nDOFR__!=self.nDOF__B:
            print('>>> Need to do Guyan Rigid body Mass, TODO TODO TODO')
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
        Ileader = self.DOF_Boundary
        Ifollow = self.DOF_Internal
        Mr, Kr, Phi_G, Phi_CB, f_G, f_CB, I1, I2 = CraigBampton(self.MM, self.KK, Ileader=Ileader, Ifollow=Ifollow, nModesCB=0)
        #! --- Insert SSI from Mass and stiffness matrix again
        #CALL InsertSoilMatrices(Init%M, Init%K, p%NodesDOFred, Init, p, ErrStat2, ErrMsg2, Substract=.False.); if(Failed()) return
        return Mr[np.ix_(I1,I1)] # return Guyan only

    
    # --------------------------------------------------------------------------------}
    # --- IO 
    # --------------------------------------------------------------------------------{


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

    if mainElementType in ['frame3d', 'timoshenko']:
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

