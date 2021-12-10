"""
A FEM Model is a special kind of Graph
   Graphs have: Nodes, Elements, Node properties, and a "Model" storing it all




"""
import numpy as np
import pandas as pd


from welib.FEM.utils import DCM, rigidTransformationMatrix
from welib.FEM.fem_core import insertFixedBCinModes

from welib.tools.clean_exceptions import *
# from welib.FEM.graph import Node as GraphNode
# from welib.FEM.graph import Element as GraphElement
# from welib.FEM.graph import NodeProperty
# from welib.FEM.graph import GraphModel
from welib.weio.tools.graph import Node as GraphNode
from welib.weio.tools.graph import Element as GraphElement
from welib.weio.tools.graph import NodeProperty
from welib.weio.tools.graph import GraphModel 

# Following the convention of SubDyn
idJointCantilever = 1
idJointUniversal  = 2
idJointPin        = 3
idJointBall       = 4
idMemberBeam      = 1
idMemberCable     = 2
idMemberRigid     = 3

idDOF_Fixed    =  0 # Fixed DOF BC
idDOF_Internal = 10 # Free/Internal DOF 
idDOF_Leader   = 20 # Leader DOF


# class MaterialProperty(NodeProperty):
#     def __init__(self):
#         Property.__init__(self)
#         pass

class FEMNode(GraphNode):
    def __init__(self, ID, x, y, z=0, Type=idJointCantilever, DOFs=[], **kwargs):
        kwargs['Type'] = Type
        kwargs['DOFs'] = DOFs
        super(FEMNode,self).__init__(ID, x, y, z, **kwargs)

    def __repr__(self):
        s='<FNode{:4d}> x:{:7.2f} y:{:7.2f} z:{:7.2f}, {:}'.format(self.ID, self.x, self.y, self.z, self.data)
        return s

class FEMElement(GraphElement):
    def __init__(self, ID, nodeIDs, nodes=None, propset=None, propIDs=None, properties=None, **kwargs):
        super(FEMElement, self).__init__(ID, nodeIDs, nodes, propset, propIDs, properties, **kwargs)
        self.data['Type'] = 'Generic'

    #@property
    def Ke(self):
        raise NotImplementedError()

    #@property
    def Me(self):
        raise NotImplementedError()

    def __repr__(self):
        s='<FElem{:4d}> NodeIDs: {} {}'.format(self.ID, self.nodeIDs, self.data)
        if self.propIDs is not None:
            s+=' {'+'propIDs:{} propset:{}'.format(self.propIDs, self.propset)+'}'
        if self.nodes is not None:
            s+=' l={:.2f}'.format(self.length)
        return s

# --------------------------------------------------------------------------------}
# --- Beam-like elements
# --------------------------------------------------------------------------------{
class Beam3dElement(FEMElement):
    def __init__(self, ID, nodeIDs, nodes=None, propset=None, propIDs=None, properties=None, **kwargs):
        super(Beam3dElement,self).__init__(ID, nodeIDs, nodes, propset, propIDs, properties, **kwargs)
        self.data['Type'] = 'Beam3d'

    def __repr__(self):
        s='<{:s}Elem{:4d}> NodeIDs: {} {}'.format(self.data['Type'],self.ID, self.nodeIDs, self.data)
        if self.propIDs is not None:
            s+=' {'+'propIDs:{} propset:{}'.format(self.propIDs, self.propset)+'}'
        if self.nodes is not None:
            s+=' l={:.2f}'.format(self.length)
        return s

    @property
    def DCM(self):
        n1, n2 = self.nodes[0], self.nodes[1]
        return DCM((n1.x,n1.y,n1.z),(n2.x, n2.y, n2.z), main_axis='z')

    @property
    def t(self):
        try:
            n1, n2 = self.nodeProps[0], self.nodeProps[1] # Using SubDyn convention
        except:
            n1, n2 = self.nodes[0], self.nodes[1]
        t1, t2 = n1.data['t'], n2.data['t']
        return (t1+t2)/2

    @property
    def r1(self):
        try:
            n1, n2 = self.nodeProps[0], self.nodeProps[1] # Using SubDyn convention
        except:
            n1, n2 = self.nodes[0], self.nodes[1]
        D1, D2 = n1.data['D'], n2.data['D']
        return (D1+D2)/4

    @property
    def r2(self):
        try:
            n1, n2 = self.nodeProps[0], self.nodeProps[1] # Using SubDyn convention
        except:
            n1, n2 = self.nodes[0], self.nodes[1]
        D1, D2 = n1.data['D'], n2.data['D']
        t  = self.t
        if t==0:
            r2 = 0
        else:
            r2 = self.r1 - t
        return r2

    @property
    def area(self):
        return np.pi*(self.r1**2- self.r2**2)

    @property
    def inertias(self):
        Ixx = 0.25*np.pi*(self.r1**4-self.r2**4)
        Iyy = Ixx
        Jzz = 2.0*Ixx
        return Ixx, Iyy, Jzz

    @property
    def D(self): 
        return 2*self.r1

    @property
    def rho(self): 
        try:
            return self.nodeProps[0].data['rho'] # same convention as SubDyn returns density of first node
        except:
            return self.data['rho'] 
    @property
    def E(self): 
        try:
            return self.nodeProps[0].data['E'] # same convention as SubDyn returns density of first node
        except:
            return self.data['E']
    @property
    def G(self): 
        try:
            return self.nodeProps[0].data['G'] # same convention as SubDyn returns density of first node
        except:
            return self.data['G'] 

    @property
    def T0(self):
        return -9.990000E+36 # pretension, only for cable


    def Fe_g(e, g, main_axis='z', local=False):
        """ Element force due to gravity """
        from .timoshenko import timoshenko_Fe_g
        R = np.eye(3) if local else e.DCM
        return timoshenko_Fe_g(e.length, e.area, e.rho, g, R=R, main_axis=main_axis)

    def Fe_o(e, main_axis='z', local=False):
        """ Element force due to other sources (e.g. pretension cable) """
        return np.zeros(12)

# --------------------------------------------------------------------------------}
# --- Frame3D
# --------------------------------------------------------------------------------{
class Frame3dElement(Beam3dElement):
    def __init__(self, ID, nodeIDs, nodes=None, propset=None, propIDs=None, properties=None, **kwargs):
        super(Frame3dElement,self).__init__(ID, nodeIDs, nodes, propset, propIDs, properties, **kwargs)
        self.data['Type'] = 'Frame3d'

    @property
    def kappa(self): return 0  # shear coefficients are zero for Euler-Bernoulli

    def Ke(e, main_axis='z', local=False):
        from .frame3d import frame3d_KeMe # TODO TODO
        from .timoshenko import timoshenko_Ke
        I = e.inertias
        R = None if local else e.DCM
        return timoshenko_Ke(e.length, e.area, I[0], I[1], I[2], e.kappa,  e.E, e.G, shear=False, R=R, main_axis=main_axis) # NOTE: shear False for 

    def Me(e, main_axis='z', local=False):
        from .frame3d import frame3d_KeMe # TODO TODO
        from .timoshenko import timoshenko_Me
        I = e.inertias
        R = None if local else e.DCM
        return timoshenko_Me(e.length, e.area, I[0], I[1], I[2], e.rho, R=R, main_axis=main_axis)

# --------------------------------------------------------------------------------}
# --- Timoshenko 3D 
# --------------------------------------------------------------------------------{
class Timoshenko3dElement(Beam3dElement):
    def __init__(self, ID, nodeIDs, nodes=None, propset=None, propIDs=None, properties=None, **kwargs):
        super(Timoshenko3dElement,self).__init__(ID, nodeIDs, nodes, propset, propIDs, properties, **kwargs)
        self.data['Type'] = 'Timoshenko3d'

    @property
    def kappa(self): 
        try:
            n1, n2 = self.nodeProps[0], self.nodeProps[1] # Using SubDyn convention
        except:
            n1, n2 = self.nodes[0], self.nodes[1]
        D1 = n1.data['D']
        D2 = n2.data['D']
        t1 = n1.data['t']
        t2 = n2.data['t']
        r1 = (D1+D2)/4
        t  = (t1+t2)/2
        nu      = self.E/(2*self.G) - 1
        D_outer = 2 * r1              # average (outer) diameter
        D_inner = D_outer - 2*t       # remove 2x thickness to get inner diameter
        ratioSq = ( D_inner / D_outer)**2
        kappa =   ( 6 * (1 + nu) **2 * (1 + ratioSq)**2 )/( ( 1 + ratioSq )**2 * ( 7 + 14*nu + 8*nu**2 ) + 4 * ratioSq * ( 5 + 10*nu + 4 *nu**2 ) )
        return kappa

    def Ke(e, main_axis='z', local=False):
        from .timoshenko import timoshenko_Ke
        I = e.inertias
        R = None if local else e.DCM
        return timoshenko_Ke(e.length, e.area, I[0], I[1], I[2], e.kappa,  e.E, e.G, shear=True, R=R, main_axis=main_axis)

    def Me(e, main_axis='z', local=False):
        from .timoshenko import timoshenko_Me
        I = e.inertias
        R = None if local else e.DCM
        return timoshenko_Me(e.length, e.area, I[0], I[1], I[2], e.rho, R=R, main_axis=main_axis)

# --------------------------------------------------------------------------------}
# --- Cable3D 
# --------------------------------------------------------------------------------{
class Cable3dElement(Beam3dElement):
    def __init__(self, ID, nodeIDs, nodes=None, propset=None, propIDs=None, properties=None, **kwargs):
        super(Cable3dElement,self).__init__(ID, nodeIDs, nodes, propset, propIDs, properties, **kwargs)
        self.data['Type'] = 'Cable3d'

    # For compatbility with other beam-like structures
    @property
    def area(self): return 1
    @property
    def G(self): return -9.99e+36
    @property
    def inertias(self): return (-9.99e+36,-9.99e+36,-9.99e+36)
    @property
    def kappa(self): return -9.99e+36
    @property
    def E(self): return self.EA/self.area

    @property
    def EA(self):
        try:
            return self.nodeProps[0].data['EA']
        except:
            return self.data['EA']

    @property
    def L0(e): 
        """ rest length for which pretension would be zero"""
        Eps0 = e.T0/(e.EA)
        return  e.length/(1+Eps0)

    @property
    def T0(self): 
        try:
            return self.nodeProps[0].data['T0'] # same convention as SubDyn returns density of first node
        except:
            return self.data['T0'] 

    def Ke(e, main_axis='z', local=False):
        from .cable import cable3d_Ke
        I = e.inertias
        R = None if local else e.DCM


        return cable3d_Ke(e.length, e.area, e.E, e.T0, R=R, main_axis=main_axis)

    def Me(e, main_axis='z', local=False):
        from .cable import cable3d_Me
        I = e.inertias
        R = None if local else e.DCM
        return cable3d_Me(e.L0, e.area, e.rho, R=R, main_axis=main_axis) # NOTE: we use L0 for the mass

    def Fe_o(e, main_axis='z', local=False):
        from .cable import cable3d_Fe_T0
        R = None if local else e.DCM
        return cable3d_Fe_T0(e.T0, R=R, main_axis=main_axis)


# --------------------------------------------------------------------------------}
# --- Main FEM model class  
# --------------------------------------------------------------------------------{
class FEMModel(GraphModel):
    @classmethod
    def from_graph(cls, graph, ndiv=None, mainElementType='frame3d', TP=(0,0,0), gravity=9.81):
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
                    g.Elements[i] = Frame3dElement(e.ID, e.nodeIDs, None, e.propset, e.propIDs, None, **e.data)
                elif mainElementType=='timoshenko':
                    g.Elements[i] = Timoshenko3dElement(e.ID, e.nodeIDs, None, e.propset, e.propIDs, None, **e.data)
                else:
                    raise NotImplementedError()
            elif e.data['Type']=='Cable':
                if mainElementType in ['frame3d', 'timoshenko']:
                    g.Elements[i] = Cable3dElement(e.ID, e.nodeIDs, None, e.propset, e.propIDs, None, **e.data)
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
                ElemPropertySets=g.ElemPropertySets, MiscPropertySets=g.MiscPropertySets, TP=TP, gravity=gravity)
        #print(self)


        return self

    def __init__(self, Elements=[], Nodes=[], NodePropertySets=dict(), ElemPropertySets=dict(), MiscPropertySets=dict(),
            mainElementType='frame3d', gravity=9.81, main_axis='z', TP=(0,0,0)): # FEM specific
        """ 

        """

        GraphModel.__init__(self, Elements=Elements, Nodes=Nodes, NodePropertySets=NodePropertySets,
                ElemPropertySets=ElemPropertySets, MiscPropertySets=MiscPropertySets)

        self.TP              = TP
        self.gravity         = gravity
        self.main_axis       = main_axis
        self.mainElementType = mainElementType
        # Main data generated
        self.MM       = None
        self.KK       = None
        self.DD       = None
        self.FF       = None

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

        self.KK = KK
        self.MM = MM 
        self.FF = FF
        return self

    # --------------------------------------------------------------------------------}
    # --- Direct elimination of Internal Constraints (Rigid links and rotational joints)
    # --------------------------------------------------------------------------------{
    def applyInternalConstraints(self):
        """ 
        Apply internal constraints such as rigid links and rotational joints
        using a direct elminination technique

        - Tc: reduction matrix such that x= Tc.xc where xc is the reduced vector of DOF

        """
        rotJoints  = [n.ID for n in self.Nodes    if n.data['Type']!=idJointCantilever]
        rigidLinks = [e.ID for e in self.Elements if e.data['Type']==idMemberRigid    ]
        if len(rotJoints)>0 or len(rigidLinks)>0:
            print('Number of Rotational joints:',len(rotJoints))
            print('Number of Rigid Links      :',len(rigidLinks))
            from .direct_elimination import nDOF_c, buildTMatrix, rigidLinkAssemblies
            raise NotImplementedError('Direct elimination')

        else:
            self.T_c = np.eye(self.MM.shape[0])
            # Store new DOF indices
            for n in self.Nodes:
                n.data['DOFs_c'] = list(n.data['DOFs'])
            for e in self.Elements:
                e.data['DOFs_c'] =[]
                for n in e.nodes:
                    e.data['DOFs_c'] +=n.data['DOFs_c']

        self.MM_c = (self.T_c.T).dot(self.MM).dot(self.T_c)
        self.KK_c = (self.T_c.T).dot(self.KK).dot(self.T_c)

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
        # --- Count nodes per types
        nNodes    = len(self.Nodes)
        nNodes_I  = len(self.interfaceNodes)
        nNodes_C  = len(self.reactionNodes)
        nNodes_L  = len(self.internalNodes)

        # --- Partition Nodes:  Nodes_L = IAll - NodesR
        Nodes_I = [n.ID for n in self.interfaceNodes]
        Nodes_C = [n.ID for n in self.reactionNodes]
        Nodes_R = Nodes_I + Nodes_C
        Nodes_L = [n.ID for n in self.Nodes if n.ID not in Nodes_R]

        # --- Count DOFs - NOTE: we count node by node
        nDOF___  = sum([len(n.data['DOFs_c'])                                 for n in self.Nodes])
        # Interface DOFs
        nDOFI__  = sum([len(n.data['DOFs_c'])              for n in self.interfaceNodes])
        nDOFI_B = sum([sum(np.array(n.data['IBC'])==idDOF_Leader)   for n in self.interfaceNodes])
        nDOFI_F  = sum([sum(np.array(n.data['IBC'])==idDOF_Fixed )   for n in self.interfaceNodes])
        if nDOFI__!=nDOFI_B+nDOFI_F: raise Exception('Wrong distribution of interface DOFs')
        # DOFs of reaction nodes
        nDOFC__ = sum([len(n.data['DOFs_c'])              for n in self.reactionNodes]) 
        nDOFC_B = sum([sum(np.array(n.data['RBC'])==idDOF_Leader)   for n in self.reactionNodes])
        nDOFC_F = sum([sum(np.array(n.data['RBC'])==idDOF_Fixed)    for n in self.reactionNodes])
        nDOFC_L = sum([sum(np.array(n.data['RBC'])==idDOF_Internal) for n in self.reactionNodes])
        if nDOFC__!=nDOFC_B+nDOFC_F+nDOFC_L: raise Exception('Wrong distribution of reaction DOFs')
        # DOFs of reaction + interface nodes
        nDOFR__ = nDOFI__ + nDOFC__ # Total number, used to be called "nDOFR"
        # DOFs of internal nodes
        nDOFL_L  = sum([len(n.data['DOFs_c']) for n in self.internalNodes])
        if nDOFL_L!=nDOF___-nDOFR__: raise Exception('Wrong distribution of internal DOF')
        # Total number of DOFs in each category:
        nDOF__B = nDOFC_B + nDOFI_B
        nDOF__F = nDOFC_F + nDOFI_F          
        nDOF__L = nDOFC_L           + nDOFL_L 
        self.nDOFR__ = nDOFR__
        self.nDOF__B = nDOF__B

        # --- Distibutes the I, L, C nodal DOFs into  B, F, L sub-categories 
        # NOTE: order is importatn for compatibility with SubDyn
        IDI__ = []
        IDI_B = []
        IDI_F = []
        for n in self.interfaceNodes:
            IDI__ += n.data['DOFs_c'] # NOTE: respects order
            IDI_B += [dof for i,dof in enumerate(n.data['DOFs_c']) if n.data['IBC'][i]==idDOF_Leader]
            IDI_F += [dof for i,dof in enumerate(n.data['DOFs_c']) if n.data['IBC'][i]==idDOF_Fixed ]
        IDI__ = IDI_B+IDI_F
        IDC__ = []
        IDC_B = []
        IDC_L = []
        IDC_F = []
        for n in self.reactionNodes:
            IDC__ += n.data['DOFs_c'] # NOTE: respects order
            IDC_B += [dof for i,dof in enumerate(n.data['DOFs_c']) if n.data['RBC'][i]==idDOF_Leader  ]
            IDC_L += [dof for i,dof in enumerate(n.data['DOFs_c']) if n.data['RBC'][i]==idDOF_Internal]
            IDC_F += [dof for i,dof in enumerate(n.data['DOFs_c']) if n.data['RBC'][i]==idDOF_Fixed   ]
        IDR=IDC__+IDI__
        IDL_L = []
        for n in self.internalNodes:
            IDL_L += n.data['DOFs_c']

        # --- Total indices per partition B, F, L
        self.DOF_Leader   =         IDC_B + IDI_B  # boundary/retained/leader DOFs
        self.DOF_Fixed    =         IDC_F + IDC_B  # Fixed DOFs
        self.DOF_Follower = IDL_L + IDC_L          # internal DOFs

        self.DOF_Boundary   = IDR    # I+C Boundary nodes for rigid body equivalent
        self.DOF_Internal   = IDL_L  # L   Internal nodes for rigid body equivalent
        self.DOF_Interface  = IDI__  # I   Interface


        ID_ALL = self.DOF_Leader + self.DOF_Fixed + self.DOF_Follower
        for i in np.arange(nDOF___):
            if i not in ID_ALL:
                raise Exception('DOF {} not found in DOF list after partition')

        # Storing variables similar to SubDyn
        SD_IO_Vars={}
        SD_IO_Vars['nDOF___']=nDOF___;
        SD_IO_Vars['nDOFI__']=nDOFI__; SD_IO_Vars['nDOFI_B']=nDOFI_B; SD_IO_Vars['nDOFI_F']=nDOFI_F;
        SD_IO_Vars['nDOFC__']=nDOFC__; SD_IO_Vars['nDOFC_B']=nDOFC_B; SD_IO_Vars['nDOFC_F']=nDOFC_F; SD_IO_Vars['nDOFC_L']=nDOFC_L;
        SD_IO_Vars['nDOFR__']=nDOFR__; SD_IO_Vars['nDOFL_L']=nDOFL_L;
        SD_IO_Vars['nDOF__B']=nDOF__B; SD_IO_Vars['nDOF__F']=nDOF__F; SD_IO_Vars['nDOF__L']=nDOF__L;
        SD_IO_Vars['IDC__']=IDC__;
        SD_IO_Vars['IDC_B']=IDC_B;
        SD_IO_Vars['IDC_F']=IDC_F;
        SD_IO_Vars['IDC_L']=IDC_L;
        SD_IO_Vars['IDI__']=IDI__;
        SD_IO_Vars['IDI_B']=IDI_B;
        SD_IO_Vars['IDI_F']=IDI_F;
        SD_IO_Vars['ID__B']=self.DOF_Leader
        SD_IO_Vars['ID__F']=self.DOF_Fixed
        SD_IO_Vars['ID__L']=self.DOF_Follower
        self.SD_IO_Vars=SD_IO_Vars



    # --------------------------------------------------------------------------------}
    # ---  
    # --------------------------------------------------------------------------------{
    def applyFixedBC(self, IFixed=None):
        """ 
        Apply boundary conditions. (Fixed boundary conditions only)
        """
        # NOTE: we use the matrices where internal constraints have been eliminated 
        MM = self.MM_c
        KK = self.KK_c
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

        self.MM_cr=Mr
        self.KK_cr=Kr
        self.T_r =Tr
        #
        self.DOF_Leader_r   = [IFull2BC[i] for i in self.DOF_Leader]
        self.DOF_Follower_r = [IFull2BC[i] for i in self.DOF_Follower]
        self.DOF_Fixed_r    = [IFull2BC[i] for i in self.DOF_Fixed]

        return Mr, Kr, Tr, IFull2BC, IBC2Full


    def eig(self, normQ='byMax'):
        from welib.system.eva import eig
        KK = self.KK_cr
        MM = self.MM_cr

        # --- Compute modes and frequencies
        [Q, freq]= eig(KK, MM, freq_out=True, normQ=normQ)


        self.Q_cr    = Q
        self.freq_cr = freq

        self.Q_c = insertFixedBCinModes(Q, self.T_r)


    # --------------------------------------------------------------------------------}
    # --- Reference point (TP for now)
    # --------------------------------------------------------------------------------{
    @property
    def TI(self):
        return rigidTransformationMatrix(self.DOF_Interface, self.TP, self.DOF_c2Nodes, self.points)
    # --------------------------------------------------------------------------------}
    # --- General FEM Utils
    # --------------------------------------------------------------------------------{
    def setFullMatrices(self,MM,KK,DD=None):
        self.MM=MM
        self.KK=KK
        if DD is not None:
            self.DD=DD

    def CraigBampton(self, nModesCB=None, zeta=None, BC_before_CB=True):
        """
        Perform Craig-Bampton reduction

        nModesCB: number of CB modes to retain
        zeta :  damping ratios for CB modes 
        BC_before_CB: if true, using the matrices where the fixed BC have been applied 
        """
        from welib.FEM.reduction import CraigBampton
        if BC_before_CB:
            M, K = self.MM_cr, self.KK_cr
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
            if zeta is not None:
                if not hasattr(zeta,'__len__'):
                    zeta = [zeta]*nModesCB

                self.CC_MM_CB = 2*np.array(zeta) * omega_CB
                # TODO CC_CB combination of Guyan + CB damping matrix
            self.nModesCB= nModesCB
        else:
            raise NotImplementedError()
        #if Ifixed is not None:
        #    M,K = self.applyFixBC()
        #else:
        #    M,K = self.MM, self.KK

    def rigidBodyEquivalent(self):
        """ 
        Compute Rigid body equivalent
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
        self.M_G  = rigidBodyMassMatrixAtP(mass, J_G, (0,0,0))
        self.M_TP = rigidBodyMassMatrixAtP(mass, J_G, -np.array(self.TP)+np.array(Ref2COG))

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
        Mr, Kr, Phi_G, Phi_CB, f_G, f_CB, I1, I2 = CraigBampton(self.MM_c, self.KK_c, Ileader=Ileader, Ifollow=Ifollow, nModesCB=0)
        #! --- Insert SSI from Mass and stiffness matrix again
        #CALL InsertSoilMatrices(Init%M, Init%K, p%NodesDOFred, Init, p, ErrStat2, ErrMsg2, Substract=.False.); if(Failed()) return
        return Mr[np.ix_(I1,I1)] # return Guyan only

    
    # --------------------------------------------------------------------------------}
    # --- IO 
    # --------------------------------------------------------------------------------{
    def toYAML(self, filename):
        """ returns a yaml file similar to the YAML file of SubDyn"""
        from .fem_model_io import toYAML
        toYAML(self, filename)

    def nodeID_py2SD(self, nodeID):
        if hasattr(nodeID,'__len__'):
            return [self.Nodes.index(self.getNode(n))+1 for n in nodeID]
        else:
            return self.Nodes.index(self.getNode(nodeID))+1

    def elemID_py2SD(self, elemID):
        #e=self.getElement(elemID)
        for ie,e in enumerate(self.Elements):
            if e.ID==elemID:
                return ie+1
    def elemType_py2SD(self, elemType):
        return {'Beam3d':idMemberBeam, 'Beam':idMemberBeam, 'Frame3d':idMemberBeam,
                'Timoshenko3d':idMemberBeam,
                'Cable3d':idMemberCable, 'Cable':idMemberCable,
                'Rigid':idMemberRigid}[elemType]

    def propID_py2SD(self, propID, propset):
        prop = self.NodePropertySets[propset]
        for ip, p in enumerate(prop):
            if p.ID == propID:
                return ip+1
        



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

