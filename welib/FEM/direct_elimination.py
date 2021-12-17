""" 
Tools to perform direct elimination, similar to SubDyn implementation

"""

import numpy as np
from welib.FEM.fem_elements import *
from welib.FEM.utils import rigidTransformationTwoPoints


def nDOF_c(self, RA):
    """ Returns number of DOF after constraint reduction (via the matrix T)
    TODO
    """
    nDOF_c = 0
    # Rigid assemblies contribution
    nDOF_c = nDOF_c + 6*len(RA)
    # Contribution from all the other joints
    for iNode, node in enumerate(self.Nodes):
        nodeID = node.ID
        #m = Init%NodesConnE(iNode,1) ! Col1: number of elements connected to this joint
        NodeType = node.data['Type']
        if   NodeType == idJointPin  :
            nDOF_c += 5 + 1*m
            print('   NodeID',nodeID, ' is a pin joint, number of members involved: ',m)
        elif NodeType == idJointUniversal :
            nDOF_c = nDOF_c + 4 + 2*m
            print('   NodeID',nodeID, ' is an universal joint, number of members involved: ',m)
        elif NodeType == idJointBall :
            nDOF_c = nDOF_c + 3 + 3*m
            print('   NodeID',nodeID, ' is a ball joint, number of members involved: ',m)
        elif NodeType == idJointCantilever:
            hasRigid, er = nodeHasRigidElem(self, node)
            if hasRigid:
                # This joint is involved in a rigid link assembly, we skip it (accounted for above)
                print('   NodeID',nodeID, ' is involved in a Rigid assembly')
            else:
                # That's a regular Cantilever joint
                nDOF_c = nDOF_c + 6
    return nDOF_c

def nodeHasRigidElem(self, node):
    elems = self.node2Elements(node)
    for e in elems:
        if e.data['TypeID']==idMemberRigid:
            return True, e
    return False, None

def rigidLinkAssemblies(self):
    """ Setup a list of rigid link assemblies (RA)
    !! Variables created by this routine:
    !! - RA(ia)= [e1,..,en]  list of elements forming each rigid link assembly "ia".
    !!                       Needed for BuildTMatrix
    !! - RAm1(e)=(RA^-1(e)= a) : for a given element give the index of a rigid assembly. 
    !!                       Needed for BuildTMatrix
    type(IList), dimension(:)   :: RA   !< RA(a) = [e1,..,en]  list of elements forming a rigid link assembly
    integer(IntKi), dimension(:):: RAm1 !< RA^-1(e) = a , for a given element give the index of a rigid assembly
    """
#         allocate(RAm1(1:Init%NElem)) ! NOTE: do not deallocate, this is an "output" of this function
#         RAm1(1:Init%NElem) = -1
# 
    # --- Establish a list of rigid link elements
    Er = [e for e in self.Elements if e.data['TypeID']==idMemberRigid ]
    EIDr = [e.ID for e in self.Elements if e.data['TypeID']==idMemberRigid ]
    print('   Rigid Elements ',EIDr)
    RA  = []
    while len(EIDr)>0:
        # Creating List Ea of elements of a given assembly
        EIDa =[]
        id0 = EIDr.pop()
        EIDa.append(id0) 
        e0 = self.getElement(id0)
        addNeighbors(self, e0, EIDr, EIDa)
        print('   Rigid Assembly, element IDs: ', EIDa)
        RA.append(EIDa) 
    return RA
   

def addNeighbors(self, e0, EIDr, EIDa) :
    """ The neighbor-elements of element e0 (that are found within the list Er) are added to the list Ea  
        - e0:  Given element
        - EIDr:  List of rigid elements ids available
        - EIDa:  List of elements ids in current rigid assembly
    """
    #print('----------------------------------------------------')
    #print('>>> Looking for neighbors of ',e0.ID, 'within',EIDr)
    if len(EIDr)==0:
        return EIDa

    EIDn =[] # List of neighbors of e0
    # Loop through all elements, setup list of e0-neighbors, add them to Ea, remove them from Er
    commonEIDs=[]
    for idk in EIDr:
        ek = self.getElement(idk)
        #print('Looking if ', e0.ID, 'and', idk, 'are directly connected.')
        if self.areElementsConnected(e0, ek):
            #print('   YES, ElementID',idk, 'is connected to elementID', e0.ID)
            commonEIDs.append(idk)
        #else:
        #    print('   NO')
    # Remove element from Er (a rigid element can belong to only one assembly)
    #print('EIDr is ',EIDr)
    for idk in commonEIDs:
        EIDr.remove(idk)
        #print('removing ',idk,'EIDr is ',EIDr)
        EIDn.append(idk) # adding to neighbors
        EIDa.append(idk) # adding to assembly
    #print('>>> The neighbors of ',e0.ID, 'are', EIDn)
    #print('')
    # Loop through neighbors and recursively add neighbors of neighbors
    if len(EIDr)>0:
        for idk in EIDn:
            ek = self.getElement(idk)
            EIDa = addNeighbors(self, ek, EIDr, EIDa)
    return EIDa

def RAElimination(self, RA):
    """
    !------------------------------------------------------------------------------------------------------
    !> Returns constraint matrix Tc for a rigid assembly (RA) formed by a set of elements. 
    !!   x_c = Tc.x_c_tilde  
    !! where x_c are all the DOF of the rigid assembly, and x_c_tilde are the 6 reduced DOF (leader DOF)
    """
    Elements = [self.getElement(eid) for eid in RA]
    # --- List of nodes stored first 
    #print('>>> Elements',Elements)
    Nodes = self.elements2nodes(Elements)
    INodesID = [n.ID for n in Nodes]
    print('   Nodes involved in assembly (unique)', INodesID)
    #--- Look for potential interface node
    NodesInterf = []
    for  iNodeID, (node,nodeID) in enumerate(zip(Nodes, INodesID)):
        if 'IBC' in node.data.keys():
            print('   Node',nodeID, ' is an interface node, selecting it for the rigid assembly')
            NodesInterf.append(nodeID)
    # --- Decide which node will be the main node of the rigid assembly
    if      (len(NodesInterf)==0):
        iiMainNode = 0 # By default we select the first node
    elif (len(NodesInterf)==1):
        # Finding the index of the interface node
        idMainNode = NodesInterf[0]
        iiMainNode = INodesID.index(idMainNode)
    else:
        raise Exception('Cannot have several interface nodes linked within a same rigid assembly')
    print('   Selecting node ID ',INodesID[iiMainNode], 'to be the main node for the rigid assembly')
    # --- Order list of joints with main node first (swapping iMainNode with INodes(1))
    iTmp  = INodesID[0]
    INodesID[0] = INodesID[iiMainNode]
    INodesID[iiMainNode] = iTmp
    print('   Nodes involved in assembly (after select):',INodesID)
    # --- Building Transformation matrix
    nNodes = len(INodesID)
    Tc = np.zeros((6*nNodes,6))
    # I6 for first node since it's the "leader"
    Tc[:6,:6]=np.eye(6)
    # Rigid transformation matrix for the other nodes 
    P1 = self.getNode(INodesID[0]).point # reference node coordinates
    for i,nID in enumerate(INodesID[1:]):
        Pi = self.getNode(nID).point # follower node coordinates
        T_r = rigidTransformationTwoPoints(P1, Pi)
        Tc[ (i+1)*6:(i+2)*6, :] = T_r
        #print('Rigid transformation from ref point',P1,' to ',Pi, T_r)
    return Tc, INodesID

def buildTMatrix(self, RA):
    """ 
     Build transformation matrix T, such that x= T.x~ where x~ is the reduced vector of DOF
     Variables set by this routine
    - p%NodesDOFred(iNode)=[list of DOF]: Created for each node, the list of DOF of this node in the 
            reduced system. 
            NOTE: follower nodes in rigid assembly have no DOFred (convention)
    - p%nDOF_red: number of DOF in reduced system (<= nDOF)
    - p%reduced: true if a reduction is needed, i.e. a T matrix is needed, and nDOF_red<nDOF
    
    Variables returned:
    - T_red: retuction matrix such that x= T_red.x~ where x~ is the reduced vector of DOF
    """
#    ! --- Misc inits
    nDOF_red = nDOF_c(self, RA)
    nDOF     = self.nDOF
    print('   Number of reduced DOF',nDOF_red, '/',nDOF)
    T_c = np.zeros((nDOF, nDOF_red)) 
    IRA = list(np.arange(len(RA)))
    # --- For each node:
    #  - create list of indices I      in the assembled vector of DOF
    #  - create list of indices Itilde in the reduced vector of DOF
    #  - increment iPrev by the number of DOF of Itilde
    iPrev =0 
    for iNode, node in enumerate(self.Nodes):
        idNodeSel = node.ID # Unless changed by Rigid assembly, using this index
        JType = node.data['Type']
        if JType == idJointCantilever: 
            hasRigid,er =  nodeHasRigidElem(self, node)
            if hasRigid:
                # --- The joint is involved in a rigid link assembly
                aID = -1
                for iRA, RA0 in enumerate(RA):
                    if er.ID in RA0:
                        aID=iRA
                        break
                if aID==-1:
                    raise Exception()
                    import pdb; pdb.set_trace()
                if aID not in IRA:
                    #print('NID',idNodeSel, 'SKIPPED, the RA',aID, 'has already been processed')
                    continue # We pass to the next joint, important so that:
                    #         - we don't increase iPrev
                    #         - we don't set Tc
                    #         - p%NodesDOFred is not set (assuming it has already been done)
                else:
                    # --- Proceessing the rigid assembly
                    # Returns TC and INodesID, do not change other variables
                    Tc, INodesID = RAElimination(self, RA[aID])
                    #print('Tc\n',Tc)
                    # The rigid assembly has been processed, delete index
                    IRA.remove(aID)
                    nj = len(INodesID) # Number of nodes in this rigid assembly
                    IDOFOld = []
                    for nID in INodesID:
                        IDOFOld += self.getNode(nID).data['DOFs']
                    # Storing DOF list for this RA (Note: same as NodesDOFred below, only for debug)
                    nc = Tc.shape[0] # Should be 6 
                    # --- Processing trigger for leader/follower Nodes
                    idNodeSel = INodesID[0]  # The first index returned is the leader of the assembly, we use this from now on
                    for nID in INodesID[1:]: # start at 2 because 1 is always the leader
                       # NEW: this node has no DOFs, so we set an empty list of DOFred for this node
                       node = self.getNode(nID)
                       node.data['DOFs_c'] = []
                       node.data['ID_link'] = idNodeSel
                       #print('Node', nID,' has no reduced DOF since its the follower of leader node ',idNodeSel)
            else:
                Tc      = np.eye(6)
                IDOFOld = node.data['DOFs']
        else:
            raise NotImplementedError('Rotational joints')
            # --- Ball/Pin/Universal joint
            # allocate(IDOFOld(1:len(p%NodesDOF(iNodeSel))))
            # IDOFOld(:) = p%NodesDOF(iNodeSel)%List(:)
            # phat = Init%Nodes(iNodeSel, iJointDir:iJointDir+2)
            # ! Return Tc, do not change other variable
            # call JointElimination(Init%NodesConnE(iNodeSel,:), JType, phat, p, Tc, ErrStat2, ErrMsg2); if(Failed()) return
        # Assemble T_red based on T_c
        nc = Tc.shape[1]
        node = self.getNode(idNodeSel) # NOTE: might be different from the one at the top of the loop
        node.data['DOFs_c'] = list( iPrev + np.arange(nc))
        # KEEP ME, VERY USEFUL
        #print('NID',idNodeSel,'I ',node.data['DOFs'])
        #print('NID',idNodeSel,'It',node.data['DOFs_c'])
        #print('NID',idNodeSel,'Ia',IDOFOld)
        T_c[np.ix_(IDOFOld, node.data['DOFs_c'])] = Tc
        iPrev = iPrev + nc
    #print('--- End of BuildTMatrix')
    #print('   - p%nDOF_red', nDOF_red)
    #print('   - p%NodesDOFred: (list of reduced DOF indices per node) ')
    #for node in self.Nodes:
    #    print('NID',node.ID, 'It', node.data['DOFs_c'])
    # --- Safety checks
    if (len(IRA)>0):
        raise Exception('Not all rigid assemblies were processed')
    if iPrev != nDOF_red :
        raise Exception('Inconsistency in number of reduced DOF')
    return T_c
