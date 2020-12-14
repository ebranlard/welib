import numpy as np
import pandas as pd

from welib.tools.clean_exceptions import *
from welib.system.eva import eig
from welib.FEM.reduction import CraigBampton
from welib.FEM.graph import Node as GraphNode
from welib.FEM.graph import Element as GraphElement
from welib.FEM.graph import NodeProperty
from welib.FEM.graph import GraphModel


class MaterialProperty(NodeProperty):
    def __init__(self):
        Property.__init__(self)
        pass

class FEMNode(GraphNode):
    def __init__(self, ID, x, y, z=0, Type=None, DOFs=[]):
        GraphNode.__init__(self, ID, x, y, z)
        self.DOFs = DOFs

    def __repr__(self):
        s='<Node{:4d}> x:{:7.2f} y:{:7.2f} z:{:7.2f}, DOFs: {}'.format(self.ID, self.x, self.y, self.z, self.DOFs)
        return s

class FEMElement(GraphElement):
    def __init__(self, ID, nodeIDs, nodes=None, properties=None):
        GraphElement.__init__(self, ID, nodeIDs, nodes, properties)
        self.Ce=[]
        self.Ke=[]
        self.Me=[]

    def __repr__(self):
        s='<Elem{:4d}> NodeIDs: {}'.format(self.ID, self.nodeIDs)
        if self.nodes is not None:
            s+=' l={:.2f}'.format(self.length)
        return s

class BeamElement(FEMElement):
    def __init__(self, ID, nodeIDs, nodes, properties=None):
        super(BeamElement,self).__init__(ID, nodeIDs, nodes=nodes, properties=properties)

class FEMModel(GraphModel):
    def __init__(self):
        GraphModel.__init__(self)
        self.MM       = None
        self.KK       = None
        self.DD       = None
        self.nDOF     = None

    def setFullMatrices(self,MM,KK,DD=None):
        self.MM=MM
        self.KK=KK
        if DD is not None:
            self.DD=DD

    def CraigBampton(self, Ileader, Ifollow=None, Ifixed=None):
        """ """
        if Ifixed is not None:
            M,K = self.applyFixBC()
        else:
            M,K = self.MM, self.KK

        return CraigBampton(M, K, Ileader, Ifollow=Ifollow)

    def DOF2Nodes(self):
        DOF2Nodes=np.zeros((self.nDOF,4),int)
        for iN,node in enumerate(self.Nodes):
            for iiDOF,iDOF in enumerate(node.DOFs):
                DOF2Nodes[iDOF,0] = iDOF
                DOF2Nodes[iDOF,1] = iN
                DOF2Nodes[iDOF,2] = len(node.DOFs)
                DOF2Nodes[iDOF,3] = iiDOF+1
        return DOF2Nodes

if __name__=='__main__':
    np.set_printoptions(linewidth=500)
    mdl=SubDynModel()
    mdl.fromSummaryFile('../../data/Monopile/Pendulum.SD.sum.yaml')

