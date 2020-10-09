import numpy as np
try:
    from pybra.clean_exceptions import *
except:
    pass

from mech_system.eva import eig


try:
    from .reduction import CraigBampton
except:
    from reduction import CraigBampton


class Node():
    def __init__(self, ID, x, y, z=0, Type=None, DOFs=[]):
        self.ID   = int(ID)
        self.x    = x
        self.y    = y
        self.z    = z
        self.DOFs = DOFs

    def __repr__(self):
        s='<Node{:4d}> x:{:7.2f} y:{:7.2f} z:{:7.2f}, DOFs: {}'.format(self.ID, self.x, self.y, self.z, self.DOFs)
        return s

class MaterialProperties():
    def __init__(self):
        pass

class Element():
    def __init__(self, ID, nodeIDs, nodes=None, properties=None):
        self.nodeIDs=nodeIDs
        self.nodes=nodes
        self.prop=properties
        self.Ce=[]
        self.Me=[]
        self.ID=ID

    @property
    def length(self):
        n1=self.nodes[0]
        n2=self.nodes[1]
        return np.sqrt((n1.x-n2.x)**2+(n1.y-n2.y)**2+(n1.z-n2.z)**2)

    def __repr__(self):
        s='<Elem{:4d}> NodeIDs: {}'.format(self.ID, self.nodeIDs)
        if self.nodes is not None:
            s+=' l={:.2f}'.format(self.length)
        return s

class BeamElement(Element):
    def __init__(self, ID, nodeIDs, nodes, properties=None):
        super(BeamElement,self).__init__(ID, nodeIDs, nodes=nodes, properties=properties)

class FEMModel():
    def __init__(self):
        self.Elements=[]
        self.Nodes=[]
        self.MM=None
        self.KK=None
        self.DD=None
        self.nDOF=None

    def addNode(self,node):
        if node.ID!=len(self.Nodes):
            raise Exception('For now, assuming that Node ID is same as index')
        self.Nodes.append(node)

    def addElement(self,elem):
        if elem.ID!=len(self.Elements):
            raise Exception('For now, assuming that elem ID ({}) is same as index ({})'.format(elem.ID,len(self.Elements)))
        self.Elements.append(elem)

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

    def extent(self):
        xmax=np.max([node.x for node in self.Nodes])
        ymax=np.max([node.y for node in self.Nodes])
        zmax=np.max([node.z for node in self.Nodes])
        xmin=np.min([node.x for node in self.Nodes])
        ymin=np.min([node.y for node in self.Nodes])
        zmin=np.min([node.z for node in self.Nodes])
        return [xmin,ymin,zmin],[xmax,ymax,zmax],[xmax-xmin,ymax-ymin,zmax-zmin]

    def maxDimension(self):
        _,_,D=self.extent()
        return np.max(D)




class SubDynModel(FEMModel):


    def fromSummaryFile(self, filename):
        # --- Read summary file
        import yaml
        with open(filename, 'r', errors="surrogateescape") as fid:
            data=yaml.load(fid)
            for k,v in data.items():
                if isinstance(v,list):
                    data[k]=np.array(v)
        DOF2Nodes = data['DOF2Nodes']
        self.nDOF  = data['nDOF_red']
        PhiB      = data['PhiM']
        PhiR      = data['PhiR']
        MM        = data['M']
        KK        = data['K']
        Nodes     = data['Nodes']
        Elements  = data['Elements']
        DOFs= np.arange(self.nDOF)
        # Reindexing with -1
        DOF_L     = data['DOF___L'].ravel()-1 # internal DOFs
        DOF_B     = data['DOF___B'].ravel()-1 # internal
        if DOF2Nodes.shape[1]==3:
            DOF2Nodes=np.column_stack((DOFs,DOF2Nodes))
        else:
            DOF2Nodes[:,0]-=1
        DOF2Nodes[:,1]-=1
        Elements[:,0]-=1
        Elements[:,1]-=1
        Elements[:,2]-=1

        # Nodes and DOFs
        for iNode,N in enumerate(Nodes):
            if len(N)==9: # Temporary fix
                N[4]=np.float(N[4].split()[0])
                N=N.astype(np.float32)
            ID = int(N[0])-1
            nodeDOFs=DOF2Nodes[(DOF2Nodes[:,1]==ID),0]
            node = Node(ID=ID, x=N[1], y=N[2], z=N[3], Type=int(N[4]), DOFs=nodeDOFs)
            print(node)
            self.addNode(node)

        # Elements
        for ie,E in enumerate(Elements):
            nodeIDs=[int(E[1]),int(E[2])]
            N1 = self.Nodes[nodeIDs[0]]
            N2 = self.Nodes[nodeIDs[1]]
            elem= BeamElement(int(E[0]), nodeIDs, nodes=[N1, N2], properties=None)
            print(elem)
            self.addElement(elem)

        # Full system matrices
        self.setFullMatrices(MM,KK)
        
        print(self.extent())
        print(self.maxDimension())

        DOF_K = np.concatenate((DOF_B,DOF_L))
        maxDisp = self.maxDimension()*0.1

        # CB modes
        Phi_CB = np.vstack((np.zeros((len(DOF_B),PhiB.shape[1])),PhiB))
        dispCB, posCB, INodes = self.NodesDisp(DOF_K, Phi_CB, maxDisp=maxDisp, sortDim=2)

        # Guyan modes
        Phi_Guyan = np.vstack((np.eye(len(DOF_B)),PhiR))
        dispGy, posGy, INodesGy = self.NodesDisp(DOF_K, Phi_Guyan, maxDisp=maxDisp, sortDim=2)

        import pandas as pd
        def toDF(pos,disp):
            columns=['z_[m]','x_[m]','y_[m]']
            dataZXY=np.column_stack((pos[:,2],pos[:,0],pos[:,1]))
            disptot=disp.copy()
            for ishape in np.arange(disp.shape[2]):
                disptot[:,:,ishape]= pos + disp[:,:,ishape]
                sMode='Mode{:d}'.format(ishape+1)
                columns+=[sMode+'x_[m]',sMode+'y_[m]',sMode+'z_[m]']
            disptot= np.moveaxis(disptot,2,1).reshape(disptot.shape[0],disptot.shape[1]*disptot.shape[2])
            disp   = np.moveaxis(disp,2,1).reshape(disp.shape[0],disp.shape[1]*disp.shape[2])
            print(columns)
            data = np.column_stack((dataZXY,disptot))
            df= pd.DataFrame(data = data ,columns = columns)

            disp[np.isnan(disp)]=0
            data = np.column_stack((dataZXY,disp))
            dfDisp= pd.DataFrame(data = data ,columns = columns)
            # remove zero 
            df2= df.loc[:, (dfDisp !=0).any(axis=0)]
            return df

        df=toDF(posGy, dispGy)


        import pdb; pdb.set_trace()


        self.toDataFrame()

    def toDataFrame(self):
        dfs={}

#         Mr, Kr, Phi_G, Phi_CB, f_CB = self.CraigBampton(Ileader=DOF_B, Ifollow=DOF_L)
#         print(data['PhiL'])

    def NodesDisp(self, IDOF, UDOF, maxDisp=None, sortDim=None):
        DOF2Nodes = self.DOF2Nodes()
        INodes = list(np.sort(np.unique(DOF2Nodes[IDOF,1]))) # NOTE: sorted
        print(INodes)
        nShapes = UDOF.shape[1]
        disp=np.empty((len(INodes),3,nShapes))
        disp.fill(np.nan)
        pos=np.empty((len(INodes),3))
        pos.fill(np.nan)
        for i,iDOF in enumerate(IDOF):
            iNode       = DOF2Nodes[iDOF,1]
            nDOFPerNode = DOF2Nodes[iDOF,2]
            nodeDOF     = DOF2Nodes[iDOF,3]
            node        = self.Nodes[iNode]
            iiNode      = INodes.index(iNode)
            if nodeDOF<=3:
                pos[iiNode, 0]=node.x
                pos[iiNode, 1]=node.y
                pos[iiNode, 2]=node.z
                for iShape in np.arange(nShapes):
                    disp[iiNode, nodeDOF-1, iShape] = UDOF[i, iShape]
        # Scaling 
        if maxDisp is not None:
            for iShape in np.arange(nShapes):
                mD=np.nanmax(np.abs(disp[:, :, iShape]))
                if mD>1e-5:
                    disp[:, :, iShape] *= maxDisp/mD

        # Sorting according to a dimension
        if sortDim is not None: 
            I=np.argsort(pos[:,sortDim])
            INodes = np.array(INodes)[I]
            disp   = disp[I,:,:]
            pos    = pos[I,:]

        return disp, pos, INodes


if __name__=='__main__':
    np.set_printoptions(linewidth=500)
    mdl=SubDynModel()
    mdl.fromSummaryFile('../Pendulum.SD.sum.yaml')



