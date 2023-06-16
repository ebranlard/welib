""" 
Convert a SubDyn summary file to a FEM model

"""
import numpy as np
import pandas as pd

from welib.FEM.graph import Node as GraphNode
from welib.FEM.graph import Element as GraphElement
from welib.FEM.graph import NodeProperty
from welib.FEM.graph import GraphModel

from welib.FEM.fem_model import FEMModel


# --------------------------------------------------------------------------------}
# ---  
# --------------------------------------------------------------------------------{
class SubDynModel(FEMModel):

    def fromSummaryFile(self, filename):
        # --- Read summary file
        from welib.weio.fast_summary_file import FASTSummaryFile
        data = FASTSummaryFile(filename) 
        #import yaml
        #print('Reading yaml...', end='')
        #with open(filename, 'r', errors="surrogateescape") as fid:
        #    data=yaml.load(fid, Loader=yaml.SafeLoader)
        #print('Done')
        #print(data.keys())
        #for k,v in data.items():
        #    if isinstance(v,list):
        #        data[k]=np.array(v)
        DOF2Nodes = data['DOF2Nodes']
        self.nDOF_red = data['nDOF_red']
        PhiM      = data['PhiM']
        PhiR      = data['PhiR']
        Nodes     = data['Nodes']
        Elements  = data['Elements']
        DOFs= np.arange(self.nDOF_red)
        # Reindexing with -1
        DOF_L = data['DOF___L'].ravel()-1 # internal DOFs
        DOF_B = data['DOF___B'].ravel()-1 # internal
        DOF_F = data['DOF___F'].ravel()

        # --- TODO merge with subdynSumToGraph(data) in fast_input_file_graph...
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
                #N[4]=np.float(N[4].split()[0])
                N=N.astype(np.float32)
            ID = int(N[0])-1
            nodeDOFs=DOF2Nodes[(DOF2Nodes[:,1]==ID),0]
            node = GraphNode(ID=ID, x=N[1], y=N[2], z=N[3], Type=int(N[4]), DOFs=nodeDOFs)
            self.addNode(node)

        # Elements
        for ie,E in enumerate(Elements):
            nodeIDs=[int(E[1]),int(E[2])]
            N1 = self.Nodes[nodeIDs[0]]
            N2 = self.Nodes[nodeIDs[1]]
            elem= GraphElement(int(E[0]), nodeIDs, nodes=[N1, N2])
            self.addElement(elem)

        #print(self.extent)
        #print(self.maxDimension)

        DOF_K = np.concatenate((DOF_B,DOF_L))
        maxDisp = self.maxDimension*0.1

        # CB modes
        if PhiM is not None:
            Q_CB = np.vstack((np.zeros((len(DOF_B),PhiM.shape[1])),PhiM, np.zeros((len(DOF_F),PhiM.shape[1]))))
            dispCB, posCB, INodesCB = self.NodesDisp(DOF_K, Q_CB, maxDisp=maxDisp)

        # Guyan modes
        Q_GY = np.vstack((np.eye(len(DOF_B)),PhiR))
        dispGY, posGY, INodesGY = self.NodesDisp(DOF_K, Q_GY, maxDisp=maxDisp, sortDim=2)


        # --- FEM Model
        self.f_G    = data['GY_frequencies'].flatten()
        self.f_CB   = data['CB_frequencies'].flatten()
        if PhiM is not None:
            self.Phi_CB = data['PhiM'] # Not Extended with BC, might need checking
            self.Q_CB   = Q_CB         # full
        self.Phi_G = data['PhiR'] # Not Extended, might need checking
        self.Q_GY  = Q_GY         # full

        # --- Add modes
        # TODO this is more or less FEMModel.setModes
        for iMode in range(dispGY.shape[2]):
            try:
                self.addMode(displ=dispGY[:,:,iMode], name='GY{:d}'.format(iMode+1), freq=self.f_G[iMode], group='GY')
            except:
                raise 
        for iMode in range(len(self.f_CB)):
            self.addMode(displ=dispCB[:,:,iMode], name='CB{:d}'.format(iMode+1), freq=self.f_CB[iMode], group='CB') 

        # --- Store
        if PhiM is not None:
            self.dispCB   = dispCB
            self.posCB    = posCB
            self.INodesCB = INodesCB
            self.dfCB = self.toDF(posCB, dispCB)
        else:
            self.dfCB = None

        self.dispGY   = dispGY
        self.posGY    = posGY
        self.INodesGY = INodesGY
        self.dfGY     = self.toDF(posGY, dispGY)

        try:
            # Full system matrices
            MM        = data['M']
            KK        = data['K']
            self.setFullMatrices(MM,KK)
        except:
            print('[INFO] Full mass and stiffness matrices not included in summary file')
            pass

    def toDF(self, pos,disp):
        columns=['z_[m]','x_[m]','y_[m]']
        dataZXY=np.column_stack((pos[:,2],pos[:,0],pos[:,1]))
        disptot=disp.copy()
        for ishape in np.arange(disp.shape[2]):
            disptot[:,:,ishape]= pos + disp[:,:,ishape]
            sMode='Mode{:d}'.format(ishape+1)
            columns+=[sMode+'x_[m]',sMode+'y_[m]',sMode+'z_[m]']
        disptot= np.moveaxis(disptot,2,1).reshape(disptot.shape[0],disptot.shape[1]*disptot.shape[2])
        disp   = np.moveaxis(disp,2,1).reshape(disp.shape[0],disp.shape[1]*disp.shape[2])
        data = np.column_stack((dataZXY,disptot))
        df= pd.DataFrame(data = data ,columns = columns)

        disp[np.isnan(disp)]=0
        data = np.column_stack((dataZXY,disp))
        dfDisp= pd.DataFrame(data = data ,columns = columns)
        # remove zero 
        df2= df.loc[:, (dfDisp !=0).any(axis=0)]
        return df

    def toDataFrame(self):
        return {'GY':self.dfGY,'CB':self.dfCB}
#         Mr, Kr, Phi_G, Phi_CB, f_CB, f_G = self.CraigBampton(Ileader=DOF_B, Ifollow=DOF_L)
#         print(data['PhiL'])

    def NodesDisp(self, IDOF, UDOF, maxDisp=None, sortDim=None):
        DOF2Nodes = self.DOF2Nodes
        INodes = list(np.sort(np.unique(DOF2Nodes[IDOF,1]))) # NOTE: sorted
        nShapes = UDOF.shape[1]
        disp=np.empty((len(INodes),3,nShapes)); disp.fill(np.nan)
        pos=np.empty((len(INodes),3))         ; pos.fill(np.nan)
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
    from welib.fast.subdyn_sum import SubDynModel
    filename = 'SybDyn-sum-MC-04JAN23.yaml'
    filename = '../../data/example_files/FASTSum_Pendulum.SD.sum.yaml'

    np.set_printoptions(linewidth=500)
    mdl=SubDynModel()
    mdl.fromSummaryFile()
    mdl.toJSON('_OUT.json')

    # --- Lower level functions
    #from welib.weio.fast_summary_file import FASTSummaryFile
    #sds = FASTSummaryFile(filename) 
    #df = sds.toDataFrame(sortDim=2)
    #print(df)
    #sds.toJSON()
    #sds.toGraph().toJSON('_OUT3.json')
    #df = sds.toDataFrame(sortDim=None)

    #import matplotlib.pyplot as plt
    #plt.plot(df['z_[m]'], df['GuyanMode1x_[m]'])
    #plt.plot(df['z_[m]'], df['GuyanMode5x_[m]'])
    #plt.show()
