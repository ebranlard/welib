import numpy as np
import pandas as pd
from io import open
import os
# Local 
from .mini_yaml import yaml_read


try:
    from .file import File, EmptyFileError
except:
    EmptyFileError = type('EmptyFileError', (Exception,),{})
    File=dict

# --------------------------------------------------------------------------------}
# --- Main Class 
# --------------------------------------------------------------------------------{
class FASTSummaryFile(File):
    """ 
    Read an OpenFAST summary file (.sum, .yaml). The object behaves as a dictionary.
    NOTE: open new subdyn format supported.

    Main methods
    ------------
    - read, toDataFrame

    Examples
    --------

        # read a subdyn summary file
        sum = FASTSummaryFile('5MW.SD.sum.yaml')
        print(sum['module']) # SubDyn
        M = sum['M'] # Mass matrix
        K = sum['K'] # stiffness matrix

    """

    @staticmethod
    def defaultExtensions():
        return ['.sum','.yaml']

    @staticmethod
    def formatName():
        return 'FAST summary file'

    def __init__(self,filename=None, **kwargs):
        self.filename = None
        if filename:
            self.read(filename, **kwargs)

    def read(self, filename=None, header_only=False):
        """ """
        if filename:
            self.filename = filename
        if not self.filename:
            raise Exception('No filename provided')
        if not os.path.isfile(self.filename):
            raise OSError(2,'File not found:',self.filename)
        if os.stat(self.filename).st_size == 0:
            raise EmptyFileError('File is empty:',self.filename)

        with open(self.filename, 'r', errors="surrogateescape") as fid:
            header= readFirstLines(fid, 4)
        if any(['subdyn' in s.lower() for s in header]):
            self['module']='SubDyn'
            data=readSubDynSum(self.filename)
            for k,v in data.items():
                self[k]=v
        else:
            raise NotImplementedError('This summary file format is not yet supported')

    def toDataFrame(self):
        if 'module' not in self.keys():
            raise Exception('');

        if self['module']=='SubDyn':
            dfs=subDynToDataFrame(self)

        return dfs


# --------------------------------------------------------------------------------}
# --- Helper functions 
# --------------------------------------------------------------------------------{
def readFirstLines(fid, nLines):
    lines=[]
    for i, line in enumerate(fid):
        lines.append(line.strip())
        if i==nLines:
            break
    return lines

# --------------------------------------------------------------------------------}
# --- Sub reader for different summary files
# --------------------------------------------------------------------------------{
def readSubDynSum(filename):
    #T=yaml.load(fid, Loader=yaml.SafeLoader)
    T=yaml_read(filename)
    # Convert lists to np arrays for  convenience
    #for k,v in T.items():
    #    if isinstance(v,list):
    #        T[k]=np.array(v)
    return T

def subDynToDataFrame(data):
    def NodesDisp( IDOF, UDOF, maxDisp=None, sortDim=None):
        INodes = list(np.sort(np.unique(DOF2Nodes[IDOF,1]))) # NOTE: sorted
        nShapes = UDOF.shape[1]
        disp=np.empty((len(INodes),3,nShapes)); disp.fill(np.nan)
        pos=np.empty((len(INodes),3))         ; pos.fill(np.nan)
        for i,iDOF in enumerate(IDOF):
            iNode       = DOF2Nodes[iDOF,1]
            nDOFPerNode = DOF2Nodes[iDOF,2]
            nodeDOF     = DOF2Nodes[iDOF,3]
            iiNode      = INodes.index(iNode)
            if nodeDOF<=3:
                pos[iiNode, 0]=X[iNode-1]
                pos[iiNode, 1]=Y[iNode-1]
                pos[iiNode, 2]=Z[iNode-1]
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

    DOF2Nodes = data['DOF2Nodes']
    PhiM      = data['PhiM']
    PhiR      = data['PhiR']
    Nodes     = data['Nodes']
    if DOF2Nodes.shape[1]==3:
        DOF2Nodes=np.column_stack((np.arange(DOF2Nodes.shape[0]),DOF2Nodes))
    else:
        DOF2Nodes[:,0]-=1
    X,Y,Z=Nodes[:,1].astype(float),Nodes[:,2].astype(float),Nodes[:,3].astype(float)
    dx,dy,dz=np.max(X)-np.min(X), np.max(Y)-np.min(Y),np.max(Z)-np.min(Z)
    maxDisp = np.max([dx,dy,dz])*0.1

    DOF_B=data['DOF___B'].ravel()
    DOF_F=data['DOF___F'].ravel()
    DOF_K = (np.concatenate((DOF_B,data['DOF___L'].ravel(), DOF_F))-1).astype(int)

    # CB modes
    Phi_CB = np.vstack((np.zeros((len(DOF_B),PhiM.shape[1])),PhiM, np.zeros((len(DOF_F),PhiM.shape[1]))))
    dispCB, posCB, INodes = NodesDisp(DOF_K, Phi_CB, maxDisp=maxDisp, sortDim=2)

    # Guyan modes
    Phi_Guyan = np.vstack((np.eye(len(DOF_B)),PhiR, np.zeros((len(DOF_F),PhiR.shape[1]))))
    dispGy, posGy, INodesGy = NodesDisp(DOF_K, Phi_Guyan, maxDisp=maxDisp, sortDim=2)

    def toDF(pos,disp,preffix=''):
        disp[np.isnan(disp)]=0
        disptot=disp.copy()
        columns=[]
        for ishape in np.arange(disp.shape[2]):
            disptot[:,:,ishape]= pos + disp[:,:,ishape]
            sMode=preffix+'Mode{:d}'.format(ishape+1)
            columns+=[sMode+'x_[m]',sMode+'y_[m]',sMode+'z_[m]']
        disptot= np.moveaxis(disptot,2,1).reshape(disptot.shape[0],disptot.shape[1]*disptot.shape[2])
        disp   = np.moveaxis(disp,2,1).reshape(disp.shape[0],disp.shape[1]*disp.shape[2])
        df= pd.DataFrame(data = disptot ,columns = columns)
        # remove zero 
        dfDisp= pd.DataFrame(data = disp ,columns = columns)
        df     = df.loc[:,     (dfDisp != 0).any(axis=0)]
        dfDisp = dfDisp.loc[:, (dfDisp != 0).any(axis=0)]
        dfDisp.columns = [c.replace('Mode','Disp') for c in dfDisp.columns.values]
        return df, dfDisp
    columns = ['z_[m]','x_[m]','y_[m]']
    dataZXY = np.column_stack((posGy[:,2],posGy[:,0],posGy[:,1]))
    dfZXY   = pd.DataFrame(data = dataZXY, columns=columns)
    df1, df1d = toDF(posGy, dispGy,'Guyan')
    df2, df2d = toDF(posCB, dispCB,'CB')
    df = pd.concat((dfZXY, df1, df2, df1d, df2d), axis=1)
    return df 





if __name__=='__main__':
    T=FASTSummaryFile('../Pendulum.SD.sum.yaml')
    df=T.toDataFrame()
    print(df)
