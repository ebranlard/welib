import numpy as np
import pickle
import glob
import os
from welib.weio.fast_linearization_file import FASTLinearizationFile
import pandas as pd


class FASTLinPeriodicOP(object):
    """ Class for a set of *.lin files, all assumed to be for the same periodic operating point
    e.g. 
       ws05mps.1.lin
              [...]
       ws05mps.36.lin

    """
    def __init__(self,prefix,nLin=None):
        if nLin is None:
            linfiles= glob.glob(prefix + '.*.lin') # TODO we want a more rigorous regexp
            self.nLinTimes = len(linfiles)
        else:
            self.nLinTimes = nLin
        #print(prefix, self.nLinTimes)

        self.prefix   = prefix
        self.Data     = []
        self.vAzim    = []
        self.vWS       = []
        self.vPitch    = []
        self.vRotSpeed = []
        self.vBu = []
        for i in np.arange(self.nLinTimes):
            linfilename= prefix+'.'+str(i+1)+'.lin'
            print(linfilename)
            if not os.path.exists(linfilename):
                print('Linearization file missing: ',linfilename)
            linfile=FASTLinearizationFile(linfilename)
            df=linfile.toDataFrame()
            self.Data.append(linfile)
            #self.A=lin['A']
            #B=linfile['B']
            #u=linfile['u']
            #self.C=lin['C']
            #self.D=lin['D']
            if linfile['WindSpeed'] is not None:
                self.vWS.append(linfile['WindSpeed'])
            else:
                try:
                    self.vWS.append(df['u']['WS_[m/s]'][0])
                except:
                    print('Wind speed not found in input, assuming 0m/s')
                    self.vWS.append(0)
            self.vRotSpeed.append(linfile['RotSpeed'])
            self.vAzim.append(linfile['Azimuth'])
            self.vPitch.append(df['u']['B1pitch_[rad]'][0]*180/np.pi)

        self.WS       = np.mean(self.vWS)
        self.Pitch    = np.mean(self.vPitch)
        self.RotSpeed = np.mean(self.vRotSpeed)

        self.x = df['x']
        self.y = df['y']
        self.u = df['u']
        try:
            self.EDdescr = linfile['EDDOF']
        except:
            self.EDdescr = None



class FASTLin(object):
    """ Class for linearization data for different operating points (typically Campbell) """
    def __init__(self, linfiles=None, folder='./', prefix='', nLin=None):
        """ 
        Init with a list of linfiles, or a folder and prefix
        """
        linfiles = [] if linfiles is None else linfiles

        if not isinstance(linfiles, list):
            linfiles=[linfiles]
        if len(linfiles)>0:
            exts =[os.path.splitext(f)[1] for f in linfiles]
            extsOK =[e.lower()=='.lin' for e in exts]
            if not all(extsOK):
                raise Exception('Not all inputs have the .lin extension. Provide a list of .lin files, or a folder and a prefix')
        else:
            linfiles= list(glob.glob(folder + prefix + '*.*.lin')) # TODO we want a more rigorous regexp
            linfiles.sort()

        Sim_Prefix=np.unique(['.'.join(f.split('.')[:-2]) for f in linfiles])
        nSim      = len(Sim_Prefix)
        # --- Read period operating points
        print('Reading linearizations for {} operating points'.format(nSim))
        self.OP_Data=[FASTLinPeriodicOP(pref,nLin=nLin) for pref in Sim_Prefix]
        # --- Sort by wind speed
        Isort = np.argsort(self.WS)
        self.OP_Data  = [self.OP_Data[i] for i in Isort]

        if self.MaxNLinTimes>1:
            IBad = [i for i in np.arange(nSim) if self.nLinTimes[i]<self.MaxNLinTimes and self.OP_Data[i].WS>0]
            if len(IBad)>0: 
                print('>>> The following simulations have insufficient number of data points:')
                for i in IBad:
                    print(self.OP_Data[i].prefix, self.OP_Data[i].nLinTimes)
            self.OP_Data = [self.OP_Data[i] for i in np.arange(nSim) if i not in IBad]

    @property
    def WS(self):
        return np.array([sim.WS for sim in self.OP_Data])

    @property
    def nLinTimes(self):
        return np.array([sim.nLinTimes for sim in self.OP_Data])

    @property
    def MaxNLinTimes(self):
        return np.max(self.nLinTimes)

    @property
    def nOP(self):
        return len(self.OP_Data)

    @property
    def xdescr(self):
        return self.OP_Data[0].x.columns.values
    @property
    def ydescr(self):
        return self.OP_Data[0].y.columns.values
    @property
    def EDdescr(self):
        return self.OP_Data[0].EDdescr
    @property
    def udescr(self):
        return self.OP_Data[0].u.columns.values
    @property
    def xop_mean(self):
        return np.mean(np.abs(np.array([op.x.values for op in self.OP_Data])),axis=0)
    @property
    def uop_mean(self):
        return np.mean(np.abs(np.array([op.u.values for op in self.OP_Data])),axis=0)
    @property
    def uop_mean(self):
        return np.mean(np.abs(np.array([op.u.values for op in self.OP_Data])),axis=0)

    @property
    def yop_mean(self):
        return np.mean(np.abs(np.array([op.y.values for op in self.OP_Data])),axis=0)

    def stats(self,matName,WS=None):
        if WS is None:
            WS = self.WS
            nOP=self.nOP
        else:
            nOP=len(WS)
        print('Returning stats for WS:',WS)
        M_mean=[]

        shape= self.OP_Data[0].Data[0][matName].shape

        M_all       = np.zeros( (nOP, self.MaxNLinTimes, shape[0],shape[1]))
        M_mean_perWS= np.zeros( (nOP, shape[0],shape[1]))
        M_std_perWS = np.zeros( (nOP, shape[0],shape[1]))

        # loop on operating points (e.g. WS)
        ii=0
        for iop, op in enumerate(self.OP_Data):
            if op.WS in WS:
                # Loop on linearization times (e.g. Azimuth)
                for iTimes in np.arange(self.MaxNLinTimes):
                    if op.nLinTimes==1:
                        M_all[ii,iTimes,:,:]=op.Data[0][matName]
                    else:
                        M_all[ii,iTimes,:,:]=op.Data[iTimes][matName]

                M_mean_perWS[ii,:,:] = np.mean(M_all[ii,:,:,:],axis=0)
                M_std_perWS [ii,:,:]  = np.std(M_all[ii,:,:,:],axis=0)
                ii+=1

        M_mean    = np.mean( M_mean_perWS, axis=0 )
        M_stdWS   = np.std ( M_mean_perWS, axis=0 ) # How much elements vary with wind speed
        M_stdAzim = np.mean( M_std_perWS , axis=0)  # How much elements vary due to azimuth

        return M_mean, M_mean_perWS, M_stdAzim, M_stdWS, M_all


    def average(self, WS=None):
        A_mean = self.stats('A',WS=WS)[0]
        B_mean = self.stats('B',WS=WS)[0]
        C_mean = self.stats('C',WS=WS)[0]
        D_mean = self.stats('D',WS=WS)[0]
        #self.M_mean = self.stats('M',WS=WS)[0]
        return A_mean, B_mean, C_mean, D_mean

    def average_subset(self, sX_sel=None, sU_sel=None, sY_sel=None, sE_sel=None, WS=None, exportFile=None, baseDict=None):
        """ 
        Average state spaces based on WS, then extract a subset based on sensor names
        """
        sX, sU, sY, sED = self.xdescr, self.udescr, self.ydescr, self.EDdescr

        if sX_sel is None:
            sX_sel=sX
        if sU_sel is None:
            sU_sel=sU
        if sY_sel is None:
            sY_sel=sY

        # Average
        A,B,C,D = self.average(WS=WS)

        # Indices
        try:
            IDOFX = np.array([list(sX).index(s) for s in sX_sel])
        except:
            print(sX)
            raise
        IDOFU = np.array([list(sU).index(s) for s in sU_sel])
        IDOFY = np.array([list(sY).index(s) for s in sY_sel])
        if sE_sel is not None:
            IDOFE = np.array([list(sE).index(s) for s in sE_sel])

        # Subset. TODO use DataFrame directly and subMat in linmodel
        Ar = A[np.ix_(IDOFX,IDOFX)]
        Br = B[np.ix_(IDOFX,IDOFU)]
        Cr = C[np.ix_(IDOFY,IDOFX)]
        Dr = D[np.ix_(IDOFY,IDOFU)]

        # Outputs
        Ar = pd.DataFrame(data = Ar, index=sX_sel, columns=sX_sel)
        Br = pd.DataFrame(data = Br, index=sX_sel, columns=sU_sel)
        Cr = pd.DataFrame(data = Cr, index=sY_sel, columns=sX_sel)
        Dr = pd.DataFrame(data = Dr, index=sY_sel, columns=sU_sel)
        if baseDict is None:
            outDict={}
        else:
            outDict=baseDict.copy()
        outDict['A']=Ar
        outDict['B']=Br
        outDict['C']=Cr
        outDict['D']=Dr
        if sE_sel is not None:
            Mr = M[np.ix_(IDOFE,IDOFE)]
            Mr = pd.DataFrame(data = Mr, index=sED_sel, columns=sED_sel)
            outDict['M']=Mr

        if exportFile is not None:
            import pickle
            with open(exportFile,'wb') as f:
                pickle.dump(outDict,f)

        if sE_sel is not None:
            return Ar, Br, Cr, Dr, Mr
        else:
            return Ar, Br, Cr, Dr

    def save(self,filename):
        with open(filename,'wb') as f:
            pickle.dump(self,f)

