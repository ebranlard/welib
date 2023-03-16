""" 
Tools to handle a set of linearization files

"""

import numpy as np
import pickle
import glob
import os
import re
from welib.weio.fast_linearization_file import FASTLinearizationFile
import pandas as pd


class FASTLinPeriodicOP(object):
    """ Class for a set of *.lin files, all assumed to be for the same periodic operating point
    e.g. 
       ws05mps.1.lin
              [...]
       ws05mps.36.lin

    """
    def __init__(self, prefix=None, nLin=None, linFiles=None):

        # --- Init data
        self.linFiles  = []
        self.prefix    = None
        self.Data      = []     # List of linFile as returned by weio
        self.vAzim     = []
        self.vWS       = []
        self.vPitch    = []
        self.vRotSpeed = []

        # --- Figure out linFiles
        def glob_re(pattern_glob, pattern_re):
            """ 
            glob_re(r'base.([0-9]*).py', 'base.*.py') 
            """
            files = glob.glob(pattern_glob)
            files = [s.replace('\\','/') for s in files]
            return list(filter(re.compile(pattern_re).match, files))
        if linFiles is None:
            if nLin is None:
                prefix = prefix.replace('\\','/')
                linFiles = glob_re(prefix + '.*.lin', prefix + r'.([0-9]*).lin',)
                self.nLinTimes = len(linFiles)
            else:
                self.nLinTimes = nLin

            linFiles = [prefix+'.'+str(i+1)+'.lin' for i in np.arange(self.nLinTimes)]
            if len(linFiles)==0:
                raise Exception('No Lin Files found with prefix: {}'.format(prefix))
        else:
            self.nLinTimes = len(linFiles)
            prefix = None # TODO

        self.linFiles  = linFiles
        self.prefix    = prefix
        self.Data      = []     # List of linFile as returned by weio
        self.vAzim     = []
        self.vWS       = []
        self.vPitch    = []
        self.vRotSpeed = []
        for i, linFilename in enumerate(linFiles):
            print(linFilename)
            if not os.path.exists(linFilename):
                print('Linearization file missing: ',linFilename)
            linfile = FASTLinearizationFile(linFilename)
            df      = linfile.toDataFrame()
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
            if 'u' in df.keys():
                self.vPitch.append(df['u']['B1pitch_[rad]'][0]*180/np.pi)
            else:
                self.vPitch.append(np.nan)

        self.WS       = np.mean(self.vWS)
        self.Pitch    = np.mean(self.vPitch)
        self.RotSpeed = np.mean(self.vRotSpeed)

        self.x = df['x']
        self.y = None
        self.u = None
        if 'y' in df.keys():
            self.y = df['y']
        if 'u' in df.keys():
            self.u = df['u']
        try:
            self.EDdescr = linfile['EDDOF']
        except:
            self.EDdescr = None


    def __repr__(self):
        s ='<FASTLinPeriodicOP object>\n'
        s+='Attributes:\n'
        s+=' - prefix   : {}\n'.format(self.prefix)
        s+=' - WS       : {}\n'.format(self.WS)
        s+=' - Pitch    : {}\n'.format(self.Pitch)
        s+=' - RotSpeed : {}\n'.format(self.RotSpeed)
        s+=' - vAzim    : {}\n'.format(self.vAzim)
        s+=' - vWS      : {}\n'.format(self.vWS)
        s+=' - vPitch   : {}\n'.format(self.vPitch)
        s+=' - vRotSpeed: {}\n'.format(self.vRotSpeed)
        s+=' - linFiles : {}\n'.format(self.linFiles)
        s+=' - Data     : list of lin files, size {}\n'.format(len(self.Data))
        return s

class FASTLin(object):
    """ Class to handle linearization data at different operating points 
        Typically Campbell, or average over many conditions.
        Can be used for one lin file as well.
    """
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

        simPrefix=np.unique(['.'.join(f.split('.')[:-2]) for f in linfiles])
        nSim      = len(simPrefix)
        self.simPrefix = simPrefix
        # --- Read period operating points
        print('Reading linearizations for {} operating points'.format(nSim))
        self.OP_Data=[FASTLinPeriodicOP(pref, nLin=nLin) for pref in simPrefix]

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

    def __repr__(self):
        s ='<FASTLin object>\n'
        s+='Attributes:\n'
        s+=' - OP_Data     : list of FASTLinPeriodicOP (size {})\n'.format(len(self.OP_Data))
        s+=' * nOP         : {}\n'.format(self.nOP)
        s+=' * MaxNLinTimes: {}\n'.format(self.MaxNLinTimes)
        s+=' * WS          : {}\n'.format(self.WS)
        s+=' * nLinTimes   : {}\n'.format(self.nLinTimes)
        s+=' * xdescr, udescr, ydescr\n'
        s+=' * xop_mean, uop_mean, yop_mean\n'
        s+=' - simPrefix   : {}\n'.format(self.simPrefix)
        s+='Methods:\n'
        s+=' - stats(matName, WS=None)\n'
        s+=' - average(WS=None)\n'
        s+=' - exportState(self, stateFile, stateDict)\n'
        return s

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
        if self.hasY:
            return self.OP_Data[0].y.columns.values
        else:
            return []
    @property
    def EDdescr(self):
        return self.OP_Data[0].EDdescr
    @property
    def udescr(self):
        if self.hasU:
            return self.OP_Data[0].u.columns.values
        else:
            return []
    @property
    def xop_mean(self):
        return np.mean(np.array([op.x.values for op in self.OP_Data]),axis=0)
    @property
    def uop_mean(self):
        if self.hasU:
            return np.mean(np.array([op.u.values for op in self.OP_Data]),axis=0)
        else:
            raise Exception('Linear model has no inputs')
     
    @property
    def yop_mean(self):
        if hasY:
            return np.mean(np.array([op.y.values for op in self.OP_Data]),axis=0)
        else:
            raise Exception('Linear model has no outputs')

    @property
    def hasU(self): return 'u' in self.OP_Data[0].Data[0].keys()
    @property
    def hasY(self): return 'y' in self.OP_Data[0].Data[0].keys()
    @property
    def hasB(self): return 'B' in self.OP_Data[0].Data[0].keys()
    @property
    def hasC(self): return 'C' in self.OP_Data[0].Data[0].keys()
    @property
    def hasD(self): return 'D' in self.OP_Data[0].Data[0].keys()


    def stats(self, matName, WS=None):
        """ 
        Compute statistics (mean and std) on a given matrix (A, B, C, D, M)
        
        """
        if WS is None:
            WS = self.WS
            nOP=self.nOP
        else:
            nOP=len(WS)
        M_mean=[]

        shape = self.OP_Data[0].Data[0][matName].shape

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

                M_mean_perWS[ii,:,:] = np.mean(M_all[ii,:,:,:],axis=0) # TODO what if MaxNLinTimes is not the same for all OP
                M_std_perWS [ii,:,:]  = np.std(M_all[ii,:,:,:],axis=0)
                ii+=1

        M_mean    = np.mean( M_mean_perWS, axis=0 )
        M_stdWS   = np.std ( M_mean_perWS, axis=0 ) # How much elements vary with wind speed
        M_stdAzim = np.mean( M_std_perWS , axis=0)  # How much elements vary due to azimuth

        return M_mean, M_mean_perWS, M_stdAzim, M_stdWS, M_all


    def averageOP(self, WS=None):
        """ return average operating point values for a given wind speed vector"""
        if WS is None:
            WS = self.WS
        xop = np.zeros(len(self.xdescr))
        uop = np.zeros(len(self.udescr))
        yop = np.zeros(len(self.ydescr))
        for iop, op in enumerate(self.OP_Data):
            if self.WS[iop] in WS:
                xop+=op.x.values.flatten()
                if self.hasU:
                    uop+=op.u.values.flatten()
                if self.hasY:
                    yop+=op.y.values.flatten()
        xop /= len(WS)
        if self.hasU:
            uop /= len(WS)
        else:
            uop=None
        if self.hasY:
            yop /= len(WS)
        else:
            yop=None
        return xop, uop, yop

    def average(self, WS=None):
        A_mean = self.stats('A',WS=WS)[0]

        B_mean = None
        C_mean = None
        D_mean = None
        if self.hasB:
            B_mean = self.stats('B',WS=WS)[0]
        if self.hasC:
            C_mean = self.stats('C',WS=WS)[0]
        if self.hasD:
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
            self.exportState(exportFile, outDict)

        if sE_sel is not None:
            return Ar, Br, Cr, Dr, Mr
        else:
            return Ar, Br, Cr, Dr

    def exportState(self, stateFile, stateDict):
        #if any(['A','B','C','D'])

        import pickle
        with open(exportFile,'wb') as f:
            pickle.dump(stateDict,f)

    def save(self,filename):
        with open(filename,'wb') as f:
            pickle.dump(self,f)



