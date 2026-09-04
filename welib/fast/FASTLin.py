""" 
Tools to handle a set of linearization files

"""

import numpy as np
import pickle
import glob
import os
import re
from welib.weio.fast_linearization_file import FASTLinearizationFile, LinHasNAError
from welib.fast.tools.lin import subMat
from welib.tools.clean_exceptions import *
from welib.tools.strings import FAIL, INFO, WARN
import pandas as pd


class FASTLinPeriodicOP(object):
    """ Class for a set of *.lin files, all assumed to be for the same periodic operating point
    e.g. 
       ws05mps.1.lin
              [...]
       ws05mps.36.lin

    """
    def __init__(self, prefix=None, nLin=None, linFiles=None, verbose=False):

        # --- Init data
        self.linFiles  = []
        self.prefix    = None
        self.Data      = []     # List of linFile as returned by weio
        self.vAzim     = []
        self.vWS       = []
        self.vPitch    = []
        self.vRotSpeed = []
        self.x       = None
        self.y       = None
        self.u       = None
        self.EDdescr = None

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
        df = None
        for i, linFilename in enumerate(linFiles):
            if not os.path.exists(linFilename):
                FAIL('Linearization file missing: ',linFilename)
                continue
            try:
                linfile = FASTLinearizationFile(linFilename)
            except LinHasNAError:
                FAIL('Linearization file has NaN: ',linFilename)
                continue
                
            print(linFilename, f'nx:{linfile.nx} ny:{linfile.ny} nu:{linfile.nu}')
            df = linfile.toDataFrame()
            self.Data.append(linfile)

            if linfile['WindSpeed'] is not None:
                self.vWS.append(linfile['WindSpeed'])
            else:
                try:
                    self.vWS.append(df['u']['WS_[m/s]'][0])
                except:
                    FAIL('Wind speed not found in input, assuming 0m/s')
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

        if df is None:
            FAIL('All lin files for this OP are problematic')
        else:
            self.x = df['x']
            if 'y' in df.keys():
                self.y = df['y']
            if 'u' in df.keys():
                self.u = df['u']
            try:
                self.EDdescr = linfile['EDDOF']
            except:
                pass

    def __repr__(self):
        s ='<FASTLinPeriodicOP object>\n'
        s+='Attributes:\n'
        s+=' - prefix   : {}\n'.format(self.prefix)
        s+=' - WS       : {}\n'.format(self.WS)
        s+=' - Pitch    : {}\n'.format(self.Pitch)
        s+=' - RotSpeed : {}\n'.format(self.RotSpeed)
        s+=' * nx       : {}\n'.format(self.nx)
        s+=' * ny       : {}\n'.format(self.ny)
        s+=' * nu       : {}\n'.format(self.nu)
        s+=' - x        : DataFrame, len: {}\n'.format(len(self.x))
        s+=' - y        : DataFrame, len: {}\n'.format(len(self.y))
        s+=' - u        : DataFrame, len: {}\n'.format(len(self.u))
        s+=' - vAzim    : {}\n'.format(self.vAzim)
        s+=' - vWS      : {}\n'.format(self.vWS)
        s+=' - vPitch   : {}\n'.format(self.vPitch)
        s+=' - vRotSpeed: {}\n'.format(self.vRotSpeed)
        s+=' - linFiles : {}\n'.format(self.linFiles)
        s+=' - Data     : list of weio lin files, size {}\n'.format(len(self.Data))
        return s


    @property
    def nx(self):
        return len(self.x.values.flatten()) if self.x is not None else 0

    @property
    def nu(self):
        return len(self.u.values.flatten()) if self.u is not None else 0

    @property
    def ny(self):
        return len(self.y.values.flatten()) if self.y is not None else 0



class FASTLin(object):
    """ Class to handle linearization data at different operating points 
        Typically Campbell, or average over many conditions.
        Can be used for one lin file as well.
    """
    def __init__(self, linfiles=None, folder='./', prefix='', nLin=None, verbose=False):
        """ 
        Init with a list of linfiles, or a folder and prefix
        """
        # Data init
        self.OP_Data = []
        self.simPrefixes = []
        # Stats Data
        self.A_mean, self.A_mean_perWS, self.A_stdAzim, self.A_stdWS = None, None, None, None
        self.B_mean, self.B_mean_perWS, self.B_stdBzim, self.B_stdWS = None, None, None, None
        self.C_mean, self.C_mean_perWS, self.C_stdCzim, self.C_stdWS = None, None, None, None
        self.D_mean, self.D_mean_perWS, self.D_stdDzim, self.D_stdWS = None, None, None, None
        self.M_mean, self.M_mean_perWS, self.M_stdMzim, self.M_stdWS = None, None, None, None
        # 
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

        _simPrefixes = np.unique(['.'.join(f.split('.')[:-2]) for f in linfiles])
        nSim         = len(_simPrefixes)
        if verbose:
            print(f'nFiles: {nSim}, prefixes: {_simPrefixes[0]}.., nLin={nLin}')
        # --- Read period operating points
        print('Reading linearizations for {} operating points'.format(nSim))
        for iOP, _prefix in enumerate(_simPrefixes):
            pOP = FASTLinPeriodicOP(_prefix, nLin=nLin)
            if len(pOP.Data)==0:
                FAIL(f'FASTLin: No Data present, skipping: {_prefix}')
                continue
            if iOP>0:
                if self.nx!=pOP.nx:
                    FAIL(f'FASTLin: Different number of states {self.nx} (first) /= {pOP.nx} (current: {_prefix})')
                    continue
            self.OP_Data.append(pOP)
            self.simPrefixes.append(_prefix)

        # --- Sort by wind speed
        Isort = np.argsort(self.WS)
        self.OP_Data  = [self.OP_Data[i] for i in Isort]

        if self.MaxNLinTimes>1:
            IBad = [i for i in np.arange(self.nOP) if self.nLinTimes[i]<self.MaxNLinTimes and self.OP_Data[i].WS>0]
            if len(IBad)>0: 
                FAIL('>>> The following simulations have insufficient number of data points:')
                for i in IBad:
                    print(self.OP_Data[i].prefix, self.OP_Data[i].nLinTimes)
            self.OP_Data = [self.OP_Data[i] for i in np.arange(self.nOP) if i not in IBad]

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
        s+=' - simPrefixes : {}\n'.format(self.simPrefixes)
        s+=' - A_mean, B_mean (after calling averate) \n'
        s+='Methods:\n'
        s+=' - stats(matName, WS=None)\n'
        s+=' - average(WS=None) (and store stats)\n'
        s+=' - exportState(self, stateFile, stateDict)\n'
        s+=' - save(picklefile)\n'
        s+=' - from_pickle(picklefile)\n'
        return s

    @property
    def WS(self): return np.array([sim.WS for sim in self.OP_Data])

    @property
    def nLinTimes(self): return np.array([sim.nLinTimes for sim in self.OP_Data])

    @property
    def MaxNLinTimes(self): return np.max(self.nLinTimes)

    @property
    def nOP(self): return len(self.OP_Data)


    @property
    def nx(self):
        return len(self.xdescr) if self.nOP>0 else 0

    @property
    def ny(self):
        return len(self.ydescr) if self.nOP>0 else 0

    @property
    def nu(self):
        return len(self.udescr) if self.nOP>0 else 0


    @property
    def xdescr(self): return self.OP_Data[0].x.columns.values

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
            for ws in WS:
                if ws not in self.WS:
                    raise Exception(f'FASTLin: Cannot compute stats for WS={WS}, it is not in the list of operating point WS: {self.WS}')
        M_mean=[]

        if matName not in self.OP_Data[0].Data[0]:
            raise KeyError(f'Column {matName} nor present in dataframe')
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

    def average(self, WS=None, return_dataframes=True):
        if len(self.OP_Data)==0:
            raise Exception('FASTLin: No Operating point data, cannot compute stats')
        if len(self.OP_Data[0].Data)==0:
            raise Exception('FASTLin: No Operating point data, cannot compute stats')
        if 'A' in self.OP_Data[0].Data[0]:
            self.A_mean, self.A_mean_perWS, self.A_stdAzim, self.A_stdWS, _ = self.stats('A', WS=WS)
        if 'B' in self.OP_Data[0].Data[0]:
            self.B_mean, self.B_mean_perWS, self.B_stdBzim, self.B_stdWS, _ = self.stats('B', WS=WS)
        if 'C' in self.OP_Data[0].Data[0]:
            self.C_mean, self.C_mean_perWS, self.C_stdCzim, self.C_stdWS, _ = self.stats('C', WS=WS)
        if 'D' in self.OP_Data[0].Data[0]:
            self.D_mean, self.D_mean_perWS, self.D_stdDzim, self.D_stdWS, _ = self.stats('D', WS=WS)
        if 'M' in self.OP_Data[0].Data[0]:
            self.M_mean, self.M_mean_perWS, self.M_stdMzim, self.M_stdWS, _ = self.stats('M', WS=WS)

        if return_dataframes:
            self.A_mean = pd.DataFrame(data = self.A_mean, index=self.xdescr, columns=self.xdescr)
            self.B_mean = pd.DataFrame(data = self.B_mean, index=self.xdescr, columns=self.udescr)
            self.C_mean = pd.DataFrame(data = self.C_mean, index=self.ydescr, columns=self.xdescr)
            self.D_mean = pd.DataFrame(data = self.D_mean, index=self.ydescr, columns=self.udescr)
            if 'M' in self.OP_Data[0].Data[0]:
                self.M_mean = pd.DataFrame(data = self.M_mean, index=self.EDdescr, columns=self.EDdescr)

        return self.A_mean, self.B_mean, self.C_mean, self.D_mean


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
        if sE_sel is None:
            sE_sel=sED

        # Average
        A,B,C,D = self.average(WS=WS, return_dataframes=True)

        Ar = subMat(A, rows=sX_sel, cols=sX_sel, check=True, name='A', removeDuplicates=True)
        Br = subMat(B, rows=sX_sel, cols=sU_sel, check=True, name='B', removeDuplicates=True)
        Cr = subMat(C, rows=sY_sel, cols=sX_sel, check=True, name='C', removeDuplicates=True)
        Dr = subMat(D, rows=sY_sel, cols=sU_sel, check=True, name='D', removeDuplicates=True)
        if self.M_mean is not None:
            Mr = subMat(self.M_mean, rows=sE_sel, cols=sE_sel, check=True, name='M', removeDuplicates=True)
        else:
            Mr = None

        if baseDict is None:
            outDict={}
        else:
            outDict=baseDict.copy()
        outDict['A'] = Ar
        outDict['B'] = Br
        outDict['C'] = Cr
        outDict['D'] = Dr
        outDict['M'] = Mr

        if exportFile is not None:
            self.exportState(exportFile, outDict)

        return Ar, Br, Cr, Dr, Mr



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



    def exportState(self, stateFile, stateDict):
        dirname = os.path.dirname(stateFile)
        if dirname:
            os.makedirs(dirname, exist_ok=True)
            
        with open(stateFile, 'wb') as f:
            pickle.dump(stateDict, f)
        INFO(f'Written StateFile : {stateFile}')

    def save(self, filename, verbose=True):
        dirname = os.path.dirname(filename)
        if dirname:
            os.makedirs(dirname, exist_ok=True)
            
        if verbose:
            INFO(f'Writting FASTLin Dump: {filename}')
        with open(filename, 'wb') as f:
            pickle.dump(self, f)

    @classmethod
    def from_pickle(cls, filename, verbose=False):
        if verbose:
            INFO(f'Loading FASTLin Dump: {filename}')
        with open(filename, "rb") as f:
            obj = pickle.load(f)
        if not isinstance(obj, cls):
            raise TypeError(f"Expected instance of {cls.__name__}, got {type(obj).__name__}")
        return obj


