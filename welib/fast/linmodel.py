##
import numpy as np
import copy
import os
import welib.weio as weio

# --------------------------------------------------------------------------------}
# --- Creating a TNSB model from a FAST model
# --------------------------------------------------------------------------------{
class FASTLinModel():
    def __init__(self,ED_or_FST_file, StateFile=None, nShapes_twr=1, nShapes_bld=0, DEBUG=False):

        # --- Input data from fst and ED file
        ext=os.path.splitext(ED_or_FST_file)[1]
        if ext.lower()=='.fst':
            FST=weio.read(ED_or_FST_file)
            rootdir = os.path.dirname(ED_or_FST_file)
            EDfile = os.path.join(rootdir,FST['EDFile'].strip('"')).replace('\\','/')
        else:
            EDfile=ED_or_FST_file
        self.ED= weio.read(EDfile)

        # --- Loading linear model
        if StateFile is not None:
            self.A,self.B,self.C,self.D,self.M = loadLinStateMatModel(StateFile)
        else:
            raise NotImplementedError()
        self.sX = self.A.columns

        self.nGear = self.ED['GBRatio']
        self.theta_tilt=-self.ED['ShftTilt']*np.pi/180 # NOTE: tilt has wrong orientation in FAST
        # --- Initial conditions
        omega_init = self.ED['RotSpeed']*2*np.pi/60 # rad/s
        psi_init   = self.ED['Azimuth']*np.pi/180 # rad
        FA_init    = self.ED['TTDspFA']
        iPsi     = list(self.sX).index('psi_rot_[rad]')
        nDOFMech = int(len(self.A)/2)
        q_init   = np.zeros(2*nDOFMech) # x2, state space

        if nShapes_twr>0:
            q_init[0] = FA_init

        q_init[iPsi]          = psi_init
        q_init[nDOFMech+iPsi] = omega_init

        self.q_init = q_init

    def __repr__(self):
        def pretty_PrintMat(M,fmt='{:11.3e}',fmt_int='    {:4d}   ',sindent='   '):
            s=sindent
            for iline,line in enumerate(M):
                s+=''.join([(fmt.format(v) if int(v)!=v else fmt_int.format(int(v))) for v in line ])
                s+='\n'+sindent
            return s
        s=''
        s+='<FASTLinModel object>\n'
        s+='Attributes:\n'
        s+=' - A: State-State Matrix  \n'
        s+=pretty_PrintMat(self.A.values)+'\n'
        s+=' - B: State-Input Matrix  \n'
        s+=pretty_PrintMat(self.B.values)+'\n'
        s+=' - q_init: Initial conditions (state) \n'
        s+=str(self.q_init)+'\n'
        return s

# --------------------------------------------------------------------------------}
# --- Helper functions 
# --------------------------------------------------------------------------------{
def unit(s):
    iu=s.rfind('[')
    if iu>1:
        return s[iu+1:].replace(']','')
    else:
        return ''
def no_unit(s):
    s=s.replace('_[',' [')
    iu=s.rfind(' [')
    if iu>1:
        return s[:iu]
    else:
        return s

def matToSIunits(Mat, name='', verbose=False):
    """ 
    Scale a matrix (pandas dataframe) such that is has standard units

    """
    # TODO columns
    for irow,row in enumerate(Mat.index.values):
        # Scaling based on units
        u  = unit(row).lower()
        nu = no_unit(row)
        if u in['deg','deg/s','deg/s^2']:
            if verbose:
                print('Mat {} - scaling deg2rad   for row {}'.format(name,row))
            Mat.iloc[irow,:] /=180/np.pi # deg 2 rad
            Mat.index.values[irow]=Mat.index.values[irow].replace('rad','deg')
        elif u=='rpm':
            if verbose:
                print('Mat {} - scaling rpm2rad/s for row {}'.format(name,row))
            Mat.iloc[irow,:] /=60/(2*np.pi) # rpm 2 rad/s
            Mat.index.values[irow]=nu+'_[rad/s]'
        elif u=='knm':
            if verbose:
                print('Mat {} - scaling kNm to Nm for row {}'.format(name,row))
            Mat.iloc[irow,:] /=1000 # to Nm
            Mat.index.values[irow]=nu+'_[Nm]'
    return Mat


def subMat(df, rows=None, cols=None, check=True):
    """ Extract relevant part from a dataframe, perform a safety check """
    if rows is None:
        rows=df.index
    if cols is None:
        cols=df.columns
    if check:
        missingRows = [l for l in rows if l not in df.index]
        missingCols = [c for c in cols  if c not in df.columns]
        if len(missingRows)>0:
            raise Exception('The following rows are missing from outputs: {}'.format(missingRows))
        if len(missingCols)>0:
            raise Exception('The following columns are missing from inputs: {}'.format(missingCols))
    return df[cols].loc[rows].copy()

def matLabelNoUnit(df):
    """ remove unit from index and columns of matrix"""
    rows = [no_unit(row) for row in df.index.values]
    cols = [no_unit(col) for col in df.columns.values]
    df.index =  rows
    df.columns = cols
    return df

def matLabelReplace(df, s1, s2):
    """ remove unit from index and columns of matrix"""
    rows = [row.replace(s1,s2) for row in df.index.values]
    cols = [col.replace(s1,s2) for col in df.columns.values]
    df.index =  rows
    df.columns = cols
    return df


def matSimpleStateLabels(df, inplace=False):
    RenameMap={
    'PtfmSurge_[m]'       : 'x',
    'PtfmSway_[m]'        : 'y',
    'PtfmHeave_[m]'       : 'z',
    'PtfmRoll_[rad]'      : 'phi_x',
    'PtfmPitch_[rad]'     : 'phi_y',
    'PtfmYaw_[rad]'       : 'phi_z',
    'psi_rot_[rad]'       : 'psi',
    'qt1FA_[m]'           : 'q_FA1',
    'd_PtfmSurge_[m/s]'   : 'dx',
    'd_PtfmSway_[m/s]'    : 'dy',
    'd_PtfmHeave_[m/s]'   : 'dz',
    'd_PtfmRoll_[rad/s]'  : 'dphi_x',
    'd_PtfmPitch_[rad/s]' : 'dphi_y',
    'd_PtfmYaw_[rad/s]'   : 'dphi_z',
    'd_psi_rot_[rad/s]'   : 'dpsi',
    'd_qt1FA_[m/s]'       : 'dq_FA1',
    }
    return df.rename(columns=RenameMap, index=RenameMap, inplace=inplace)


def loadLinStateMatModel(StateFile, ScaleUnits=True, Adapt=True, ExtraZeros=False, nameMap={'SvDGenTq_[kNm]':'Qgen_[kNm]'}, ):
    """ 
    Load a pickle file contains A,B,C,D matrices either as sequence or dictionary.
    Specific treatments are possible:
       - ScaleUnits: convert to SI units deg->rad, rpm-> rad/s, kN->N
       - Adapt: 
       - Adapt: 
       - nameMap: rename columns and indices

    If a "model" is given specific treatments can be done
    
    """
    import pickle
    # --- Subfunctions
    def load(filename):
        with open(filename,'rb') as f:
            dat=pickle.load(f)
        return dat

    # --- Load model
    dat = load(StateFile)
    if isinstance(dat,dict):
        A=dat['A']
        B=dat['B']
        C=dat['C']
        D=dat['D']
        M=None
        model =dat['model']
    else:
        model='TNSB'
        if len(dat)==4:
            M=None
            (A,B,C,D) = dat
        else:
            (A,B,C,D,M) = dat

    # --- Renaming
    for S,Mat in zip(['A','B','C','D'],[A,B,C,D]):
        for irow,row in enumerate(Mat.index.values):
            # Changing names
            if row=='SvDGenTq_[kNm]':
                Mat.index.values[irow]='Qgen_[kNm]'
                row='Qgen_[kNm]'



    # --- Scale units
    if ScaleUnits:
        # Changing rows
        for S,Mat in zip(['A','B','C','D'],[A,B,C,D]):
            Mat = matToSIunits(Mat, name=S, verbose=True)
    # --- ColMap
    if nameMap is not None:
        for S,Mat in zip(['A','B','C','D'],[A,B,C,D]):
            Mat.rename(nameMap, axis='columns', inplace=True)
            Mat.rename(nameMap, axis='index', inplace=True)

    # --- Numerics, 0
    for S,Mat in zip(['A','B','C','D'],[A,B,C,D]):
        Mat[np.abs(Mat)<1e-14]=0


    if model=='FNS' and A.shape[0]==6:
        pass
        #print(C)
        #print(D)
    elif model=='F1NS' and A.shape[0]==4:
        pass
    elif model=='F010000NS' and A.shape[0]==4:
        pass
    elif model=='F010010NS' and A.shape[0]==6:
        pass
    elif model=='F011010NS' and A.shape[0]==6:
        pass

    elif model=='FN' and A.shape[0]==4:
        pass

        
    elif model=='TNSB' and A.shape[0]==4:
        if Adapt==True:
            A.iloc[3,:]=0 # No state influence of ddpsi ! <<<< Important
            A.iloc[2,1]=0 # No psi influence of  ddqt
            A.iloc[2,3]=0 # No psi_dot influence of ddqt
            if ExtraZeros:
                B.iloc[0,:]=0 # No thrust influence on dqt
                B.iloc[1,:]=0 # No thrust influence on dpsi
            B.iloc[:,2]=0 # no pitch influence on states ! <<<< Important since value may only be valid around a given pitch
            if ExtraZeros:
                B.iloc[2,1]=0 # No Qgen influence on qtdot
                B.iloc[3,0]=0 # No thrust influence on psi
                D.iloc[0,1]=0  # No Qgen influence on IMU
            D.iloc[0,2]=0  # No pitch influences on IMU

            C.iloc[3,:]=0 # No states influence pitch
            C.iloc[2,3]=0 # No influence of psi on Qgen !<<< Important
    else:
        raise NotImplementedError('Model {} shape {}'.format(model,A.shape))

    # ---
    try:
        D['Qgen_[Nm]']['Qgen_[Nm]']=1
    except:
        pass

    return A,B,C,D,M
