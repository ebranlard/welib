##
import numpy as np
import copy
import os
import weio

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


# Temporary hack, loading model from state file
def loadLinStateMatModel(StateFile,nDOF=2, Adapt=True, ExtraZeros=False):
    import pickle
    def load(filename):
        with open(filename,'rb') as f:
            dat=pickle.load(f)
        return dat

    try:
        (A,B,C,D,M) = load(StateFile)
    except:
        M=None
        (A,B,C,D) = load(StateFile)
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

        #A.index.values[1]='psi_rot_[rad]'
        #A.columns.values[1]='psi_rot_[rad]'
        #A.index.values[3]='d_psi_rot_[rad]'
        #A.columns.values[3]='d_psi_rot_[rad]'
        #C.columns.values[1]='psi_rot_[rad]'
        #C.columns.values[3]='d_psi_rot_[rad]'

        C.iloc[3,:]=0 # No states influence pitch
        C.iloc[2,3]=0 # No influence of psi on Qgen !<<< Important
        C.index.values[1]='RotSpeed_[rad/s]'
        D.index.values[1]='RotSpeed_[rad/s]'
        C.index.values[2]='Qgen_[Nm]'
        D.index.values[2]='Qgen_[Nm]'
        C.index.values[3]='BPitch1_[rad]'
        D.index.values[3]='BPitch1_[rad]'
        C.iloc[1,:]/=60/(2*np.pi) # RotSpeed output in radian/s
        D.iloc[1,:]/=60/(2*np.pi) # RotSpeed output in radian/s
        C.iloc[2,:]*=1000 # GenTq in Nm
        D.iloc[2,:]*=1000 # GenTq in Nm
        D['Qgen_[Nm]']['Qgen_[Nm]']=1
        C.iloc[3,:]/=180/np.pi # Pitch output in radian
        D.iloc[3,:]/=180/np.pi # Pitch output in radian

    return A,B,C,D,M
