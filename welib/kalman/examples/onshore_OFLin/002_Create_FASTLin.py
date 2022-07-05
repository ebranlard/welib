""" 
Creates a fast linearized model from a set of "OpenFAST .lin files"

"""

import numpy as np
import pandas as pd
import os
import welib.weio as weio # https://github.com/ebranlard/weio

# Local
# from helper_functions import *
from wtDigiTwin.fast.fastlinfiles import FASTLin
from wtDigiTwin.tools.clean_exceptions import *

# --- Script parameters
bLoadAll=True
# bLoadAll=False
bExport=True
Method='mean'
Method='ws'
WS=1000 # Hack to avoid overriding files from repo


nDOF=2
sWT='NREL5MW'
LinDump   = '{:s}_Lin_2DOF_FL_dump.dat'.format(sWT)
LinFolder = '../../_data/{:s}_Linearizations/'.format(sWT)
OutFile   = '../../_data/{:s}_{:d}DOF_ABCD_{:s}.dat'.format(sWT,nDOF,Method)

# --- Helper functions
import pickle
def save(dat,filename):
    with open(filename,'wb') as f:
        pickle.dump(dat,f)
def load(filename):
    with open(filename,'rb') as f:
        dat=pickle.load(f)
    return dat
def findsignal(sL, pat):
    I=[]
    for i,s in enumerate(sL):
        if s.find(pat)>=0:
            print(i,s)
            I.append(i)
    return I

# --- Read all the linearizations from FAST
if bLoadAll:
    FL=FASTLin(LinFolder, prefix='ws_')
    FL.save(LinDump)
else:
    FL = load(LinDump.replace('.dat','.dat'))

# FL= FASTLin(LinFolder, prefix='ws_4', nLin=1)

# --- Mean matrices
if Method=='mean':
    A,B,C,D,M = FL.stats('A')[0], FL.stats('B')[0], FL.stats('C')[0], FL.stats('D')[0],FL.stats('M')[0]
elif Method=='ws':
    OutFile   = '../../_data/{:s}_{:d}DOF_ABCD_{:s}{}.dat'.format(sWT,nDOF,Method,WS)
    try:
        A,B,C,D,M = FL.stats('A',WS=[WS])[0], FL.stats('B',WS=[WS])[0], FL.stats('C',WS=[WS])[0], FL.stats('D',WS=[WS])[0],FL.stats('M',WS=[WS])[0]
    except:
        A,B,C,D = FL.stats('A',WS=[WS])[0], FL.stats('B',WS=[WS])[0], FL.stats('C',WS=[WS])[0], FL.stats('D',WS=[WS])[0]
    pass
else:
    raise NotImplementedError()

# A=FL.OP_Data[0].Data[0]['A']
# B=FL.OP_Data[0].Data[0]['B']
# C=FL.OP_Data[0].Data[0]['C']
# D=FL.OP_Data[0].Data[0]['D']


# --- Extracting relevant model
if nDOF==2:
    sX_DOF  = ['qt1FA_[m]', 'psi_gen_[rad]', 'd_qt1FA_[m/s]', 'd_psi_gen_[rad/s]']
    sED_DOF = ['7_TwFADOF1' ,'13_GeAz']
    sU_DOF  = ['HubFxN1_[N]','Qgen_[Nm]','PitchColl_[rad]']
    sY_DOF  = ['NcIMUTAxs_[m/s^2]', 'RotSpeed_[rpm]','SvDGenTq_[kNm]', 'BPitch1_[deg]']
#     sY_DOF  = ['NcIMUTAxs_[m/s^2]', 'RotSpeed_[rpm]','SvDGenTrq_[Nm]', 'BPitch1_[deg]']

# SvDGenTq_[kNm]
# SvDGenPwr_[kW]
# RotSpeed_[rpm]
# BPitch1_[deg]
# >>> findsignal(sY,'Speed')

sX, sU, sY = FL.xdescr,FL.udescr,FL.ydescr
try:
    sED_DOF = FL.EDdescr
except:
    sED = sED_DOF
# Renaming signals
I=findsignal(sY, 'SvDBlPitchCom_[rad]')
if len(I)==3:
    sY[I[0]]='SvDBlPitchCom1_[rad]'
    sY[I[1]]='SvDBlPitchCom2_[rad]'
    sY[I[2]]='SvDBlPitchCom3_[rad]'
# findsignal(sY,'Speed')

# IDOFM = np.array([i for i,s in enumerate(sED) if s in sED_DOF])
IDOFU = np.array([list(sU).index(s) for s in sU_DOF])
IDOFY = np.array([list(sY).index(s) for s in sY_DOF])
try:
    IDOFE = np.array([list(sED).index(s) for s in sED_DOF])
except:
    pass

# ---
Ar = A
Br = B[:,IDOFU]
Cr = C[IDOFY,:]
Dr = D[np.ix_(IDOFY,IDOFU)]

try:
    Mr = M[np.ix_(IDOFE,IDOFE)]
except:
    Mr=np.zeros((2,2))
# Br[2,1]=0
# Br[3,0]=0
# Br[:,2]=0 # No pitch influence

Ar = pd.DataFrame(data = Ar, index=sX, columns=sX)
Br = pd.DataFrame(data = Br, index=sX, columns=sU_DOF)
Cr = pd.DataFrame(data = Cr, index=sY_DOF, columns=sX)
Dr = pd.DataFrame(data = Dr, index=sY_DOF, columns=sU_DOF)
Mr = pd.DataFrame(data = Mr, index=sED_DOF, columns=sED_DOF)

print('A\n',Ar)
print('B\n',Br)
print('C\n',Cr)
print('D\n',Dr)
print('M\n',Mr)

print('>>> Saving to', OutFile)
save((Ar,Br,Cr,Dr,Mr), OutFile)
