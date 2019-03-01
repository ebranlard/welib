import os
import glob
import subprocess
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys

try:
    import weio
except:
    import sys
    sys.path.append('../../_libs/pyDatView/')
    import weio
    

FAST_EXE='openfast'
bDEBIAN=False


def readFASTOutAscii(filename):
    # Read with panda
    f = weio.FASTOutFile(filename)
    df = f.toDataFrame()
    #df=pd.read_csv(filename, sep='\t', skiprows=[0,1,2,3,4,5,7])
    #df.rename(columns=lambda x: x.strip(),inplace=True)
    return df


def createStepWind(filename,WSstep=1,WSmin=3,WSmax=25,tstep=100,dt=0.5,tmin=0,tmax=999):
    f = weio.FASTWndFile()
    Steps= np.arange(WSmin,WSmax+WSstep,WSstep)
    print(Steps)
    nCol = len(f.colNames)
    nRow = len(Steps)*2
    M = np.zeros((nRow,nCol));
    M[0,0] = tmin
    M[0,1] = WSmin
    for i,s in enumerate(Steps[:-1]):
        M[2*i+1,0] = tmin + (i+1)*tstep-dt 
        M[2*i+2,0] = tmin + (i+1)*tstep
        M[2*i+1,1] = Steps[i]
        if i<len(Steps)-1:
            M[2*i+2,1] = Steps[i+1]
        else:
            M[2*i+2,1] = Steps[-1]
    M[-1,0]= max(tmax, (len(Steps)+1)*tstep)
    M[-1,1]= WSmax
    f.data=pd.DataFrame(data=M,columns=f.colNames)
    #
    print(f.data)
    f.write(filename)
    #plt.plot(M[:,0],M[:,1])
    #plt.show()

    #print(f.toDataFrame())
    #pass
#createStepWind('test.wnd',tstep=200,WSmax=28)
# createStepWind('test.wnd',tstep=200,WSmin=5,WSmax=7,WSstep=2)

def run_fastfiles(fastfiles,parallel=True,Outputs=True,nCores=4):
    ps=[]
    iProcess=0
    for i,f in enumerate(fastfiles):
        print('Process {}/{}: {}'.format(i+1,len(fastfiles),f))
        ps.append(run_fast(f,wait=not parallel,Outputs=Outputs))
        iProcess += 1
        # waiting once we've filled the number of cores
        if parallel:
            if iProcess==nCores:
                for p in ps:
                    p.wait()
                ps=[]
                iProcess=0
    for p in ps:
        p.wait()

def run_fast(input_file,wait=True,Outputs=False):
    if not os.path.isabs(input_file):
        input_file=os.path.abspath(input_file)
    workdir  = os.path.dirname(input_file)
    basename = os.path.basename(input_file)
    if bDEBIAN:
        input_file=input_file.replace('C:','/mnt/c')
        input_file=input_file.replace('Work','work')
        input_file=input_file.replace('\\','/')
    if os.path.isabs(FAST_EXE):
        fast_exe = FAST_EXE
    else:
        if bDEBIAN:
            fast_exe = './'+FAST_EXE.strip()
        else:
            fast_exe = '.\\'+FAST_EXE.strip()
    if bDEBIAN:
        args = ['debian', 'run', 'cd '+workdir+' && '+ fast_exe +' '+basename]
    else:
#         args = ['cmd','/k \"cd '+workdir+' && '+ fast_exe +' '+basename+'\"']
        args = 'cd '+workdir+' && '+ fast_exe +' '+basename
    #print(args)
    FNULL = open(os.devnull, 'w')
    if wait:
        if Outputs:
            if bDEBIAN:
                p=subprocess.call(args)
            else:
                p=subprocess.call(args, shell=True)
        else:
            if bDEBIAN:
                p=subprocess.call(args, stdout=FNULL, stderr=subprocess.STDOUT)
            else:
                p=subprocess.call(args, stdout=FNULL, stderr=subprocess.STDOUT, shell=True)
    else:
        if Outputs:
            p=subprocess.Popen(args, shell=True)
        else:
            p=subprocess.Popen(args, stdout=FNULL, stderr=subprocess.STDOUT, shell=True)
    return p


def removeFASTOuputs(workdir):
    # Cleaning folder, just in case
    for f in glob.glob(os.path.join(workdir,'*.out')):
        os.remove(f)
    for f in glob.glob(os.path.join(workdir,'*.outb')):
        os.remove(f)
    for f in glob.glob(os.path.join(workdir,'*.ech')):
        os.remove(f)


# def detectFastFiles(workdir):
#     FstFiles=glob.glob(os.path.join(workdir,'*.fst'))+glob.glob(os.path.join(workdir,'*.FST'))
#     DatFiles=glob.glob(os.path.join(workdir,'*.dat'))+glob.glob(os.path.join(workdir,'*.DAT'))
#     Files=dict()
#     Files['Main']      = FstFiles
#     Files['Inflow']    = None
#     Files['Aero']      = None
#     Files['Tower']     = None
#     Files['Blade']     = None
#     Files['AeroBlade'] = None
#     Files['ServoDyn']  = None
#     for f in DatFiles:
#         b = os.path.basename(f).lower()
#         if b.find('inflow'):
#             Files['Inflow'] = f
#     windfile_ref = 'InflowWind.dat';
#     fastfile_ref = 'Turbine.fst';
#     elasfile_ref = 'ElastoDyn.dat';
#         remove
   
