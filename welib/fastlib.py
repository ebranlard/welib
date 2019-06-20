
# --- For cmd.py
import os
import subprocess
import multiprocessing

import collections
import glob
import pandas as pd
import numpy as np
import distutils.dir_util
from shutil import copytree, ignore_patterns, rmtree, copyfile

# --- External library for io
import weio
    

FAST_EXE='openfast'

# --------------------------------------------------------------------------------}
# ---  
# --------------------------------------------------------------------------------{
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


# --------------------------------------------------------------------------------}
# --- Tools for running FAST
# --------------------------------------------------------------------------------{
# --- START cmd.py
def run_cmds(inputfiles, exe, parallel=True, ShowOutputs=True, nCores=None, ShowCommand=True): 
    """ Run a set of simple commands of the form `exe input_file`
        By default, the commands are run in parallel
    """
    ps=[]
    iProcess=0
    if nCores is None:
        nCores=multiprocessing.cpu_count()
    if nCores<0:
        nCores=len(inputfiles)+1
    for i,f in enumerate(inputfiles):
        #print('Process {}/{}: {}'.format(i+1,len(inputfiles),f))
        ps.append(run_cmd(f, exe, wait=(not parallel), ShowOutputs=ShowOutputs, ShowCommand=ShowCommand))
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
def run_cmd(input_file, exe, wait=True, ShowOutputs=False, ShowCommand=True):
    """ Run a simple command of the form `exe input_file`  """
    # TODO TODO TODO capture Exit code!
    if not os.path.isabs(input_file):
        input_file=os.path.abspath(input_file)
    if not os.path.exists(exe):
        raise Exception('FAST Executable not found: {}'.format(exe))
    args= [exe,input_file]
    #args = 'cd '+workdir+' && '+ exe +' '+basename
    shell=False
    if ShowOutputs:
        STDOut= None
    else:
        STDOut= open(os.devnull, 'w') 
    if ShowCommand:
        print('Running: '+' '.join(args))
    if wait:
        p=subprocess.call(args , stdout=STDOut, stderr=subprocess.STDOUT, shell=shell)
    else:
        p=subprocess.Popen(args, stdout=STDOut, stderr=subprocess.STDOUT, shell=shell)
    return p
# --- END cmd.py

def run_fastfiles(fastfiles, fastExe=None, parallel=True, ShowOutputs=True, nCores=None, ShowCommand=True):
    if fastExe is None:
        fastExe=FAST_EXE
    run_cmds(fastfiles, fastExe, parallel=parallel, ShowOutputs=ShowOutputs, nCores=nCores, ShowCommand=ShowCommand)

def run_fast(input_file, fastExe=None, wait=True, ShowOutputs=False, ShowCommand=True):
    if fastExe is None:
        fastExe=FAST_EXE
    return run_cmd(input_file, fastExe, wait=wait, ShowOutputs=ShowOutputs, ShowCommand=ShowCommand)

def removeFASTOuputs(workdir):
    # Cleaning folder
    for f in glob.glob(os.path.join(workdir,'*.out')):
        os.remove(f)
    for f in glob.glob(os.path.join(workdir,'*.outb')):
        os.remove(f)
    for f in glob.glob(os.path.join(workdir,'*.ech')):
        os.remove(f)
    for f in glob.glob(os.path.join(workdir,'*.sum')):
        os.remove(f)


# --------------------------------------------------------------------------------}
# --- Template replace 
# --------------------------------------------------------------------------------{
def templateReplace(template_dir, PARAMS, workdir=None, main_file=None, name_function=None, RemoveAllowed=False, RemoveRefSubFiles=False):
    """ Replace parameters in a fast folder using a list of dictionaries where the keys are for instance:
        'FAST|DT', 'EDFile|GBRatio', 'ServoFile|GenEff'
    """
    def fileID(s):
        return s.split('|')[0]
    def basename(s):
        return os.path.splitext(os.path.basename(s))[0]
    def rebase(s,sid):
        split = os.path.splitext(os.path.basename(s))
        return os.path.join(workdir,split[0]+sid+split[1])
    def rebase_rel(s,sid):
        split = os.path.splitext(s)
        return os.path.join(workdir,split[0]+sid+split[1])

    # Default value of workdir if not provided
    if template_dir[-1]=='/'  or template_dir[-1]=='\\' :
        template_dir=template_dir[0:-1]
    if workdir is None:
        workdir=template_dir+'_Parametric'

    # Copying template folder to workdir
    if os.path.exists(workdir) and RemoveAllowed:
        rmtree(workdir)
    distutils.dir_util.copy_tree(template_dir, workdir)
    #copytree(template_dir, workdir, ignore=ignore_patterns('.git'))
    if RemoveAllowed:
        removeFASTOuputs(workdir)

    # --- Fast main file use as "master"
    if main_file is None:
        FstFiles=set(glob.glob(os.path.join(template_dir,'*.fst'))+glob.glob(os.path.join(template_dir,'*.FST')))
        if len(FstFiles)>1:
            print(FstFiles)
            raise Exception('More than one fst file found in template folder, provide `main_file` or ensure there is only one `.fst` file') 
        main_file=rebase(FstFiles.pop(),'')
    else:
        main_file=os.path.join(workdir, os.path.basename(main_file))

    # Params need to be a list
    if not isinstance(PARAMS,list):
        PARAMS=[PARAMS]

    fastfiles=[]
    # TODO: Recursive loop splitting at the pipes '|', for now only 1 level supported...
    for p in PARAMS:
        if name_function is None:
            raise NotImplementedError('')
        strID =name_function(p)
        FileTypes = set([fileID(k) for k in list(p.keys()) if k!='__index__'])

        # ---Copying main file and reading it
        fst_full = rebase(main_file,strID)
        copyfile(main_file, fst_full )
        Files=dict()
        Files['FAST']=weio.FASTInFile(fst_full)
        # --- Looping through required files and opening them
        for t in FileTypes: 
            # Doing a naive if
            # The reason is that we want to account for more complex file types in the future
            if t=='FAST':
                continue
            org_filename   = Files['FAST'][t].strip('"')
            org_filename_full =os.path.join(workdir,org_filename)
            new_filename_full = rebase_rel(org_filename,strID)
            new_filename      = os.path.relpath(new_filename_full,workdir)
            copyfile(org_filename_full, new_filename_full)
            Files['FAST'][t] = '"'+new_filename+'"'
            # Reading files
            Files[t]=weio.FASTInFile(new_filename_full)
        # --- Replacing in files
        for k,v in p.items():
            if k =='__index__':
                continue
            t,kk=k.split('|')
            Files[t][kk]=v
            #print(t+'|'+kk+'=',v)
        # --- Rewritting all files
        for t in FileTypes:
            Files[t].write()

        fastfiles.append(fst_full)
    # --- Remove extra files at the end
    if RemoveRefSubFiles:
        FST = weio.FASTInFile(main_file)
        for t in FileTypes:
            if t=='FAST':
                continue
            filename   = FST[t].strip('"')
            #fullname   = rebase(filename,'')
            fullname   = os.path.join(workdir,filename)
            os.remove(fullname)
    os.remove(main_file)

    return fastfiles


# --------------------------------------------------------------------------------}
# --- Tools for template replacement 
# --------------------------------------------------------------------------------{
def paramsSteadyAero(p=dict()):
    p['AeroFile|AFAeroMod']=1 # remove dynamic effects dynamic
    return p

def paramsNoGen(p=dict()):
    p['EDFile|GenDOF' ]  = 'False'
    return p

def paramsNoController(p=dict()):
    p['ServoFile|PCMode']   = 0;
    p['ServoFile|VSContrl'] = 0;
    p['ServoFile|YCMode']   = 0;
    return p

def paramsStiff(p=dict()):
    p['EDFile|FlapDOF1']  = 'False'
    p['EDFile|FlapDOF2']  = 'False'
    p['EDFile|EdgeDOF' ]  = 'False'
    p['EDFile|TeetDOF' ]  = 'False'
    p['EDFile|DrTrDOF' ]  = 'False'
    p['EDFile|YawDOF'  ]  = 'False'
    p['EDFile|TwFADOF1']  = 'False'
    p['EDFile|TwFADOF2']  = 'False'
    p['EDFile|TwSSDOF1']  = 'False'
    p['EDFile|TwSSDOF2']  = 'False'
    p['EDFile|PtfmSgDOF'] = 'False'
    p['EDFile|PtfmSwDOF'] = 'False'
    p['EDFile|PtfmHvDOF'] = 'False'
    p['EDFile|PtfmRDOF']  = 'False'
    p['EDFile|PtfmPDOF']  = 'False'
    p['EDFile|PtfmYDOF']  = 'False'
    return p

def paramsWS_RPM_Pitch(WS,RPM,Pitch,BaseDict=None,FlatInputs=False):
    """ """
    # --- Naming function appropriate for such parametric study
    def default_naming(p): # TODO TODO CHANGE ME
        return '_{:03d}_ws{:04.1f}_pt{:04.2f}_om{:04.2f}'.format(p['__index__'],p['InflowFile|HWindSpeed'],p['EDFile|BlPitch(1)'],p['EDFile|RotSpeed'])

    # --- Ensuring everythin is an iterator
    def iterify(x):
        if not isinstance(x, collections.Iterable): x = [x]
        return x
    WS    = iterify(WS)
    RPM   = iterify(RPM)
    Pitch = iterify(Pitch)
    # --- If inputs are not flat but different vectors to length through, we flatten them (TODO: meshgrid and ravel?)
    if not FlatInputs :
        WS_flat    = []
        Pitch_flat = []
        RPM_flat   = []
        for pitch in Pitch:
            for rpm in RPM:
                for ws in WS:
                    WS_flat.append(ws)
                    RPM_flat.append(rpm)
                    Pitch_flat.append(pitch)
    else:
        WS_flat, Pitch_flat, RPM_flat = WS, Pitch, RPM

    # --- Defining the parametric study 
    PARAMS=[]
    i=0
    for ws,rpm,pitch in zip(WS_flat,RPM_flat,Pitch_flat):
        if BaseDict is None:
            p=dict()
        else:
            p = BaseDict.copy()
        p['EDFile|RotSpeed']       = rpm
        p['InflowFile|HWindSpeed'] = ws
        p['InflowFile|WindType']   = 1 # Setting steady wind
        p['EDFile|BlPitch(1)']     = pitch
        p['EDFile|BlPitch(2)']     = pitch
        p['EDFile|BlPitch(3)']     = pitch

        p['__index__']  = i
        i=i+1
        PARAMS.append(p)
    return PARAMS, default_naming


# --------------------------------------------------------------------------------}
# --- Tools for PostProcessing several simulations
# --------------------------------------------------------------------------------{
def averagePostPro(outFiles,TimeAvgWindow=10,ColMap=None,ColKeep=None,ColSort=None):
    """ Opens a list of FAST output files, perform average of its signals and return a panda dataframe
    For now, the scripts only computes the mean between the End of the simulation time and End-TimeAvgWindow . In the future more options will be provided.
    The script only computes the mean for now. Other stats will be added

    `TimeAvgWindow`: Time used to perform an average of the value. 
                     The average is done between the last time step tEnd, and tEnd-TimeAvgWindow
    `ColMap` :  dictionary where the key is the new column name, and v the old column name.
                Default: None, output is not sorted
                NOTE: the mapping is done before sorting and `ColKeep` is applied
                ColMap = {'WS':Wind1VelX_[m/s], 'RPM': 'RotSpeed_[rpm]'}
    `ColKeep` : List of strings corresponding to the signals to analyse. 
                Default: None, all columns are analysed
                Example: ColKeep=['RotSpeed_[rpm]','BldPitch1_[deg]','RtAeroCp_[-]']
                     or: ColKeep=list(ColMap.keys())
    `ColSort` : string 
                Default: None, output is not sorted
                Example:  ColSort='RotSpeed_[rpm]'
    """
    def renameCol(x):
        for k,v in ColMap.items():
            if x==v:
                return k
        return x
    result=None
    for i,f in enumerate(outFiles):
        df=weio.FASTOutFile(f).toDataFrame()
        # Before doing the colomn map we store the time
        time = df['Time_[s]'].values
        # Column mapping
        if ColMap is not None:
            df.rename(columns=renameCol,inplace=True)
        if ColKeep is not None:
            df=df[ColKeep]
        ## Defining a window for stats (start time and end time)
        tEnd       = time[-1]
        tStart     = tEnd-TimeAvgWindow
        IWindow    = np.where((time>=tStart) & (time<=tEnd))[0]
        iStart = IWindow[0]
        iEnd   = IWindow[-1]
        ## Absolute and relative differences at window extremities
        #DeltaValuesAbs=(df.iloc[iEnd]-df.iloc[iStart]).abs()
        #DeltaValuesRel=(df.iloc[iEnd]-df.iloc[iStart]).abs()/df.iloc[iEnd]
        #EndValues=df.iloc[iEnd]
        ## Stats values during window
        MeanValues = pd.DataFrame(df.iloc[IWindow].mean()).transpose()
        #StdValues  = df.iloc[IWindow].std()
        if i==0:
            result = MeanValues.copy()
        else:
            result=result.append(MeanValues, ignore_index=True)
    if ColSort is not None:
        # Sorting 
        result.sort_values([ColSort],inplace=True,ascending=True)
        result.reset_index(drop=True,inplace=True) 
    return result 

# --------------------------------------------------------------------------------}
# --- Tools for typical wind turbine study 
# --------------------------------------------------------------------------------{
def CPCT_LambdaPitch(refdir,main_fastfile,Lambda=None,Pitch=np.linspace(-10,40,5),WS=None,Omega=None, # operating conditions
          TMax=20,bStiff=True,bNoGen=True,bSteadyAero=True, # simulation options
          fastExe=None,ShowOutputs=True,nCores=4): # execution options
    """ Computes CP and CT as function of tip speed ratio (lambda) and pitch.
    There are two main ways to define the inputs:
      - Option 1: provide Lambda and Pitch (deg)
      - Option 2: provide WS (m/s), Omega (in rpm) and Pitch (deg), in which case len(WS)==len(Omega)
    """

    WS_default=5 # If user does not provide a wind speed vector, wind speed used

    # --- Reading main fast file to get rotor radius 
    fst = weio.FASTInFile(os.path.join(refdir,main_fastfile))
    ed  = weio.FASTInFile(os.path.join(refdir,fst['EDFile'].replace('"','')))
    R = ed['TipRad']

    # --- Making sure we have 
    if (Omega is not None):
        if (Lambda is not None):
            WS = np.ones(Omega.shape)*WS_default
        elif (WS is not None):
            if len(WS)!=len(Omega):
                raise Exception('When providing Omega and WS, both vectors should have the same dimension')
        else:
            WS = np.ones(Omega.shape)*WS_default
    else:
        Omega = WS_default * Lambda/R*60/(2*np.pi) # TODO, use more realistic combinations of WS and Omega
        WS    = np.ones(Omega.shape)*WS_default


    # --- Defining flat vectors of operating conditions
    WS_flat    = []
    RPM_flat   = []
    Pitch_flat = []
    for pitch in Pitch:
        for (rpm,ws) in zip(Omega,WS):
            WS_flat.append(ws)
            RPM_flat.append(rpm)
            Pitch_flat.append(pitch)
    # --- Setting up default options
    BaseDict={'FAST|TMax': TMax, 'FAST|DT': 0.01, 'FAST|DT_Out': 0.1}
    BaseDict = paramsNoController(BaseDict)
    if bStiff:
        BaseDict = paramsStiff(BaseDict)
    if bNoGen:
        BaseDict = paramsNoGen(BaseDict)
    if bSteadyAero:
        BaseDict = paramsSteadyAero(BaseDict)

    # --- Creating set of parameters to be changed
    # TODO: verify that RtAeroCp and RtAeroCt are present in AeroDyn outlist
    PARAMS,naming = paramsWS_RPM_Pitch(WS_flat,RPM_flat,Pitch_flat,BaseDict=BaseDict, FlatInputs=True)

    # --- Generating all files in a workdir
    workdir = refdir.strip('/').strip('\\')+'_CPLambdaPitch'
    print('>>> Generating inputs files in {}'.format(workdir))
    fastFiles=templateReplace(refdir,PARAMS,workdir=workdir,name_function=naming,RemoveRefSubFiles=True,RemoveAllowed=True,main_file=main_fastfile)

    # --- Running fast simulations
    print('>>> Running {} simulations...'.format(len(fastFiles)))
    run_fastfiles(fastFiles, ShowOutputs=ShowOutputs, fastExe=fastExe, nCores=nCores)

    # --- Postpro - Computing averages at the end of the simluation
    print('>>> Postprocessing...')
    outFiles = [os.path.splitext(f)[0]+'.outb' for f in fastFiles]
    # outFiles = glob.glob(os.path.join(workdir,'*.outb'))
    ColKeepStats  = ['RotSpeed_[rpm]','BldPitch1_[deg]','RtAeroCp_[-]','RtAeroCt_[-]','Wind1VelX_[m/s]']
    result = averagePostPro(outFiles,TimeAvgWindow=5,ColKeep=ColKeepStats,ColSort='RotSpeed_[rpm]')
    # print(result)        

    # --- Adding lambda, sorting and keeping only few columns
    result['lambda_[-]'] = result['RotSpeed_[rpm]']*R*2*np.pi/60/result['Wind1VelX_[m/s]']
    result.sort_values(['lambda_[-]','BldPitch1_[deg]'],ascending=[True,True],inplace=True)
    ColKeepFinal=['lambda_[-]','BldPitch1_[deg]','RtAeroCp_[-]','RtAeroCt_[-]']
    result=result[ColKeepFinal]
    print('>>> Done')

    #  --- Converting to a matrices
    CP = result['RtAeroCp_[-]'].values
    CT = result['RtAeroCt_[-]'].values
    MCP =CP.reshape((len(Lambda),len(Pitch)))
    MCT =CT.reshape((len(Lambda),len(Pitch)))
    LAMBDA, PITCH = np.meshgrid(Lambda, Pitch)
    #  --- CP max
    i,j = np.unravel_index(MCP.argmax(), MCP.shape)
    MaxVal={'CP_max':MCP[i,j],'lambda_opt':LAMBDA[j,i],'pitch_opt':PITCH[j,i]}

    return  MCP,MCT,Lambda,Pitch,MaxVal,result


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
   


if __name__=='__main__':
    pass
    # --- Test of templateReplace
    def naming(p):
        return '_ws_'+str(p['InflowFile|HWindSpeed'])
    PARAMS                          = {}
    PARAMS['FAST|TMax']             = 10
    PARAMS['FAST|DT']               = 0.01
    PARAMS['FAST|DT_Out']           = 0.1
    PARAMS['EDFile|RotSpeed']       = 100
    PARAMS['EDFile|BlPitch(1)']     = 1
    PARAMS['EDFile|GBoxEff']        = 0.92
    PARAMS['ServoFile|VS_Rgn2K']    = 0.00038245
    PARAMS['ServoFile|GenEff']      = 0.95
    PARAMS['InflowFile|HWindSpeed'] = 8
    templateReplace(ref_dir,PARAMS,name_function=naming,RemoveRefSubFiles=True)

