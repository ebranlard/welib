
import os
import pandas as pd
import numpy as np
from pandas import ExcelFile

try:
    from . import fastlib
except:
    from welib.fastlib import fastlib

def campbell(Cases, mainFst, tStart, nPerPeriod, workDir, toolboxDir, matlabExe, fastExe,  
             baseDict=None, generateInputs=True, runFast=True, runMBC=True, prefix='', suffix='', sortedSuffix=None, ylim=None):
    """ 
    Wrapper function to perform a Campbell diagram study
       see: writeFASTLinInputs, matlabMBC, and postproMBC for more description of inputs.

    INPUTS:
      - Cases        pandas DataFrame with WS, RPM, Pitch, and potentially gen torque and tower top FA
      - mainFst      Main file, used as a template
      - tStart       Time after which linearization is done (need to reach periodic steady state)
      - nPerPeriod   Number of linearizations per revolution
      - workDir      Output folder for FAST input files and linearization (will be created)
      - fastExe      Full path to a FAST exe (and dll)
      - toolboxDir   path to matlab-toolbox
      - matlabExe    path the matlab or octave exe
      - prefix       prefix used for filenames, only needed to distinguish different sets of simulations
      - suffix       suffix used for filenames, idem
      - sortedSuffix use a separate file where IDs have been sorted
      - runFast      Logical to specify whether to run the simulations or not
    """
    # --- Generating input files
    if generateInputs:
        WS, RPM, Pitch = Cases['WindSpeed_[m/s]'], Cases['Omega_[RPM]'], Cases['Pitch_[deg]']
        try:
            TT = Cases['TTDispFA_[m]'] # Optional
        except:
            TT = None
        #GenTorq = Cases['GenTrq_[kNm]'] # TODO
        # Generate input files
        fastfiles= writeFASTLinInputs(mainFst, workDir, WS, RPM, Pitch, TT=TT, nPerPeriod=nPerPeriod, baseDict=baseDict, tStart=tStart, prefix=prefix, suffix=suffix)
        # Create a batch script (optional)
        fastlib.writeBatch(os.path.join(workDir,'_RUN_ALL.bat'),fastfiles,fastExe=fastExe)

    # --- Run the simulations
    if runFast:
        fastlib.run_fastfiles(fastfiles, fastExe=fastExe, parallel=True, showOutputs=True, nCores=3)

    # --- Postprocess linearization outputs (MBC + modes ID)
    if runMBC:
        outBase = matlabMBC(workDir, toolboxDir, matlabExe, prefix=prefix, suffix=suffix)
    OP, Freq, Damp, UnMapped, _ = postproMBC(csvBase=os.path.join(workDir,prefix+suffix), sortedSuffix=sortedSuffix)

    # ---  Plot Campbell
    fig, axes = plotCampbell(OP, Freq, Damp, sx='WS_[m/s]', UnMapped=UnMapped, ylim=ylim)
    #fig, axes = fastlib.plotCampbell(Freq, Damp, sx='RPM_[rpm]', UnMapped=UnMapped)
    #  fig.savefig('{:s}CampbellWS.png'.format(suffix))
    return OP, Freq, Damp, fig


def writeFASTLinInputs(main_fst, workDir, WS, RPM, Pitch, TT=None, 
        nPerPeriod=36, baseDict=None, tStart=100,
        LinInputs=0, LinOutputs=0,
        prefix='', suffix=''):
    """
    Write FAST inputs files for linearization, to a given directory `workDir`.


    INPUTS:
      - main_fst: path to an existing  .fst file
                  This file (and the ones it refers to) will be used as templates.
                  Values of the templates can be modified using `baseDict`.
                  The parent directory of the main fst file will be copied to `workDir`.
      - workDir: directory (will be created) where the simulation files will be generated

      - operating conditions: (1d arrays, all of the same size)
        WS:    wind speed [m/s]
        RPM:   rotational speed [rpm]
        Pitch: pitch angle [deg]
        TT:    tower tower position [m]. If not provided, 0s are assumed
      - nPerPeriod : number of linearization points per rotation (usually 12 or 36)
      - baseDict : a dictionary of inputs files keys to be applied to all simulations
                   Ignored if not provided.
                   e.g. baseDict={'DT':0.01, 'EDFile|ShftTilt':-5, 'InflowFile|PLexp':0.0}
                   see templateReplaceGeneral.
      - tStart: time at which the linearization will start. 
                When triming option is not available, this needs to be sufficiently large for
                the rotor to reach an equilibrium
      - LinInputs:  linearize wrt. inputs (see OpenFAST documentation). {0,1,2, default:1}
      - LinOutputs: linearize wrt. outputs (see OpenFAST documentation). {0,1, default:0}

    OUTPUTS:
       - list of fst files created
    """
    # --- Optional values
    if TT is None:
        TT= [0]*len(WS)
    if baseDict is None:
        baseDict=dict()

    # --- Checking main fst file
    fst = fastlib.FASTInFile(main_fst)
    hasTrim = 'TrimCase' in fst.keys()
    if fst['CompServo']==1:
        print('[WARN] For now linearization is done without controller.')
        baseDict['CompServo']=0

    # --- Generating list of parameters that vary based on the operating conditions provided
    PARAMS     = []
    for i,(ws,rpm,pitch,tt) in enumerate(zip(WS,RPM,Pitch,TT)):
        # Determine linearization times based on RPM and nPerPeriod
        Omega = rpm/60*2*np.pi
        if abs(Omega)<0.001:
            LinTimes = [tStart]
            Tmax     = tStart+1
        else:
            T = 2*np.pi/Omega
            LinTimes = np.linspace(tStart,tStart+T,nPerPeriod+1)[:-1]
            Tmax       = tStart+1.01*T
        # --- Creating "linDict", dictionary of changes to fast input files for linearization
        linDict=dict()
        linDict['__name__']     = prefix+'ws{:04.1f}_rpm{:05.2f}'.format(ws,rpm)+suffix
        # --- Main fst options
        linDict['TMax']         = Tmax
        linDict['TStart']       = 0
        if abs(ws)<0.001:
            linDict['CompAero']    = 0
            linDict['CompInflow']  = 0
        # --- Linearization options
        linDict['Linearize']    = True
        linDict['NLinTimes']    = len(LinTimes)
        linDict['LinTimes']     = list(LinTimes)
        linDict['OutFmt']       = '"ES16.9E2"'  # Important for decent resolution
        linDict['LinInputs']    = LinInputs     # 0: none, 1: standard, 2: to get full linearizations
        linDict['LinOutputs']   = LinOutputs    # 0: none, 1: based on outlist 
        # --- New Linearization options
        # TrimCase - Controller parameter to be trimmed {1:yaw; 2:torque; 3:pitch} [used only when CalcSteady=True]
        # TrimTol - Tolerance for the rotational speed convergence [>eps] [used only when CalcSteady=True]
        # TrimGain - Proportional gain for the rotational speed error (rad/(rad/s) or Nm/(rad/s)) [>0] [used only when CalcSteady=True]
        # Twr_Kdmp - Damping factor for the tower (N/(m/s)) [>=0] [used only when CalcSteady=True]
        # Bld_Kdmp - Damping factor for the blade (N/(m/s)) [>=0] [used only when CalcSteady=True]
        # --- Mode shape vizualization options
        if hasTrim:
            linDict['WrVTK']        = 3
            linDict['VTK_type']     = 1
            linDict['VTK_fields']   = True
            linDict['VTK_fps']      = 30
        else:
            linDict['WrVTK']        = 0
        # --- Aero options
        linDict['AeroFile|WakeMod']   = 1 # Needed for linearization
        linDict['AeroFile|AFAeroMod'] = 1 # Needed for linearization
        linDict['AeroFile|FrozenWake'] = True # Needed for linearization
        # --- Inflow options
        linDict['InflowFile|WindType'] = 1
        linDict['InflowFile|HWindSpeed'] = ws
        # --- ElastoDyn options
        linDict['EDFile|BlPitch(1)'] = pitch
        linDict['EDFile|BlPitch(2)'] = pitch
        linDict['EDFile|BlPitch(3)'] = pitch
        linDict['EDFile|RotSpeed']   = rpm
        linDict['EDFile|TTDspFA']    = tt
        # --- Servo options

        # --- Merging linDict dictionary with user override inputs
        for k,v in baseDict.items():
            if k in linDict and v != linDict[k]:
                print('Overriding key {} with value {} (previous value {})'.format(k,v,linDict[k]))
            linDict[k]=v

        PARAMS.append(linDict)

    # --- Generating all files in a workDir
    refDir   = os.path.dirname(main_fst)
    main_file = os.path.basename(main_fst)
    fastfiles = fastlib.templateReplace(PARAMS,refDir,workDir=workDir,RemoveRefSubFiles=True,main_file=main_file)

    return fastfiles

def matlabMBC(workDir, toolboxDir, matlabExe, prefix='', suffix='', outputFormat='csv'):
    """ 
    Run the matlab script `campbellFolderPostPro` located in the matlab-toolbox.
    This matlab script generates a set of csv files, or one excel file.
    Lin files are looked for with the pattern: workDir+prefix+'*'+suffix+'*.lin'

    INPUTS:
      - workDir   : directory where .lin files are to be found
      - toolboxDir: directory where the matlab toolbox is located
      - matlabExe : path to MATLAB or OCTAVE executable
      - prefix/suffix: strings such that the files looked for are: [folder prefix '*' suffix '*.lin']
      - ouputFormat : 'csv' or 'xls', choose between output as CSV files or one Excel file
    """
    workDir    = os.path.abspath(workDir).replace('\\','/')
    toolboxDir = os.path.abspath(toolboxDir).replace('\\','/')

    if matlabExe.lower().find('matlab')>=0:
        args=['-nodisplay','-nojvm']
        matlabRun = '{} -nodisplay -nojvm -r'.format(matlabExe)
    else: # here we assume Octave
        args=['--eval']
        matlabRun = '{} --eval '.format(matlabExe)
    matlabCmd = 'addpath(genpath("{}")); campbellFolderPostPro("{}","{}","{}", "{}"); quit'.format(toolboxDir, workDir, prefix, suffix, outputFormat)
    args+=[matlabCmd]

    # system call to matlab
    p=fastlib.run_cmd(args, matlabExe, wait=True, showOutputs=True, showCommand=True)
    if p.returncode==1:
        raise Exception('Matlab command failed to run, non 0 exit code return')
    # checks
    outBase = os.path.join(workDir,prefix+suffix)
    if outputFormat.lower()=='csv':
        outName = outBase+'CampbellModesID.csv'
    else:
        outName = outBase+'CampbellDataSummary.xlsx'
    if not os.path.exists(outName):
        raise Exception('Matlab command failed to run, main xlsx or csv file was not created')
    print('[ OK ] Matlab command ran successfully')
    return outBase

def postproMBC(xlsFile=None,csvBase=None, sortedSuffix=None):
    """ 
    Generate Cambell diagram data from an xls file, or a set of csv files
    INPUTS:
      - xlsFile: path to an excel file, or, basename for a set of csv files generated by campbellFolderPostPro
      - csvBase: basename for a set of csv files generated by campbellFolderPostPro
            the csv files are asssumed to be named:
             - csvBase + 'CampbellOP.csv'
             - csvBase + 'CampbellModesID.csv'
             - csvBase + 'CampbellPointI.csv' I=1..n_Op
      - sortedSuffix: when present, the following file will be read instead of the 'ModesID' mentioned above
             - csvBase + 'CampbellModesID' + sortedSuffix + '.csv'
    OUTPUTS:
        - Freq: dataframe with columns [WS, RPM, Freq_Mode1,.., Freq_ModeN] for the N identified modes
        - Damp: dataframe with columns [WS, RPM, Damp_Mode1,.., Damp_ModeN] (damping ratios), for the N identified modes
        - UnMapped: dataframe with columns [WS, RPM, Freq, Damp] for all unidenfied modes
        - ModesData: low-level data, dictionaries for each OP with frequencies and damping
    """
    rpmSweep=False
    if sortedSuffix is None:
        sortedSuffix=''

    if xlsFile is not None:
        # --- Excel file reading
        OPFileName=xlsFile;
        IDFileName=xlsFile;
        sheets=dict()
        # Reading all sheets
        xls = pd.ExcelFile(xlsFile)
        dfs = {}
        for sheet_name in xls.sheet_names:
            # Reading sheet
            df = xls.parse(sheet_name, header=None)
            if df.shape[0]>0:
                sheets[sheet_name]=df
        OP   = sheets['OP']
        WS   = OP.iloc[1,1:].values
        RPM  = OP.iloc[2,1:].values
        ID   = None
        if 'WS_ModesID' in sheets.keys():
            ID = sheets['WS_ModesID']
            rpmSweep  = False
            sweepVar  = WS
            sweepUnit = 'mps'
        elif 'ModesID' in sheets.keys():
            ID = sheets['ModesID']
            rpmSweep  = True
            sweepVar  = RPM
            sweepUnit = 'rpm'
        else:
            raise Exception('Mode identification sheet not found in excel file')
        # Storing data for each points, we try a bunch of keys since matlab script uses a loose num2str for now
        Points=dict()
        for i,v in enumerate(sweepVar):
            keys=[('{:.'+str(ires)+'f} {:s}').format(v,sweepUnit) for ires in [0,1,2,3,4]]
            for k in keys:
                try:
                    Points[i] = sheets[k]
                except:
                    pass
            if i not in Points.keys():
                raise Exception('Couldnf find sheet for operating point {:d}'.format(i))
    else:
        # --- csv file reading
        OPFileName=csvBase+'CampbellOP.csv'
        IDFileName=csvBase+'CampbellModesID'+sortedSuffix+'.csv'

        OP      = pd.read_csv(OPFileName, sep = ',')
        ID      = pd.read_csv(IDFileName, sep = ',',header=None)
        nCol    = OP.shape[1]-1
        naCount = OP.isna().sum(axis=1)
        WS     = OP.iloc[0,1:].values
        RPM    = OP.iloc[1,1:].values
        del OP
        # Storing data for each points into a dict
        Points=dict()
        for i,v in enumerate(WS):
            Points[i] = pd.read_csv(csvBase+'CampbellPoint{:d}.csv'.format(i+1), sep = ',', header=None)
    # --- Mode Identification
    ID.iloc[:,0].fillna('Unknown', inplace=True) # replace nan
    ModeNames = ID.iloc[2: ,0].values
    ModeIDs   = ID.iloc[2: ,1:].values
    nModesIDd = len(ModeNames) 

    if ModeIDs.shape[1]!=len(WS):
        print('OP Windspeed:',WS)
        raise Exception('Inconsistent number of operating points between OP ({} points) and ID ({} points) data.\nOP filename: {}\nID filename: {}\n'.format(ModeIDs.shape[1], len(WS), OPFileName, IDFileName))

    # --- Extract Frequencies and Damping from Point table
    ModeData=[]
    ioff=0
    coff=0
    for i,ws in enumerate(WS):
        P = Points[i]
        opData = dict()
        opData['Fnat']  = P.iloc[1+ioff, 1::5+coff].values[:].astype(float) # natural frequencies
        opData['Fdmp']  = P.iloc[2+ioff, 1::5+coff].values[:].astype(float) # damped frequencies
        opData['Damps'] = P.iloc[3+ioff, 1::5+coff].values[:].astype(float) # Damping values
        ModeData.append(opData)

    # --- Creating a cleaner table of operating points
    OP = pd.DataFrame(np.nan, index=np.arange(len(WS)), columns=['WS_[m/s]', 'RotSpeed_[rpm]'])
    OP['WS_[m/s]']        = WS
    OP['RotSpeed_[rpm]'] = RPM
    print(WS)
    print(RPM)
    print(OP.shape)

    UnMapped_WS   = []
    UnMapped_RPM  = []
    UnMapped_Freq = []
    UnMapped_Damp = []
    # --- Unidentified modes, before "nModes"
    for iOP,(ws,rpm) in enumerate(zip(WS,RPM)):
        m = ModeData[iOP]
        nModesMax = len(m['Fnat'])
        nModes = min(nModesMax, nModesIDd) # somehow sometimes we have 15 modes IDd but only 14 in the ModeData..
        Indices = (np.asarray(ModeIDs[:,iOP])-1).astype(int)
        IndicesMissing = [i for i in np.arange(nModes) if i not in Indices]
        f   = np.asarray([m['Fnat'] [iiMode] for iiMode in IndicesMissing])
        d   = np.asarray([m['Damps'][iiMode] for iiMode in IndicesMissing])
        ws  = np.asarray([ws]*len(f))
        rpm = np.asarray([rpm]*len(f))
        UnMapped_Freq = np.concatenate((UnMapped_Freq, f))
        UnMapped_Damp = np.concatenate((UnMapped_Damp, d))
        UnMapped_WS   = np.concatenate((UnMapped_WS, ws))
        UnMapped_RPM  = np.concatenate((UnMapped_RPM, rpm))
    # --- Unidentified modes, beyond "nModes"
    for m,ws,rpm in zip(ModeData,WS,RPM):
        f   = m['Fnat'][nModesIDd:]
        d   = m['Damps'][nModesIDd:]
        ws  = np.asarray([ws]*len(f))
        rpm = np.asarray([rpm]*len(f))
        UnMapped_Freq = np.concatenate((UnMapped_Freq, f))
        UnMapped_Damp = np.concatenate((UnMapped_Damp, d))
        UnMapped_WS   = np.concatenate((UnMapped_WS, ws))
        UnMapped_RPM  = np.concatenate((UnMapped_RPM, rpm))

    # --- Put identified modes into a more convenient form
    cols=[ m.split('-')[0].strip().replace(' ','_') for m in ModeNames]
    cols=[v + str(cols[:i].count(v) + 1) if cols.count(v) > 1 else v for i, v in enumerate(cols)]
    Freq = pd.DataFrame(np.nan, index=np.arange(len(WS)), columns=cols)
    Damp = pd.DataFrame(np.nan, index=np.arange(len(WS)), columns=cols)
    for iMode in np.arange(nModesIDd):
        ModeIndices = (np.asarray(ModeIDs[iMode,:])-1).astype(int)
        ModeName= ModeNames[iMode].replace('_',' ')
        if ModeName.find('(not shown)')>0:
            f   = np.asarray([m['Fnat'][iiMode]  for m,iiMode in zip(ModeData,ModeIndices) if iiMode>=0 ])
            d   = np.asarray([m['Damps'][iiMode] for m,iiMode in zip(ModeData,ModeIndices) if iiMode>=0 ])
            ws  = np.asarray([ws  for ws,iiMode in zip(WS,ModeIndices)   if iiMode>=0 ])
            rpm = np.asarray([rpm for rpm,iiMode in zip(RPM,ModeIndices) if iiMode>=0 ])
            UnMapped_Freq = np.concatenate((UnMapped_Freq, f))
            UnMapped_Damp = np.concatenate((UnMapped_Damp, d))
            UnMapped_WS   = np.concatenate((UnMapped_WS, ws))
            UnMapped_RPM  = np.concatenate((UnMapped_RPM, rpm))
        else:
            if all(ModeIndices==-1):
                print('Skipping mode number ',iMode)
            else:
                f=np.asarray([m['Fnat'] [iiMode] if iiMode>=0 else np.nan for m,iiMode in zip(ModeData,ModeIndices)])
                d=np.asarray([m['Damps'][iiMode] if iiMode>=0 else np.nan for m,iiMode in zip(ModeData,ModeIndices)])
                Freq.iloc[:, iMode]=f
                Damp.iloc[:, iMode]=d
    #  Removing modes that are full nan (not_shown ones)
    # NOTE: damgerous since OP is not part of it anymore
    # Freq.dropna(how='all',axis=0,inplace=True)
    # Freq.dropna(how='all',axis=1,inplace=True)
    # Damp.dropna(how='all',axis=0,inplace=True)
    # Damp.dropna(how='all',axis=1,inplace=True)

    # --- UnMapped modes into a dataframe
    M = np.column_stack((UnMapped_WS, UnMapped_RPM, UnMapped_Freq, UnMapped_Damp))
    UnMapped = pd.DataFrame(data=M, columns=['WS_[m/s]','RPM_[rpm]','Freq_[Hz]','Damping_[-]'])

    return OP, Freq, Damp, UnMapped, ModeData

def plotCampbell(OP, Freq, Damp, sx='WS_[m/s]', UnMapped=None, fig=None, axes=None, ylim=None):
    """ Plot Campbell data as returned by postproMBC """
    import matplotlib.pyplot as plt
    FullLineStyles = [':', '-', '-+', '-o', '-^', '-s', '--x', '--d', '-.', '-v', '-+', ':o', ':^', ':s', ':x', ':d', ':.', '--','--+','--o','--^','--s','--x','--d','--.'];
    Markers    = ['', '+', 'o', '^', 's', 'd', 'x', '.']
    LineStyles = ['-', ':', '-.', '--'];
    Colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    MW_Light_Blue    = np.array([114,147,203])/255.
    MW_Light_Orange  = np.array([225,151,76])/255.
    MW_Light_Green   = np.array([132,186,91])/255.
    MW_Light_Red     = np.array([211,94,96])/255.
    MW_Light_Gray    = np.array([128,133,133])/255.
    MW_Light_Purple  = np.array([144,103,167])/255.
    MW_Light_DarkRed = np.array([171,104,87])/255.
    MW_Light_Kaki    = np.array([204,194,16])/255.
    MW_Blue     =     np.array([57,106,177])/255.
    MW_Orange   =     np.array([218,124,48])/255.
    MW_Green    =     np.array([62,150,81])/255.
    MW_Red      =     np.array([204,37,41])/255.
    MW_Gray     =     np.array([83,81,84])/255.
    MW_Purple   =     np.array([107,76,154])/255.
    MW_DarkRed  =     np.array([146,36,40])/255.
    MW_Kaki     =     np.array([148,139,61])/255.

    def modeStyle(i, lbl):
        lbl=lbl.lower().replace('_',' ')
        ms = 4
        c  = Colors[np.mod(i,len(Colors))]
        ls = LineStyles[np.mod(int(i/len(Markers)),len(LineStyles))]
        mk = Markers[np.mod(i,len(Markers))]
        # Color
        if any([s in lbl for s in ['1st tower']]):
            c=MW_Blue
        elif any([s in lbl for s in ['2nd tower']]):
            c=MW_Light_Blue
        elif any([s in lbl for s in ['1st blade edge','drivetrain']]):
            c=MW_Red
        elif any([s in lbl for s in ['1st blade flap']]):
            c=MW_Green
        elif any([s in lbl for s in ['2nd blade flap']]):
            c=MW_Light_Green
        elif any([s in lbl for s in ['2nd blade edge']]):
            c=MW_Light_Red
        # Line style
        if any([s in lbl for s in ['tower fa','collective','drivetrain']]):
            ls='-'
        elif any([s in lbl for s in ['tower ss','regressive']]):
            ls='--'
        elif any([s in lbl for s in ['tower ss','progressive']]):
            ls='-.'
        # Marker
        if any([s in lbl for s in ['collective']]):
            mk='2'; ms=8
        elif any([s in lbl for s in ['blade','tower','drivetrain']]):
            mk=''; 
        return c, ls, ms, mk

    # Init figure
    if fig is None:
        fig,axes_ = plt.subplots(1,2)
        fig.set_size_inches(13,7.0,forward=True) # default is (6.4,4.8)
        fig.subplots_adjust(top=0.78,bottom=0.11,left=0.04,right=0.98,hspace=0.06,wspace=0.16)
    if axes is None:
        axes=axes_

    # Estimating figure range
    FreqRange = [0                         , np.nanmax(Freq.iloc[:,:])*1.01]
    DampRange = [np.nanmin(Damp.iloc[:,2:]), np.nanmax(Damp.iloc[:,:])*1.01]
    if ylim is not None:
        FreqRange=ylim
    if DampRange[0]>0:
        DampRange[0]=0

    # Plot "background"


    # Plot mapped modes
    iModeValid=0
    xPlot=[]; yPlot=[]
    for iMode,lbl in enumerate(Freq.columns.values):
        if lbl.find('not_shown')>0:
            # TODO ADD TO UNMAPPED
            continue
        iModeValid+=1
        c, ls, ms, mk = modeStyle(iModeValid, lbl)
        axes[0].plot(OP[sx].values, Freq[lbl].values, ls, marker=mk, label=lbl.replace('_',' '), markersize=ms, color=c)
        axes[1].plot(OP[sx].values, Damp[lbl].values, ls, marker=mk                            , markersize=ms, color=c)
        xPlot=np.concatenate((xPlot, OP[sx].values))
        yPlot=np.concatenate((yPlot, Freq[lbl].values))

    # Unmapped modes (NOTE: plotted after to over-plot)
    if UnMapped is not None:
        axes[0].plot(UnMapped[sx].values, UnMapped['Freq_[Hz]'  ].values, '.', markersize=6, color=[0.5,0.5,0.5])
        axes[1].plot(UnMapped[sx].values, UnMapped['Damping_[-]'].values, '.', markersize=1, color=[0.5,0.5,0.5])
    # Highligh duplicates (also after)
    Points=[(x,y) for x,y in zip(xPlot,yPlot)]
    Plot = pd.Series(Points)
    for xDupl,yDupl in Plot[Plot.duplicated()]:
        axes[0].plot(xDupl,yDupl, 'o',color='r')

    axes[0].set_xlabel(sx.replace('_',' '))
    axes[1].set_xlabel(sx.replace('_',' '))
    axes[0].set_ylabel('Frequencies [Hz]')
    axes[1].set_ylabel('Damping ratios [-]')
    axes[0].legend(bbox_to_anchor=(0., 1.02, 2.16, .802), loc='lower left', ncol=4, mode="expand", borderaxespad=0.)
    axes[0].set_ylim(FreqRange)
    
    XLIM=axes[1].get_xlim()
    axes[1].plot(XLIM, [0,0],'-', color='k', lw=0.5)
    axes[1].set_xlim(XLIM)
    axes[1].set_ylim(DampRange)
    return fig, axes
