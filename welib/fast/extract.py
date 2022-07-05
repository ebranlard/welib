""" 

"""
import os
import numpy as np
import pandas as pd
import tempfile
from numpy.linalg import inv
# Local
import welib.weio.fast_input_file as fi
import welib.weio.fast_input_deck as fd
import welib.weio.fast_linearization_file as fl
import welib.fast.runner as runner
import welib.fast.postpro as postpro
from welib.fast.linmodel import matToSIunits, subMat
from welib.yams.utils import identifyRigidBodyMM, translateInertiaMatrixFromCOG



def extractRNAInertia(fstFile, fastExe=None, workDir=None, method='fastlin', nAzim=1, cleanUp=True, 
        includeYawBearingMass=True, showOutputs=False, verbose=False, EDDict=None):
    """ 
    Compute the RNA inertia from an openfast model

    INPUTS:
      - fstFile: fast file for as template. Must point to existing elastodyn files
      - fastExe: path to an OpenFAST executable  compatible with the fstFile provided.
                 The executable is used to extract the RNA information
    OPTIONAL INPUTS:
      - workdir
      - nAzim :  Number of azimuthal steps to use. If >1, the average inertia over a revolution is computed
      - cleanUp: if true, delete file created
      - includeYawBearingMass: if true, the yaw bearing mass will be included in the RNA mass.
      - showOutputs: show Outputs when running the fast executable, useful to debug
      - verbose: display some status on commands being run
      - EDDict: Additional dictionary to set values of ElastoDyn, e.g.
              EDDict={'PreCone(1)':-4, 'HubMass':0}

    OUTPUT: dictionary with keys:
       - mass: mass
       - J_TT: Inertia at tower top
       - J_G:  Inertia at RNA COG
       - CM:   Distance from tower top to RNA COG
       - azim_*: azimuthal variation of the variables above (useful if nAzim>1)



        Extra RNA from linearization with massless tower, rigid structure, free platform DOF
    TO obtain RNA inertia information using FAST linearization:
         - FST file
             - TMax=0
             - Linearize=True, CalcSteady=False, LinTime=0, NLinTimes=1, LinOutputs=1, LinInputs=2
         - ED file: 
              - Set Ptfm DOF to True, rest to false
              - Set to 0: PtfmMass PtfmRIner PtfmPIner PtfmYIner
              - Maybe: set TowerBsHt=0 PtfmRefzt=0
                  - Add to ElastoDyn output list: PtfmTAxt PtfmTAyt PtfmTAzt PtfmRAxt PtfmRAyt PtfmRAzt
         - Set tower massden to 0 in tower file



    """
    deck = fd.FASTInputDeck(fstFile, readlist=['ED','EDtwr','EDbld','BD','BDbld']) # TODO BeamDyn

    # 
    if workDir is None:
        wd=tempfile.TemporaryDirectory()
        workDir = wd.name
    else:
        wd = None
        os.makedirs(workDir, exist_ok=True) 
    if verbose:
        print('ExtractRNA directory:',workDir)

    if method =='fastlin':

        fst = deck.fst
        ED = deck.ED
        TwrHeight = ED['TowerHt'] - ED['TowerBsHt']
        PtfmRef2TwrTop = ED['TowerHt'] - ED['PtfmRefzt']
        if ED is None:
            raise Exception('ElastoDyn File was not found')

        if nAzim>1:
            #vpsi = np.linspace(0,360/ED['NumBl'],nAzim+1)[:-1]
            vpsi = np.linspace(0,360,nAzim+1)[:-1]
        else:
            vpsi = [ED['Azimuth']]



        fstFiles=[]
        createdFiles=[]
        for psi in vpsi:
            # Filenames
            fstFile2    = os.path.join(workDir, 'ExtractRNA_psi{:03.0f}.fst'.format(psi))
            fstFiles.append(fstFile2)

            # --- Changes to fst file
            fst['TMax']       = 0
            if fst['CompElast']==0:
                fst['CompElast'] = 1
            fst['CompInflow']  = 0
            fst['CompAero']    = 0
            fst['CompServo']   = 0
            fst['CompHydro']   = 0
            fst['CompSub']     = 0
            fst['CompMooring'] = 0
            fst['CompIce']    = 0
            fst['Linearize']  = True
            fst['CalcSteady'] = False
            fst['LinTimes']   = 0
            fst['NLinTimes']  = 1
            fst['LinInputs']  = 2
            fst['LinOutputs'] = 1
            fst['OutFileFmt'] = 2
            fst['OutFmt']     = "ES20.11E2"
            fst['SumPrint']   = False
            fst['Echo']       = False

            # --- Changes to ED file
            ED['FlapDOF1']  = False
            ED['FlapDOF2']  = False
            ED['EdgeDOF']   = False
            ED['TeetDOF']   = False
            ED['DrTrDOF']   = False
            ED['GenDOF']    = False
            ED['YawDOF']    = False
            ED['TwFADOF1']  = False
            ED['TwFADOF2']  = False
            ED['TwSSDOF1']  = False
            ED['TwSSDOF2']  = False
            ED['PtfmSgDOF'] = True
            ED['PtfmSwDOF'] = True
            ED['PtfmHvDOF'] = True
            ED['PtfmRDOF']  = True
            ED['PtfmPDOF']  = True
            ED['PtfmYDOF']  = True
            ED['Azimuth'] = psi
            ED['TwrNodes']   = 3    # No need for resolution here
            ED['NTwGages']   = 0
            ED['NBlGages']   = 0
            ED['SumPrint']   = True
            ED['Echo']       = False
            if not includeYawBearingMass:
                ED['YawBrMass'] = 0
            ED['OutList']  = ['','PtfmTAxt', 'PtfmTAyt', 'PtfmTAzt', 'PtfmRAxt', 'PtfmRAyt', 'PtfmRAzt']
            # Override ED
            if EDDict is not None:
                for k,v in EDDict.items():
                    if k not in ED.keys():
                        raise Exception('Key `{}` not present in ElastoDyn file'.format(k))
                    ED[k] = v

            # --- Changes to ED twr
            twr = deck.fst_vt['ElastoDynTower'] 
            if twr is None:
                raise Exception('ElastoDyn Tower File was not found')
            twrProp = twr['TowProp']
            twrProp[:,1] = 1e-4/TwrHeight
            twrProp[:,2] = 1e13
            twrProp[:,3] = 1e13
            twr['TowProp'] = twrProp
            twr['AdjTwMa'] = 1.0

            deck.write(fstFile2, suffix='_psi{:03.0f}'.format(psi), prefix='ExtractRNA_')
            createdFiles+=deck.inputFiles

        # --- Run FAST
        success, fail = runner.run_fastfiles(fstFiles, fastExe=fastExe, parallel=True, showOutputs=showOutputs, nCores=None, showCommand=verbose, reRun=True, verbose=verbose)
        if not success:
            raise Exception('Some simulation failed. Run with `showOutputs=True` to debug.')

        # --- Postprocess lin
        vmass = []
        vJ_TT = []
        vJ_G  = []
        vCM   = []
        for f in fstFiles:
            linFile = os.path.splitext(f)[0]+'.1.lin'
            mass, J_G, CM = extractPtfmInertiaFromLinFile(linFile)
            TT_2_CM = np.array([CM[0], CM[1], CM[2]-PtfmRef2TwrTop])
            # --- Translate inertia matrix to the tower top
            J_TT = translateInertiaMatrixFromCOG(J_G, mass, r_GP=-TT_2_CM)

            vmass.append(mass)
            vJ_G.append(J_G)
            vJ_TT.append(J_TT)
            vCM.append(CM)
        out_avg={'mass':None, 'CM':None,'J_G':None,'J_TT':None}
        out_azim={'mass':None,'CM':None,'J_G':None,'J_TT':None}
        out_azim['psi']  = np.asarray(vpsi)
        out_azim['mass'] = np.asarray(vmass)
        out_azim['J_G']  = np.asarray(vJ_G)
        out_azim['J_TT'] = np.asarray(vJ_TT)
        out_azim['CM']   = np.asarray(vCM)
        out_avg['mass']  = np.mean(out_azim['mass'], axis=0)
        out_avg['J_G']   = np.mean(out_azim['J_G'] , axis=0)
        out_avg['J_TT']  = np.mean(out_azim['J_TT'], axis=0)
        if nAzim==1:
            out_avg['CM']    = out_azim['CM'][0]
        else:
            out_avg['CM']    = np.mean(out_azim['CM'], axis=0)

        if cleanUp:
            def rm(filename):
                try: 
                    os.remove(filename)
                except:
                    pass

            if wd is None:
                for f in fstFiles:
                    rm(os.path.splitext(f)[0]+'.1.lin')
                    rm(os.path.splitext(f)[0]+'.outb')
                    rm(os.path.splitext(f)[0]+'.ED.sum')
                    rm(f)
                for f in createdFiles:
                    rm(f)
            else:
                wd.cleanup()

        return out_avg, out_azim


    else:
        raise NotImplementedError()



def extractRigidBodyInertia(fstFile, fastExe=None, workDir=None, method='fastlin', nAzim=1, cleanUp=True, 
        zeroInit=True,
        includeYawBearingMass=True, showOutputs=False, verbose=False, EDDict=None):
    """ 
    Compute the RigidBody inertia of an openfast model

    INPUTS:
      - fstFile: fast file for as template. Must point to existing elastodyn files
      - fastExe: path to an OpenFAST executable  compatible with the fstFile provided.
                 The executable is used to extract the RNA information
    OPTIONAL INPUTS:
      - workdir
      - nAzim :  Number of azimuthal steps to use. If >1, the average inertia over a revolution is computed
      - cleanUp: if true, delete file created
      - includeYawBearingMass: if true, the yaw bearing mass will be included in the RNA mass.
      - showOutputs: show Outputs when running the fast executable, useful to debug
      - verbose: display some status on commands being run
      - zeroInit: zero the initial conditions
      - EDDict: Additional dictionary to set values of ElastoDyn, e.g.
              EDDict={'PreCone(1)':-4, 'HubMass':0}

    OUTPUT: dictionary with keys:
       - mass: mass
       - J_Ref: Inertia at PtfmRefz
       - J_G:  Inertia at RNA COG
       - CM:   Distance from PtfmRefz to structure COG
       - azim_*: azimuthal variation of the variables above (useful if nAzim>1)
    """
    deck = fd.FASTInputDeck(fstFile, readlist=['ED','EDtwr','EDbld','SD']) # TODO BeamDyn

    # 
    if workDir is None:
        wd=tempfile.TemporaryDirectory()
        workDir = wd.name
    else:
        wd = None
        os.makedirs(workDir, exist_ok=True) 

    fst = deck.fst
    ED = deck.ED

    if nAzim>1:
        #vpsi = np.linspace(0,360/ED['NumBl'],nAzim+1)[:-1]
        vpsi = np.linspace(0,360,nAzim+1)[:-1]
    else:
        vpsi = [ED['Azimuth']]


    fstFiles=[]
    createdFiles=[]
    for psi in vpsi:
        # Filenames
        fstFile2    = os.path.join(workDir, 'ExtractRB_psi{:03.0f}.fst'.format(psi))
        fstFiles.append(fstFile2)

        # --- Changes to fst file
        fst['TMax']       = 0
        if fst['CompElast']==0:
            fst['CompElast'] = 1
        fst['CompInflow']  = 0
        fst['CompAero']    = 0
        fst['CompServo']   = 0
        fst['CompHydro']   = 0
        #fst['CompSub']     = 0 # <<< We keep SubDyn
        fst['CompMooring'] = 0
        fst['CompIce']    = 0
        fst['Linearize']  = True
        fst['CalcSteady'] = False
        fst['LinTimes']   = 0
        fst['NLinTimes']  = 1
        fst['LinInputs']  = 2
        fst['LinOutputs'] = 1
        fst['OutFileFmt'] = 2
        fst['OutFmt']     = "ES20.11E2"
        fst['SumPrint']   = False
        fst['Echo']       = False

        # --- Changes to ED file
        ED['FlapDOF1']  = False
        ED['FlapDOF2']  = False
        ED['EdgeDOF']   = False
        ED['TeetDOF']   = False
        ED['DrTrDOF']   = False
        ED['GenDOF']    = False
        ED['YawDOF']    = False
        ED['TwFADOF1']  = False
        ED['TwFADOF2']  = False
        ED['TwSSDOF1']  = False
        ED['TwSSDOF2']  = False
        ED['PtfmSgDOF'] = True
        ED['PtfmSwDOF'] = True
        ED['PtfmHvDOF'] = True
        ED['PtfmRDOF']  = True
        ED['PtfmPDOF']  = True
        ED['PtfmYDOF']  = True
        ED['Azimuth'] = psi
        ED['NTwGages']   = 0
        ED['NBlGages']   = 0
        ED['SumPrint']   = True
        ED['Echo']       = False
        if not includeYawBearingMass:
            ED['YawBrMass'] = 0
        ED['OutList']  = ['','PtfmTAxt', 'PtfmTAyt', 'PtfmTAzt', 'PtfmRAxt', 'PtfmRAyt', 'PtfmRAzt']
        # Remove initial conditions
        if zeroInit:
            ED['OoPDefl']     = 0
            ED['IPDefl']      = 0
            ED['BlPitch(1)']  = 0
            ED['BlPitch(2)']  = 0
            ED['BlPitch(3)']  = 0
            ED['TeetDefl']    = 0
            ED['Azimuth']     = 0
            ED['RotSpeed']    = 0
            ED['NacYaw']      = 0
            ED['TTDspFA']     = 0
            ED['TTDspSS']     = 0
            ED['PtfmSurge']   = 0
            ED['PtfmSway']    = 0
            ED['PtfmHeave']   = 0
            ED['PtfmRoll']    = 0
            ED['PtfmPitch']   = 0
            ED['PtfmYaw']     = 0

        # Override ED
        if EDDict is not None:
            for k,v in EDDict.items():
                if k not in ED.keys():
                    raise Exception('Key `{}` not present in ElastoDyn file'.format(k))
                ED[k] = v

        deck.write(fstFile2, suffix='_psi{:03.0f}'.format(psi), prefix='ExtractRB_')
        createdFiles+=deck.inputFiles

    # --- Run FAST
    success, fail = runner.run_fastfiles(fstFiles, fastExe=fastExe, parallel=True, showOutputs=showOutputs, nCores=None, showCommand=verbose, reRun=True, verbose=verbose)
    if not success:
        raise Exception('Some simulation failed. Run with `showOutputs=True` to debug.')
    # --- Postprocess lin
    vmass  = []
    vJ_Ref = []
    vJ_G   = []
    vCM    = []
    for f in fstFiles:
        linFile = os.path.splitext(f)[0]+'.1.lin'
        mass, J_G, CM = extractPtfmInertiaFromLinFile(linFile)
        Ref_2_CM = np.array([CM[0], CM[1], CM[2]])
        # --- Translate inertia matrix to the Ref point top
        J_Ref = translateInertiaMatrixFromCOG(J_G, mass, r_GP=-Ref_2_CM)

        vmass.append(mass)
        vJ_G.append(J_G)
        vJ_Ref.append(J_Ref)
        vCM.append(CM)
    out_avg={'mass':None, 'CM':None,'J_G':None,'J_Ref':None}
    out_azim={'mass':None,'CM':None,'J_G':None,'J_Ref':None}
    out_azim['psi']  = np.asarray(vpsi)
    out_azim['mass'] = np.asarray(vmass)
    out_azim['J_G']  = np.asarray(vJ_G)
    out_azim['J_Ref'] = np.asarray(vJ_Ref)
    out_azim['CM']   = np.asarray(vCM)
    out_avg['mass']  = np.mean(out_azim['mass'], axis=0)
    out_avg['J_G']   = np.mean(out_azim['J_G'] , axis=0)
    out_avg['J_Ref']  = np.mean(out_azim['J_Ref'], axis=0)
    if nAzim==1:
        out_avg['CM']    = out_azim['CM'][0]
    else:
        out_avg['CM']    = np.mean(out_azim['CM'], axis=0)

    if cleanUp:
        def rm(filename):
            try: 
                os.remove(filename)
            except:
                pass

        if wd is None:
            for f in fstFiles:
                rm(os.path.splitext(f)[0]+'.1.lin')
                rm(os.path.splitext(f)[0]+'.outb')
                rm(os.path.splitext(f)[0]+'.ED.sum')
                rm(f)
            for f in createdFiles:
                rm(f)
        else:
            wd.cleanup()

    return out_avg, out_azim



def extractPtfmInertiaFromLinFile(linFile):
    """ 
    Use lin file to obtain Ptfm inertia
    Based on force at platofmr and acceleration at platform

    """
    # np.set_printoptions(precision=5)
    # pd.set_option("precision", 5)

    # --- Read the linearization file
    lin = fl.FASTLinearizationFile(linFile)
    #ED  = weio.read(EDFile)
    #TwrHeight = ED['TowerHt']-ED['TowerBsHt']

    # --- Read lin file and extract relevant part of the D matrix (Minv)
    dfs = lin.toDataFrame()
    # print(dfs.keys())
    # print(lin.keys())
    D   = lin['D']
    dfD = dfs['D']
    # --- Extract the relevant 6x6 matrix
    Cols  = ['PtfmFxN1_[N]', 'PtfmFyN1_[N]', 'PtfmFzN1_[N]', 'PtfmMxN1_[Nm]', 'PtfmMyN1_[Nm]', 'PtfmMzN1_[Nm]']
    Lines = ['PtfmTAxt_[m/s^2]', 'PtfmTAyt_[m/s^2]', 'PtfmTAzt_[m/s^2]', 'PtfmRAxt_[deg/s^2]', 'PtfmRAyt_[deg/s^2]', 'PtfmRAzt_[deg/s^2]']
    Minv = subMat(dfD, rows=Lines, cols=Cols, check=True)
    #print(Minv)

    # --- Convert the units to SI units
    #print('------- Convert to SI units')
    Minv=matToSIunits(Minv, 'D')
    #print(Minv)

    # ---  Inverse the matrix
    #print('------- Inverse the matrix')
    M=inv(Minv)

    # --- Identify mass, inertia and position of center of mass, based on mass matrix 
    m=M[0,0]
    M[np.abs(M)<m*0.001]=0
    np.set_printoptions(precision=4)
    mass, J_G, CM = identifyRigidBodyMM(M)
    return mass, J_G, CM

    #print(np.around(M,1))
    #print('M\n',M)
    #print('PtfmRef2TwrTop\n',PtfmRef2TwrTop)
    #print('J_G')
    #print(J_G)
    #print('CM',CM)
    #print('Mass',mass)
    #print('xCM',CM[0])
    #print('yCM',CM[1])
    #print('zCM',CM[2]-PtfmRef2TwrTop) # NOTE: the center of mass was calculated wrt tower base, we want it wrt tower top



def mainLinInputs(all=False, hub=1, nac=1, ptfm=1, gen=1, pitch=0):
    if all:
        cols=None
    else:
        cols=[]
        if pitch==1:
            cols+=['PitchColl_[rad]']
        if gen==1:
            cols+=['Qgen_[Nm]']
        if hub>=1:
            cols+=['HubFxN1_[N]'] 
            cols+=['HubFyN1_[N]'] 
            cols+=['HubFzN1_[N]'] 
        if hub>=2:
            cols+=['HubMxN1_[Nm]'] 
            cols+=['HubMyN1_[Nm]'] 
            cols+=['HubMzN1_[Nm]'] 
        if nac>=1:
            cols+=['NacFxN1_[N]'] 
            cols+=['NacFyN1_[N]'] 
            cols+=['NacFzN1_[N]'] 
        if nac>=2:
            cols+=['NacMxN1_[Nm]'] 
            cols+=['NacMyN1_[Nm]'] 
            cols+=['NacMzN1_[Nm]'] 
        if ptfm>=1:
            cols+=['PtfmFxN1_[N]'] 
            cols+=['PtfmFyN1_[N]'] 
            cols+=['PtfmFzN1_[N]'] 
        if ptfm>=2:
            cols+=['PtfmMxN1_[Nm]'] 
            cols+=['PtfmMyN1_[Nm]'] 
            cols+=['PtfmMzN1_[Nm]'] 
    return cols



def extractIMUAccFromLinFile(linFile, all=False, hub=1, nac=1, ptfm=1, gen=1, pitch=0):
    """ 
    Extract C and D matrices relevant for IMU
    """

    colIMU=['NcIMUTAxs_[m/s^2]', 'NcIMUTAys_[m/s^2]', 'NcIMUTAzs_[m/s^2]']

    # --- Read the linearization file
    if not isinstance(linFile, fl.FASTLinearizationFile):
        lin = fl.FASTLinearizationFile(linFile)
    else:
        lin=linFile
    dfs = lin.toDataFrame()
    colAugForce = mainLinInputs(all=all, hub=hub, nac=nac, ptfm=ptfm, gen=gen, pitch=pitch)

    # --- Extract the relevant 6x6 matrix
    C = subMat(dfs['C'], rows=colIMU, cols=None       , check=True)
    D = subMat(dfs['D'], rows=colIMU, cols=colAugForce, check=True)
    C=matToSIunits(C, 'C')
    D=matToSIunits(D, 'D')
    return C, D













if __name__ == '__main__':
    scriptDir = os.path.dirname(__file__)
    fstFile = os.path.join(scriptDir, '../../data/NREL5MW/Main_Onshore.fst')
    workDir = '_extractRNA'
    workDir = None
    fastExe =  os.path.join(scriptDir, '../../data/openfast_x64s_dev_2021-12-15.exe')
    RNA, avg = extractRNAInertia(fstFile, workDir=workDir, fastExe=fastExe, nAzim=1) 
    # ---



    print(RNA)
# {'mass': 7815719.392152023, 'CM': array([-9.0578e-03,  0.0000e+00, -8.3596e+01]), 'J_G': array([[6.9453e+09, 0.0000e+00, 1.1100e+07],
#        [0.0000e+00, 6.9303e+09, 0.0000e+00],
#        [1.1100e+07, 0.0000e+00, 1.8958e+08]]), 'J_TT': array([[6.1564e+10, 0.0000e+00, 5.1821e+06],
#        [0.0000e+00, 6.1549e+10, 0.0000e+00],
#        [5.1821e+06, 0.0000e+00, 1.8958e+08]])}
# {'mass': 7815719.392152023, 'CM': array([-9.0578e-03,  0.0000e+00, -8.3596e+01]), 'J_G': array([[6.9453e+09, 0.0000e+00, 1.1100e+07],
#        [0.0000e+00, 6.9303e+09, 0.0000e+00],
#        [1.1100e+07, 0.0000e+00, 1.8958e+08]]), 'J_TT': array([[6.1564e+10, 0.0000e+00, 5.1821e+06],
#        [0.0000e+00, 6.1549e+10, 0.0000e+00],
#        [5.1821e+06, 0.0000e+00, 1.8958e+08]])}

# 
