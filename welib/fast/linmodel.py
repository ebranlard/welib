""" 
Tools to extract a linear model from OpenFAST


"""

##
import numpy as np
import pandas as pd
import copy
import os
import matplotlib.pyplot as plt

# Local
import welib.weio as weio
from welib.weio.fast_output_file import writeDataFrame

from welib.fast.tools.lin import * # backward compatibility
from welib.fast.tools.lin import matToSIunits, renameList, dfToSIunits

from welib.system.statespacelinear import LinearStateSpace
from welib.yams.windturbine import FASTWindTurbine


DEFAULT_COL_MAP_LIN ={
  'psi_rot_[rad]'       : 'psi'      ,
  'Azimuth_[rad]'       : 'psi'      ,
  'RotSpeed_[rad/s]'    : 'dpsi'     ,
  'd_psi_rot_[rad/s]'   : 'dpsi'     ,
  'qt1FA_[m]'           : 'q_FA1'    ,
  'd_qt1FA_[m/s]'       : 'dq_FA1'   ,
  'd_PtfmSurge_[m/s]'   : 'dx'       ,
  'd_PtfmSway_[m/s]'    : 'dy'       ,
  'd_PtfmHeave_[m/s]'   : 'dz'       ,
  'd_PtfmRoll_[rad/s]'  : 'dphi_x'   ,
  'd_PtfmPitch_[rad/s]' : 'dphi_y'   ,
  'd_PtfmYaw_[rad/s]'   : 'dphi_z'   ,
  'PtfmSurge_[m]'       : 'x'        ,
  'PtfmSway_[m]'        : 'y'        ,
  'PtfmHeave_[m]'       : 'z'        ,
  'PtfmRoll_[rad]'      : 'phi_x'    ,
  'PtfmPitch_[rad]'     : 'phi_y'    ,
  'PtfmYaw_[rad]'       : 'phi_z'    ,
  'NcIMUTAxs_[m/s^2]'   : 'NcIMUAx'   ,
  'NcIMUTAys_[m/s^2]'   : 'NcIMUAy'   ,
  'NcIMUTAzs_[m/s^2]'   : 'NcIMUAz'   ,
  'NcIMUTVxs_[m/s]'     : 'NcIMUVx'   ,
  'NcIMUTVys_[m/s]'     : 'NcIMUVy'   ,
  'NcIMUTVzs_[m/s]'     : 'NcIMUVz'   ,
  'BPitch1_[rad]'       : 'pitchB1'    , # Also B1Pitch_[rad]
  'PitchColl_[rad]'     : 'pitch'    , # Also B1Pitch_[rad]
  'Qgen_[Nm]'           : 'Qgen'     ,
  'HubFxN1_[N]'         : 'Thrust'   ,  # A bit controversial to call those "Aero", they come from ED
  'HubFyN1_[N]'         : 'fay'   ,   
  'HubFzN1_[N]'         : 'faz'   ,
  'HubMxN1_[Nm]'        : 'Qaero'    , 
  'HubMyN1_[Nm]'        : 'may'      , 
  'HubMzN1_[Nm]'        : 'maz'      , 
  'PtfmFxN1_[N]'        : 'fhx'      ,
  'PtfmFyN1_[N]'        : 'fhy'      ,
  'PtfmFzN1_[N]'        : 'fhz'      ,
  'PtfmMxN1_[Nm]'       : 'mhx'      ,
  'PtfmMyN1_[Nm]'       : 'mhy'      ,
  'PtfmMzN1_[Nm]'       : 'mhz'      ,
  'Q_Sg_[m]'        : 'x',
  'Q_Sw_[m]'        : 'y',
  'Q_Hv_[m]'        : 'z',
  'Q_R_[rad]'       : 'phi_x',
  'Q_P_[rad]'       : 'phi_y',
  'Q_Y_[rad]'       : 'phi_z',
  'QD_Sg_[m/s]'     : 'dx',
  'QD_Sw_[m/s]'     : 'dy',
  'QD_Hv_[m/s]'     : 'dz',
  'QD_R_[rad/s]'    : 'dphi_x'   ,
  'QD_P_[rad/s]'    : 'dphi_y'   ,
  'QD_Y_[rad/s]'    : 'dphi_z'   ,
  'QD2_Sg_[m/s^2]'  : 'ddx',
  'QD2_Sw_[m/s^2]'  : 'ddy',
  'QD2_Hv_[m/s^2]'  : 'ddz',
  'QD2_R_[rad/s^2]' : 'ddphi_x'   ,
  'QD2_P_[rad/s^2]' : 'ddphi_y'   ,
  'QD2_Y_[rad/s^2]' : 'ddphi_z'   ,
  'NacYaw_[rad]'    : 'yaw',
#           'NacFxN1_[N]'       : 'fnx'   ,
#           'NacMxN1_[N]'       : 'mnx'    , 
#    Qgen_[Nm]  HubFxN1_[N]  HubFyN1_[N]  HubFzN1_[N] 
#               NacFxN1_[N]  NacFyN1_[N]  NacFzN1_[N]
}


DEFAULT_COL_MAP_OF ={
#   'NcIMUTAxs_[m/s^2]'   : 'TTIMUx'   ,
#   'NcIMUTAys_[m/s^2]'   : 'TTIMUy'   ,
#   'NcIMUTAzs_[m/s^2]'   : 'TTIMUz'   ,
  'BldPitch1_[rad]'     : 'pitch'    , # Also B1Pitch_[rad]
#   'Qgen_[Nm]'           : 'Qgen'     ,
  'RtAeroFxh_[N]'         : 'Thrust'   ,
  'RtFldFxh_[N]'          : 'Thrust'   ,
  'RtAeroMxh_[N-m]'       : 'Qaero'   ,
  'RtFldMxh_[N-m]'        : 'Qaero'   ,
#   'HubFyN1_[N]'         : 'fay'   ,
#   'HubFzN1_[N]'         : 'faz'   ,
#   'HubMxN1_[Nm]'        : 'Qaero'    , 
#   'HubMyN1_[Nm]'        : 'may'      , 
#   'HubMzN1_[Nm]'        : 'maz'      , 
  'FxhO_[N]'        : 'fhx'      ,
  'FyhO_[N]'        : 'fhy'      ,
  'FzhO_[N]'        : 'fhz'      ,
  'MxhO_[Nm]'       : 'mhx'      ,
  'MyhO_[Nm]'       : 'mhy'      ,
  'MzhO_[Nm]'       : 'mhz'      ,
  'PtfmSurge_[m]'       : 'x'        ,
  'PtfmSway_[m]'        : 'y'        ,
  'PtfmHeave_[m]'       : 'z'        ,
  'PtfmRoll_[rad]'      : 'phi_x'    ,
  'PtfmPitch_[rad]'     : 'phi_y'    ,
  'PtfmYaw_[rad]'       : 'phi_z'    ,
  'Q_Sg_[m]'        : 'x',
  'Q_Sw_[m]'        : 'y',
  'Q_Hv_[m]'        : 'z',
  'Q_R_[rad]'       : 'phi_x',
  'Q_P_[rad]'       : 'phi_y',
  'Q_Y_[rad]'       : 'phi_z',
  'Q_TFA1_[m]'      : 'q_FA1',
  'Q_TSS1_[m]'      : 'q_FA1',
  'Q_TFA2_[m]'      : 'q_FA2',
  'Q_TSS2_[m]'      : 'q_SS2',
  'Q_Yaw_[m]'       : 'yaw',
  'NacYaw_[rad]'    : 'yaw',
  'Azimuth_[rad]'   : 'psi'      ,
  'Q_DrTr_[rad]'    : 'nu'   ,
  'QD_Sg_[m/s]'     : 'dx',
  'QD_Sw_[m/s]'     : 'dy',
  'QD_Hv_[m/s]'     : 'dz',
  'QD_R_[rad/s]'    : 'dphi_x'   ,
  'QD_P_[rad/s]'    : 'dphi_y'   ,
  'QD_Y_[rad/s]'    : 'dphi_z'   ,
  'QD_TFA1_[m/s]'   : 'dq_FA1',
  'QD_TSS1_[m/s]'   : 'dq_FA1',
  'QD_TFA2_[m/s]'   : 'dq_FA2',
  'QD_TSS2_[m/s]'   : 'dq_SS2',
  'QD_Yaw_[rad/s]'  : 'dyaw',
  'RotSpeed_[rad/s]': 'dpsi',
  'Q_DrTr_[rad/s]'  : 'dnu'   ,
  'QD_GeAz_[rad/s]' : 'dpsi',
  'QD2_Sg_[m/s^2]'  : 'ddx',
  'QD2_Sw_[m/s^2]'  : 'ddy',
  'QD2_Hv_[m/s^2]'  : 'ddz',
  'QD2_R_[rad/s^2]' : 'ddphi_x'   ,
  'QD2_P_[rad/s^2]' : 'ddphi_y'   ,
  'QD2_Y_[rad/s^2]' : 'ddphi_z'   ,
  'QD2_TFA1_[m/s^2]': 'ddq_FA1',
  'QD2_TSS1_[m/s^2]': 'ddq_FA1',
  'QD2_TFA2_[m/s^2]': 'ddq_FA2',
  'QD2_TSS2_[m/s^2]': 'ddq_SS2',
  'QD2_Yaw_[rad/s^2]': 'ddyaw',
  'QD2_GeAz_[rad/s^2]': 'ddpsi',
  'QD2_DrTr_[rad/s^2]': 'ddnu',
  'Q_B1F1_[m]': 'q_B1Fl1',
  'Q_B1F2_[m]': 'q_B1Fl2',
  'Q_B1E1_[m]': 'q_B1Ed1',
}



def _loadOFOut(filename, tMax=None, tRange=None, zRef=None):
    """ 
    see also welib.yams.model.simulator 
    """
    ext = os.path.splitext(filename)[1].lower()
    if ext=='.fst':
        if os.path.exists(filename.replace('.fst','.outb')): 
            outfile=filename.replace('.fst','.outb')
        elif os.path.exists(filename.replace('.fst','.out')): 
            outfile=filename.replace('.fst','.out')
        else:
            raise Exception('Cannot find an OpenFAST output file near: {}'.format(filename))
    else:
        outfile=filename
    print('FASTLinModel: loading OF :', outfile)
    dfFS = weio.read(outfile).toDataFrame()
    if tMax is not None:
        dfFS=dfFS[dfFS['Time_[s]']<tMax]
    if tRange is not None:
        dfFS = dfFS[np.logical_and(dfFS['Time_[s]']>=tRange[0],dfFS['Time_[s]']<=tRange[1])]
    time =dfFS['Time_[s]'].values
    dfFS.reset_index(inplace=True)

    # Remove duplicate
    dfFS = dfFS.loc[:,~dfFS.columns.duplicated()].copy()

    # --- Convert hydro loads to loads at zref
    #if zRef is not None:
    #    from welib.FEM.utils import transferRigidLoads
    #    from welib.yams.utils import transferLoadsZPoint
    #    P_HDRef = np.array((0,0,0))
    #    P_EDRef = np.array((0,0,zRef))
    #    # Input loads are at the body origin (ED ref point)
    #    cols = ['HydroFxi_[N]', 'HydroFyi_[N]', 'HydroFzi_[N]', 'HydroMxi_[N-m]', 'HydroMyi_[N-m]', 'HydroMzi_[N-m]']
    #    if 'Q_R_[rad]' in dfFS.columns:
    #        vphi_x = dfFS['Q_R_[rad]']
    #    else:
    #        vphi_x = dfFS['PtfmRoll_[deg]'].values*np.pi/180
    #    if 'Q_P_[rad]' in dfFS.columns:
    #        vphi_y = dfFS['Q_P_[rad]']
    #    else:
    #        vphi_y = dfFS['PtfmPitch_[deg]'].values*np.pi/180
    #    if 'Q_Y_[rad]' in dfFS.columns:
    #        vphi_z = dfFS['Q_Y_[rad]']
    #    else:
    #        vphi_z = dfFS['PtfmYaw_[deg]'].values*np.pi/180
    #    M = dfFS[cols].values
    #    MT = transferLoadsZPoint(M.T, zRef, vphi_x, vphi_y, vphi_z, rot_type='default').T
    #    cols = ['FxhO_[N]', 'FyhO_[N]', 'FzhO_[N]', 'MxhO_[Nm]', 'MyhO_[Nm]', 'MzhO_[Nm]']
    #    dfHydro = pd.DataFrame(data=MT, columns=cols)
    #    dfFS = pd.concat((dfFS, dfHydro), axis=1)

    return dfFS, time

# --------------------------------------------------------------------------------}
# --- Class to handle a linear model from OpenFAST
# --------------------------------------------------------------------------------{
class FASTLinModel(LinearStateSpace):

    def __init__(self, fstFilename=None, linFiles=None, pickleFile=None, usePickle=False):
        # Init parent class
        LinearStateSpace.__init__(self)
        # --- DATA
        self.WT          = None
        self.fstFilename = fstFilename
        self.pickleFile  = None
        # --- DATA for time simulation
        self.WT_sim = None
        self.sys_sim = None
        self.qop_sim = None
        self.uop_sim = None
        self.yop_sim = None
        self.fstFilename_sim = None
        self.time        = None
        self.dfFS        = None
        self.df          = None

        if usePickle:
            if fstFilename is None:
                raise Exception('Provide an fstFilename to figure out the pickle file')
            self.fstFilename = fstFilename
            if os.path.exists(self.defaultPickleFile):
                self.load()
                return

        if pickleFile is not None:
            # Load pickle File
            self.load(pickleFile)

        elif linFiles is not None:
            # Load all the lin File
            A, B, C, D, xop, uop, yop, sX, sU, sY = self.loadLinFiles(linFiles)
            LinearStateSpace.__init__(self, A=A, B=B, C=C, D=D, 
                    qop=xop, uop=uop, yop=yop,
                    sX=sX, sU=sU, sY=sY,
                    verbose=False)

        elif fstFilename is not None:
            # 
            linFiles = [os.path.splitext(fstFilename)[0]+'.1.lin']
            A, B, C, D, xop, uop, yop, sX, sU, sY = self.loadLinFiles(linFiles)
            LinearStateSpace.__init__(self, A=A, B=B, C=C, D=D, 
                    qop=xop, uop=uop, yop=yop,
                    sX=sX, sU=sU, sY=sY,
                    verbose=False)
        else:
            raise Exception('Input some files')

        if fstFilename is not None:
            print('FASTLinModel: loading WT :',fstFilename)
            self.WT = FASTWindTurbine(fstFilename, algo='OpenFAST')
            self.fstFilename     = fstFilename

        # Set A, B, C, D to SI units
        # will scale the matrices
        self.toSI(verbose=False)


        if usePickle:
            if self.pickleFile is None:
                self.save() # No pickle file was found (.lin was read, we save the pickle to speedup)


    def loadLinFiles(self, linFiles):
        from welib.fast.FASTLin import FASTLin # TODO rename me
        FL = FASTLin(linfiles=linFiles)
        A, B, C, D    = FL.average(WS = None)
        xop, uop, yop = FL.averageOP(WS = None)
        sX, sU, sY = FL.xdescr, FL.udescr, FL.ydescr
        return A, B, C, D, xop, uop, yop, sX, sU, sY

    def rename(self, colMap=None, verbose=False):
        """ Rename labels """
        if colMap is None:
            colMap = DEFAULT_COL_MAP_LIN
        LinearStateSpace.rename(self, colMap=colMap, verbose=verbose)

    # --------------------------------------------------------------------------------}
    # --- SETUP FROM OF 
    # --------------------------------------------------------------------------------{
    def setupSimFromOF(self, outFile=None, fstFilename=None, tRange=None, renameFS=True, colMap=None, inPlace=True, 
            qopMethod='lin', 
            uopMethod='zero', uMethod='zero',
            yopMethod='mean',
            **kwargs):
        """ 
        Set simulation input and model based on OpenFAST simulation.
        INPUTS:
         - tRange : if provided, limit the time vector to tRange
         - rename: if True, rename out file columns based on colMap or on DEFAULT_COL_MAP_OF
         - inPlace: if True, will potentially reduce the dimension of the A, B, C, D of the current lin model
        """
        # TODO: Harmonize with yams.models.simulator

        # --- Load turbine config
        if fstFilename is not None:
            print('FASTLinModel: loading WT :',fstFilename)
            self.WT_sim = FASTWindTurbine(fstFilename, algo='OpenFAST')
            self.fstFilename_sim = fstFilename 
        else:
            self.WT_sim = self.WT
            self.fstFilename_sim = self.fstFilename 

        # --- Load Reference simulation
        #zRef =  -self.p['z_OT'] 
        zRef = - self.WT_sim.twr.pos_global[2]  
        if outFile is None:
            self.dfFS, self.time = _loadOFOut(self.fstFilename_sim, tRange=tRange, zRef=zRef)
        else:
            self.dfFS, self.time = _loadOFOut(outFile, tRange=tRange, zRef=zRef)

        # --- Scale to SI
        self.dfFS_raw = self.dfFS.copy()
        self.dfFS = dfToSIunits(self.dfFS, 'dfOF', verbose=False)

        # --- Rename OpenFAST dataframe (inputs/outputs)
        if renameFS:
            if colMap is None:
                colMap = DEFAULT_COL_MAP_OF
            self.dfFS.columns = renameList(list(self.dfFS.columns), colMap, False)
            # Remove duplicate
            self.dfFS = self.dfFS.loc[:,~self.dfFS.columns.duplicated()].copy()


        # --- Create a linear model for this simulation
        ## DOF names
        sX = self.WT_sim.DOFname
        # sanity check that sX is fully included in self.sX
        sMissingInSelf = [s for s in sX if s not in self.sX]
        if len(sMissingInSelf)>0:
            print('[WARN] The simulation has the following DOFs which are not in the lin model: {}'.format(sMissingInSelf))
            sX = [s for s in sX if s in self.sX]
            print('       The intersection of the two will be used: {}'.format(sX))
        sX = list(sX) + ['d'+s for s in sX]
        sY = self.sY # TODO?
        sU = self.sU # TODO?
        if inPlace:
            self.extract(sX=sX, sU=sU, sY=sY, check=True, inPlace=True) # NOTE: self.s*, A, B, etc. will be updated
            self.sys_sim = None
            self.q0_NL = self.WT_sim.z0[self.sX]
            if qopMethod=='lin':
                # We use the linear operating point from the linear model
                self.qop_sim = self.qop
            elif qopMethod=='zero':
                self.qop_sim = self.qop*0
            elif qopMethod=='mean':
                self.qop_sim = self.qop*0
                for s in self.sX:
                    if s in self.dfFS.columns:
                        self.qop_sim[s] = self.dfFS[s].mean()
                    else:
                        print('[WARN] column missing in OpenFAST DataFrame {}'.format(s))
            else:
                self.qop_sim = self.qop*0
                if len(qopMethod)!=len(self.q_op_sim):
                    raise Exception()
                for i,s in enumerate(self.sX):
                    self.qop_sim[s] = qopMethod[i]

            ## --- qop mean
            self.setStateInitialConditions(self.q0_NL-self.qop_sim)

        else:
            raise NotImplementedError()
            ## Initial conditions (Non-linear)!
            q0_NL = self.WT_sim.z0 #  WT_sim
            q0_lin = q0_NL - self.q0_ # Linear initial conditions are wrt to lin model op
            A, B ,C, D  = self.extract(sX=sX, sU=sU, sY=sY, check=True, inPlace=False)
            #A, B, C, D, xop, uop, yop, sX, sU, sY = self.loadLinFiles(linFiles)
            self.sys_sim = LinearStateSpace(A=A, B=B, C=C, D=D, q0=q0_lin,
                    sX=sX, sU=sU, sY=sY,
                    verbose=False)

        # --- Inputs
        self.setupInputs(uMethod, uopMethod, df=self.dfFS)

        # --- Outputs..
        self.setupOutputs(yopMethod, df=self.dfFS)

        # --- Initial parameters
        #if self.modelName[0]=='B':
        #    self.p = self.WT.yams_parameters(flavor='onebody',**kwargs)
        #else:
        #    self.p = self.WT.yams_parameters(**kwargs)

        return self.time, self.dfFS #, self.p


    def setupInputs(self, uMethod, uopMethod, df=None):
        """ 
            uopMethod='zero', uMethod='zero',
        """

        # --- Non linear inputs
        nu = len(self.sU)
        u  = np.zeros((nu,len(self.time)))
        if uMethod=='zero':
            self._zeroInputs() #uopMethod=uopMethod)

        elif uMethod=='DF':
            for iu, su in enumerate(self.sU):
                if su not in ['fhx','fhy','fhz','mhx','mhy','mhz']: # NOTE: for hydro, we would double count the stiffness
                    print('>>> Getting input from OpenFAST: ',su)
                    vu = df[su].values
                    bNaN = np.isnan(vu)
                    if sum(bNaN)>0:
                        print('[WARN] INPUT {} has {} NaN values, replaced by 0'.format(su, sum(bNaN)))
                    vu[bNaN] = 0
                    u[iu,:] = vu
        else:
            raise NotImplementedError('uMethod ',uMethod)

        # --- Non linear operating point
        if uopMethod=='zero':
            # --- uop 0
            self.uop_sim  = np.zeros(nu)
        elif uopMethod=='mean':
            # --- uop mean
            self.uop_sim  = np.zeros(nu)
            for iu, su in enumerate(self.sU):
                if su not in ['fhx','fhy','fhz','mhx','mhy','mhz']: # NOTE: for hydro, we would double count the stiffness
                    self.uop_sim[iu] = np.mean(u[iu,:])
        elif uopMethod=='lin':
            self.uop_sim  = self.uop.values
        else:
            raise NotImplementedError('uopMethod ',uopMethod)

        # This is the "true" inputs for a linear model
        du = (u.T-self.uop_sim).T
        self.setInputTimeSeries(self.time, du)

            #  u[0,:] = dfOF['HydroFzi_[N]'].values - model.WT.mass*model.WT.gravity
            #  u[0,:] = dfOF['HydroFzi_[N]'].values  
            #  # --- Linear model input operating point
            #  uop['F_hx'] = np.mean(dfFS['HydroFxi_[N]'].values)  *0
            #  uop['F_hy'] = np.mean(dfFS['HydroFyi_[N]'].values)  *0
            #  #uop['F_hz'] = np.mean(dfFS['HydroFzi_[N]'].values)
            #  uop['F_hz'] = p['M_B']*p['g']
            #  uop['M_hx'] = np.mean(dfFS['HydroMxi_[N-m]'].values)*0
            #  if meanMhy:
            #      uop['M_hy'] = np.mean(dfFS['HydroMyi_[N-m]'].values)
            #  else:
            #      uop['M_hy'] = np.mean(dfFS['HydroMyi_[N-m]'].values)*0
            #  #uop['M_hy'] = dfFS['HydroMyi_[N-m]'].values[0]
            #  uop['M_hz'] = np.mean(dfFS['HydroMzi_[N-m]'].values)*0
            #  # --- Linear pertubation inputs
            #  for i,su in enumerate(self.info['su']):
            #      if su=='F_hz': du[i,:] = dfFS['HydroFzi_[N]'].values     - uop[su]  #- p['M_B']*p['g']
            #      if su=='M_hz': du[i,:] = dfFS['HydroMzi_[N-m]'].values   - uop[su]

    def setupOutputs(self, yopMethod, df=None):
        # --- Outputs..
        ny = len(self.sY)
        self.yop_sim  = np.zeros(ny)
        if yopMethod == 'zero':
            pass
        elif yopMethod == 'mean':
            for iy, sy in enumerate(self.sY):
                print('>>> Getting output from DataFrame: ',sy)
                yi = df[sy].values 
                self.yop_sim[iy] = np.mean(yi)
        elif yopMethod == 'lin':
            self.yop_sim = self.yop.values
        else:
            raise Exception('yopMethod')

    def resFromOF(self, removeOP=True):
        """ """
        from scipy.optimize import OptimizeResult as OdeResultsClass 
        #self.WT_sim = self.WT
        #self.fstFilename_sim = self.fstFilename 
        #self.dfFS, self.time = _loadOFOut(outFile, tRange=tRange, zRef=zRef)
        #self.dfFS_raw = self.dfFS.copy()
        #self.dfFS = dfToSIunits(self.dfFS, 'dfOF', verbose=False)
        state_ts = np.zeros((len(self.sX),len(self.time)))
        self.qop_sim
        for i,s in enumerate(self.sX):
            if s not in self.dfFS.keys():
                print('[WARN] fastlinModel: resFromOF: state missing from dataframe: ', s)
            else:
                if removeOP:
                    state_ts[i,:] = self.dfFS[s].values - self.qop_sim[s]
                else:
                    state_ts[i,:] = self.dfFS[s].values
#                 if s!='q_FA1':
#                     state_ts[i,:] = 0
#                 else:
#                     state_ts[i,:] *= 3.32
#                 if s=='x':
#                     state_ts[i,:] = 0
#                 if s=='dx':
#                     state_ts[i,:] = 0
#                 if s=='phi_y':
#                     state_ts[i,:] = 0
#                 if s=='dphi_y':
#                     state_ts[i,:] = 0

        res = OdeResultsClass(t=self.time, y=state_ts) # To mimic result class of solve_ivp
        return res

    def accFromOF(self):
        """ """
        sq = self.WT_sim.DOFname
        cols = ['dd'+ s for s in sq]
        data    = np.full((len(self.time), len(cols)), np.nan)
        df      = pd.DataFrame(columns = cols, data = data)
        for s in cols:
            if s not in self.dfFS.keys():
                print('[WARN] fastlinModel: accFromOF: acceleration missing from dataframe: ', s)
            else:
                df[s] = self.dfFS[s].values
        return df

    # --------------------------------------------------------------------------------}
    # --- INPUTS
    # --------------------------------------------------------------------------------{
    # From welib.system.statespacelinear.py
    #def setInputTimeSeries(self, vTime, vU):
    # From welib.system.statespace.py
    #def setInputFunction(self, fn, signature_u=None):
    def _zeroInputs(self, uopMethod='zero'):
        """ 
        u:   dictionary of functions of time
        uop: dictionary
        du : nu x nt array, time series of time
        """
        # TODO: Harmonize with yams.models.simulator
        # Examples of su: T_a, M_y_a M_z_a F_B

        nu = len(self.sU)
        # --- linear inputs "u" is a "du"
        #u=dict()
        #for su in self.sU:
        #    u[su] = lambda t, q=None: 0  # NOTE: qd not supported yet
        #    #u[su] = lambda t, q=None, qd=None: 0  # Setting inputs as zero as funciton of time

        uop_sim=dict() # Inputs at operating points
        if uopMethod=='zero':
            for su in self.sU:
                uop_sim[su] = 0  # Setting inputs as zero as function of time
        else:
            raise NotImplementedError('uopMethod {}'.format(uopMethod))

        u = np.zeros((nu, len(self.time))) # Zero for all time

        # --- Steady State states
        #qop  = None
        qdop  = None

        # --- Store in class
        # equivalent to: self.setInputTimeSeries(self.time, u)
        self.setInputTimeSeries(self.time, u)

        #self.du  = du
        self.uop_sim = uop_sim
        #self.qop = self.qop_default
        #self.qdop = qdop


    # --------------------------------------------------------------------------------}
    # --- Time integration 
    # --------------------------------------------------------------------------------{
    # From yams.system.statespacelinear :
    #def integrate(self, t_eval, method='RK45', y0=None, u=None, calc='', xoffset=None, uoffset=None, **options):
    # From yams.system.statespace :
    #def integrate(self, t_eval, method='RK45', y0=None, p=None, u=None, calc='u,y,qd', xoffset=None, **options):
    #def res2DataFrame(self, res=None, calc='u,y,xd', sStates=None, xoffset=None, uoffset=None, yoffset=None):
    #def store_states(self, res, sStates=None, xoffset=None, Factors=None):
    #def calc_outputs(self, res=None, insertTime=True, dataFrame=True, yoffset=None):
    #def calc_inputs(self, time=None, res=None, insertTime=True, dataFrame=True, uoffset=None):

    def simulate_fake(self, calc='u,y'):
        """ Use OpenFAST states instead of performing time integration """
        if self.dfFS is None:
            raise Exception('Call `setupSimFromOF` first')
        # Get fake "res" from OF dataframe 
        res  = self.resFromOF()
        # Store states in dataframe, compute inputs & outputs
        #    dfLIX = self.store_states(res=res)
        #    dfLIY = self.calc_outputs(res=res)
        #    dfLII = self.calc_inputs(res=res)
        dfLI  = self.res2DataFrame(res=res, calc=calc, xoffset=self.qop_sim, uoffset=self.uop_sim, yoffset=self.yop_sim)
        # Store OF accelerations TODO Check for duplicate if one day we add acceleration to calc (we should)
        dfLIA = self.accFromOF()
        dfLI = pd.concat((dfLI, dfLIA), axis=1)
        # Store (watch out for that...)
        self.df  = dfLI
        self.res = res
        return dfLI


    def simulate(self, sys_sim=None, out=False, prefix='', calc='', renameStates=False, **options):
        """
        Simulate based on a 

        - calc: string comma separated, e.g. 'u,y,qd' : compute inputs, outputs, and "accelerations"
        """
        # TODO: harmonize with yams.models.simulator
        if sys_sim is None:
            if self.sys_sim is None:
                sys_sim = self
            else:
                sys_sim = self.sys_sim

        # TODO TODO TODO WHY U IS NONE?????
        self.res, self.df =  sys_sim.integrate(self.time, method='RK45', y0=sys_sim.q0, u=None, calc=calc, xoffset=self.qop_sim, uoffset=self.uop_sim, yoffset=self.yop_sim, **options)
        #self.WT_sim.FASTDOFScales, 
        # NOTE: this is a duplication, res2DataFrame is already called in integrate...
        # --- Using OpenFAST DOF names
        if renameStates:
            # TODO TODO TODO THIS IS NASTY
            self.df = self.res2DataFrame(self.res, sStates=self.WT_sim.channels, calc=calc, xoffset=self.qop_sim, uoffset=self.uop_sim, yoffset=self.yop_sim, Factors=self.WT_sim.FASTDOFScales) #x0=qop, xd0=qdop, 
#         if calcOutput:
#             for it, t in enumerate(time):
#                 # Reference point motion
#                 q   = dfOF[qCol].iloc[it].values
#                 qd  = dfOF[qdCol].iloc[it].values
#                 qdd = dfOF[qddCol].iloc[it].values
#                 fh[it,:] = -M.dot(qdd) - C.dot(qd) - K.dot(q)
#                 pass
#         self.resLI, self.sysLI, self.dfLI = resLI, sysLI, dfLI
# 
        if out:
            outFile = self.fstFilename_sim.replace('.fst', '{}_FASTLin.outb'.format(prefix))
            print('FASTLinModel: writing {}'.format(outFile))
            writeDataFrame(self.df, outFile)

        return self.df

    def plotCompare(self, export=False, nPlotCols=2, prefix='', fig=None, figSize=(12,10), title='', useRenamedFS=True,
            columns=None
            ):
        """ 
        NOTE: taken from simulator. TODO harmonization
        """
        from welib.tools.colors import python_colors
        from welib.tools.stats import comparison_stats
        # --- Simple Plot
        dfLI = self.df
        if useRenamedFS:
            dfFS = self.dfFS # <<<
        else:
            dfFS = self.dfFS_raw # <<<

        if dfLI is None and dfFS is None:
            df = dfFS
        else:
            df = dfLI

        if columns is None:
            columns=df.columns[1:]

        if fig is None:
            fig,axes = plt.subplots(int(np.ceil((len(columns))/nPlotCols)), nPlotCols, sharey=False, sharex=True, figsize=figSize)
        else:
            axes=fig.axes
            assert(len(axes)>0)
        if nPlotCols==2:
            fig.subplots_adjust(left=0.07, right=0.98, top=0.955, bottom=0.05, hspace=0.20, wspace=0.20)
        else:
            fig.subplots_adjust(left=0.07, right=0.98, top=0.955, bottom=0.05, hspace=0.20, wspace=0.33)
        for i,ax in enumerate((np.asarray(axes).T).ravel()):
            if i>=len(columns):
                continue
            t2=None; y2=None
            chan=columns[i]
            if dfLI is not None:
                if chan in dfLI.columns:
                    t2=dfLI['Time_[s]']; y2=dfLI[chan]
                    ax.plot(t2, y2, '--' , label='linear', c=python_colors(1))
                else:
                    print('Missing column in Lin: ',chan)
            if dfFS is not None:
                if chan in dfFS.columns:
                    t1=dfFS['Time_[s]']; y1=dfFS[chan]
                    ax.plot(t1, y1, 'k:' , label='OpenFAST')
                    if t2 is not None:
                        stats2, sStats2 =  comparison_stats(t1, y1, t2, y2, stats='eps,R2', method='mean')
                        print(chan, sStats2)
                else:
                    print('Missing column in OpenFAST: ',chan)
            ax.set_xlabel('Time [s]')
            ax.set_ylabel(chan)
            ax.tick_params(direction='in')
            if i==0:
                ax.legend()

        # Scale axes if tiny range
        for i,ax in enumerate((np.asarray(axes).T).ravel()):
            mi, mx = ax.get_ylim()
            mn = (mx+mi)/2
            if np.abs(mx-mn)<1e-6:
                ax.set_ylim(mn-1e-5, mn+1e-5)
        fig.suptitle(title)

        if export:
            fig.savefig(self.fstFilename.replace('.fst','{}_linModel.png'.format(prefix)))

        return fig

    # --------------------------------------------------------------------------------}
    # --- Matrix manipulation
    # --------------------------------------------------------------------------------{
    # From welib.system.statespacelinear
    #def rename(self, colMap, verbose=False):
    #def extract(self, sX=None, sU=None, sY=None, verbose=False, check=True, inPlace=True):
    #def toSI(self, verbose=False):
    #def fromDataFrames(self, A, B=None, C=None, D=None):
    #def toDataFrames(self):


    # --------------------------------------------------------------------------------}
    # ---  IO functions for printing
    # --------------------------------------------------------------------------------{
    # From welib.system.statespacelinear
    #def __repr__(self):

    def printOP(self):
        print('>>> q0_NL  ', dict(self.q0_NL))
        print('>>> qop    ', dict(self.qop))
        print('>>> q0     ', dict(self.q0))
        print('>>> qop_sim', dict(self.qop_sim))
        print('>>> uop_sim', self.uop_sim)
        print('>>> yop_sim', self.yop_sim)


    def picklable(self):
        """ Make the object picklable..."""
        if self.WT:
            self.WT.picklable()
        if self.WT_sim:
            self.WT_sim.picklable()
#         def noneIfLambda(obj):
#             if callable(obj) and obj.__name__ == "<lambda>":
#                 obj=None

    @property
    def defaultPickleFile(self):
        if self.fstFilename is None:
            raise NotImplementedError('Default pickle with no fstFilename')
        return self.fstFilename.replace('.fst','_linModel.pkl')

    def load(self, pickleFile=None):
        if pickleFile is None:
            pickleFile = self.defaultPickleFile
        print('FASTLinModel: loading PKL:',pickleFile)

        d = LinearStateSpace.load(self, pickleFile)
        self.WT          = d['WT']
        self.fstFilename = d['fstFilename']
        self.pickleFile = pickleFile


    def save(self, pickleFile=None):
        if pickleFile is None:
            pickleFile = self.defaultPickleFile
        # Remove MAP dll (problematic in pickle file)
        if self.WT.MAP is not None:
            self.WT.MAP.lib=None
            self.WT.MAP = None
        d = {'fstFilename':self.fstFilename, 'WT':self.WT}
        LinearStateSpace.save(self, pickleFile, d)
        print('FASTLinModel: writing PKL: ', pickleFile)
        self.pickleFile = pickleFile


# --------------------------------------------------------------------------------}
# --- Class to handle a FTNSB model from OpenFAST
# --------------------------------------------------------------------------------{
class FASTLinModelFTNSB(FASTLinModel):

    def __init__(self, fstFilename=None, linFiles=None, pickleFile=None, usePickle=False):
        """ 
        inputFile: a lin file or a fst file
        """
        FASTLinModel.__init__(self, fstFilename, linFiles, pickleFile, usePickle=usePickle)





# --------------------------------------------------------------------------------}
# --- Creating a TNSB model from a FAST model. Used to be called FASTLinModel
# --------------------------------------------------------------------------------{
class FASTLinModelTNSB():
    def __init__(self, ED_or_FST_file, StateFile=None, nShapes_twr=1, nShapes_bld=0, DEBUG=False):

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
        # TODO use printMat from welib.tools.strings
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


def loadLinStateMatModel(StateFile, ScaleUnits=True, Adapt=True, ExtraZeros=False, nameMap={'SvDGenTq_[kNm]':'Qgen_[kNm]'}, ):
    """ 


    Load a pickle file contains A,B,C,D matrices either as sequence or dictionary.
    Specific treatments are possible:
       - ScaleUnits: convert to SI units deg->rad, rpm-> rad/s, kN->N
       - Adapt: 
       - Adapt: 
       - nameMap: rename columns and indices

    If a "model" is given specific treatments can be done

    NOTE:  the A, B, C, D matrices and state file were likely created by 
        FASTLin.average_subset()

    
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
