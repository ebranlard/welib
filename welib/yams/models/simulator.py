""" 

"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# Local
import welib.weio as weio
from welib.weio.fast_output_file import writeDataFrame
from welib.yams.windturbine import FASTWindTurbine
from welib.yams.models.packman import loadPackage


# --------------------------------------------------------------------------------}
# --- Helper functions 
# --------------------------------------------------------------------------------{


def _loadOFOut(filename, tMax=None):
    ext = os.path.splitext(filename)[1].lower()
    if ext=='.fst':
        if os.path.exists(filename.replace('.fst','.outb')): 
            dfFS = weio.read(filename.replace('.fst','.outb')).toDataFrame()
        elif os.path.exists(filename.replace('.fst','.out')): 
            dfFS = weio.read(filename.replace('.fst','.out')).toDataFrame()
        else:
            raise Exception('Cannot find an OpenFAST output file near: {}'.format(filename))
    else:
        dfFS = weio.read(filename).toDataFrame()
    if tMax is not None:
        dfFS=dfFS[dfFS['Time_[s]']<tMax]
    time =dfFS['Time_[s]'].values
    return dfFS, time


def hydroMatToSysMat(M, su, sq=None):
    """ 
    Returns a dataframe with row "su" and columns "sq", filling in the 6x6 hydro matrix where necessary
    """
    # Create an empty dataframe with relevant columns/index
    if sq is None:
        Mout = pd.DataFrame(data=np.zeros(len(su)), index=su)
        # Transform matrix to a dataframe
        M = pd.DataFrame(data=M, index=['F_hx','F_hy','F_hz','M_hx','M_hy','M_hz'])
    else:
        Mout = pd.DataFrame(data=np.zeros((len(su),len(sq))), index=su, columns=sq)
        # Transform matrix to a dataframe
        M = pd.DataFrame(data=M, columns=['x','y','z','phi_x','phi_y','phi_z'], index=['F_hx','F_hy','F_hz','M_hx','M_hy','M_hz'])

    # Columns and rows that are part of the 6x6
    if sq is not None:
        sq_h = [s for s in sq if s in M.columns]
    su_h = [s for s in su if s in M.index]

    # Extract relevant hydro matrix
    if sq is None:
        Mh = M.loc[su_h]
        Mout.loc[su_h] = Mh
    else:
        Mh = M.loc[su_h,sq_h]
        Mout.loc[su_h,sq_h] = Mh
    return Mout
    
def moorMatToSysMat(M, sq=None):
    """ 
    Returns a dataframe with row "sq" and columns "sq", filling in the 6x6 mooring matrix where necessary
    """
    # Create an empty dataframe with relevant columns/index
    sq6 = ['x','y','z','phi_x','phi_y','phi_z']
    if sq is None:
        Mout = pd.DataFrame(data=np.zeros((len(sq),len(sq))), index=sq, columns=sq)
        # Transform matrix to a dataframe
        M = pd.DataFrame(data=M, index=sq6, columns=sq6)
    else:
        Mout = pd.DataFrame(data=np.zeros((len(sq),len(sq))), index=sq, columns=sq)
        # Transform matrix to a dataframe
        M = pd.DataFrame(data=M, columns=sq6, index=sq6)

    # Columns and rows that are part of the 6x6
    if sq is not None:
        sq_h = [s for s in sq if s in M.columns]

    # Extract relevant sub matrix
    Mh = M.loc[sq_h,sq_h]
    Mout.loc[sq_h,sq_h] = Mh
    return Mout


# --------------------------------------------------------------------------------}
# --- SimulatorFromOF 
# --------------------------------------------------------------------------------{
class SimulatorFromOF():
    def __init__(self, WT=None, fstFilename=None, modelName=None, packageDir=''):

        # --- Load the wind turbine model, NOTE relevant parameters "p" are in WT.yams_parameters()
        if WT is not None:
            self.WT = WT
            self.fstFilename = WT.FST.filename
        else:
            self.WT = FASTWindTurbine(fstFilename, twrShapes=[0,2], nSpanTwr=50)  # TODO
            self.fstFilename = fstFilename

        self.dfNL = None
        self.dfLI = None
        # --- Import the python module that was generated
        self.modelName=modelName
        if modelName is not None:
            self.loadPackage(modelName, packageDir)
            self.WT.checkPackage(self.pkg)

    def loadPackage(self, modelName=None, packageDir='', packagePath=None):
        pkg, packagePath=loadPackage(modelName=modelName, packageDir=packageDir, packagePath=packagePath)
        self.pkg          = pkg
        self._packagePath = packagePath
        self.info         = self.pkg.info()

    def unloadPackage(self):
        self.pkg=None
        self.info=None

    def reloadPackage(self):
        self.unloadPackage()
        self.loadPackage(packagePath=self._packagePath)

    def setupSim(self, outFile=None, tMax=None, **kwargs):
        # --- Load Reference simulation
        if outFile is None:
            self.dfFS, self.time = _loadOFOut(self.fstFilename, tMax)
        else:
            self.dfFS, self.time = _loadOFOut(outFile, tMax)

        # --- Initial inputs to zero
        self._zeroInputs()

        # --- Initial parameters
        #p = WT.yams_parameters(flavor='onebody', J_at_Origin=True) # TODO TODO Change for B or F
        if self.modelName[0]=='B':
            self.p = self.WT.yams_parameters(flavor='onebody',**kwargs)
        else:
            self.p = self.WT.yams_parameters(**kwargs)

        return self.time, self.dfFS, self.p


    def _zeroInputs(self):
        """ 
        u:   dictionary of functions of time
        uop: dictionary
        du : nu x nt array, time series of time
        """
        # Examples of su: T_a, M_y_a M_z_a F_B

        nu = len(self.info['su'])
        # --- Non-linear Inputs
        u=dict()
        for su in self.info['su']:
            u[su] = lambda t, q=None, qd=None: 0  # Setting inputs as zero as funciton of time

        # --- Linear inputs
        uop=dict() # Inputs at operating points
        for su in self.info['su']:
            uop[su] = 0  # Setting inputs as zero as function of time

        du = np.zeros((nu, len(self.time))) # Zero for all time

        # --- Steady State states
        qop  = None
        qdop  = None

        # --- Store in class
        self.u   = u
        self.du  = du
        self.uop = uop
        self.qop = qop
        self.qdop = qdop


    def setInputs(self, u, du, uop, qop=None, qdop=None):
        """ 
        u:   dictionary of functions of time
        uop: dictionary
        du : nu x nt array, time series of time
        """
        if not all(k in u.keys() for k in self.info['su']):
            print('u keys()     :',u.keys())
            print('self.u keys():',self.info['su'])
            raise Exception('Some u keys are missing')
        if not all(k in uop.keys() for k in self.info['su']):
            print('uop keys()     :',uop.keys())
            print('self.uop keys():',self.info['su'])
            raise Exception('Some uop keys are missing')

        self.u   = u
        self.du  = du
        self.uop = uop
        self.qop = qop
        self.qdop = qdop



    def linmodel(self, MCKextra=None, MCKu=None, noBlin=False):
        # --- Time Integration
        sysLI = self.WT.py_lin(self.pkg, self.p, self.time, uop=self.uop, du=self.du, qop=self.qop, qdop=self.qdop, MCKextra=MCKextra, MCKu=MCKu, noBlin=noBlin)
        return sysLI


    def simulate(self, out=False, prefix='', NL=True, Lin=True, MCKextra=None, MCKu=None, calcOutput=True, noBlin=False):
        dfNL = None
        dfLI = None
        if calcOutput:
            acc=True
            forcing=False
        # --- Time Integration
        if Lin:
            resLI, sysLI, dfLI = self.WT.simulate_py_lin(self.pkg, self.p, self.time, uop=self.uop, du=self.du, qop=self.qop, qdop=self.qdop, MCKextra=MCKextra, MCKu=MCKu, acc=acc, forcing=forcing, noBlin=noBlin)
        if NL:
            resNL, sysNL, dfNL = self.WT.simulate_py    (self.pkg, self.p, self.time, u=self.u, acc=acc, forcing=forcing)

        if calcOutput:
            if Lin:
                pass
                #uop=self.uop, du=self.du, qop=self.qop, qdop=self.qdop, MCKextra=MCKextra, MCKu=MCKu
#             for it, t in enumerate(time):
#                 # Reference point motion
#                 q   = dfOF[qCol].iloc[it].values
#                 qd  = dfOF[qdCol].iloc[it].values
#                 qdd = dfOF[qddCol].iloc[it].values
# 
#                 fh[it,:] = -M.dot(qdd) - C.dot(qd) - K.dot(q)
#                 pass
        # Store in object
        self.MCKu=MCKu
        if Lin:
            self.resLI, self.sysLI, self.dfLI = resLI, sysLI, dfLI
        if NL:
            self.resNL, self.sysNL, self.dfNL = resNL, sysNL, dfNL

        if out:
            if Lin:
                writeDataFrame(dfLI, self.fstFilename.replace('.fst', '{}_yamsSim_Lin.outb'.format(prefix)))
            if NL:
                writeDataFrame(dfNL, self.fstFilename.replace('.fst', '{}_yamsSim_NL.outb'.format(prefix)))

        return dfNL, dfLI

    # --------------------------------------------------------------------------------}
    # --- property
    # --------------------------------------------------------------------------------{
    # Linear properties
    @property
    def M_lin(self): return self.sysLI.M
    @property
    def K_lin(self): return self.sysLI.K
    @property
    def C_lin(self): return self.sysLI.C
    @property
    def B_lin(self): return self.sysLI._B
    @property
    def q0_lin(self): return self.sysLI.q0
    @property
    def forcing0_lin(self): return self.sysLI._forcing0
    # Non-linear properties
    @property
    def M0(self): return self.sysNL._M0
    @property
    def forcing0(self): return self.sysNL._forcing0
    @property
    def q0(self): return self.sysNL.q0


    def plot(self, export=False, nPlotCols=2, prefix='', fig=None, figSize=(12,10), title=''):
        from welib.tools.colors import python_colors
        # --- Simple Plot
        dfNL = self.dfNL
        dfLI = self.dfLI
        dfFS = self.dfFS
        if dfLI is None and dfNL is None:
            df = dfFS
        elif dfLI is None:
            df = dfNL
        else:
            df = dfLI

        if fig is None:
            fig,axes = plt.subplots(int(np.ceil((len(df.columns)-1)/nPlotCols)), nPlotCols, sharey=False, sharex=True, figsize=figSize)
        else:
            axes=fig.axes
            assert(len(axes)>0)
        if nPlotCols==2:
            fig.subplots_adjust(left=0.07, right=0.98, top=0.955, bottom=0.05, hspace=0.20, wspace=0.20)
        else:
            fig.subplots_adjust(left=0.07, right=0.98, top=0.955, bottom=0.05, hspace=0.20, wspace=0.33)
        for i,ax in enumerate((np.asarray(axes).T).ravel()):
            if i+1>=len(df.columns):
                continue
            chan=df.columns[i+1]
            if dfNL is not None:
                if chan in dfNL.columns:
                    ax.plot(dfNL['Time_[s]'], dfNL[chan], '-'  , label='non-linear', c=python_colors(0))
                else:
                    print('Missing column in NL: ',chan)
            if dfLI is not None:
                if chan in dfLI.columns:
                    ax.plot(dfLI['Time_[s]'], dfLI[chan], '--' , label='linear', c=python_colors(1))
                else:
                    print('Missing column in Lin: ',chan)
            if dfFS is not None:
                if chan in dfFS.columns:
                    ax.plot(dfFS['Time_[s]'], dfFS[chan], 'k:' , label='OpenFAST')
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
            fig.savefig(self.fstFilename.replace('.fst','{}_yamsSim.png'.format(prefix)))

        return fig

    def save(self, filename=None, prefix=''):
        if filename is None:
            filename = self.fstFilename.replace('.fst','{}_yamsSim.pkl'.format(prefix))

        import pickle
        D={'dfNL':self.dfNL, 'dfLI':self.dfLI, 'dfFS':self.dfFS, 'p':self.p}
        print('>>> Export to:',filename)
        pickle.dump(D,  open(filename,'wb'))

    # --------------------------------------------------------------------------------}
    # --- Model specific
    # --------------------------------------------------------------------------------{
    def setPrescribedHydroInputs(self, zRef=None, meanMhy=False):
        """ Set inputs based on OpenFAST"""
        dfFS    = self.dfFS
        p       = self.p
        time    = self.time
        u       = self.u
        uop     = self.uop
        du      = self.du
        if zRef is None:
            if self.modelName[0]=='B':
                zRef =  self.p['z_B0']
            else:
                zRef =  -self.p['z_OT'] 
        P_HDRef = np.array((0,0,0))
        P_EDRef = np.array((0,0,zRef))

        if 'hydro0' in self.modelName:
            print('>>> Precribed Hydro Loads at ',(0,0,0))
        elif 'hydroO' in self.modelName:
            print('>>> Precribed Hydro Loads at ',(0,0,zRef), 'NOTE: WEIRD SIGN<<<<<< TODO TODO TODO TODO')
            from welib.FEM.utils import transferRigidLoads
            from welib.yams.utils import transferLoadsZPoint
            # Input loads are at the body origin (ED ref point)
            cols = ['HydroFxi_[N]', 'HydroFyi_[N]', 'HydroFzi_[N]', 'HydroMxi_[N-m]', 'HydroMyi_[N-m]', 'HydroMzi_[N-m]']
            if 'Q_R_[rad]' in dfFS.columns:
                vphi_x = dfFS['Q_R_[rad]']
            else:
                vphi_x = dfFS['PtfmRoll_[deg]'].values*np.pi/180
            if 'Q_P_[rad]' in dfFS.columns:
                vphi_y = dfFS['Q_P_[rad]']
            else:
                vphi_y = dfFS['PtfmPitch_[deg]'].values*np.pi/180
            M = dfFS[cols].values
            #MT = transferRigidLoads(M.T, P_HDRef, P_EDRef).T
            MT = transferLoadsZPoint(M.T, zRef, vphi_x, vphi_y).T
            dfFS = pd.DataFrame(data=MT, columns=cols)
        else:
            raise NotImplementedError()

        # Input loads are at the "0" hydro point
        u['F_hx'] = lambda t,q=None,qd=None: np.interp(t, time, dfFS['HydroFxi_[N]'])
        u['F_hy'] = lambda t,q=None,qd=None: np.interp(t, time, dfFS['HydroFyi_[N]'])
        u['F_hz'] = lambda t,q=None,qd=None: np.interp(t, time, dfFS['HydroFzi_[N]'])
        u['M_hx'] = lambda t,q=None,qd=None: np.interp(t, time, dfFS['HydroMxi_[N-m]'])
        u['M_hy'] = lambda t,q=None,qd=None: np.interp(t, time, dfFS['HydroMyi_[N-m]'])
        u['M_hz'] = lambda t,q=None,qd=None: np.interp(t, time, dfFS['HydroMzi_[N-m]'])
        # --- Linear model input operating point
        uop['F_hx'] = np.mean(dfFS['HydroFxi_[N]'].values)  *0
        uop['F_hy'] = np.mean(dfFS['HydroFyi_[N]'].values)  *0
        #uop['F_hz'] = np.mean(dfFS['HydroFzi_[N]'].values)
        uop['F_hz'] = p['M_B']*p['g']
        uop['M_hx'] = np.mean(dfFS['HydroMxi_[N-m]'].values)*0
        if meanMhy:
            uop['M_hy'] = np.mean(dfFS['HydroMyi_[N-m]'].values)
        else:
            uop['M_hy'] = np.mean(dfFS['HydroMyi_[N-m]'].values)*0
        #uop['M_hy'] = dfFS['HydroMyi_[N-m]'].values[0]
        uop['M_hz'] = np.mean(dfFS['HydroMzi_[N-m]'].values)*0
        # --- Linear pertubation inputs
        for i,su in enumerate(self.info['su']):
            if su=='F_hx': du[i,:] = dfFS['HydroFxi_[N]'].values     - uop[su]
            if su=='F_hy': du[i,:] = dfFS['HydroFyi_[N]'].values     - uop[su]
            if su=='F_hz': du[i,:] = dfFS['HydroFzi_[N]'].values     - uop[su]  #- p['M_B']*p['g']
            if su=='M_hx': du[i,:] = dfFS['HydroMxi_[N-m]'].values   - uop[su]
            if su=='M_hy': du[i,:] = dfFS['HydroMyi_[N-m]'].values   - uop[su]
            if su=='M_hz': du[i,:] = dfFS['HydroMzi_[N-m]'].values   - uop[su]

if __name__ == '__main__':
    pass
