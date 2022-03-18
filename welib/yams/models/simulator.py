""" 

"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import importlib
# Local
from welib.yams.windturbine import FASTWindTurbine
import welib.weio as weio


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
            raise NotImplementedError()
    else:
        dfFS = weio.read(filename).toDataFrame()
    if tMax is not None:
        dfFS=dfFS[dfFS['Time_[s]']<tMax]
    time =dfFS['Time_[s]'].values
    return dfFS, time
    


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

        # --- Import the python module that was generated
        self.modelName=modelName
        if modelName is not None:
            self.loadPackage(modelName, packageDir)

    def loadPackage(self, modelName=None, packageDir='', packagePath=None):
        if modelName is not None:
            packagePath = os.path.join(packageDir, modelName)
        packageString = packagePath.replace('/','.').replace('\\','.')
        self.pkg = importlib.import_module(packageString)
        self.info = self.pkg.info()

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
            u[su] = lambda t, q, qd: 0  # Setting inputs as zero as funciton of time

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

    def simulate(self):
        # --- Time Integration
        resLI, sysLI, dfLI = self.WT.simulate_py_lin(self.pkg, self.p, self.time, uop=self.uop, du=self.du, qop=self.qop)
        resNL, sysNL, dfNL = self.WT.simulate_py    (self.pkg, self.p, self.time, u=self.u)
        # Store in object
        self.resLI, self.sysLI, self.dfLI = resLI, sysLI, dfLI
        self.resNL, self.sysNL, self.dfNL = resNL, sysNL, dfNL

        return dfNL, dfLI


    def plot(self, export=True, nPlotCols=2):
        # --- Simple Plot
        dfNL=self.dfNL
        dfLI=self.dfLI
        dfFS=self.dfFS
        fig,axes = plt.subplots(int(np.ceil((len(dfNL.columns)-1)/nPlotCols)), nPlotCols, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
        fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
        for i,axi in enumerate((np.asarray(axes).T).ravel()):
            if i+1<len(dfNL.columns):
                ax=axi
                chan=dfNL.columns[i+1]
                print(chan)
                ax.plot(dfNL['Time_[s]'], dfNL[chan], '-'  , label='non-linear')
                ax.plot(dfLI['Time_[s]'], dfLI[chan], '--' , label='linear')
                if dfFS is not None:
                    ax.plot(dfFS['Time_[s]'], dfFS[chan], 'k:' , label='OpenFAST')
                ax.set_xlabel('Time [s]')
                ax.set_ylabel(chan)
                ax.tick_params(direction='in')
        ax.legend()

        if export:
            fig.savefig(self.fstFilename.replace('.fst','_yamsSim.png'))

        return fig


    # --------------------------------------------------------------------------------}
    # --- Model specific
    # --------------------------------------------------------------------------------{
    def setPrescribedHydroInputs(self, modelName):
        pass
        # TODO TODO TODO
        # B MODEL
        #if noHydro:
        #    pass
        #    #u['F_hz'] = lambda t,q=None,qd=None: 0
        #    #u['M_hy'] = lambda t,q=None,qd=None: 0
        #else:
        #    u['F_hx'] = lambda t,q=None,qd=None: np.interp(t, time, dfFS['HydroFxi_[N]'])
        #    u['F_hz'] = lambda t,q=None,qd=None: np.interp(t, time, dfFS['HydroFzi_[N]'])
        #    u['M_hy'] = lambda t,q=None,qd=None: np.interp(t, time, dfFS['HydroMyi_[N-m]'])
        # B model
        #uop=dict()
        #uop['F_hx'] = np.mean(dfFS['HydroFxi_[N]'].values)*0
        ##uop['F_hz'] = np.mean(dfFS['HydroFzi_[N]'].values)   
        #uop['F_hz'] = p['M_B']*p['g']
        #uop['M_hy'] = np.mean(dfFS['HydroMyi_[N-m]'].values)*0
        #du = np.zeros((nu, len(time)))
        #du[0,:] = dfFS['HydroFxi_[N]'].values     - uop['F_hx']
        #du[1,:] = dfFS['HydroFzi_[N]'].values     - uop['F_hz']  #- p['M_B']*p['g']
        #du[2,:] = dfFS['HydroMyi_[N-m]'].values   - uop['M_hy']
        #du[0,:] *= 0
        #du[2,:] *= 0

if __name__ == '__main__':
    pass
