import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from welib.essentials import *
from welib.yams.windturbine import FASTWindTurbine
from welib.tools.stats import comparison_stats
from welib.weio.dataframe import WEIODataFrame
from welib.tools.strings import latexStrip




class YAMSSectionLoadCalculator():
    def __init__(self, fstFile=None, WT=None, HD_compFile=None):
        self.fstFile = fstFile
        self.FAST = None
        self.dfRef = None

        if WT is None:
            self.WT = FASTWindTurbine(fstFile, algo='OpenFAST', HD_compFile=HD_compFile).WT
        else:
            self.WT = WT

    def emptyInputDF(self, nt, inputFrame='R_xs'):
        # TODO use OpenFAST Units
        DOFNames = ['Sg', 'Sw', 'Hv' ,'R', 'P', 'Y', 'TFA1', 'TSS1', 'Yaw']
        sq   = ['Q_'+s for s in DOFNames]
        sqd  = ['QD_'+s for s in DOFNames]
        sqdd = ['QD2_'+s for s in DOFNames]
        if inputFrame == 'R_xs':
            # Input at point R in coordinate system xs
            sLoads = ['Fadd_R_xs', 'Madd_R_xs']
        else:
            raise NotImplementedError()
#             sqdd = ['FN', 'Madd_R_xs']
        cols = ['Time']+sq+sqd+sqdd +sLoads
        data = np.zeros( (nt, len(cols)))
        #df   = pd.DataFrame(columns=cols, data=data)
        df   = pd.DataFrame()
        return df

    def fromDF(self, df, useTopLoadsFromDF=False, useInterfaceLoadsFromDF=False, noAcc=False, dt_resample=None, tRange=None):
        """Compute OpenFAST-like section-load outputs in a single call.

        Typical usage:
            YSL = YAMSSectionLoadCalculator(fstFile=fstFile)
            dfOut = YSL.calc(dfIn)

        For monopile/foundation section outputs:
            dfOut, sections = YSL.fromSD(dfIn)
        """
        from scipy.interpolate import interp1d


        if tRange is not None:
            df=df[df['Time_[s]']>=tRange[0]]
            df=df[df['Time_[s]']<=tRange[1]]
            df.reset_index(inplace=True)

        if dt_resample is not None:
            INFO(f'Resampling at {dt_resample}')
            t_new = np.arange(df['Time_[s]'].min(), df['Time_[s]'].max(), dt_resample )
            f = interp1d(df['Time_[s]'], df.drop(columns=['Time_[s]']), axis=0, kind='linear', fill_value='extrapolate')
            df = pd.DataFrame(f(t_new), columns=df.columns.drop('Time_[s]'))
            df.insert(0, 'Time_[s]', t_new)

        df = WEIODataFrame(df) # Case incensitive
        self.dfRef = df


        dfOut, sections = self.WT.calcOutputsFromDF(df, useTopLoadsFromDF=useTopLoadsFromDF, useInterfaceLoadsFromDF=useInterfaceLoadsFromDF, noAcc=noAcc)

        self.dfOut = dfOut
        self.sec   = sections

#         dfOut.export(outBase + '.YAMS.outb')

        zDepth   = None
        if self.WT.pHD is not None:
            zDepth   = self.WT.pHD['zDepth']

        return dfOut, sections


    # --------------------------------------------------------------------------------}
    # --- Plots 
    # --------------------------------------------------------------------------------{
    def plot_tower_section_loads(self, IsecTwr=None, component='MLyt', figFilename=None, printStats=False):
        dfRef = self.dfRef
        dfOut = self.dfOut

        if IsecTwr is None:
            IsecTwr = [0,4,8] # Section indices

        fig,ax = plt.subplots(1, 1, sharey=False, figsize=(12.4,12.8)) # (6.4,4.8)
        fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
        for iiED,iED in enumerate(IsecTwr):
            sT = f'TwHt{iED+1}{component}_[kN-m]'
            t1 = dfRef['Time_[s]'].values;
            y1 = dfRef[sT].values
            t2 = dfOut['Time_[s]'].values
            y2 = dfOut[sT].values
            ax.plot(t1, y1, 'k-'                        , label='OpenFAST' if iiED==0 else None)
            ax.plot(t2, y2,  '--' , color=fColrs(iiED)  , label='Sim Ht{}'.format(iED+1))
            stats, sStats =  comparison_stats(t1,y1,t2,y2, stats='sigRatio,eps,R2', method='meanabs'); 
            addStats(ax, sT, sStats, printStats=printStats, factY=0.8)
        ax.set_xlabel('Time [s]')
        ax.set_ylabel('')
        ax.legend()
        if figFilename is not None:
            fig.savefig(figFilename)
        return fig

    def plot_tower_accelerations(self, IsecTwr=None, component='x', figFilename=None, printStats=False):
        dfRef = self.dfRef
        dfOut = self.dfOut

        if IsecTwr is None:
            IsecTwr = [0,4,8] # Section indices

        n = len(IsecTwr)
        if f'NcIMUTA{component}s_[m/s^2]' in dfRef:
            n+=1

        j=-1
        fig,axes = plt.subplots(n, 1, sharey=False, figsize=(6.4,5.8)) # (6.4,4.8)
        fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)

        if f'NcIMUTA{component}s_[m/s^2]' in dfRef:
            j=j+1;sig = 'NcIMUTAxs_[m/s^2]'; t1, y1, t2, y2 = dfRef['Time_[s]'].values ,dfRef[sig].values ,dfOut['Time_[s]'].values ,dfOut[sig].values; 
            axes[j].plot(t1, y1, 'k-'); axes[j].plot(t1, y2, '--'); axes[j].set_ylabel(sig); 
            stats, sStats =  comparison_stats(t1,y1,t2,y2, stats='sigRatio,eps,R2', method='meanabs'); 
            addStats(ax, sig, sStats, printStats=printStats, factY=0.8)
            #print(sig, stats) # TODO
        for iiED,iED in enumerate(IsecTwr):
            j=j+1; sig = f'TwHt{iED+1}AL{component}t_[m/s^2]'; t1, y1, t2, y2 = dfRef['Time_[s]'].values ,dfRef[sig].values ,dfOut['Time_[s]'].values ,dfOut[sig].values; 
            axes[j].plot(t1, y1, 'k-'); axes[j].plot(t1, y2, '--'); axes[j].set_ylabel(sig); 
            stats, sStats =  comparison_stats(t1,y1,t2,y2, stats='sigRatio,eps,R2', method='meanabs'); 
            addStats(ax, sig, sStats, printStats=printStats, factY=0.8)
            #print(sig, stats)
        if figFilename is not None:
            fig.savefig(figFilename)
        return fig



    def plot_monopile_section_loads_stats(self, IZ=None, component=0, figFilename=None, tRange=None, stat='meanabs'):
        F_sec    = self.sec['monopile']['F_sec']
        F_secRef = self.sec['monopile']['F_secRef']
        zDepth   = self.sec['monopile']['z']
        vTime    = self.dfRef['Time_[s]']
        IZ = [int(3*len(zDepth)/6)-5, int(2*len(zDepth)/6), int(1*len(zDepth)/6), 0]
        fig = sec_plotFM(vTime, zDepth, F_sec, F_secRef, stat=stat, tRange=tRange)
        if figFilename is not None:
            fig.savefig(figFilename)

    
    def plot_monopile_section_loads(self, IZ=None, component=0, figFilename=None, tRange=None, printStats=False):
        """ 
         - component: 0 Fx, 4:My
        """

        F_sec    = self.sec['monopile']['F_sec']
        F_secRef = self.sec['monopile']['F_secRef']
        zDepth   = self.sec['monopile']['z']
        vTime    = self.dfRef['Time_[s]']


        fig, axes = plt.subplots(len(IZ), 1, sharey=True, sharex=True, figsize=(8.4,7.5))
        fig.subplots_adjust(left=0.12, right=0.96, top=0.95, bottom=0.05, hspace=0.07, wspace=0.20)
        for ii, iz in enumerate(IZ):
            time_plot (vTime, F_secRef[0, iz, :]/1e6,  F_sec[0, iz, :]/1e6, f'z={zDepth[iz]:.0f}m', tRange=tRange, ax=axes[ii], fig=fig)
            stats, sStats =  comparison_stats(vTime, F_secRef[0,iz,:]/1e6, vTime, F_sec[0,iz,:]/1e6, stats='sigRatio,eps,R2', method='meanabs')
            addStats(ax, 'Fsec'+str(component), sStats, printStats=printStats, factY=0.8)
            #print(f'z {zDepth[iz]:5.0f}: ', stats)
        axes[-1].set_xlabel('Time [s]')
        fig.suptitle('FxSec [MN]')
        if figFilename is not None:
            fig.savefig(figFilename)

        return fig


    def plot_comp(self, sig, ylabel=None, figFilename=None, tRange=None, scale=1, ax=None, factY=0.8, printStats=True):
        if ylabel is None:
            ylabel=sig
        vTime = self.dfRef['Time_[s]']
        fig, ax = time_plot (vTime, self.dfRef[sig]/scale, self.dfOut[sig]/scale, ylabel, tRange=tRange, ax=ax)
        stats, sStats =  comparison_stats(vTime, self.dfRef[sig]/scale, vTime, self.dfOut[sig]/scale, stats='sigRatio,eps,R2', method='1-2')
        addStats(ax, sig, sStats, printStats=printStats, factY=0.8)

        if figFilename is not None:
            fig.savefig(figFilename)
        return fig, stats

    def plot_section_loads_stats(self, IZ=None, component=0, figFilename=None, tRange=None, stat='meanabs'):
        F_sec    = self.sec['combined']['F_sec']
        F_secRef = self.sec['combined']['F_secRef']
        zDepth   = self.sec['combined']['z']
        zRef     = self.sec['combined']['zRef']
        vTime    = self.dfRef['Time_[s]']
        fig = sec_plotFM(vTime, zDepth, F_sec, F_secRef, zRef=zRef, stat=stat, tRange=tRange)
        if figFilename is not None:
            fig.savefig(figFilename)
        return fig



# --------------------------------------------------------------------------------}
# --- Helper functions TODO
# --------------------------------------------------------------------------------{
from welib.tools.colors import MW_Orange, fColrs
import matplotlib
COLRS   = [fColrs(1), MW_Orange ,fColrs(2), fColrs(3)]
try:
    cmap = matplotlib.colormaps['viridis']
except:
    cmap = matplotlib.cm.get_cmap('viridis')
CMAP = [(cmap(v)[0],cmap(v)[1],cmap(v)[2],0.9) for v in np.linspace(0,1,5+1)]
COLRS=[CMAP[1], CMAP[4]]
colRef=COLRS[0]
colSim=COLRS[1]
LWRef=2.4
LWSim=1.5



def addStats(ax, sig, sStats, printStats=False, factY=0.8):
    Ylim = ax.get_ylim(); Xlim = ax.get_xlim()
    ax.text(Xlim[0]+(Xlim[1]-Xlim[0])/1000 ,Ylim[0]+(Ylim[1]-Ylim[0])*factY, sStats, fontsize=10)
    if printStats:
        print(f"{sig:10s} "+latexStrip(sStats))


def meanabs(x, **kwargs):
    return np.mean(np.abs(x), **kwargs)


def time_plot(t, ref=None, sim=None, label='', other=None, ax=None, tRange=None, refLab='OpenFAST', otherLab='Other'):
    figNotProvided = ax is None
    if ax is None:
        fig=plt.figure()
        fig.subplots_adjust(left=0.18, right=0.94, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
        ax = fig.add_subplot(111)
    else:
        fig = ax.figure
    if ref is not None:
        ax.plot(t, ref,'-'  , c=colRef, lw=LWRef, label=refLab)
    if sim is not None:
        ax.plot(t, sim,'--' , c=colSim, lw=LWSim, label='YAMS')
    if other is not None:
        ax.plot(t, other,':' , label=otherLab)
    if tRange is not None:
        ax.axvline(x = tRange[0] , color = 'k', ls='--')
        ax.axvline(x = tRange[1], color = 'k', ls='--')
    if figNotProvided:
        ax.legend()
        ax.set_xlabel('Time [s]')
    ax.set_ylabel(label)
    ax.tick_params(direction='in')
    return fig, ax


def sec_plotFM(vTime, zDepth, F_sec, F_secRef=None, zRef=None, stat='mean', tRange=None, label='Section Force', other=None, otherLab='Other'):
    if tRange is None:
        IT = np.arange(0,F_sec.shape[1])
    else:
        IT = np.logical_and(vTime>tRange[0], vTime<tRange[1])
        if len(IT)==0:
            IT = np.arange(0, sec.shape[1])
    if zRef is None:
        zRef = zDepth

    fstat = {'mean':np.mean, 'max':np.max, 'std':np.std, 'meanabs':meanabs}[stat]
    if F_sec is not None:
        Fx_sec0 = fstat(F_sec   [0,:,IT],axis=0) # 3, nSpan nT
        My_sec0 = fstat(F_sec   [4,:,IT],axis=0)

    if F_secRef is not None:
        Fx_sec1 = fstat(F_secRef[0,:,IT],axis=0) # 3, nSpan nT
        My_sec1 = fstat(F_secRef[4,:,IT],axis=0)
        bNaN = np.isnan(Fx_sec1) # if sections are missing
        Fx_sec1 = Fx_sec1[~bNaN] 
        My_sec1 = My_sec1[~bNaN] 
        zRef    = zRef[~bNaN] 
        
    if other is not None:
        Fx_sec2 = fstat(other[0,:,IT],axis=0) # 3, nSpan nT
        My_sec2 = fstat(other[4,:,IT],axis=0)

    fig,axes = plt.subplots(1, 2, sharey=True, figsize=(6.4,4.8)) # (6.4,4.8)
    fig.subplots_adjust(left=0.14, right=0.95, top=0.95, bottom=0.12, hspace=0.20, wspace=0.20)
    ax=axes[0]
    if F_secRef is not None:
        ax.plot(Fx_sec1/1e6, zRef,    '-' , c=colRef, lw=LWRef, label='OpenFAST')
    if F_sec is not None:
        ax.plot(Fx_sec0/1e6, zDepth,  '--', c=colSim, lw=LWSim, label='YAMS')
    if other is not None:
        ax.plot(Fx_sec2/1e6, zDepth,  ':' ,label=otherLab)

    ax.tick_params(direction='in')

    ax.set_xlabel(label)
    ax.set_ylabel('Vertical position [m]')

    ax=axes[1]
    if F_secRef is not None:
        ax.plot(My_sec1/1e6, zRef , '-',  c=colRef, lw=LWRef, label='OpenFAST')
    if F_sec is not None:
        ax.plot(My_sec0/1e6, zDepth, '--', c=colSim, lw=LWSim, label='YAMS')
    if other is not None:
        ax.plot(My_sec2/1e6, zDepth,  ':' ,label=otherLab)
    ax.set_xlabel('Section Moment My [MNm]')
    ax.tick_params(direction='in')
    if F_sec is not None:
        ax.legend()
    fig.suptitle(stat)
    return fig

def sec_plot(vTime, zDepth, sec, ref=None, other=None, stat='mean', tRange=None, label='Section Force', fig=None, refLab='OpenFAST', otherLab='Other'):
    if tRange is None:
        IT = np.arange(0, sec.shape[1])
    else:
        IT = np.logical_and(vTime>tRange[0], vTime<tRange[1])
        if len(IT)==0:
            IT = np.arange(0, sec.shape[1])
    fstat = {'mean':np.mean, 'max':np.max, 'std':np.std, 'meanabs':meanabs}[stat]
    sec0 = fstat(sec[:,IT], axis=1)
    if ref is not None:
        sec1 = fstat(ref[:,IT], axis=1)
    if other is not None:
        sec2 = fstat(other[:,IT], axis=1)

    if fig is None:
        fig, ax = plt.subplots(1, 1, sharey=True, figsize=(6.4,4.8)) # (6.4,4.8)
        fig.subplots_adjust(left=0.14, right=0.95, top=0.95, bottom=0.12, hspace=0.20, wspace=0.20)
    else:
        ax=fig.axes[0]

    if ref is not None:
        ax.plot(sec1, zDepth,  '-' , c=colRef, lw=LWRef, label=refLab)
    ax.plot(sec0,     zDepth,  '--', c=colSim, lw=LWSim, label='YAMS')
    if other is not None:
        ax.plot(sec2, zDepth,  '-' , c=colRef, lw=LWRef, label=refLab)
    ax.tick_params(direction='in')
    ax.legend()

    ax.set_xlabel(label)
    ax.set_ylabel('Vertical position [m]')
    fig.suptitle(stat)
    return fig


