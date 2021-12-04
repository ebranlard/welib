import numpy as np
import pandas as pd
import scipy.interpolate as si
from scipy.optimize import minimize_scalar

# ---
def interp2d_pairs(*args,**kwargs):
    """ Same interface as interp2d but the returned interpolant will evaluate its inputs as pairs of values.
    """
    # Internal function, that evaluates pairs of values, output has the same shape as input
    def interpolant(x,y,f):
        x,y = np.asarray(x), np.asarray(y)
        return (si.dfitpack.bispeu(f.tck[0], f.tck[1], f.tck[2], f.tck[3], f.tck[4], x.ravel(), y.ravel())[0]).reshape(x.shape)
    # Wrapping the scipy interp2 function to call out interpolant instead
    return lambda x,y: interpolant(x,y,si.interp2d(*args,**kwargs))


def Paero(WS, Pitch, Omega, R, rho_air, fCP):
    """ Taero returns the aerodynamic power
         -lambda (tip speed ratio)
         - pitch 
         - R : the blade radius
         - fCP : an interpolant for CP(Pitch,lambda)
         - rho_air : the air density
    """
    Lambda = Omega * R / WS
    CP     = fCP(Pitch,Lambda)
    P      = 1/2*rho_air*np.pi*R**2*WS**3*CP
    return P


def Qaero(WS, Pitch, Omega, R, rho_air, fCP):
    """ Qaero returns the aerodynamic torque
         - Omega [rad/s]
         - pitch [deg]
         - R : the blade radius
         - fCP : an interpolant for CP(Pitch,lambda)
         - rho_air : the air density
    """
    Lambda = Omega * R / WS
    CP = fCP(Pitch,Lambda)
    Q = 1/2*rho_air*np.pi*R**2*WS**3/Omega*CP
    return Q

def Taero(WS, Pitch, Omega, R, rho_air, fCT):
    """ Taero returns the aerodynamic thrust of a given turbine
         -lambda (tip speed ratio)
         - pitch 
         - R : the blade radius
         - fCP : an interpolant for CP(Pitch,lambda)
         - rho_air : the air density
    """
    Lambda = Omega * R / WS
    CT = fCT(Pitch,Lambda)
    T = 1/2*rho_air*np.pi*R**2*WS**2*CT
    return T



class TabulatedWSEstimator():

    def __init__(self, R=None, rho_air=1.225, fst_file=''):
        if len(fst_file)>0:
            import welib.weio as weio
            import weio.fast_input_deck
            fst=weio.fast_input_deck.FASTInputDeck(fst_file, )
            R       = fst.ED['TipRad']
            if fst.AD is None:
                raise Exception('AeroDyn file not read but needed for wind speed estimator, while reading {}'.format(fst_file))
            rho_air = fst.AD['AirDens']
        self.R       = R
        self.rho_air = rho_air
        pass

    def load_files(self,LambdaFile=None,PitchFile=None,CPFile=None,CTFile=None,OperFile=None,base=None,suffix=''):
        if base is not None:
            LambdaFile = base+'_Lambda'+suffix+'.csv'
            PitchFile  = base+'_Pitch'+suffix+'.csv'
            CPFile     = base+'_CP'+suffix+'.csv'
            CTFile     = base+'_CT'+suffix+'.csv'
            OperFile   = base+'_Oper'+suffix+'.csv'
        self.PITCH  = pd.read_csv(PitchFile ,header = None).values.ravel()
        self.LAMBDA = pd.read_csv(LambdaFile,header = None).values.ravel()
        self.CP     = pd.read_csv(CPFile,header     = None).values
        self.CP[self.CP<=0]=0
        if CTFile is not None:
            self.CT     = pd.read_csv(CTFile,header     = None).values
            self.CT[self.CT<=0]=0
        else:
            self.CT=None
        if OperFile is not None:
            import welib.weio as weio
            self.Oper = weio.read(OperFile).toDataFrame()
            #print(self.Oper)
            self.WS   =self.Oper['WS_[m/s]'].values
            self.Omega=self.Oper['RotSpeed_[rpm]'].values*2*np.pi/60
            self.RtAeroMxh=self.Oper['RtAeroMxh_[kN-m]'].values*1000
            self.OmegaRated=np.max(self.Omega)
            self.OmegaLow  =0.4*self.OmegaRated
            self.WSRated=np.interp(self.OmegaRated*0.98,self.Omega,self.WS)
            self.WSCutOff=28
        else:
            self.Oper=None

        self.computeWeights()

    def computeWeights(self):
        self.fCP = interp2d_pairs(self.PITCH,self.LAMBDA,self.CP,kind='cubic')
        if self.CT is not None:
            self.fCT = interp2d_pairs(self.PITCH,self.LAMBDA,self.CT,kind='cubic')
        else:
            self.fCT = None

    def Power(self,WS,Pitch,Omega):
        return Paero(WS, Pitch, Omega, self.R, self.rho_air, self.fCP)

    def Thrust(self,WS,Pitch,Omega):
        return Taero(WS, Pitch, Omega, self.R, self.rho_air, self.fCT)

    def Torque(self,WS,Pitch,Omega):
        return Qaero(WS, Pitch, Omega, self.R, self.rho_air, self.fCP)


    def estimate(self,Qaero_hat,pitch,omega, WS0, relaxation=0, WSavg=None): # TODO compute rolling average on the fly
        """
         - omega [rad/s]
         - pitch [deg]
        """
        def estim(WS0,delta, maxiter=50, tol=0.1):
            z = lambda WS : abs(Qaero_hat - Qaero(WS,pitch, omega, self.R, self.rho_air, self.fCP))
            res = minimize_scalar(z,bounds=[max(1,WS0-delta),WS0+delta],method='bounded', options={'xatol': tol, 'maxiter': maxiter})
            residual=Qaero_hat - Qaero(res.x,pitch, omega, self.R, self.rho_air, self.fCP)
            return res.x, residual

#         if WSavg is not None:
#             WS0=(WS0+WSavg)/2
        if omega<self.OmegaLow:
            ws_guess=np.interp(omega,self.Omega,self.WS)
            WS1,residual = estim(WS0,2)
            WS=(4*WS1+ws_guess)/5
        else:
            WS,residual = estim(WS0,1)
        WS= WS0*relaxation + (1-relaxation)*WS

#         if omega<self.OmegaRated*0.95:
#             #  Below rated, we have the quasi steady ws as functin of omega as a best guess
#             ws_qs=np.interp(omega,self.Omega,self.WS)
#             if omega<self.OmegaLow:
#                 WS1,residual = estim(WS0,3)
#                 WS=(4*WS1+ws_qs)/5
#             else:
#                 WS,residual = estim(WS0,3)
# 
#             if np.abs(WS-ws_qs)>4:
#                 WS,residual = estim(ws_qs,3, maxiter=1000)
# 
#             if np.abs(residual)/Qaero_hat>0.1:
#                 WS,residual = estim(ws_qs, 17, maxiter=1000)
# 
#             if np.abs(residual)/Qaero_hat>0.1:
#                 WS,residual = estim(ws_qs+10, 10, maxiter=1000)
# 
#             if np.abs(residual)/Qaero_hat>0.1:
#                 if WSavg is not None:
#                     print('NOT GOOD 1 - WS={:.1f} WSqs={:.1f} WSavg={:.1f} - om={:.2f} pitch={:.2f}'.format(WS,ws_qs,WSavg,omega,pitch))
#                 else:
#                     print('NOT GOOD 1 - WS={:.1f} WSqs={:.1f} - om={:.2f} pitch={:.2f}'.format(WS,ws_qs,omega,pitch))
#         else:
#             # above omega rated, we are between WSrated-3 and WSCutoff
#             WS,residual = estim(WS0,3)
#             if WS<self.WSRated:
#                 WSmid=(self.WSCutOff+self.WSRated)/2
#                 WS,residual = estim(WSmid, 16, maxiter=1000)
# 
# #                 if WSavg is not None:
# #                     if np.abs(WS-WSavg)>4:
# #                         WS,residual = estim(WSavg, 6, maxiter=1000)
# 
#             if np.abs(residual)/Qaero_hat>0.1:
#                 print('NOT GOOD 2 - WS={:.1f} WS0={:.1f} - om={:.2f} pitch={:.2f}'.format(WS,WS0,omega,pitch))

#         WS= WS0*relaxation + (1-relaxation)*WS


        return WS



    def __repr__(self):
        s=''
        s+='<ws_estimator.TabulatedWSEstimator object> \n'
        s+='  Lambda : [min={:8.3f}, max={:8.3f}, delta={:8.4f}, n={}]  \n'.format(np.min(self.LAMBDA),np.max(self.LAMBDA),self.LAMBDA[1]-self.LAMBDA[0], len(self.LAMBDA))
        s+='  Pitch  : [min={:8.3f}, max={:8.3f}, delta={:8.4f}, n={}]  \n'.format(np.min(self.PITCH) ,np.max(self.PITCH) ,self.PITCH[1]-self.PITCH[0]  , len(self.PITCH))
        s+='  CP     : [min={:8.3f}, max={:8.3f}, n={}x{}]  \n'.format(np.min(self.CP),np.max(self.CP),self.CP.shape[0],self.CP.shape[1])
        if self.CT is not None:
            s+='  CT     : [min={:8.3f}, max={:8.3f}, n={}x{}]  \n'.format(np.min(self.CT),np.max(self.CT),self.CT.shape[0],self.CT.shape[1])
        s+='  R      : {}  \n'.format(self.R)
        s+='  rho    : {}  \n'.format(self.rho_air)
        return s


if __name__=='__main__':
    import pandas as pd
    import matplotlib.pyplot as plt
    from spectral import fft_wrap
    import welib.weio as weio
    # --- Parameters
    # InputFile = 'GeneratorDynamics.outb'
    InputFile = 'DLC120_ws13_ye000_s1_r1.outb'

    # --- Turbine data
    turbine = dict()
    turbine['R']         = 63
    turbine['rho_air']   = 1.225
    g=9.81

    # --- Reading aerodynamic data for the turbine
    PITCH  = pd.read_csv('PITCH_data.csv',header  = -1).values
    LAMBDA = pd.read_csv('LAMBDA_data.csv',header = -1).values
    CP     = pd.read_csv('CP_data.csv',header     = -1).values
    CT     = pd.read_csv('CT_data.csv',header     = -1).values
    # Create the interpolant for CP and CT, CP(pitch,lambda) (same interface interp2d) 
    turbine['fCP'] = interp2d_pairs(PITCH,LAMBDA,CP,kind='cubic')
    turbine['fCT'] = interp2d_pairs(PITCH,LAMBDA,CT,kind='cubic')



    # --- Reading in some "measuremensts" or "simulation"
    # --- Removing units from columns
    df=weio.read(InputFile).toDataFrame()
    df.columns = [  v.split('_[')[0] for v in df.columns.values] 
    time      = df['Time'].values
    genspeed  = df['GenSpeed'].values * 2*np.pi/60 # converted to rad/s
    rotspeed  = df['RotSpeed'].values * 2*np.pi/60 # converted to rad/s
    thrust    = df['RotThrust']*1000
    gentq     = df['GenTq']*1000*97                # Convert to rot torque
    azimuth   = df['Azimuth']                      # deg
    windspeed = df['Wind1VelX']
    pitch     = df['BldPitch3']
    rottorq   = df['RotTorq']*1000
    rottorq2  = df['RtAeroMxh']
    thrust2   = df['RtAeroFxh']




    # --- Evaluate the interpolant on each pairs of x and y values
    F = Taero(windspeed, pitch, rotspeed, turbine['R'], turbine['rho_air'], turbine['fCT'])
    
    Q = Qaero(windspeed, pitch, rotspeed, turbine['R'], turbine['rho_air'], turbine['fCP'])
    
    # --- normalize F, Q, rottorq2 and thrust2 to imrpove spectra analysis
    Fnorm = F - np.average(F)
    Qnorm = Q - np.average(Q)
    thrust2norm = thrust2 - np.average(thrust2)
    rottorq2norm = rottorq2 - np.average(rottorq2)
    # --- Import functio for calculating spectra from Q and F data
#    f1, S1, Info  = fft_wrap(time,Fnorm,output_type='amplitude',averaging='Welch')
#    f2, S2, Info  = fft_wrap(time,rottorq2norm,output_type='amplitude',averaging='Welch')
    f1, S1, Info  = fft_wrap(time,Fnorm,output_type='PSD',averaging='Welch')
    f2, S2, Info  = fft_wrap(time,thrust2norm,output_type='PSD',averaging='Welch')
    
    
    
    # --- Figures
    
    # -- Figure 1: Qaero
    plt.figure(num=1, figsize=(11, 6), dpi=80, facecolor='w', edgecolor='k')
    plt.plot(time,Q/1000,'r--',label='Qaero calc')
    plt.plot(time,rottorq2/1000,'k',label = 'Qaero fast')
    plt.gca().legend()
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.xlabel('time (s)')
    plt.ylabel('Qaeo (kN.m)')
    
    # --- Figure 2: Faero
    plt.figure(num=2, figsize=(11, 6), dpi=80, facecolor='w', edgecolor='k')
    plt.plot(time,F      ,'c--',label='Faero calc')
    plt.plot(time,thrust2,'m'  ,label='Faero (FAST)')
    plt.plot(time,thrust ,'b:' ,label='Faero (FAST Fstruct - Fweight)')
    plt.gca().legend()
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.xlabel('time (s)')
    plt.ylabel('Faero (kN)')

    # --- Figure 3: Spectra between Faero calc and Faero FAST
    plt.figure(num=3, figsize=(11, 6), dpi=80, facecolor='w', edgecolor='k')
    plt.plot(f1,S1,'r:',label='Faero calc')
    plt.plot(f2,S2,'k--',label='Faero (FAST)')
    plt.gca().legend()
    plt.axvline(x=.201)
    plt.axvline(x=3*.201)
    plt.axvline(x=9*.201)
    
    plt.xlim([-.0001, 5])
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Power Spectral Density (Welch Avg.)') 
    plt.yscale('log')
