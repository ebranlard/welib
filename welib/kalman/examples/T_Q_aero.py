import numpy as np
import scipy.interpolate as si

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

def Qaero(WS, Pitch, Omega, R, rho_air, fCP):
    """ Taero returns the aerodynamic thrust of a given turbine
         -lambda (tip speed ratio)
         - pitch 
         - R : the blade radius
         - fCP : an interpolant for CP(Pitch,lambda)
         - rho_air : the air density
    """
    Lambda = Omega * R / WS
    CP = fCP(Pitch,Lambda)
    Q = 1/2*rho_air*np.pi*R**2*WS**3/Omega*CP
    return Q



if __name__=='__main__':
    import pandas as pd
    import matplotlib.pyplot as plt
    from spectral import fft_wrap
    import weio
    # --- Parameters
    # InputFile = 'GeneratorDynamics.outb'
    InputFile = 'DLC120_ws13_ye000_s1_r1.outb'

    # --- Turbine data
    turbine = dict()
    turbine['R']         = 63
    turbine['tilt']      = 5
    turbine['RotorMass'] = 109389;
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
