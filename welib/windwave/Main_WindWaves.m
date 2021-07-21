clear all;
close all;
clc;


h=30; % water depth % positive !!!! NOTE
g=9.81; % gravity
zBot=-5; % depth of floater negative!!!! NOTE
bWheeler=0; % if true (1) use wheeler stretching
bFloating=0; % if true (1) the moment is computed at z=0 and only from zBot and not from -h
bVzTo0=1; % if true (1) vz goes to 0

%%  Jonswap for 50 year sea state
disp('Jonswap spectrum')
Hs=8.1;
Tp=12.70;
SigmaApprox=Hs/4;

nt=3601;
vt=linspace(0,3600,nt); % time vector [s]
dt = vt(2)-vt(1);          % Sample time
Fs = 1/dt;          % Sampling frequency [Hz}
df=1/max(vt);      % smallest frequency 
fHighCut=Fs/2;     % Nyquist 
[ vf_Jonswap,S_Jonswap ] = fJonswap( Hs,Tp,df,fHighCut );

% integration of the spectrum
Area=trapz(vf_Jonswap,S_Jonswap);
normalization_factor=(Hs/4)^2/Area;
S_Jonswap_=normalization_factor*S_Jonswap;
Area_=trapz(vf_Jonswap,S_Jonswap_);
figure(222)
plot(vf_Jonswap,S_Jonswap,'k-'),hold all
plot(vf_Jonswap,S_Jonswap_,'b.')
grid on
xlim([0 0.3])
legend('Function','Normalized')
xlabel('Frequencies [Hz]')
ylabel('JONSWAP Spectral density S [m2.s]')


%% Summation of various wave components using Jonswap spectrum
disp('Generation of Wave')
N       = length(vf_Jonswap)      ;  % number of component for the sum that will be used to compute eta(t)  
vphases_Jonswap = rand(1,N)*2*pi ;  % random phases between [0 2*pi]
vA_Jonswap      = sqrt(2*S_Jonswap*df)   ;  % scaled amplitudes according to Jonswap spectrum
% Dispersion for all frequencies
[ vk_Jonswap ] = fgetDispersion( vf_Jonswap,h,g );


vf      = vf_Jonswap      ;  % vector of frequencies
vphases = vphases_Jonswap ;  % random phases between [0 2*pi]
vA      = vA_Jonswap      ; 
vk      = vk_Jonswap      ; 
% param
nz=1;
x=0;
[ eta Stored_u Stored_dudt]=  fWaveKin(vt,vf,vk,vA,vphases,x,h,zBot,nz,bWheeler,bFloating,bVzTo0);

figure,hold all,box on, grid on
plot(vt,eta),xlim([0 3600])
xlabel('Time t [s]'), ylabel('Elevation [m]');
title('ElevationTime')


%% Verification of spectrum 
disp('Verification of spectrum')
% Fs = 1/T;          % Sampling frequency [Hz}
% df=1/max(vt);      % smallest frequency 
% fHighCut=Fs/2;     % Nyquist 
% Matlab method
nt=length(vt);                % Length of signal
NFFT = 2^nextpow2(nt); % Next power of 2 from length of y
Y = fft(eta,NFFT)/nt;
Ab=2*abs(Y(1:NFFT/2+1)); % Single sided Amplitude spectrum
fb = Fs/2*linspace(0,1,NFFT/2+1);
psb=Ab.^2/2/df;

% Method2
fRef=[1:nt]*df-df;
aRef=abs(fft(eta))/nt; % double sided amplitude spectrum
psEtaRef=aRef.^2/df;psEta(1)=0;

figure,hold all
plot(fb,psb)
plot(fRef,psEtaRef*2)
plot(vf_Jonswap,S_Jonswap,'k')
legend('Matlab','Ref','Jonswap S')
xlim([0 0.3]);
ylabel('PSD of \eta [m^2/Hz]'),
xlabel('Hz');


%% Plot single-sided amplitude spectrum.
figure,grid on,hold all,box on
hold all
plot(fb,Ab) 
plot(fRef,aRef*2) 
plot(vf_Jonswap,vA_Jonswap,'k','LineWidth',2) 
xlim([0 0.5])
ylabel('Single Sided Amplitude Spectrum of eta(t) - a')
xlabel('Frequency (Hz)')
legend('From Time series','From JONSWAP spectrum')
title('SpectrumComparison')






%% Wind and wave spectra
clear all;
disp('Wind and Wave spectra');
% ------------------- Jonswap spectrum
h=30; % water depth % positive !!!! NOTE
g=9.81; % gravity
nt=3601;
vt=linspace(0,3600,nt); % time vector [s]
T = vt(2)-vt(1);          % Sample time
Fs = 1/T;          % Sampling frequency [Hz}
df=1/max(vt);      % smallest frequency 
fHighCut=Fs/2;     % Nyquist 
Hs=6;
Tp=10;
SigmaApprox=Hs/4;
[ vf_Jonswap,S_Jonswap ] = fJonswap( Hs,Tp,df,fHighCut ); % returns vectos of wave frequencies and amplitudes
N       = length(vf_Jonswap)      ;  % number of component for the sum that will be used to compute eta(t)  see slides 02 page 10
vphases_Jonswap = rand(1,N)*2*pi ;  % random phases between [0 2*pi]
vA_Jonswap      = sqrt(2*S_Jonswap*df)   ;  % scaled amplitudes according to Jonswap spectrum
% Dispersion for all frequencies
[ vk ] = fgetDispersion( vf_Jonswap,h,g );

% ------------------- Kaimal spectrum 
[ vUw_Kaimal, fraw, Sraw,fKailman,SKailman,fBin,SBin ]  = fStocWind(2*vt(end));
vUw_Kaimal=vUw_Kaimal(1:length(vt)) ;


% --------- Stochastic Times series
eta=zeros(1,length(vt));
for it=1:length(vt)
    eta_Jonswap(it) =  sum(vA_Jonswap.*cos((2*pi*vf_Jonswap)*vt(it) - vk*0 +vphases_Jonswap)); % water elevation, summation of waves
end
figure,plot(vt,vUw_Kaimal ,'Color',[0 0 0.7]),box on, grid on,title('Wind Time Series'),xlabel('t [s]'),ylabel('Wind speed [m/s]');
figure,plot(vt,eta_Jonswap,'Color',[0 0 0.7]),box on, grid on,title('Wave Time Series'),xlabel('t [s]'),ylabel('Water elevation [m]');

% [ S,f ] = fSpectrum( vt(Iselect),vUw_Kaimal(Iselect,1),1,{'Kaimal wind','','Uw'});
[ SBin_eta,fBin_eta,Sraw_eta,fraw_eta ] = fSpectrum( vt,eta_Jonswap,0,{'Jonswap','','Eta'});

% --------- Spectra of Stochastic Times series
figure,
loglog(fraw,Sraw,'-','Color',[0.66 0.66 0.66]) 
hold on, % stupid matlab
% loglog(fWelch,SWelch,'-','Color',[0.35 0.35 0.35])
loglog(fKailman,SKailman,'k-','linewidth',2)
loglog(fBin,SBin,'r-')
title('GeneratedTimeSeriesKaimalSpectrum')
xlabel('Frequency [Hz]')
ylabel('PSD of Uw [m^2/s^2 s]')
grid on
lg=legend('Raw spectrum','Kaimal spectrum','Bin-averaged spectrum');
set(lg,'location','southwest')

figure,
Ikp=S_Jonswap>10^-5;
loglog(fraw_eta,Sraw_eta,'-','Color',[0.66 0.66 0.66]) 
hold on, % stupid matlab
% loglog(fWelch,SWelch,'-','Color',[0.35 0.35 0.35])
loglog(vf_Jonswap(Ikp),S_Jonswap(Ikp),'k','linewidth',2)
loglog(fBin_eta,SBin_eta,'r-')
title('GeneratedTimeSeriesJonswapSpectrum')
xlabel('Frequency [Hz]')
ylabel('PSD of eta [m^2 s]')
grid on
lg=legend('Raw spectrum','Jonswap spectrum','Bin-average spectrum');
set(lg,'location','southwest')% axis tight



%% Spectra of time series
% [ S,f ] = fSpectrum( tnwsw(Iselect),qnwsw(Iselect,1),1,{' No wave',' stoc. wind','x'});
