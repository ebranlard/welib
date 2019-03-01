
% InitClear
% setFigurePath('./report/figs/');
% setMatFigurePath('./report/matfigs/');
% setFigureLatex(0);
% setFigureHeight('7');
% 


%% Call for Monopile, one wave
disp('Monopile hydro force calculation - One Wave')
% ------- Position
q(1)=0;% x0
q(2)=0;% theta
q(3)=0;
q(4)=0;
% ------- Component related inputs
vf= 1/12; %: vector of frequencies of the different wave component
vA=3 ; % wave amplitudes ( the wave goes from +-A) !!! NOTE
vphases=0 ;% wave phases
% ------- Environment
g=9.81; % gravity
rhow = 1025; % kg/m3 %rho water!!! NOTE
% ------- Cylinder properties
h=30; % water depth % positive !!!! NOTE
D=6  ;% cylinder diameter
zBot=-5; % depth of floater negative!!!! NOTE
Cm= 1.0; %  !!!!!!!!!!!!!!!!!!!!!!!!!!!! This is Cm not CM, CM=2 Cm=1 NOTE
CD= 1.0; % Drag coeff
% ------- Algo
nz=50;  % number of computation points over the depth
bWheeler=0; % if true (1) use wheeler stretching
bFloating=0; % if true (1) the moment is computed at z=0 and only from zBot and not from -h
bVzTo0=1; % if true (1) vz goes to 0

vk= fgetDispersion(vf,h,g); % vector of wave number as returned by fgetDisperion
% ------ Time
vt=linspace(0,1/vf,100) ; % either scalar or a vector.. but watch out q for now is scalar

[ eta Ftot Mtot Fdrag Finertia vz ut dudt dFtot dMtot dFdrag dFinertia]=  fHydroCalcFinal(vt,q,vf,vk,vA,vphases,g,rhow,h,zBot,D,Cm,CD,nz,bWheeler,bFloating,bVzTo0);

figure,hold all, grid on,box on;
plot(vt,eta)
plot(vt,ut)
plot(vt,dudt)
xlabel('Time [s]')
ylabel('')
lg=legend('$\eta$','u','du/dt');
set(lg,'location','northwest');
title('ElevationVelocityAcceleration')

figure;
subplot(1,3,1),box on,grid on
plot(vt,eta)
xlabel('Time [s]')
ylabel('Eta')
subplot(1,3,2),box on,grid on
hold all
plot(vt,Ftot)
plot(vt,Fdrag)
plot(vt,Finertia)
xlabel('Time [s]')
ylabel('Inline Force [N] ')
legend('Tot','Drag','Inertia')
subplot(1,3,3),box on,grid on
hold all
plot(vt,Mtot)
xlabel('Time [s]')
ylabel('Moment [N] ')


%% Plots of force as function of height
T=12;
bVzTo0=0; % if true (1) vz goes to 0
vt=linspace(0,T,9); % time vector
for it=1:length(vt)
    t=vt(it);
    [ eta Ftot Mtot Fdrag Finertia vz Stored_u Stored_dudt dFtot dMtot dFdrag dFinertia]=  fHydroCalcFinal(t,q,vf,vk,vA,vphases,g,rhow,h,zBot,D,Cm,CD,nz,bWheeler,bFloating,bVzTo0);

    figure(100+it),clf,hold all,grid on,box on;
    plot(dFtot,vz,'k')
    plot(dFinertia,vz,'b.')
    plot(dFdrag,vz,'r')
    if(it==1)
        lg=legend('f_{tot}','f_{inertia}','f_{drag}');
    end
    LIM=[-1 1]*85000;
    plot(LIM,LIM*0,'k'), plot(LIM,LIM*0-3,'k--'), plot(LIM,LIM*0+3,'k--')
    xlim(LIM)
    ylim([-30 4])
    xlabel('Inline Force f [N/m]');
    ylabel('Depth z [m]');
    title(sprintf('ForceWithDepthT%.2f',t/T))
%     set(gca);

end



%%  Jonswap for 50 year sea state
disp('Jonswap spectrum')
Hs=8.1;
Tp=12.70;
SigmaApprox=Hs/4;

nt=3601;
vt=linspace(0,3600,nt); % time vector [s]
T = vt(2)-vt(1);          % Sample time
Fs = 1/T;          % Sampling frequency [Hz}
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
nz=10;
% ------- Position
q(1)=0;% x0
q(2)=0;% theta
q(3)=0;
q(4)=0;
[ eta Ftot Mtot Fdrag Finertia vz Stored_u Stored_dudt dFtot dMtot dFdrag dFinertia]=  fHydroCalcFinal(vt,q,vf,vk,vA,vphases,g,rhow,h,zBot,D,Cm,CD,nz,bWheeler,bFloating,bVzTo0);

%% Verification of spectra 
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
% Areab=trapz(fb,Sb);
% Area=trapz(f,S);

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


%%
% Plot single-sided amplitude spectrum.
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


%% trying to get right power spectrum of Force 
% Ref method
fRef=[1:nt]*df-df;
aRef=abs(fft(Ftot))/nt; % double sided amplitude spectrum
psEtaRef=aRef.^2/df;psEta(1)=0;

figure,hold all
plot(fRef,psEtaRef)
xlim([0 0.3]);
ylabel('PSD of \eta [m^2/Hz]'),
xlabel('Hz');












%% Plots
figure,hold all,box on, grid on
plot(vt,eta),xlim([0 3600])
xlabel('Time t [s]'), ylabel('Elevation [m]');
title('ElevationTime')
figure,hold on,grid on,box on
plot(vt,Mtot(:)),xlim([0 3600])
xlabel('Time t [s]') , ylabel('Integrated Moment [Nm]');
title('IntegratedMoment')

figure,hold on,grid on,box on
plot(vt,Ftot(:)),xlim([0 3600])
xlabel('Time t [s]'), ylabel('Integrated Force [N]');
title('IntegratedForce')



