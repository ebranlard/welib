function [ Uw, fraw, Sraw,fKailman,SKailman,fBin,SBin,fWelch,SWelch ] = fStocWind(Tend )
addpath(genpath('fft'))
%% Get theoretical Kaimal spectrum
% inputs to the spectrum
V10=8;
I=0.14;
l=340.2;

% Tmanu=0:0.1:500
T=Tend; % duration of time serie in sec
% df=1/dt;
fHighCut= 10; % cut off frequency in Hz
df=1/T;
f = 0:df:fHighCut; % freq. vector [Hz]

% kaimal spectrum
Sw = 4*I^2*V10*l ./ (1+(6*f*l)/V10).^(5/3);

%% Generate time serie
% Create two random normal distributed variables
N=T*fHighCut+1;
a_m=0;
a_std=1;
a=a_m+a_std*randn(floor(N/2)+1,1);

b_m=0;
b_std=1;
b=b_m+b_std*randn(floor(N/2)+1,1);

% Create first part of the complex vector
xlhalf1 =a+b*j;
xlhalf1(1)=0;

for i=1:N/2
  wj=2*pi*f(i)/2; % omega
  S_Kaimal_j = 4*I^2*V10*l ./ (1+(6*wj*l)/V10).^(5/3); % theoretical spectrum S(omega)
  sigma_K_j=sqrt(T/(2*pi)*S_Kaimal_j);
  xlhalf1(i)=xlhalf1(i)*sigma_K_j;
end

% Create second part of the complex vector
xlhalf2=flipud(xlhalf1(2:end));

% Total complex vector
xl=[xlhalf1;xlhalf2];

% Fourier inverse
Uw=ifft(xl);

% Remove zero imaginairy part and normalize
Uw=sqrt(2)*pi*fHighCut*real(Uw);

%% Plot spectrum and verify agreement
[Sm,F1] = fSpectrumAvg(Uw,fHighCut,1,2*pi*fHighCut);
NN = length(Sm);
[Smav,F1av]=fSpectrumCalc(Sm(2:round(NN/2)),F1(2:round(NN/2)));

A=trapz(F1av,Smav)*2*pi*2;
B=var(Uw);
C=mean(Uw);
verif=A-B; % check diff integrate with var


% other method using pwelch 
fN= length(Uw)/(2*T); % frequency
nw = 16; % Number of windows (adjust until noise is removed sufficiently)
nfft = round(length(Uw)./nw); % Choose window length
% [palt,falt] = pwelch(Uw,nfft); % Normalized results based on Welch's method
% nfac = fN./pi; falt = falt.*nfac; palt = palt./nfac; % De-normalize results

% figure
% loglog(F1(2:floor(NN/2)),2*pi*2*Sm(2:floor(NN/2)),'-','Color',[0.6 0.6 0.6])  % used to be y
% hold on
% loglog(falt,palt,'-','Color',[0.25 0.25 0.25]) % used to be g
% loglog(f,Sw,'k-','linewidth',2) % used to be b
% loglog(F1av,2*pi*2*Smav,'-','linewidth',1.5,'Color',fColrs(2)) % used to be k
% title('GeneratedTimeSeriesKaimalSpectrum')
% xlabel('Frequency [Hz]')
% ylabel('PSD')
% grid on
% lg=legend('Raw spectrum','pWelch spectrum','Kaimal spectrum','Bin-averaged spectrum');
% set(lg,'location','southwest')
% % axis tight


%% Preparing outputs
fraw=F1(2:floor(NN/2));
Sraw=2*pi*2*Sm(2:floor(NN/2));
fKailman=f;
SKailman=Sw;
fBin=F1av;
SBin=2*pi*2*Smav;
% fWelch=falt;
% SWelch=palt;
fWelch=0;
SWelch=0;


end

