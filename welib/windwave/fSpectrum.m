function [ SBin,fBin,Sraw,fraw ] = fSpectrum( t,q,bPlot,ttlstr )
if t(1)~=0
    t=t-t(1);
end

N_samples= length(q); % number of samples
dt= t(end); % time increment of 10 minutes (600 seconds)
Fs= N_samples/(dt); % frequencyy
% for the annual pattern
K=1; % splitting factor

omega=2*pi*Fs;

[S_,F_] = fSpectrumAvg(q,Fs,K,omega);
NN = length(S_);

[S,f]=fSpectrumCalc(S_(2:floor(NN/2)),F_(2:floor(NN/2)));


eigenfreq = [8.0403e-003 38.1133e-003];

%%
if bPlot==1
%     figure('name',['PSD with',ttlstr{1},' and',ttlstr{2}])
    figure()
%     loglog(F_(2:floor(NN/2)),2*pi*2*S_(2:floor(NN/2)),'-','Color',[0.66 0.66 0.66])
    loglog(F_(2:floor(NN/2)),2*pi*2*S_(2:floor(NN/2)),'-','Color',[0.66 0.66 0.66])
    hold on
%     loglog(f,2*pi*2*S,'k-','linewidth',2)
    loglog(f,2*pi*2*S,'k-','linewidth',2)
    text(eigenfreq(1),2*pi*2*S(whichvalue(f,eigenfreq(1)))*1.1,'f1','HorizontalAlignment','center','VerticalAlignment','Bottom')
    text(eigenfreq(2),2*pi*2*S(whichvalue(f,eigenfreq(2)))*1.1,'f2','HorizontalAlignment','center','VerticalAlignment','Bottom')
%     title(['PSD with',ttlstr{1},' and',ttlstr{2},' for ',ttlstr{3}])
    xlabel('Frequency [Hz]')
%     ylabel(['PSD of ' ttlstr{3} ''])
%     title(['PSD with',ttlstr{1},' and',ttlstr{2},' for',ttlstr{3}])
%     xlabel('Frequency [Hz]')
%     ylabel(['PSD of ',ttlstr{3}],'interpreter','latex','fontsize',16)
    grid on
    legend('Raw PSD','Average PSD')
end

% Outputs
SBin=2*pi*2*S;
fBin=f;
Sraw=2*pi*2*S_(2:floor(NN/2));
fraw=F_(2:floor(NN/2));


% 
% sig=q;
% figure;subplot(2,1,1)
% plot(t,sig,'r')
% % ylabel('\eta [m]'); legend('\eta [m]')
% xlabel('t [s]')
% % Fast Fourier Transform
% dfFFT=1/t(end);
% fFFT=[1:length(t)]*dfFFT-dfFFT;
% a=abs(fft(sig))/length(t);
% psSig=a.^2/dfFFT;psSig(1)=0;
% subplot(2,1,2),hold all
% plot(fFFT,psSig)
% plot(F_(2:floor(NN/2)),2*pi*2*S_(2:floor(NN/2)),'m-')
% set(gca,'xlim',[0 0.3])
% % ylabel('PSD, \eta [m^2/Hz]'),
% % legend('PSD, \eta [m^2/Hz]')
% xlabel('Hz');


