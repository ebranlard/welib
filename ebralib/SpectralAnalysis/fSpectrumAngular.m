function [Sp,vomega]=fSpectrumAngular(y,N,fs)
% Returns the Spectrum of the signal y, in terms of angular frequency omega
%  If requires in terms of cyclic frequency f, see fSpectrum
%
% 
y=y(:)';

vomega=2*pi*(0:floor(N/2))*fs/N;
Sp=abs(fft(y)).^2 / (2*pi*N*fs);
Sp=Sp(1:floor(N/2)+1);
% figure
% loglog(vf(2:end),Sp(2:end))

%% Manual version
% i=sqrt(-1);
% Sp=zeros(1,floor(N/2)+1);
% for l=0:floor(N/2)
%     Sp(l+1)=0;
%     for k=0:N-1
%         Sp(l+1)=Sp(l+1)+ y(k+1)*exp(-2*i*pi*l*k/N);
%     end
%     Sp(l+1)=abs(Sp(l+1))^2;
% end
% Sp=Sp/(2*pi*N*fs);


%% Vectorized Manual version
% i=sqrt(-1);
% Sp=zeros(1,floor(N/2)+1);
% vk=0:N-1;
% for l=0:floor(N/2)
%     Sp(l+1)=abs(sum(y(1:N).*exp(-2*i*pi*l*vk/N))).^2;
% end
% Sp=Sp/(2*pi*N*fs);
