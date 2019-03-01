function [S,F1] = fSpectrumAvg(u,fs,K,omega)


    n=floor(length(u)/K); % compute the size of the sample
    j=1;
    F1 = [0:n-1]/n*fs; % Nyquist frequency
    
    %-----------Splitting the time series---------------------
    for i=1:K
        x(1:length(u(j:j+n-1)),i)=u(j:j+n-1)';
        j=j+n;
    end
    
    %-----------------calculating the spectrum-----------------
    for i=1:length(x(1,:))
        X1(1:length(x(:,i)),i)=(abs((fft(x(:,i))')).^2)/(omega*length(x(:,1)));
    end
    
    %--------------- mean of X1 ----------------------------------
    for i=2:length(x(:,1))
        S(i)=nanmean(X1(i,:));
    end
end



