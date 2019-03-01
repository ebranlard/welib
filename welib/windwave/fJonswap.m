function [ f,a ] = fJonswap( Hs,Tp,df,fHighCut )
%FJONSWAP Summary of this function goes here
%   Detailed explanation goes here


fp=1/Tp; % from previous question

%% get gamma value

if (Tp/sqrt(Hs)<=3.6)==1
    gamma=5;
elseif (Tp/sqrt(Hs)>=3.6)==1 && (Tp/sqrt(Hs)<=5)==1
    gamma=exp(5.75-1.15*(Tp/sqrt(Hs)));
elseif  (Tp/sqrt(Hs)>3.6)==1
    gamma=1;
else
    disp('Something is wrong with your inputs, model not valid')
end
gamma=3.3;


%% sigma
vFreq=df:df:fHighCut;
for iFreq=1:length(vFreq)
    freq=vFreq(iFreq);
    if (freq<=fp)==1
        sigma=0.07;
    elseif (freq>fp)==1
        sigma=0.09;
    end
    S(iFreq)=0.3125*Hs^2*Tp*((freq/fp)^(-5))*exp(-1.25*(freq/fp)^(-4))*(1-0.287*log(gamma))*(gamma)^(exp(-0.5*(((freq/fp)-1)*(1/sigma))^2));  
end

f=vFreq; 
a=S;










 
