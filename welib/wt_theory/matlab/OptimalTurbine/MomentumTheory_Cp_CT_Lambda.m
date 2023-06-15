%% Documentation   
% Contact: E. Branlard 
% 
% Used to generate figures 9.8 of [1]

% Reference:
%  [1]  Branlard 2017 Wind turbine aerodynamics and vorticity-based methods
% 
%% Initialization
clear all; close all; clc; % addpath()
restoredefaultpath;
addpath(genpath('C:/Config/path/MatlabPath/libs/')) % http://github.com/ebranlard/matlab-path.git

%% Parameters

% setFigurePath('/work/publications/book/figsdump/')
% setFigurePath('/home/manu/Dropbox/springer-book/figs_ad/')
% setFigureFont('14')
% I explain this a bit in my book section Maximum Power extraction 1D momentum+rotation


%
R=100;
nr=200;

n2=10000; % for second method
n3=100000; % for third method

vlambdaRef=[0.5 1 1.5 2 2.5 5 7.5 10];
CPRef=[0.288 0.416 0.480 0.512 0.532 0.570 0.582 0.593]; %WilsonLissaman 1974 p57
CPMartin=[0.486 0.703 0.811 0.865 0.899 0.963 0.983 0.987]*16/27;% martin p39


options = optimset('TolX',1e-11);
%% Overview of the function whose zeroes are sought
% figure, hold all
% for x=vlambdaRef
%     va=linspace(0.25,1/3,100);
%     plot(va,16*va.^3-24*va.^2+va.*(9-3*x^2)-1+x^2)
% end
% plot(va,va*0,'k')
% ylim([-0.01 0.01])
% 

%% That's one way of doing it
tic()
disp('First method...')
vlambda=unique(sort([vlambdaRef 10.^-(9:-1:1) 0.1:0.5:1 1:15]));
vr=linspace(0,R,nr);
a=zeros(1,nr);
CP=zeros(1,length(vlambda));
CT=zeros(1,length(vlambda));
for il=1:length(vlambda)    
    lambda=vlambda(il);
    lambda_r=vr/R*lambda;  
    for e = 1:nr
        x=lambda_r(e);
        a(e)=fzero(@(a) 16*a^3-24*a^2+a*(9-3*x^2)-1+x^2,0.3,options); %martin 8.1 combined from 4.32 and 4.38
    end
    aprime=(1-3*a)./(4*a-1);
    CP(il)=8/lambda^2*trapz(lambda_r,aprime.*(1-a).*lambda_r.^3);    % martin 4.30
    CT(il)=8/lambda^2*trapz(lambda_r,a.*(1-a).*lambda_r);    % See my book simplified 2D with rotation
end
toc()

%% Another way of doing it
tic()
disp('Second method...')
vlambda2=unique(sort([vlambdaRef linspace(0,max(vlambda),n2)]));
nn=length(vlambda2);
a2=zeros(1,nn);
aprime2=zeros(1,nn);
for il=1:nn
    x=vlambda2(il);
    a2(il)=fzero(@(a) 16*a^3-24*a^2+a*(9-3*x^2)-1+x^2,0.3,options);  %martin 8.1 combined from 4.32 and 4.38
end
aprime2=(1-3*a2)./(4*a2-1);
CP2=8./(vlambda2.^2).*cumtrapz(vlambda2,aprime2.*(1-a2).*vlambda2.^3);
CT2=8./(vlambda2.^2).*cumtrapz(vlambda2,a2.*(1-a2).*vlambda2);
toc()

%% and yet another way of doing it, the deterministic way
tic()
disp('Third method...')
a3=linspace(0.25+10^-9,1/3,n3); % why do i start at 0.25, because of the expression for aprime
aprime3=(1-3*a3)./(4*a3-1);
vlambda3=sqrt(a3.*(1-a3)./(aprime3.*(1+aprime3)));
[vlambda3 Isort]=sort(vlambda3);% in fact not needed...
CP3=8./(vlambda3.^2).*cumtrapz(vlambda3,aprime3(Isort).*(1-a3(Isort)).*vlambda3.^3);
CT3=8./(vlambda3.^2).*cumtrapz(vlambda3,a3(Isort).*(1-a3(Isort)).*vlambda3);
a3=a3(Isort);
aprime3=aprime3(Isort);
toc()


figure,hold all, 
plot(vlambda2,aprime2.*(1-a2).*vlambda2.^3)
plot(vlambda3,aprime3(Isort).*(1-a3(Isort)).*vlambda3.^3)

% CP
figure,hold on,grid,box
plot(vlambda,vlambda*0+16/27 ,'k--')
plot(vlambda2,CP2 ,'k-')
% plot(vlambda,CP,'b+')
% plot(vlambda3,CP3,'r--')
% plot(vlambdaRef,CPRef,'ko')
legend('Betz limit','Ideal rotor with wake rotation')
% ,'Same but different way of thinking it','Same but deterministic','Wilson Lissaman',0)
xlim([0 10])
ylim([0 0.62])
xlabel('Tip speed ratio \lambda [-]')
ylabel('C_P [-]')
title('MomentumTheoryActuatorDiskOptimalCP')



% CT
figure,hold on,grid,box
plot(vlambda,vlambda*0+8/9 ,'k--')
plot(vlambda2,CT2 ,'k-')
% plot(vlambda,CT,'b+')
% plot(vlambda3,CT3,'r--')
legend('Betz limit','Ideal rotor with wake rotation')
% ,'Same but different way of thinking it','Same but deterministic',0)
xlim([0 10])
xlabel('Tip speed ratio \lambda [-]')
ylabel('C_T [.]')
title('MomentumTheoryActuatorDiskOptimalCT')

% a
figure,hold on,grid,box
plot(vlambda,vlambda*0+1/3 ,'k--')
plot(vlambda2,a2 ,'k-')
% plot(vlambda3,a3,'r--')
legend('Betz limit','Ideal rotor with wake rotation')
% ,'Same but deterministic',0)
xlim([0 10])
ylim([0 0.35])
xlabel('Local tip speed ratio \lambda_r [-]')
ylabel('Axial induction factor a [-]')
title('MomentumTheoryActuatorDiskOptimala')

% aprime
figure,hold on,grid,box
plot(vlambda2,aprime2 ,'k-')
% plot(vlambda3,aprime3,'r--')
ylim([0 0.1])
legend('Ideal rotor with wake rotation')
% ,'Same but deterministic',0)
xlim([0 10])
title('MomentumTheoryActuatorDiskOptimalaprime')
xlabel('Local tip speed ratio \lambda_r [-]')
ylabel('Tangential induction factor a'' [-]')


%% 
I =whichvalue(vlambda2,vlambdaRef);
I2=whichvalue(vlambda,vlambdaRef);
I3=whichvalue(vlambda3,vlambdaRef);

[vlambda3(I3)' CP2(I)' CP(I2)' CP3(I3)' CPRef' CPMartin']

fprintf('\n\nlambda\t CP2\t CPref CPMartin\n');
fprintf('%.3f\t %.4f\t %.4f %.3f\n', [vlambda2(I)' CP2(I)' CPRef' CPMartin']');
%% Cp
vlambdaAD=linspace(0,30,151);
vaAD=interp1(vlambda3(1:end-2),a3(1:end-2),vlambdaAD);
vaprimeAD=interp1(vlambda3(1:end-2),aprime3(1:end-2),vlambdaAD);
vCPAD=interp1(vlambda3(1:end-2),CP3(1:end-2),vlambdaAD);
vCTAD=interp1(vlambda3(1:end-2),CT3(1:end-2),vlambdaAD);
vCPAD(1)=0;
vCTAD(1)=0;
vaAD(1)=0.25;



% save([./MomentumTheoryActuatorDisk.mat'],'vlambdaAD','vCPAD','vCTAD','vaAD','vaprimeAD');

%%
