%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% TWIST OPTIMIZATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialization for the BEM code
InitClear
InitRotor
InitSystemProperties
R=Rotor.R;
r=Rotor.r;
omega=Rotor.Omega;

%% BEM code parameters
BEMParam.nbIt=200;
BEMParam.alphaCrit=0.01;
BEMParam.relaxation=0.3;
BEMParam.correction='Glauert';

tic()
%% calulating the power for each vind speed
Vspeed=(4:25)+0.5;
CPv=zeros(1,length(Vspeed));
Pv=zeros(1,length(Vspeed));
for i=1:length(Vspeed)
    V0=Vspeed(i);
    lambda=omega*R/V0;
    % calling the BEM code 
    BEM = fBEMsteady(0,0,0);
    CPv(i)=BEM.CP;
    Pv(i)=BEM.P;
end    
toc()
%%
figure(51)
plot(Vspeed,Pv/1000)
xlabel('V [m/s]');
ylabel('P [kW]')
grid
box

figure(52)
plot(omega*R./Vspeed,CPv)
xlabel('lambda [.]');
ylabel('Cp [.]')
grid
box
%%

hpy=24*365.25;% hours per year
hoursRepartition=wblpdf(Vspeed,6,2)*hpy;
AEO=sum(hoursRepartition.*Pv)/1000/1000; % in MW.h/year
AEO

% %%
% alpha=[4.327 6.712 6.716 5.3  3.123 4.327 6.712 6.716 5.3  3.123 4.327 6.712 6.716 5.3  3.123]
% 
% plot(alpha,'+')
% xlabel('Iterations')
% ylabel('alpha [deg]')
% xlim([0 15])



%AEO beta4 = 10.149
%AEO beta1 = 10.193
%AEO beta0 = 10.