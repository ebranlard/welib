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
BEMParam.correction='Spera';

tic()
%% calculating the power for each wind speed
Vspeed=(4:1:25)+0.5;
CP=zeros(1,length(Vspeed));
P=zeros(1,length(Vspeed));
for i=1:length(Vspeed)
    Aero.Wind.V0(3)=Vspeed(i);
    lambda=omega*R/Vspeed(i);
    % calling the BEM code 
    BEM = fBEMsteady();
    CP(i)=BEM.CP;
    P(i)=BEM.Power;
end    
toc()
%%
figure(51)
plot(Vspeed,P/1000)
xlabel('V [m/s]');
ylabel('P [kW]')
grid
box

figure(52)
plot(omega*R./Vspeed,CP)
xlabel('lambda [.]');
ylabel('Cp [.]')
grid
box
%%

hpy=24*365.25;% hours per year
hoursRepartition=wblpdf(Vspeed,6,2)*hpy;
AEO=sum(hoursRepartition.*P)/1000/1000; % in MW.h/year
AEO
