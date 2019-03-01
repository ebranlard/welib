%% question 5
InitClear
InitRotor
InitSystemProperties
% 
% % geometry of the turbine for question 1 and 4
% %values for lambda=10.11  
% %question 1 : theoretical twist
% beta1=[23.6845 14.3652 4.9604 1.4195 -0.4047 -1.5122 -1.6833 -1.8419 -1.9893 -2.1266 -2.1918];
% %question 4 : twist optimized element per element
% beta4=[23.8269 14.2840 4.6998 1.3165 -0.4972 -2.1530 -2.1530 -2.6530 -3.0223 -3.4203 -3.7896];
% c=[0.3178 0.2526 0.1477 0.1017 0.0772 0.0612 0.0580 0.0539 0.0479 0.0371 0.0276];
% %values for lambda =8
% beta0 =[ 28.2268   18.5799    7.7040    3.3575    1.0798   -0.3128   -0.5285   -0.7286   -0.9146   -1.0880   -1.1703];
% c0 =[   0.4190    0.3660    0.2304    0.1616    0.1233    0.0966    0.0905    0.0829    0.0721    0.0545    0.0399];
      
    
%% BEM code parameters
BEMParam.nbIt=200;
BEMParam.alphaCrit=0.05;
BEMParam.relaxation=0.3;
BEMParam.correction='Spera';
BEMParam.BigStorage=0;
%% calulating the power for each vind speed
Vspeed=(4:25)+0.5;
CPv=zeros(1,length(Vspeed));
Pv=zeros(1,length(Vspeed));
for i=1:length(Vspeed)
    Aero.Wind.V0=[0; 0; Vspeed(i)];
    % calling the BEM code 
    BEM = fBEMsteady(0,0,0);
    %Pv(i)=CP*(0.5*rho*Vspeed(i)^3*pi*R^2);
    Pv(i)=BEM.Power;
    CPv(i)=BEM.CP;
end    

%%
figure(51)
plot(Vspeed,Pv/1000)
xlabel('V [m/s]');
ylabel('P [kW]')
grid
box

figure(52)
plot(Rotor.Omega*Rotor.R./Vspeed,CPv)
xlabel('lambda [.]');
ylabel('Cp [.]')
grid
box
%%

hpy=24*365.25;% hours per year
hoursRepartition=wblpdf(Vspeed,6,2)*hpy;
AEO=sum(hoursRepartition.*Pv)/1000/1000; % in MW.h/year
AEO

%%
% alpha=[4.327 6.712 6.716 5.3  3.123 4.327 6.712 6.716 5.3  3.123 4.327 6.712 6.716 5.3  3.123]
% 
% plot(alpha,'+')
% xlabel('Iterations')
% ylabel('alpha [deg]')
% xlim([0 15])



%AEO beta4 = 10.149
%AEO beta1 = 10.193
%AEO beta0 = 10.