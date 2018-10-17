InitClear

%%
user.MainPath = 'data/test/';
user.HtcFileName = 'HtcFindPitch.htc';
user.BladeMainBodyName = 'blade1';

par.eta=0.93;
par.Pr=500000.0;
par.Pref = par.Pr/par.eta;

par.r = 1.5; % Hub radius

% AERO Parameters
par.Omega = 27.1/60*2*pi;
par.rho = 1.225;
par.Uvec = 5:1:25;


% BEM parameters
par.k = [-0.0017077,0.251163,0.0544955,0.0892074]; % glauert correction
par.relax = 0.85;  % relaxation of solution of induction
par.Induction = 1;
par.TipLoss = 1;


% load hawc2 data
Data = fReadHtcFileAll(user);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% call main function
[par res]=fFindPitch(user,par,Data);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%

load('data/Workspace.mat')
close all
P=res.all(2,:);
CP=res.all(5,:);
U=res.all(1,:);
beta=res.all(3,:);
dpdbeta=res.all(4,:);

% compute kk factor based on frozen wake
I = find(beta > 0);
temp = polyfit(beta(I(1):end),dpdbeta(I(1):end),1);
KK = temp(2)/temp(1);

figure(1)
hold all
plot(beta(I(1):end),dpdbeta(I(1):end),'k+-','LineWidth',2)
plot(beta(I(1)-1:end),beta(I(1)-1:end)*temp(1)+temp(2),'--','LineWidth',2,'Color',[0.5 0.5 0.5])
text('Interpreter','latex','String',sprintf('$$\\frac{dP}{d\\theta} = -$$ %.2f [kW / deg]',abs(dPdpitch_0)), 'Position',[2.5 -10], 'FontSize',14)



grid on
box on
xlabel('Pitch \theta [deg]')
ylabel('dP/d\theta (frozen wake) [kW/deg]')
title('dpdpitch')
%%
figure(2)
plot(U,beta,'k','LineWidth',2)
grid on
box on
xlabel('Wind speed U_\infty [m/s]')
ylabel('Pitch \theta [deg]')
title('pitchvsWS')



%%

% figure
% plot(U,P)
% 
% figure
% plot(U,CP)
% 
% figure
% plot(par.s,par.c)
% figure
% plot(par.s,par.thickness)
% 
% figure
% plot(par.s,par.theta*180/pi)


%%



