InitClear
%%
user.MainPath = '/work/data/WT/NTK500p/';
user.HtcFileName = 'HtcFindPitch.htc';
%  user.MainPath = '/work/data/WT/NTK500/';
%  user.HtcFileName = 'NTK500.htc';
user.BladeMainBodyName = 'blade1';

par.eta=0.93;
par.Pr=500000.0;
par.Pref = par.Pr/par.eta;
par.Pref=500000;


par.r = 1.5; % Hub radius

% AERO Parameters
par.Omega = 27.1/60*2*pi;
par.rho = 1.225;
par.Uvec = 4:2:24;


% BEM parameters
par.k = [-0.0017077,0.251163,0.0544955,0.0892074]; % glauert correction
par.relax = 0.85;  % relaxation of solution of induction
par.Induction = 1;
par.TipLoss = 1;
par.res=1E-3; % parameter I added

% load hawc2 data
Data = ReadHtcFileAll(user);

[par res Cp P]=fPower(user,par,Data);
par.P=P;

%%
figure
plot(par.Uvec,par.P/1000)
grid on
figure
plot(par.Uvec,res.beta*180/pi)
grid on
% hold on
% plot(par.Uvec,par2.P/1000)


% par2=par;
%save('NTK500_stall2.mat','par','res')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% call main function
% [par res]=fFindPitch(user,par,Data);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%
% 
% load('data/Workspace.mat')
% close all
% P=res.all(2,:);
% CP=res.all(5,:);
% U=res.all(1,:);
% beta=res.all(3,:);
% dpdbeta=res.all(4,:);
% 
% % compute kk factor based on frozen wake
% I = find(beta > 0);
% temp = polyfit(beta(I(1):end),dpdbeta(I(1):end),1);
% KK = temp(2)/temp(1);
% 
% figure(1)
% hold all
% plot(beta(I(1):end),dpdbeta(I(1):end),'k+-','LineWidth',2)
% plot(beta(I(1)-1:end),beta(I(1)-1:end)*temp(1)+temp(2),'--','LineWidth',2,'Color',[0.5 0.5 0.5])
% text('Interpreter','latex','String',sprintf('$$\\frac{dP}{d\\theta} = -$$ %.2f [kW / deg]',abs(dPdpitch_0)), 'Position',[2.5 -10], 'FontSize',14)
% 
% 
% 
% grid on
% box on
% xlabel('Pitch \theta [deg]')
% ylabel('dP/d\theta (frozen wake) [kW/deg]')
% title('dpdpitch')
% %%
% figure(2)
% plot(U,beta,'k','LineWidth',2)
% grid on
% box on
% xlabel('Wind speed U_\infty [m/s]')
% ylabel('Pitch \theta [deg]')
% title('pitchvsWS')



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



