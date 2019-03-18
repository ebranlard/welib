
InitClear
% setFigurePath('/work/publications/phdthesis/figs_airfoil/');
% setMatFigPath('/work/publications/phdthesis/figs_airfoil/matfig/');
setFigureTitle(0);
setFigurePath('./');
setMatFigure(0);
setFigureFont('16');



%% Initialization
% require('PROFILES','v01')
addpath('../')

%% Parameters
% inputs for fully separated polar (in degrees)
alpha_merge=35;
delta_alpha=10;
dclda_fs=pi*(pi/180); % deg



%% 

% Loading profile data
% Data=load(['../data/tjaere11_ds.dat']);
% alpha=Data(:,1);
% Cl=Data(:,2);
% Cd=Data(:,3);
% Cl_inv=Data(:,5);
% Cl_fs=Data(:,6);
% f_st=Data(:,4);

Data=load(['data/FFA-W3-241-Re12e6.dat']);
alpha=Data(:,1);
Cl=Data(:,2);
Cd=Data(:,3);
Cm=Data(:,4);

%%
[Cl_inv Cl_inv_sin alpha0]      = fPolarInviscid(alpha,Cl)                                        ; 
[Cl_fs,alpha1,alpha2,alr1,alr2] = fPolarFullySeparated(alpha,Cl,dclda_fs,alpha_merge,delta_alpha) ; 
f_st=(Cl-Cl_fs)./(Cl_inv-Cl_fs);



%%
figure, hold all, box on
title('AirfoilFFAW3Clfullyseparated')
plot(alpha,Cl,'LineWidth',2)
plot(alpha,Cl_inv)
plot(alpha,Cl_fs)
vline(alpha0,'Color','k')
vline(alpha_merge,'Color','k')
vline(alpha1-alr1,'LineStyle','--','Color',fGray(0.5))
vline(alpha1+alr1,'LineStyle','--','Color',fGray(0.5))
% vline(alpha2-alr2)
% vline(alpha2+alr2)
xlabel('\alpha [^o]')
ylabel('C_l [-]')
% grid on
xlim([-50 50])
ylim([-2 4])
legend('C_l steady','C_l inviscid', 'C_l fully-separated',2)


%% 
% figure,hold all 
% plot(alpha,f_st)
% vline(alpha_fs)


