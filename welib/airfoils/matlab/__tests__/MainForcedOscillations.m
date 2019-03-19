% InitClear
% setFigurePath('/work/publications/phdthesis/figs_airfoil/');
% setMatFigPath('/work/publications/phdthesis/figs_airfoil/matfig/');
% setFigureTitle(0);
addpath(genpath('C:/Config/path/MatlabPath/libs/'))
InitClear
% setFigurePath('/work/publications/phdthesis/figs_airfoil/');
% setMatFigPath('/work/publications/phdthesis/figs_airfoil/matfig/');
setFigureTitle(0);
% setFigurePath('/home/manu/Dropbox/springer-book/figs_airfoil/');
setFigurePath('./');
setMatFigure(0);
setFigureFont('16');

%% Initialization
addpath('../')

%% Parameters
% 
dt=0.01; % time step
valpha_mean=[5 10 15 20 25 30 35]; %[deg]
XLIM=[0 40]
omega=12.57; % frequency of oscillations rad/s
t_max=1.3*(2*pi)/omega % simulation length
% tau=4 * chord / Vrel;
tau=0.08;
% tau=0.8;
bAnimate=false;



%% Loading Polar
Polar=fInitPolar('filename','data/FFA-W3-241-Re12e6.dat')
vt=0:dt:t_max;
nt=length(vt);

vCl=zeros(length(valpha_mean),nt);
if bAnimate
    figure,box on, hold on, grid on
    plot(Polar.alpha,Polar.Cl,'-')
    xlim(XLIM)
end
for ia=1:length(valpha_mean)
    valpha{ia}=valpha_mean(ia)+2*sin(omega*vt);    
    fs_prev=interp1(Polar.alpha,Polar.f_st,valpha_mean(ia)); % init with steady value
    for it=1:nt
        [vCl(ia,it), fs_prev]=fDynaStallOye(valpha{ia}(it),Polar,tau,fs_prev,dt);
        if bAnimate
            plot(valpha{ia}(it),vCl(ia,it),'k.')
            pause(0.01)
        end
    end


end
close(1)

%%
I=1:nt;
figure,box on, hold all, grid on
plot(Polar.alpha,Polar.Cl,'-','LineWidth',2)
for ia=1:length(valpha_mean)
    hvis='off';
    if ia==1; hvis='on'; end
    plot(valpha{ia}(I),vCl(ia,I),'k','handlevisibility',hvis)
end
% plot(Polar.alpha,Polar.Cl_inv,'-')
% plot(Polar.alpha,Polar.Cl_fs,'-')

xlim(XLIM)
ylim([0 2.2])
xlabel('\alpha [^o]')
ylabel('Lift coefficient C_l [-]')
title('DynaStallOyeForcedOscillations')
legend('Static airfoil data','Dynamic airfoil data','Location','South')
% legen('Static airfoil data','Dynamic airfoil data','C_l inviscid','C_l fully-separated','Location','South')
% 
