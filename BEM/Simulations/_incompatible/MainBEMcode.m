%% Initialization for the BEM code
InitClear;
InitDefault;
path(path,'./f_format')
path(path,'./f_aeroelastic')
path(path,'./f_optimal')


%
%Opts.AeCols=4;
%fReadSpec('data/NORDTANK_Spec.dat'); Algo.Ngrid=30;
%Files={'data/test/NTK500_pitch.htc'};
%fInit('hawc',Files,Opts)


% Opts.AeCols=4;
% Files={'data/NORDTANK/NTK500.htc','data/NORDTANK_Spec.dat'};
% fInit('hawc',Files,Opts)

% 
% Files={'data/Tjaere/Tjaere_BladeGeometry.dat','data/Tjaere/Tjaere_BladeProfiles.dat','data/Tjaere/Tjaere_Spec.dat'};
% fInit('flex',Files,Opts)
% 
% Files={'data/NORDTANK_BladeGeometry.dat','data/NORDTANK_BladeProfiles.dat','data/NORDTANK_Spec.dat'};
% fInit('flex',Files,Opts)

%  
%Files={'data/WTperf_aerodyn/Test03_CART3.wtp'}; OutputFile='data/WTperf_aerodyn/Test03_CART3.oup';
Files={'data/WTperf_aerodyn/Test022_AWT27.wtp'};OutputFile='data/WTperf_aerodyn/Test022_AWT27.oup';
fInit('wtperf',Files,Opts)
% 
% Files={'data/export_wt_perf/Test022_AWT27.wtp'};OutputFile='data/WTperf_aerodyn/Test022_AWT27.oup';
% fInit('wtperf',Files,Opts)

% Files={'data/export_wt_perf/Test03_CART3.wtp'}; OutputFile='data/WTperf_aerodyn/Test03_CART3.oup';
% fInit('wtperf',Files,Opts)



% Opts.AeCols=6;
% Files={'data/b45/htc','data/NORDTANK_Spec.dat'};
%  P=fReadAeFile('data/b45/BH_Aerodyn_Ae.b45',2,7);
%  P=fReadPcFile('data/b45/BH_Aerodyn_PC.b45',2);
% P=fReadPcFile('data/b45/BH_Aerodyn_PC.b45',2);
% P=fReadBtcFile('data/b45/SWT-23-93_F03_v29_1.btc');
% Files={'data/b45/SWT-23-93_F03_v29_1.btc','data/b45/B45_Spec.dat'};
% fInit('hawc',Files,Opts)
% Opts.AeCols=4;


%A=fReadXbladeParamFile('C:/work/data/Aero/B49/B49_basecase.txt')
% PCFileExtendedRead('C:/work/data/Aero/B49/B49_extendedPCFile_NACA_18_21_FFA_24_30_B49_34-59.txt');
% A=fReadPcFileExtended('C:/work/data/Aero/B49/B49_extendedPCFile_NACA_18_21_FFA_24_30_B49_34-59.txt');
%% calling the BEM code 
Algo.Steady=1;
Algo.NumSect=1;
Algo.nbIt=500;
%Algo.alphaCrit=0.001;
%Algo.relaxation=0.8;
%Algo.correction='Spera';
Algo.YawModel=0;
Simulation.Parametric.Omega(1:2)=Rotor.Omega;
Simulation.Parametric.Pitch=[-2 -2 0];
Simulation.Parametric.WS=[5 5 1];
%%

BEMSimulation


%%
fid=fopen(OutputFile);
A=cell2mat(textscan(fid,'%c %*[^\n]\n',8));
P=cell2mat(textscan(fid,'%f %f %f %f %*[^\n]\n'));
A=cell2mat(textscan(fid,'%c %*[^\n]\n',4));
CP=cell2mat(textscan(fid,'%f %f %f %f %*[^\n]\n'));
A=cell2mat(textscan(fid,'%c %*[^\n]\n',4));
Q=cell2mat(textscan(fid,'%f %f %f %f %*[^\n]\n'));
A=cell2mat(textscan(fid,'%c %*[^\n]\n',4));
Flp=cell2mat(textscan(fid,'%f %f %f %f %*[^\n]\n'));
A=cell2mat(textscan(fid,'%c %*[^\n]\n',4));
Thr=cell2mat(textscan(fid,'%f %f %f %f %*[^\n]\n'));
fclose(fid);

%%

i=1;j=1;

%%
figure
hold all
for(i=1:length(Vpitch))
    plot(Vws,squeeze(Thrust(i,j,:)))
end
plot(Thr(:,1),Thr(:,2:4),'+')
xlabel('WS [m/s]')
ylabel('Thrust [kN]')
title('1thrust')
grid on
%%

figure(2)
hold all
for(i=1:length(Vpitch))
    plot(Vws,squeeze(Flap(i,j,:))) 
end
%plot(Vws,squeeze(Edge(i,j,:)),'g')
plot(Flp(:,1),Flp(:,2:4),'+')
xlabel('WS [m/s]')
ylabel('Flap [kNn]')
title('1flap')
grid on

%%
i=1;j=1;
figure(3)
hold all
for(i=1:length(Vpitch))
    plot(Vws,squeeze(Power(i,j,:)),'o')
end
plot(P(:,1),P(:,2:4),'+')
%xlim([6 8])
xlabel('WS [m/s]')
ylabel('Power [kW]')
title('1power')
grid on



%%


