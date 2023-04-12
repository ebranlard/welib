% This script looks only at the wake expansion (it doens' look at tip-losses or BEM etc.)
% There are some nasty copy paste to look at different CT

%% 
InitClear;
PATH.EXPANSION='/work/lib/WTTheory/Expansion/';
require('THEODORSEN','v01');
% require('EXPANSION','');
require('OPTIMCIRC','v00');
require('WTlib','v03');
require('VC_LIB_MAT','v04');
require('BEM','v03');
require('VL','v03');
require('VC_LIB_C','v04');
% require('VC_LIB_C','v01_siemens');


% setFigurePath({'./' , '/work/publications/articles/2012-tiploss-theoretical/figs/'});
setFigurePath({'./' , '/work/publications/phdthesis/figsdump/'});

vCT=[0.35 0.65 0.8];
% vCT=0.85;
% vCT=0.737;
for iCT=1:length(vCT)
    CT=vCT(iCT);
    vz_bar=linspace(0,30,100); %z/R !!! not z/D

    Algo.nr=100;
    Algo.ntheta=length(vz_bar); % dangerous
    Algo.z_inf=max(vz_bar);     % dangerous
    Algo.bWT=1;
%     lambda=15
    lambda=7.06;
    nB=3;
    U0=10;
    R=1;

    [ exp_franksen ] = fExpansionFranksen(CT,0.7,2,vz_bar);
    [ exp_rathmann ] = fExpansionRathmann(CT,vz_bar);
    [ exp_theodorsen ] = fTheodorsenExpansionWrap(lambda,nB,U0,R,1/2*(1-sqrt(1-CT)),CT, Algo); % !!!! mean a, that's not really good
    [ exp_vortexrings ] = fExpansionVortexRings(CT,vz_bar); 


    %%
    figure,hold all,grid on,box on
    plot(vz_bar/2,exp_franksen,'')
    plot(vz_bar/2,exp_rathmann,'')
    plot(vz_bar/2,exp_theodorsen,'k')
    plot(vz_bar/2,exp_vortexrings,'k--')
    k=0.05;m=1;
    plot(vz_bar/2,k*(vz_bar/2).^m+1,'','Color',[0.4 0.4 0.4])
    xlim([0 6])
    xlabel('z/D [.]')
    ylabel('r/R [.]')
    title(sprintf('WakeExpansionDifferentTheoriesCT%2d',CT*100));
    legend('Franksen','Rathman','Theodorsen','Vortex Cylinder','k x^m')


    %% BEM run to then use local value for Mac's prescribed wake
    vx=linspace(0,R,100);
    [w_bar l_bar CTgold G]=fGoldsteinFarWakeParams(CT,lambda,nB);

    [ WT ]   = fInitWT( 'SB1' ,'xblade',PATH.DATA_WT);
    WT.Rotor.cone=0; WT.Nacelle.tilt=0;
    WT.Sources.Rotor.chord=WT.Sources.Rotor.chord/WT.Rotor.R;
    WT.Sources.Rotor.twist=WT.Sources.Rotor.twist*0+20;
    WT.Rotor.R=R; WT.Rotor.rhub=0.00; WT.Rotor.BladeLength=R; WT.Rotor.R_coned=R; WT.Rotor.SweptArea=pi*R^2;
    WT.Sources.Rotor.r=linspace(WT.Rotor.rhub,R,length(WT.Sources.Rotor.r));
    Scale=60/WT.Rotor.R;
    [ WT ] = fSetRotorGrid(40,WT);
    % BEM Algorithm init
    [ Algo ] = fInitAlgo();
    Algo.bPrescribedGamma=1;
    RPM=lambda*U0/WT.Rotor.R*60/2/pi;
    [ Sim ]  = fInitSim( WT , [U0 RPM -10]);
    [ Wind ] = fInitWind( Sim );
    Gamma=G/nB*w_bar*U0*2*pi*R*l_bar;
    Algo.r_PrescribedGamma=vx;
    Algo.PrescribedGamma=Gamma;
    [ BEM ] = fRunBEM(WT,Sim,Wind,Algo);
    fprintf('CT asked %.2f - CT gold %.2f - CT BEM %.2f (Gold circ)\n',CT,CTgold,BEM.CT);

    %% Mac's prescribed wake
    Algo.WakeExpansionMethod='Mac';
    Algo.WakeType='Prescribed';
    zcoord=linspace(0,max(vz_bar),10000);
    Algo.Viscous.t=zeros(1,length(zcoord)+1);
    algo.Viscous.Model=0;
    rhelix=(WT.Rotor.rfull(1:end-1)+WT.Rotor.rfull(2:end))/2;
    [ Wake ] = fGeometryPrescribedWake(lambda,nB,U0,R,rhelix, BEM.r, BEM.Gamma, BEM.CTloc , BEM.a, BEM.aprime , zcoord ,[],[],[], [],[],[], Algo );

    %% Free wake with prescribed Gamma
%     [ Sim ]  = fInitSim( WT , [U0 lambda*U0/R*60/(2*pi) 0] );
%     [ Wind ] = fInitWind( Sim );
    % Main Parameters
    Algo.VL.bCorrector=0; % !!!!
    Algo.bSteady=0;
    Algo.bPrescribedGamma=1;
    Algo.bPlots=0; 
    Algo.bExport=0; 
    Algo.bInspectResults=0;
    % Resolution parameters
    Algo.n_span=6;    % number of PANELS spanwise
    Algo.n_chord=1;    % number of PANELS chordwise
    Algo.nRev=30;
    Algo.nPhi=15;
    Algo.Viscous.Model=4;    %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Algo.Viscous.delta=1;
    Algo.Viscous.t0=0.5;
    Algo.Viscous.delta_wake=0.0001;
    Algo.Viscous.delta_bound=0.001;
    % Running with function
%     [ VL WTout ] = fRunVL(WT,Sim,Wind,Algo );



% 
    %% wake extremities
    % find wake extremities for mac's prescribed wake
    MI=getIntervals(Wake.X(:,end)>1);
    I=[];
    for i=1:size(MI,1)
        [~,ii]=max(Wake.X(MI(i,2):MI(i,3),end));
        I=[I MI(i,2)+ii-1];
    end
    %% Free wake - extremities
%     % find wake extremities for free wake, because of rollup I look at several of them
%     Ztip={};
%     Xtip={};
%     for iend=0:2
%         MI=getIntervals(VL.Wake.X(:,end-iend)>1);
%         Ifree=[];
%         for i=1:size(MI,1)
%             [~,ii]=max(VL.Wake.X(MI(i,2):MI(i,3),end-iend));
%             Ifree=[Ifree MI(i,2)+ii-1];
%         end
%         Ztip{iend+1}=VL.Wake.Z(Ifree,end-iend);
%         Xtip{iend+1}=VL.Wake.X(Ifree,end-iend);
%     end
% 

    %% plotting
    figure(12),hold all,grid on,box on
    plot(vz_bar/2,exp_theodorsen,'-','Color',fColrs(iCT))
    plot(vz_bar/2,exp_vortexrings,'--','Color',fColrs(iCT))
    plot(Wake.Z(I,end)/2,Wake.X(I,end),'+','Color',fColrs(iCT));
%     for iend=1:length(Ztip)
%         plot(Ztip{iend}/2,Xtip{iend},'o','Color',fColrs(iend))
%     end
    xlim([0 4])
    ylim([1 1.3])
    xlabel('z/D [.]')
    legend('Theodorsen','Vortex Cylinder','Prescribed Wake','Location','North')
    title('WakeExpansionTheodorsenCylinderMac')
    ylabel('r/R [.]')
    %%
end
%%
