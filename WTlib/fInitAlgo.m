function [ Algo ] = fInitAlgo( varargin  )

Algo.CodeVersion=6;
Algo.bSteady=logical(1);
Algo.bWT=logical(1); % Wind turbine or propeller

%% Convergence stuff
Algo.nbIt=200;
Algo.aTol=10^-6;
Algo.swTol=10^-6;
Algo.relaxation=0.5;



%% Aero stuff - some are obscure - Typicaly used by AeroCoeffWrap
Algo.bAIDrag=logical(1); % !!!
Algo.bTIDrag=logical(1); % !!!
Algo.bUseCm=logical(0); % !!!
Algo.bSwirl=logical(1); 
Algo.bDynaStall=logical(0);
Algo.bReInterp=logical(1);
Algo.bRoughProfiles=logical(0);
Algo.bCl2piAlpha=logical(0);
Algo.bCl2piAlphaCorr=logical(0);
Algo.bNoDrag=logical(logical(0));
Algo.ClOverCd=100000;
Algo.bSkewWake=logical(0);
Algo.bAdvanceBrakeState=logical(1);
Algo.bIndProp=logical(0);
Algo.bAlphaHacked=logical(0);

Algo.bProfilesData=logical(1);  % important parameter for VL code
Algo.bPrescribedGamma=logical(0);
Algo.r_PrescribedGamma=[];
Algo.PrescribedGamma=[];


%% resolution
Algo.n_span=15;    % number of PANELS spanwise
Algo.n_chord=1;    % number of PANELS chordwise
Algo.nRev=4;
Algo.nPhi=40;
Algo.nRev_farwake=5;
Algo.nRev_farwake_cut=20;
Algo.nRev_bigwake=2.5;

Algo.sChordMeshFunction='Cos';
Algo.sSpanMeshFunction='Cosb'; %Cosb

%% plot switches 
Algo.bPlotGeometry=0; % plot geometry
Algo.bPlots=0; % plot geometry
Algo.bMovie=0;
Algo.bSimpleMovie=0;
Algo.bComplexMovie=0;
Algo.bExport=1;
Algo.ntotMarkers=0;


%% output and storage
Algo.vr_bar_out=linspace(1.02,1.5,5);
Algo.vpsi_out=linspace(0,2*pi,90);
Algo.vz_out=[0];
Algo.bInspectResults=0;
Algo.bUnsteadyStorage=logical(0);

%% General stuff
Algo.bParallel=logical(1);
Algo.bPolarIn=logical(0);
Algo.bPolarOut=logical(0);
Algo.bGridIn=logical(0); % is the input three vectors from which the grid should be formed
Algo.bFlatOut=logical(0);
Algo.bAllRotor=logical(0);

Algo.bBoundInTipLoss=0;
Algo.bTipLossPlaneAvg=1;
Algo.bTipLossInfiniteB=0;
Algo.bWinglet=0;
Algo.WingletLength=0.02;
Algo.WingletAngle=90;

Algo.bWeight=logical(0);
Algo.bSteady=logical(1);
Algo.DOF=logical(0);




%% BEM 
if exist('fInitBEMAlgo','file')
    Algo.BEM=fInitBEMAlgo();
end


%% Vortex Cylinder
if exist('fInitVCYLAlgo','file')
    Algo.VCYL=fInitVCYLAlgo();
end


%% Vortex Filament Code
if exist('fInitVCFILAlgo','file')
    Algo.VCFIL=fInitVCFILAlgo();
end


%% Vortex Lattice code
if exist('fInitVLAlgo','file')
    Algo.VL=fInitVLAlgo();
end


%% Helix stuff
Algo.Helix.Method='line_c';
Algo.Helix.bWT=logical(1);
Algo.Helix.bInf=logical(0); %infinite helix solutions?

%% Ring Stuff
Algo.Ring.ntot=100;




Algo.WakeType='none';  %Helix, Prescribed, Free
Algo.WakeExpansionMethod='none';

%% Changing some default parameters depending on the code
if nargin>0
    sCode=varargin{1};
	switch(sCode)
		case 'VFil'
            Algo.bPolarIn=logical(1);
            Algo.bPolarOut=logical(1);
            Algo.bGridIn=logical(1);
            Algo.Viscous.Model=0;  % 0: no model ; 1: Rankine ; 2:Lamb-Olsen ; 3:Vatistas n=2 ; 4:Cut-off
            Algo.Viscous.delta= 10; %   viscous   core   diffusitivity   coeff.
            Algo.Viscous.t0   = 0.01;% time offset in the viscous core growth model
            Algo.Viscous.bDeformation=logical(0); %core deformation
            Algo.Viscous.bCoreGrowth=logical(0);
            Algo.Viscous.t=zeros(1,10);
		case 'VCYL'
            Algo.bPolarIn=logical(1);
            Algo.bPolarOut=logical(1);
            Algo.bGridIn=logical(1);
		case 'VC'
			Algo.bSteady=logical(0);
			Algo.Viscous.Model=3;  % 0: no model ; 1: Rankine ; 2:Lamb-Olsen ; 3:Vatistas n=2 ; 4:Cut-off
		case 'BEM'
            %
        otherwise
            error('Code unknown')

	end
end


Algo.Viscous.Model=3;  % 0: no model ; 1: Rankine ; 2:Lamb-Olsen ; 3:Vatistas n=2 ; 4:Cut-off
Algo.Viscous.delta  = 900; %   viscous   core   diffusitivity   coeff.
Algo.Viscous.t0      = 20;% time offset in the viscous core growth model
Algo.Viscous.bDeformation=0;
Algo.Viscous.Epsilon=0;
Algo.Viscous.alpha = 1.25643;
Algo.Viscous.t      = 0;% 


%% Specific stuff for WT
if nargin==2
    sGeometry=varargin{2};
    % Geometry
    if(isequal(sGeometry,'SB2'))
        Algo.Viscous.delta  = 500; %   viscous   core   diffusitivity   coeff.
        Algo.Viscous.t0      = 10;% time offset in the viscous core growth model
    elseif(isequal(sGeometry,'NTK500') || isequal(sGeometry,'NTK500p'))
        Algo.Viscous.delta  = 900; %   viscous   core   diffusitivity   coeff.
        Algo.Viscous.t0      = 20;% time offset in the viscous core growth model
    elseif(isequal(sGeometry,'NREL5MW'))
        Algo.Viscous.delta  = 900; %   viscous   core   diffusitivity   coeff.
        Algo.Viscous.t0      = 20;% time offset in the viscous core growth model
    elseif(isequal(sGeometry,'Riso10MW'))
        Algo.Viscous.delta  = 900; %   viscous   core   diffusitivity   coeff.
        Algo.Viscous.t0      = 20;% time offset in the viscous core growth model
    elseif(isequal(sGeometry,'SB1'))
        Algo.Viscous.delta  = 900; %   viscous   core   diffusitivity   coeff.
        Algo.Viscous.t0      = 20;% time offset in the viscous core growth model
    elseif(isequal(sGeometry,'SB1_swept'))
        Algo.Viscous.delta  = 500; %   viscous   core   diffusitivity   coeff.
        Algo.Viscous.t0      = 10;% time offset in the viscous core growth model
    elseif(isequal(sGeometry,'FlatTUD'))
        %         load('Validation/data/SantCirc1.mat');
        %         r_GammaPrescribed=SantCirc1(:,1)*R;
        %         GammaPrescribed=SantCirc1(:,2);
    elseif(isequal(sGeometry,'elliptic'))
        Algo.sChordMeshFunction='Linear';
        Algo.sSpanMeshFunction='Cos';
        Algo.c0=1;
        Algo.alpha=5.7106;
        Algo.B=1;
    elseif(isequal(sGeometry,'flatplate'))
        Algo.sChordMeshFunction='Linear';
        Algo.sSpanMeshFunction='Cos';
        Algo.opts.alpha=5;
        Algo.opts.c0=1;
    elseif(isequal(sGeometry,'trapezoid'))
        Algo.sChordMeshFunction='Linear';
        Algo.sSpanMeshFunction='Linear';
        Algo.opts.alpha=10;
        Algo.opts.lambda=45;
        Algo.opts.b=5;
        Algo.opts.c0=0.2*Algo.opts.b;
    end


end

end
