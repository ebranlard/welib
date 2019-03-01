function [ Sim ] = fSetSim(Sim,WT,varargin )

% Called with three arguments, it's considered as a combined case, ie a matrix of [WS RPM Pitch] or [WS RPM Pitch Yaw]
% Sim,WT,[WS RPM Pitch (YAW)]

% Called with five or six arguments it's considered as a parametric cases:
% Sim,WT, vWS,vRPM,vPitch,(vYaw)

% Now, WT can be empty, and then lambda is set to 0;


if(nargin==3)
    Sim.ParametricPITCH=[];
    Sim.ParametricRPM=[];
    Sim.ParametricWS=[];
    Sim.ParametricYAW=[];
    %considered as a Combined Case
    WS=varargin{1}(:,1);
    RPM=varargin{1}(:,2);
    PITCH=varargin{1}(:,3);
    if(size(varargin{1},2)==3)
        %yaw has been omitted, set it to zero
        YAW=PITCH*0;
    else
        YAW=varargin{1}(:,4);
    end
    if( length(WS)~=length(PITCH) || length(WS)~=length(RPM) || length(WS)~=length(YAW) )
        error('wrong combined case')
    end
    Sim.CombinedCases=[ WS RPM PITCH YAW ];
    Sim.Run=fSetRun(Sim.CombinedCases(1,:));
else
    Sim.CombinedCases=[];
    %considered as a parametric case
    WS=varargin{1};
    RPM=varargin{2};
    PITCH=varargin{3};
    if(nargin==5)
        YAW=0;
    else
        YAW=varargin{4};
    end
    Sim.ParametricPITCH=PITCH;
    Sim.ParametricRPM=RPM;
    Sim.ParametricWS=WS;
    Sim.ParametricYAW=YAW;

    Sim.Run=fSetRun(Sim.ParametricWS(1),Sim.ParametricRPM(1),Sim.ParametricPITCH(1),Sim.ParametricYAW(1));
end
Sim.Run.Omega=Sim.Run.RPM*2*pi/60;

if isempty(WT)
    Sim.Run.lambda=0;
else
    Sim.Run.lambda=Sim.Run.Omega*WT.Rotor.R/Sim.Run.WS;
end
end

