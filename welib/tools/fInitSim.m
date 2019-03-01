function [Sim] = fInitSim(varargin)

% Call it with:
% - 0 arguments:             (then a default empty Sim is returned, a bit pointless..)
% - 1 arguments: WT          (then WT.Spec.vSIM is used if it exists)
% - 1 arguments: sWTName     (then the WT name is given to Sim.WT and a default empty Sim is set)
% - 2 arguments: WT, vSIM    (then vSIM is used for the simulation)
%              See Documentation of fSetSim for the format of vSIM. 
%              In short: if vSIM is a matrix(nx3 or nx4), then it's considered as a combined case. [WS RPM Pitch] or [WS RPM Pitch Yaw]
%              I'm not sure parametric cases can be set with fInitSim for now. For that, use fSetSim



% A simulation is made of several runs defined by a parametric study or a combined case
% Sim.Run is the current run

% It has been extended(more hacked) for non WT sim

bWT=0;
%% Simulation
if(nargin==0)
    Sim.WT='None'
    Sim.Name='None';
    Sim.rho=1.225;
    Sim.KinVisc=15.68*10^-6;
elseif(isstr(varargin{1})) % user just provided a geometry name
    Sim.WT=varargin{1};
    Sim.Name=varargin{1};
    Sim.rho=1.225;
    Sim.KinVisc=15.68*10^-6;
else % the user hopefully provided a WT
    WT=varargin{1};
    bWT=1;% we have a WT
    Sim.WT=WT.Name;
    Sim.Name=WT.Name;
    Sim.rho=WT.DefaultEnvironment.rho;
    Sim.KinVisc=WT.DefaultEnvironment.KinVisc;
end
Sim.bKiloOut=1;
Sim.CombinedCases=[];
Sim.ParametricPITCH=[];
Sim.ParametricRPM=[];
Sim.ParametricWS=[];
Sim.ParametricYAW=[];

if(bWT)
    if(nargin==2 || size(WT.Spec.vSIM,1)~=0)
        % For compatibility with previous use of fInitSim which is followed by a fInitWind
        % We set the sim if it was provided
        % and we set the run to the first one in the simulation
        if nargin==2
            Sim=fSetSim(Sim,WT,varargin{2});
        else
            Sim=fSetSim(Sim,WT,WT.Spec.vSIM);
        end

    else
        fprintf('Remember this is just init, not set. Use Set instead\n');
    end
end
