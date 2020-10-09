function [ Algo ] = fInitVLAlgo_old( varargin )



% AeroCoeffWrap Params
Algo.bAIDrag=logical(1);
Algo.bTIDrag=logical(0);
Algo.bDynaStall=logical(0);

Algo.bReInterp=logical(1);
Algo.bRoughProfiles=logical(0);

Algo.bCl2piAlpha=logical(0);
Algo.bCl2piAlphaCorr=logical(0);
Algo.bNoDrag=logical(0);
Algo.ClOverCd=100000;

Algo.bProfilesData=logical(1);  % important parameter for VL code
Algo.bPrescribedGamma=logical(0);
Algo.r_PrescribedGamma=[];
Algo.PrescribedGamma=[];









%% VL algo switches
%switches
Algo.VL.bRollUp=1;       % rollup of the wake
Algo.VL.bFarWake=1;      % Adding an helical far wake
Algo.VL.bBigWake=0;      % We start with an attached wake
Algo.VL.bCorrector=1;    % using corrector and predictor for computation of propagating velocity
Algo.VL.bAsymmetric=1;    % Asymmetric simulation for shear yaw, turbulence, aeroelasticity....
Algo.VL.bInducedUV=1;
Algo.VL.bInducedW=1;
Algo.VL.bStreching=0;
Algo.VL.bInducedDrag=0;
Algo.VL.bStopWhenConverged=0;
Algo.VL.ConvergedCrit=10^(-6);

% Using Tabulated profiles data
Algo.VL.bFeedSolution=0;
Algo.VL.gamma_crit=0.01;
Algo.VL.gamma_relax=0.22;  % I changed this value it used to be 0.02, probably for flat plates and stuff... it makes convergence way faster at the beginning of WT simulation









