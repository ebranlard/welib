function [Parametric,CPmax,lambda_opt,pitch_opt]= fWTParametricCP(WT,vLambda,vPitch,bPlot,Algo)
% Computes a parametric study on vLambda and vPitch
% Returns the CP curve and the optimal parameters

Code='BEM';

[ Sim ]  = fInitSim(WT);
[ Wind ] = fInitWind(  );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Cp curve for lambda and pitch - optimal pitch
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Computing Cp-lambda-pitch curve...\n')

% Setting up a parametric study for various Pitch and lambda (ie RPM with WS constant)
WS=10;
vRPM     = vLambda*WS/WT.Rotor.R *60/(2*pi);
[ Sim ]  = fSetSim( Sim, WT, WS,vRPM, vPitch );
%
[ RCPlp ]= fWTSimulation(Code,WT,Sim,Wind,Algo); % Results CP Lambda Pitch RCPlp
% CP matrix of dimension nPitch x nRPM
RCPlp.CP;

RCPlp.CP(RCPlp.CP<0)=NaN;
% optimal Cp
[MC,irows]          = max(RCPlp.CP);
[CPmax,icol_lambda] = max(MC);
irow_pitch=irows(icol_lambda);
pitch_opt =RCPlp.PITCH (irow_pitch,icol_lambda); %PITCH  dimension nPitch x nRPM
lambda_opt=RCPlp.lambda(irow_pitch,icol_lambda);%lambda dimension nPitch x nRPM

fprintf('Optimal parameters : Pitch =  %.2f - Lambda = %.2f - Cp = %.2f\n',pitch_opt,lambda_opt,CPmax);


%% Returning data
Parametric.vLambda = vLambda;
Parametric.vPitch  = vPitch;
Parametric.Data    = RCPlp     ;

%% Plotting
if(bPlot)
    %colrs=fColrs(1:4);
    %sty={'-','+-','--'};
    figure, hold on, grid on, box on
    plot(vLambda,RCPlp.CP(irow_pitch,:));
    xlabel('\lambda [.]')
    ylabel('C_P [.]')
   
   
    figure, hold on,grid on, box on
    surf(RCPlp.lambda,RCPlp.PITCH,RCPlp.CP)
    plot3(lambda_opt,pitch_opt,CPmax,'ko','MarkerSize',6,'MarkerFaceColor','k');
    xlabel('\lambda [.]')
    ylabel('Pitch [deg]')
    zlabel('C_P [.]')
    view(3)
    axis ij
end

