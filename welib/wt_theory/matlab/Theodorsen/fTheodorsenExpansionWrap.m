function [ Expansion vz_bar R_over_RW_Wald] = fTheodorsenExpansionWrap(lambda,nB,U0,R0,a,CT,varargin )
% User friendly wrapper to compute expansion factor according to Theodorsen Theory
% TODO: put the iterative method in there, with a switch to chose which method to use
%       The iterative method is required since cs depend on R_w and the far wake parameters are in principle not known
        % TODO : use TheodorsenFarWakeParameters


%--- OUTPUT
% Expansion : R/R_w 
% vz_bar= z/R : downstream dimension
% R_over_Rw_Wald : expansion according to Wald's formula
%--- INPUT
% lambda: tip-speed ratio
% nB: number of blades
% U0: free stream
% a: axial induction
% CT: Thrust coefficient (at the rotor)

tic();
fprintf('Computing Theodorsen''s expansion...');
if isempty(varargin)
    % Default arguments
    Algo.nr=100;
    Algo.ntheta=3000;
    Algo.z_inf=100;
    Algo.bWT=1;
else
    Algo=varargin{1};
end
    
    %% Non iterative method
    % these formulae are propably according to one of Okulov paper of 2007
    l_bar=(1-a)/lambda;
    w_bar=2*a;

    vx=linspace(0,1,Algo.nr);  % Note: there is a singularity at x=1, see the denominator of y1 for theta=tau=0. The denaminator cancels out. Normal, vortex singularity 
    [ K kappa eps_over_k epsz cs cp e eta a0 ] = fTheodorsenAtOperatingPoint( l_bar, nB, vx ,w_bar );
    [ vz_bar ExpFactor Y Y1 Y2 kappa] = fTheodorsenExpansion(l_bar,vx,K,nB,Algo.z_inf,Algo.ntheta,Algo.bWT );
    Rw=ExpFactor(end)*CT*R0+R0;

    Expansion=1+ExpFactor*CT;
fprintf('Done. \t'); toc();


%% Wald 4.4.5 (p116)
xh=0; % hub raidus dimensionless, divided by Rotor radius R0 (not the R0 for Wald, which is Rw)
I1=trapz(vx,2*K.*vx.^3./(l_bar^2+vx.^2)*1./(vx.^2+xh^2/1-xh^2) );

R_over_RW_Wald=sqrt(  ( (1+w_bar*(1/2+eps_over_k))/(1+w_bar)  ) / ( (1-xh^2)*1/2*w_bar*(1+w_bar)/lambda^2*I1/kappa  )  ); % TODO review the formula



end

