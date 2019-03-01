function [ w_bar l_bar Rw Expansion vz_bar ] = fTheodorsenFarWakeParams(CT,lambda,nB,varargin)
% -- NOTES
% Find farwake parameters that correspond to a set of rotor parameter
% See also fGoldsteinFarWakeParams
% As far as I remember they give quite different results
% Most of the stuff commented below are from the goldstein equivalent function

% -- Usage
% See in CirculationForCanonicalWTSimulation....
% Should be as well in CompareWakeShape..

% -- INPUT/OUTPUT
% Remember vz_bar and ExpFactor are given relative to Radius!
% varargin is for Algo used in Expansion factor


% if nB==3
%     load('OptimalCTGoldsteinB3.mat');% vw vCP vCT -> values obtained by iteration such that w_bar = 2/3
% else
%     error('dsf');
% end
bGoOptimal=0;
% if(~isempty(varargin))
%     bGoOptimal=1;
% end
% CT_opt=interp1(vLambda,vCT,lambda);
% if CT > CT_opt
%     warning(sprintf('there is no hope to reach this CT ofr this lambda, the max CT is %.2f',CT_opt));
%     bGoOptimal=1;
% end
% 
% if bGoOptimal
%     w_bar=interp1(vLambda,[1 vw],lambda);
%     CP_out=interp1(vLambda,[0 vCP],lambda);
%     CT_out=CT_opt;
%     l_bar=(1-w_bar/2)/lambda;
%     vr=linspace(0,1,100);  
%     G=fGoldsteinFactor( l_bar,nB,vr);
% else
if isempty(varargin)
    Algo.nr=100;
    Algo.ntheta=3000;
    Algo.z_inf=100;
    Algo.bWT=1;
else
    Algo=varargin{1};
end

% initial guess
R=1;
Rw=1.01;


w_bar0=(1-sqrt(1-CT))*0.6;  %1D momentum result used as a first guess with a 0.6 factor as rule of thumb..
cs_new=CT*R^2/Rw^2;
cs_old=0;
tic(); 
fprintf('Iterating till converging on right expansion...\n');
cpt=0;
while abs(cs_old-cs_new)/max(cs_old,cs_new) >5*10^-3
    opts=optimset('TolX',1e-3);
%     opts=optimset('TolX',1e-3,'Display','iter');
    tic(); 
    fprintf('Iterating to find Theodorsen far wake parameters...');
    % find the w_bar which gives the right cs
    % the if below is useless if the optimization help is commented
    if(cpt==0)
        [w_bar,~,~,Output]=fzero(@(w) cs_new-fcs(w,lambda,nB,R,Rw),w_bar0,opts);
    else
        %let's try to help the optimzation
%          s=Output.message;
%          keyboard
% %          [~,xstart,xend]=strread(s,'%s[%f,%f]');
%         [xstart,xend]=strread(s,'Zero found in the interval [%f, %f]');
%          if(Rw_old<Rw)
%              %the lower bound will be smaller
% %              xstart=xstart-0.05*xstart*Rw_old/Rw;
%              xstart=xstart*Rw_old/Rw*0.99;
%          else
%              %the upper bound should probably be increase
%             xend=xend*Rw/Rw_old*1.01;
%         end
%         disp([xstart xend]);
        [w_bar,~,~,Output]=fzero(@(w) cs_new-fcs(w,lambda,nB,R,Rw),w_bar0,opts);
    end
    l_bar=1/lambda*(1-w_bar)*R/Rw;
    fprintf('Done,w_bar: %.3f  - l_bar %.3f \t',w_bar,l_bar),toc();
    % now find the expansion....
    vx=linspace(0,1,Algo.nr); 
    [ K kappa eps_over_k epsz cs cp ] = fTheodorsenAtOperatingPoint( l_bar, nB, vx ,w_bar );
    [ vz_bar ExpFactor ] = fTheodorsenExpansion(l_bar,vx,K,nB,Algo.z_inf,Algo.ntheta,Algo.bWT );
    Rw_old=Rw;
    Rw=ExpFactor(end)*cs*R+R;
    cs_new=CT*R^2/Rw^2;
    cs_old=cs;
    w_bar0=w_bar;
    
    Expansion=1+ExpFactor*CT;

    fprintf('CT design %.2f - cs found %.3f - following cs %.3f - Rw found %.3f\n',CT,cs_old,cs_new,Rw); 
    cpt=cpt+1;
end
fprintf('Done after %d iterations:\t',cpt),toc();



function [cs  l_bar G]= fcs(w_bar,lambda,nB,R,Rw)
vx=linspace(0,1,100);  
l_bar=1/lambda*(1-w_bar)*R/Rw;
[ G kappa eps_over_k epsz cs ] = fTheodorsenAtOperatingPoint( l_bar, nB, vx ,w_bar );


