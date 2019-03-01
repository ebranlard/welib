function [ w_bar l_bar CT_out G] = fGoldsteinFarWakeParams(CT,lambda,nB,varargin)
if nB==3
    load('OptimalCTGoldsteinB3.mat');% vw vCP vCT -> values obtained by iteration such that w_bar = 2/3
else
    error('dsf');
end
bGoOptimal=0;
n=100;
if(isempty(varargin))
    bGoOptimal=0;
elseif nargin==4
    bGoOptimal=varargin{1};
elseif nargin==5
    bGoOptimal=varargin{1};
    n=varargin{2};
end
if(n>3000)
    error('Mex file only for n<3000')
end


CT_opt=interp1(vLambda,vCT,lambda);
if CT > CT_opt
    warning(sprintf('there is no hope to reach this CT ofr this lambda, the max CT is %.2f',CT_opt));
    bGoOptimal=1;
end

if bGoOptimal
    w_bar=interp1(vLambda,[1 vw],lambda);
    CP_out=interp1(vLambda,[0 vCP],lambda);
    CT_out=CT_opt;
    l_bar=(1-w_bar/2)/lambda;
    vr=linspace(0,1,100);  
    G=fGoldsteinFactor( l_bar,nB,vr);
else
    w_bar0=1-sqrt(1-CT);  %1D momentum result used as a first guess
    opts=optimset('TolX',1e-3);
    tic(); 
    fprintf('Iterating to find Goldstein far wake parameters ...');
%     w_bar=fzero(@(w) CT-fCT(w,lambda,nB),[0 1.5],opts);
    w_bar=fzero(@(w) CT-fCT(w,lambda,nB,n),w_bar0,opts);
    fprintf('Done,\t'),toc();
    [CT_out CP_out l_bar G]=fCT(w_bar, lambda,nB,n);
end




function [CT CP l_bar G]= fCT(w_bar,lambda,nB,n)
vr=linspace(0,1,n);  
l_bar=(1-w_bar/2)/lambda;  % that's the link given by Okulov and Sorensen.
G=fGoldsteinFactor( l_bar,nB,vr);
I1=2*trapz(vr,G.*vr);
%I2=2*trapz(vr,G*l_bar.*vr.^2./(l_bar^2+vr.^2));
I3=2*trapz(vr,G.*vr.^3./(l_bar^2+vr.^2));
CP=2*w_bar*(1-w_bar/2)*(I1-w_bar/2*I3);
CT=2*w_bar*(I1-w_bar/2*I3);
