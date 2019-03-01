function [Expansion vz_bar]=fExpansionFranksen(CT,alpha,k,varargin)
% the variable argument provides vz_bar (distance behind the rotor in terms of rotor radius);
% Remeber dimensionsless with radius not diameters!!!!!!

% default values are alpha =0.7 and k=2

if(nargin==1)
    vz_bar=linspace(0,20,100);
else
    vz_bar=varargin{1};
end
if(vz_bar(1)<0)
    error('that should not be');
end

a=(1-sqrt(1-CT));  % check with article!!!!!! they use a weird a!!!!
beta=(1-a/2)/(1-a);

Expansion = (beta.^(k/2)+1/2*alpha*vz_bar).^(1/k);

if vz_bar(1)==0
    Expansion(1)=1;
end

