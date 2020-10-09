function [Expansion vz_bar]=fExpansionRathmann(CT,varargin)
% the variable argument provides vz_bar (distance behind the rotor in terms of rotor radius);
% Remeber dimensionsless with radius not diameters!!!!!!

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
alpha=0.7;
Expansion = max(beta^(1/2),(1/2*alpha*vz_bar).^(1/2));

if vz_bar(1)==0
    Expansion(1)=1;
end

