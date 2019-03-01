function [Expansion vz_bar ]=fExpansionVortexRings(CT,varargin)
% Uses the Solution from vortex cylinder theory to assess the expansion.
% It's called vortex rings because it turns out that if you use vortex rings with a radius changing according to the vortex cylinder centerline velocity, then you get the right constant induction 1/3 at the rotor plane. The intensity of the vortex rings has also to be changed accordingly (though I need to double check).

% -- INPUT/OUTPUT stuff
% as usual, Expansion and vz_bar are as function of Radius.
% It is possible to specify the downstream coordinate as a varargin, e.g. linspace(0,5,100); (meaning 0 to 5 Radius donwstream)
% Otherwise a weird position is used with more resolution near the "rotor"


if(nargin==1)
    x=[0 exp((log(100)/500).*(1:500))-1]; x=x./x(end); x=x.^2; x=x.*100;
else
    x=varargin{1};
end

R=1;
a=1/2*(1-sqrt(1-CT));

% xm=0.5.*(x(1:(end-1))+x(2:end)); % mid points...

% just lazyness
xm=x;
adumR=a;

Expansion=R.*sqrt((1-adumR)./(1-adumR.*(1+xm./sqrt(1+xm.^2)))); %Radius corr to axial vel in center of vortex cylinder

vz_bar=xm;  % I like to call it vz_bar

