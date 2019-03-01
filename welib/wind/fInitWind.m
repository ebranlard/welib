function [ Wind ] = fInitWind( varargin )
%% Wind
Wind.nu=0.0; % exponent for the power law
Wind.Model='Constant'; %;
Wind.V0=[0; 0; 0]; % wind speed at hub height
Wind.fV0=@(x)[0; 0; 0] ; % wind speed at hub height

if(nargin==1)
    Wind=fSetWind(Wind,varargin{1});
end




end
