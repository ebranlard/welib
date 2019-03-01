function [ G kappa eps_over_k epsz cs cp e eta a0 ] = fTheodorsenAtOperatingPoint( l_bar, nB, vx ,w_bar )
% -- NOTES
% Simple application of Theodorsen's theory - AS FUNCTION OF FAR WAKE PARAMETERS !!!
% See fTheodorsenFarWakeParams to find those, and pplication scripts such as Expansion Comparisons and FarWake stuff (though they might be moved soon)
% We need to get the derivative of the mass coefficient, for this we use two points surrounding l_bar

% -- INPUT/OUTPUT
% Usual Parameters and variable names according to Theodorsen


[ G ] = fGoldsteinFactor( l_bar,nB,vx );
kappa=2*trapz(vx,G.*vx);


%% Loss function estimate - Theodorsen p37
dl_bar=10^-3;
vl_bar=l_bar+([-dl_bar dl_bar]); % we estimate the derivative at two points surrounding the operational point l_bar
k=zeros(1,length(vl_bar));
for ilb=1:length(vl_bar)
    [ Gtmp ] = fGoldsteinFactor( vl_bar(ilb),nB,vx );
    k(ilb)=2*trapz(vx,Gtmp.*vx); % mass coefficient
end
dkdl=diff(k)./diff(vl_bar); 
k0=mean(k);
eps_over_k=1+0.5*l_bar/k0*dkdl;
epsz=eps_over_k*k0;


% Now everything is known
cs=2*kappa*w_bar*(1+w_bar*(1/2+eps_over_k));
e=2*kappa*w_bar^2*(1/2+eps_over_k*w_bar);
cp=2*kappa*w_bar*(1+w_bar)*(1+eps_over_k*w_bar);

eta=(1+w_bar*(0.5+eps_over_k))./((1+w_bar).*(1+eps_over_k*w_bar));

a0=(1/2*w_bar+eps_over_k*w_bar^2)/(1+w_bar*(1/2+eps_over_k));


end

