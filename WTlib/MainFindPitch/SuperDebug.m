InitClear;
require('WTlib','v02');
require('VC_LIB_C','v03');
require('VC_LIB_MAT','v02');
require('BEM','v02');
%%
WS=20;
RPM=27.1;
Ngrid = 30;
TipLoss=1;
pitch=0;
MXIT=200
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5%
user.MainPath = '/work/data/WT/NTK500p/';
user.HtcFileName = 'HtcFindPitch.htc';
user.MainPath = '/work/data/WT/NTK500/';
user.HtcFileName = 'NTK500.htc';
user.BladeMainBodyName = 'blade1';
par.eta=0.93;
par.Pr=500000.0;
par.Pref = par.Pr/par.eta;
par.r = 1.5; % Hub radius
% AERO Parameters
par.Omega = RPM/60*2*pi;
par.rho = 1.225;
par.Uvec = WS
% BEM parameters
par.k = [-0.0017077,0.251163,0.0544955,0.0892074]; % glauert correction
par.relax = 0.85;  % relaxation of solution of induction
par.Induction = 1;
par.TipLoss = TipLoss;
% load hawc2 data
Data = ReadHtcFileAll(user);


dbeta = 0.01; % used for computing slope

%% stations without hub length
% find thickness and profile data
if par.r > 0
    par.sfull = linspace(Data.AeData(1,1),Data.AeData(end,1),Ngrid+1)
    par.s = par.sfull(1:end-1)';
else % skip center point if no hub element
    par.sfull = linspace(Data.AeData(1,1),Data.AeData(end,1),Ngrid+2);
    par.s = par.sfull(2:end-1)';
end

for i=1:Ngrid
    par.c(i,1) = interp1(Data.AeData(:,1),Data.AeData(:,2),par.s(i));
    par.thickness(i,1) = interp1(Data.AeData(:,1),Data.AeData(:,3),par.s(i));
    par.theta(i,1) = interp1(Data.PitchAxis(:,4),Data.PitchAxis(:,5),par.s(i))*pi/180;
   profilesup=find(Data.ThicknessVec(:,1) >= par.thickness(i));
   profileinf=find(Data.ThicknessVec(:,1) <=par.thickness(i));
   profileinf=max(profileinf);
   if isempty(profilesup)
%        disp(' oué oué')
       profilesup=profileinf;
   end
   profilesup=profilesup(1);
   par.ProfileSet(i,:) = [(profileinf~=profilesup)+1,profileinf,profilesup];
end

   profilesup=find(Data.ThicknessVec(:,1) >= par.thickness(i));
   profileinf=find(Data.ThicknessVec(:,1) <=par.thickness(i));
   profileinf=max(profileinf);
   if isempty(profilesup)
%        disp(' oué oué')
       profilesup=profileinf;
   end
   profilesup=profilesup(1);
   par.ProfileSet(i,:) = [(profileinf~=profilesup)+1,profileinf,profilesup];

%% stations now have hub length
par.s = par.s+par.r;
par.sfull = par.sfull+par.r;
par.sigma = Data.Nb*par.c./(2*pi*par.s);
par.an = zeros(Ngrid,1);
par.at = par.an;
par.Uinf = WS;
i=1;
res = 1;
cpt=0; 
 while res > 1e-6;
%  while cpt<MXIT;
    cpt=cpt+1;
    % finding wind speed and angle of attack
    Ut = par.Omega*par.s.*(1+par.at);
    Un = par.Uinf*(1-par.an);
    U = sqrt(Un.^2+Ut.^2);
    psi = atan2(Un,Ut); %(this returns at positive flow angle)

    alpha = psi+par.theta-pitch;
    ClCdCm = LiftDataFun(alpha*180/pi,par,Data);
    Ct = ClCdCm(:,1).*sin(psi)-ClCdCm(:,2).*cos(psi);
    Cn = ClCdCm(:,1).*cos(psi)+ClCdCm(:,2).*sin(psi);

    % Induced velocities
    an_old = par.an;
    at_old = par.at;
    if par.TipLoss == 1 % Prandtl's tip loss, as described in M.O.L.Hansen's book
        Ftiploss = Data.Nb*(par.sfull(end)-par.s)./(2*par.s.*sin(psi));
        Ftiploss = 2./pi*acos(exp(-Ftiploss));
    else
        Ftiploss = 1;
    end
    if par.Induction == 1 % induction as done in HAWC2
        CT = (U.^2.*Cn.*par.c*Data.Nb)./(2*pi.*par.s*par.Uinf^2);
        CT = CT./Ftiploss;
        temp = par.k(4)*CT.^3+par.k(3)*CT.^2+par.k(2)*CT+par.k(1);
        par.an = par.relax*par.an+(1-par.relax)*temp;
        par.at = (U.^2.*Ct.*par.c*Data.Nb)./(8*pi.*par.s.^2.*(1-par.an)*par.Uinf*par.Omega);
    end
    res = norm(an_old-par.an)+norm(at_old-par.at);
end
cpt
Ft = 1/2*par.rho*par.c.*U.^2.*Ct;
Fn = 1/2*par.rho*par.c.*U.^2.*Cn;

if par.r > 0
    P = trapz(par.sfull',[Data.Nb*Ft*par.Omega;0].*(par.sfull'));
else % skip center point if no hub element
    P = trapz(par.sfull',[0;Data.Nb*Ft*par.Omega;0].*(par.sfull'));
end
beta(i) = pitch;
disp(['Uinf = ',num2str(par.Uvec(i)),'  P = ',num2str(P(i)),'  pitch = ',num2str(beta(i)*180/pi)])
% compute dP/dpitch, frozen wake
par.Induction = 0;
Cp = P./(0.5*par.rho*par.Uvec.^3*pi*par.sfull(end)^2);
par.P=P;







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5%

%% To be Made More general
sGeometry='NTK500'; Format='hawc'; opts.bRoughProfiles=logical(0);
[ WT ]   = fInitWT( sGeometry ,Format,PATH.DATA_WT);
WT.Nacelle.tilt=0;
[ WT ]   = fSetRotorGrid(Ngrid-1,WT);
[ WT ]   = fSetRotorGrid(par.s,WT);
[ Sim ]  = fInitSim( WT ); %load default specs as sim
[ Sim ]  = fSetSim( Sim, WT,WS ,RPM,pitch,0  ); 
[ Wind ] = fInitWind(  ); 
[ Wind ] = fSetWind(Wind,Sim  ); 

% Setting algorithm parameters
[ BEMAlgo ]   = fInitBEMAlgo();
[ VFILAlgo ]  = fInitVFilAlgo();
% BEMAlgo.Correction='none';
BEMAlgo.bHubLoss=0;
BEMAlgo.bReInterp=0;
BEMAlgo.aTol=10^-6;
BEMAlgo.relaxation=0.15;
BEMAlgo.Correction='Hawc';
BEMAlgo.bTipLoss=TipLoss;
BEMAlgo.nbIt=MXIT;

BEM = fBEMsteady(WT,Sim,Wind,BEMAlgo);





%%
% figure
% plot(par.Uvec,P/1000)
% grid on
% 
% 
% 
% 
% 
% 
% 
% % 
% 
% 
% 
% 
% 
% 
% 
% Codes={BEM}; legds={'BEM'};
% % Codes={BEM,VCFIL} ; legds={'BEM','VCFIL'};
% colrs=fColrs(1:4);
% sty={'-','+-','--'};
% R=WT.Rotor.R;
% figure, fplotCodesComparison('WS','Power',Codes,legds,colrs,sty,1,1,[],[],'','')
% load('/work/lib/WTlib/v02/MainFindPitch/NTK500p_stall.mat')
% plot(par.Uvec,par.P/1000,'k-+');
