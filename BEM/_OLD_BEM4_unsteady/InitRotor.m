global Rotor Profiles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Rotor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Params
Rotor.R = 20.54; % Blade length [m]
Rotor.rhub = 1.5; % Hub  length [m]
Rotor.Ngrid = 10;
Rotor.Omega=27.1*2*pi/60; % Nominal Rotational speed [rad/s]
Rotor.rb_center_in4=[0; 0; 0];

Rotor.M=3*1.844; % mass [kg]
Rotor.I=0.482*10^6;% [kg/m/m]??? kg m^2
%Rotor.I=3*trapz([0;Rotor.r],[0;Rotor.r.^2.*Rotor.Blade.Mass])
Rotor.cone=0;  % cone angle [deg]
Rotor.nB=3;            %number of blades

%% HAWC CONVETION - Profiles
Profiles = fReadPcFile('data/NewBlade.pc');


%% HAWC CONVETION - Geometry
Data = fReadHtcFile('data/NTK500_pitch_structure.htc','blade1');
AeData = fReadAeFile('data/NewBlade82.ae');  %chord and thickness
Stations=Data.PitchAxis;clear('Data');       %twist

% Creating linspace stations, first stations without hub length
if Rotor.rhub > 0
    Rotor.rfull = linspace(Stations(1,4),Stations(end,4),Rotor.Ngrid+1);
    Rotor.r = Rotor.rfull(1:end-1)';
else % skip center point if no hub element
    Rotor.rfull = linspace(Stations(1,4),Stations(end,4),Rotor.Ngrid+2);
    Rotor.r = Rotor.rfull(2:end-1)';
end
Rotor.chord=interp1(AeData(:,1),AeData(:,2),Rotor.r) ;
Rotor.thickness_rel = interp1(AeData(:,1),AeData(:,3),Rotor.r);
Rotor.beta = - interp1(Stations(:,4),Stations(:,5),Rotor.r); 
for i=1:Rotor.Ngrid
   profilesup=find(Profiles.ThicknessVec(:,1) >= Rotor.thickness_rel(i));
   profileinf=max(find(Profiles.ThicknessVec(:,1) <= Rotor.thickness_rel(i)));
   if isempty(profilesup)
       profilesup=profileinf;
   end
   profilesup=profilesup(1);
   Rotor.ProfileSet(i,:) = [(profileinf~=profilesup)+1,profileinf,profilesup];
end
% stations now have hub length
Rotor.r = Rotor.r+Rotor.rhub;
Rotor.rfull = Rotor.rfull+Rotor.rhub;



%% FLEX CONVENTION - Loading Geometry and Aeroelastic properties
% Data=load('../data/TjaerborgUsedInterpolatedData.mat','-ASCII');
% Rotor.Blade.Mass=[0;Data(:,9)];
% Rotor.Blade.r=[0;Data(:,1);]
% Rotor.Blade.eigen1e=load('../data/EigenMore/eigen1e.dat','-ASCII');
% Rotor.Blade.eigen1f=load('../data/EigenMore/eigen1f.dat','-ASCII');
% Rotor.Blade.eigen2f=load('../data/EigenMore/eigen2f.dat','-ASCII');
%Data=dlmread('data/TjaereborgBladeGeometry.csv',',',1,0);
% Data=dlmread('data/SmallWTBladeGeometry.csv',',',1,0);
% Rotor.r=Data(:,1);
% Rotor.chord=Data(:,2);
% Rotor.beta=Data(:,4);
% %Rotor.Blade.Mass=Data(:,9);
% Rotor.ne=length(Rotor.r);    %number of element or stations [.]
% clear('Data')

%% FLEX CONVENTION - Profiles Aero coeff
% % Loading reference file for PROFILES
% fid = fopen('data/SmallWTBladeProfiles.dat', 'r');
%      Buffer = textscan(fid, '%f %f %f %s %s');
% fclose(fid);
% Rotor.Profiles.rP=Buffer{1};           %radial positions on the blades of the profiles [m]
% Rotor.Profiles.nP=length(Buffer{1});   %number of profiles [.]
% % Rotor.chord=Buffer{2}; %chord[m]
% % Rotor.beta=Buffer{3};  %beta in [deg]
% 
% % Loading profiles
% for p=1:Rotor.Profiles.nP
%     M(:,:,p)=load(['data/' eval(Buffer{4}{p})]);
%     Rotor.Profiles.alpha(p,:)=squeeze(M(:,1,p));
%     Rotor.Profiles.Cl(p,:)=squeeze(M(:,2,p));
%     Rotor.Profiles.Cd(p,:)=squeeze(M(:,3,p));
% %     Rotor.Profiles.Cl_inv(p,:)=M(:,5,p);
% %     Rotor.Profiles.Cl_fs(p,:)=M(:,6,p);
% %     Rotor.Profiles.f_st(p,:)=M(:,4,p);
% end
% clear('Buffer');clear('M');
% 
% %% Loading modes 
% Rotor.Blade.eigen1e=zeros(ne,3);
% Rotor.Blade.eigen1f=zeros(ne,3);
% Rotor.Blade.eigen2f=zeros(ne,3);
% % 
% % m1e=load('data/eigen1e.dat','-ASCII');
% % m1f=load('data/eigen1f.dat','-ASCII');
% % m2f=load('data/eigen2f.dat','-ASCII');
% % if(length(m1e(1,:))==2 )
% %     Rotor.Blade.eigen1e(:,2:3)=m1e;
% %     Rotor.Blade.eigen1f(:,2:3)=m1f;
% %     Rotor.Blade.eigen2f(:,2:3)=m2f;
% % end    
% % Rotor.Blade.eigenF=[1.2146 2.3908 3.5659];
% 
% 
% for e=1:ne
%    [mmm Rotor.pe(e)]=min(abs(Rotor.Profiles.rP-Rotor.r(e)));
% end

clear('mmm')
clear('fid')
clear('e')
clear('p')

% Element used for the estimate of khi (at r/R=70%)
[mmm Rotor.e_ref_for_khi]=min(abs(Rotor.r/Rotor.R-0.7));
Rotor.ne=length(Rotor.r);
nB=Rotor.nB;
ne=Rotor.ne;