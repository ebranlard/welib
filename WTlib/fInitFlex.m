function WT = fInitFlex(Files,WT,Opts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Initialize properties
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Scanning files
Iprof=whichfile(Files, 'BladeProfile');
Igeom=whichfile(Files, 'BladeGeometry'); %like AE
Ispec=whichfile(Files, '(Spec\.dat)$');

%% Read spec if present
if(~isempty(Ispec))
%     fprintf('fInitFlex:\t Reading Spec file\n')
    WT=fReadSpec(WT,Files{Ispec});
end


if(isempty(Igeom))
    error('Provide BladeGeometry file');
else
    Data=dlmread(Files{Igeom},'',1,0);
    r=Data(:,1);
    chord=Data(:,2);
    beta=Data(:,3);
    t_rel=Data(:,4);
    if(max(t_rel)>10)
        t_rel=t_rel/100;
    end

end
% Loading reference file for PROFILES
if(isempty(Iprof))
    error('Provide BladeProfiles file');
else
    MainPath=dirname(Files{Iprof});
    fid = fopen(Files{Iprof}, 'r');
    Buffer = textscan(fid, '%f %f %f %s %s');
    fclose(fid);
    WT.Profiles.r=Buffer{1};          %radial positions on the blades of the profiles [m]
    WT.Profiles.n=length(Buffer{1});   %number of profiles [.]
    %chord=Buffer{2}; %chord[m]
    %beta=Buffer{3};  %beta in [deg]
end

%% Loading profiles
 for p=1:WT.Profiles.n
     M=load([MainPath,eval(Buffer{4}{p})]);
     WT.Profiles.alpha(:,p)=M(:,1);
     WT.Profiles.Cl(:,p)=M(:,2);
     WT.Profiles.Cd(:,p)=M(:,3);
     if(length(M(1,:))==6)
        WT.Profiles.Cl_inv(:,p)=M(:,5);
        WT.Profiles.Cl_fs(:,p)=M(:,6);
        WT.Profiles.f_st(:,p)=M(:,4);
     end
end   
clear('Buffer');clear('M');   
    
%% Stations and Profiles Sets interpolated for Ngrid data
%!!! FLEX CONVENTION Put Hub length in R since the beginning 
WT.Sources.Rotor.r             = r;
WT.Sources.Rotor.chord         = chord ;
WT.Sources.Rotor.thickness_rel = t_rel;
if(mean(beta)<0) sign=-1; else sign =1; end
WT.Sources.Rotor.twist        = sign*beta;




    %% Loading modes 
%     Rotor.Blade.eigen1e=zeros(ne,3);
%     Rotor.Blade.eigen1f=zeros(ne,3);
%     Rotor.Blade.eigen2f=zeros(ne,3);


%% FLEX CONVENTION - Loading Geometry and Aeroelastic properties
% Data=load('../data/TjaerborgUsedInterpolatedData.mat','-ASCII');
% Rotor.Blade.Mass=[0;Data(:,9)];
% Rotor.Blade.r=[0;Data(:,1);]
% Rotor.Blade.eigen1e=load('../data/EigenMore/eigen1e.dat','-ASCII');
% Rotor.Blade.eigen1f=load('../data/EigenMore/eigen1f.dat','-ASCII');
% Rotor.Blade.eigen2f=load('../data/EigenMore/eigen2f.dat','-ASCII');
%Data=dlmread('data/TjaereborgBladeGeometry.csv',',',1,0);


%% FLEX CONVENTION - Profiles Aero coeff
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

