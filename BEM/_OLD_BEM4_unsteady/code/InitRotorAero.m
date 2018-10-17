%% Rotor
Rotor.R=30.56;        % Blade length [m]

%Loading reference file for geometry
fid = fopen('../data/tjaereblade.dat', 'r');
    Buffer = textscan(fid, '%f %f %f %s %s');
fclose(fid);

Rotor.r=Buffer{1};     %radial positions on the blades [m]
Rotor.chord=Buffer{2}; %chord[m]
Rotor.beta=Buffer{3};  %beta in [deg]
Rotor.ne=length(Rotor.r);    %number of element [.]
Rotor.nB=3;            %number of blades

% Loading profiles
for e=1:Rotor.ne
    M(:,:,e)=load(['../data/' eval(Buffer{4}{e})]);
    Profiles.alpha(e,:)=M(:,1,e);
    Profiles.Cl(e,:)=M(:,2,e);
    Profiles.Cd(e,:)=M(:,3,e);
    Profiles.Cl_inv(e,:)=M(:,5,e);
    Profiles.Cl_fs(e,:)=M(:,6,e);
    Profiles.f_st(e,:)=M(:,4,e);
end
Rotor.Profiles=Profiles;
% Element used for the estimate of khi (at r/R=70%)
[m Rotor.e_ref_for_khi]=min(abs(Rotor.r/Rotor.R-0.7));

