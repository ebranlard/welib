function R= fLoadHawc2(SIM_DIR,sfilter,Blade)
% Open Hawc2 result file in a given folder and according to a given filter (should be such that only one sel file is returned)

% Looking for proper file
if(SIM_DIR(end)=='/'), SIM_DIR=SIM_DIR(1:end-1); end
Files=dir([SIM_DIR '/*' sfilter '*.sel']);
if length(Files)==0
    fprintf('! Can''t find a file *%s*.sel in folder: %s \n',sfilter,SIM_DIR);
    error('')
    return
elseif length(Files)>1
    fprintf('! Pattern *%s*.sel not unique to select in folder: %s \n',sfilter,SIM_DIR);
    error('')
    return
else
    [~,basefile]=fileparts(Files(1).name);
    resfile=sprintf('%s/%s',SIM_DIR,basefile);
end


[Results,Freq,Time,Flag,Binary,labels,units,descriptions] = fReadHawc2(resfile); 
I=1:length(labels);
iPower=I(cellfun(@(x)~isempty(x),strfind(labels,'power')));
iThrust=I(cellfun(@(x)~isempty(x),strfind(labels,'thrust')));
R.dt=1/Freq;
R.vt=R.dt:R.dt:Time;


%% Scanning labels to detect radial positions
A1=cellfun(@(x)sscanf(x,'%*s R= %f'),labels,'UniformOutput',false);
A2=cellfun(@(x)sscanf(x,'%*s%*s%*s R= %f'),labels,'UniformOutput',false);

A=cellfun(@(x)regexp(x{1},'(?<=radius[ ]*)(\d+.\d+)','match','once'),descriptions,'UniformOutput',false);
R.vr= sort(unique(str2num(cell2mat(A))));
R.vr_lowres= sort(unique([cell2mat(A1);cell2mat(A2)]));

fprintf('Hawc2 res file has data at the following r values:\n  ');
fprintf('%.2f ',R.vr);
fprintf('\n  ');
fprintf('%.2f ',R.vr_lowres);
fprintf('\n');

if(length(R.vr)~=length(R.vr_lowres))
    warning('Hawc2 radius in output are not unique due to simple format')
    R.vr=R.vr_lowres;
end

R.nr=length(R.vr);

% Interpolating blade data to hawc2 data
twist=interp1(Blade.r,Blade.twist,R.vr); % deg

    


%% Looping on radial positions and trying to see what we got. If data are not present we insert NaN values
Iprocessed=[];
for ir=1:R.nr
    r=R.vr_lowres(ir); % USING VR LOW RES TO BE SURE TO MATCH LABEL
    % the function below is defined in
    [R.Cl(:,ir)    ,ires]= fhawc2_signal(labels,Results,'Cl',r)              ; Iprocessed=[Iprocessed(:);ires];
    [R.Cd(:,ir)    ,ires]= fhawc2_signal(labels,Results,'Cd',r)              ; Iprocessed=[Iprocessed(:);ires];
    [R.Cm(:,ir)    ,ires]= fhawc2_signal(labels,Results,'Cm',r)              ; Iprocessed=[Iprocessed(:);ires];
    [R.L(:,ir)    ,ires]= fhawc2_signal(labels,Results,'Lift',r)          ; Iprocessed=[Iprocessed(:);ires];
    [R.D(:,ir)    ,ires]= fhawc2_signal(labels,Results,'Drag',r)          ; Iprocessed=[Iprocessed(:);ires];
    [R.alpha(:,ir) ,ires]= fhawc2_signal(labels,Results,'Alfa',r)            ; Iprocessed=[Iprocessed(:);ires];
    [R.phi(:,ir)   ,ires]= fhawc2_signal(labels,Results,'Inflow ang',r)      ; Iprocessed=[Iprocessed(:);ires];

    % Hawc2: secfrc Fx:   section%sl%secfrc(1)=F*(-section%D*dcos(section%alfa)+section%L*dsin(section%alfa))
    % Hawc2: secfrc Fy:	  section%sl%secfrc(2)=F*(section%L*dcos(section%alfa)+section%D*dsin(section%alfa))
    % Hawc2: secfrc Fz:	  section%sl%secfrc(3)=0.0
    % Forces normal and parallel to the chord !!!!! not to the rotor
    [R.Fnc(:,ir)    ,ires]= fhawc2_signal(labels,Results,'secfrc Fy',r)       ; Iprocessed=[Iprocessed(:);ires];
    [R.Ftc(:,ir)    ,ires]= fhawc2_signal(labels,Results,'secfrc Fx',r)       ; Iprocessed=[Iprocessed(:);ires];

    R.Ft(:,ir)  = R.Ftc(:,ir)*cosd(twist(ir))+R.Fnc(:,ir)*sind(twist(ir));
    R.Fn(:,ir)  =-R.Ftc(:,ir)*sind(twist(ir))+R.Fnc(:,ir)*cosd(twist(ir));
    
    [R.Uin(:,ir)    ,ires]= fhawc2_signal(labels,Results,'Induc. Vy, rpco',r) ; Iprocessed=[Iprocessed(:);ires];
    [R.Uit(:,ir)    ,ires]= fhawc2_signal(labels,Results,'Induc. Vx, rpco',r) ; Iprocessed=[Iprocessed(:);ires];
    [R.Uiy(:,ir)    ,ires]= fhawc2_signal(labels,Results,'Induc. Vy, glco',r) ; Iprocessed=[Iprocessed(:);ires];
    [R.Uix(:,ir)    ,ires]= fhawc2_signal(labels,Results,'Induc. Vx, glco',r) ; Iprocessed=[Iprocessed(:);ires];
    [R.Uiz(:,ir)    ,ires]= fhawc2_signal(labels,Results,'Induc. Vz, glco',r) ; Iprocessed=[Iprocessed(:);ires];
    [R.Vrel(:,ir)    ,ires]= fhawc2_signal(labels,Results,'Vrel',r) ; Iprocessed=[Iprocessed(:);ires];
    [R.Vrel3D(:,ir)    ,ires]= fhawc2_signal(labels,Results,'Vrel 3D',r) ; Iprocessed=[Iprocessed(:);ires];
    
    % Body velocity
    [R.Velx_rp(:,ir)    ,ires]= fhawc2_signal(labels,Results,'Ae vel Vx, rpco',r) ; Iprocessed=[Iprocessed(:);ires];
    [R.Vely_rp(:,ir)    ,ires]= fhawc2_signal(labels,Results,'Ae vel Vy, rpco',r) ; Iprocessed=[Iprocessed(:);ires];

    [R.U0x_rp(:,ir)    ,ires]= fhawc2_signal(labels,Results,'WSP Vx, rpco',r) ; Iprocessed=[Iprocessed(:);ires];
    [R.U0y_rp(:,ir)    ,ires]= fhawc2_signal(labels,Results,'WSP Vy, rpco',r) ; Iprocessed=[Iprocessed(:);ires];
    [R.U0z_rp(:,ir)    ,ires]= fhawc2_signal(labels,Results,'WSP Vz, rpco',r) ; Iprocessed=[Iprocessed(:);ires];
    [R.U0x_gl(:,ir)    ,ires]= fhawc2_signal(labels,Results,'WSP Vx, glco',r) ; Iprocessed=[Iprocessed(:);ires];
    [R.U0y_gl(:,ir)    ,ires]= fhawc2_signal(labels,Results,'WSP Vy, glco',r) ; Iprocessed=[Iprocessed(:);ires];
    [R.U0z_gl(:,ir)    ,ires]= fhawc2_signal(labels,Results,'WSP Vz, glco',r) ; Iprocessed=[Iprocessed(:);ires];

    [R.Gamma(:,ir)    ,ires]= fhawc2_signal(labels,Results,'Gamma',r) ; Iprocessed=[Iprocessed(:);ires];
end


R.L(:,:)=R.L(:,:)*1000;
R.D(:,:)=R.D(:,:)*1000;
R.Fn(:,:)=R.Fn(:,:)*1000;
R.Ft(:,:)=R.Ft(:,:)*1000;
R.Fnc(:,:)=R.Fnc(:,:)*1000;
R.Ftc(:,:)=R.Ftc(:,:)*1000;
R.Un(:,:)= R.U0y_rp+R.Uin;
R.Ut(:,:)= R.U0x_rp+R.Uit;

% Signals that are not radial position dependent
[R.vt       ,ires]= fhawc2_signal(labels,Results,'Time') ; Iprocessed=[Iprocessed(:);ires];
[vpsi_hawc  ,ires]= fhawc2_signal(labels,Results,'Azi') ; Iprocessed=[Iprocessed(:);ires]; % Hopefully blade 1
[R.Thrust   ,ires]= fhawc2_signal(labels,Results,'Ae rot. thrust') ; Iprocessed=[Iprocessed(:);ires];
[R.Power    ,ires]= fhawc2_signal(labels,Results,'Ae rot. power') ; Iprocessed=[Iprocessed(:);ires];
[R.Omega    ,ires]= fhawc2_signal(labels,Results,'Omega') ; Iprocessed=[Iprocessed(:);ires];
[R.U0x      ,ires]= fhawc2_signal(labels,Results,'WSP gl. coo.,Vx') ; Iprocessed=[Iprocessed(:);ires];
[R.U0y      ,ires]= fhawc2_signal(labels,Results,'WSP gl. coo.,Vy') ; Iprocessed=[Iprocessed(:);ires];
[R.U0z      ,ires]= fhawc2_signal(labels,Results,'WSP gl. coo.,Vz') ; Iprocessed=[Iprocessed(:);ires];


R.Power=R.Power*1000;
R.Thrust=R.Thrust*1000;
%
I=1:length(labels);
Iprocessed=sort(unique(Iprocessed));
Idiff=setdiff(I,Iprocessed);
if ~isempty(Idiff)
    fprintf('Label not processed: %s\n',labels{Idiff})
end





%% Preparation for azimuthal plot 
% Hawc convention : psi=0 blade downward
% Measurement convention : psi=0 blade upward
R.vpsi=mod(vpsi_hawc+180,360);
%
% determining number of periods by looking at intervals where 0<psi<90
MI=getIntervals(R.vpsi>=0 & R.vpsi<=90);
R.nPeriods=size(MI,1)-1; % we don't know if the last perdio is complete

R.IT=MI(end-1,2):(MI(end,2)-1); % Selecting the last period
R.psiT=mod(R.vpsi(R.IT),360); % azimuth vector for the last period








% Finding positions where the blade is at NSEW
R.ipsi0   = whichvalue(R.psiT,0) ; 
R.ipsi90  = whichvalue(R.psiT,90) ; 
R.ipsi180 = whichvalue(R.psiT,180) ; 
R.ipsi270 = whichvalue(R.psiT,270) ; 
R.ipsiT=[R.ipsi0 R.ipsi90  R.ipsi180 R.ipsi270 ];


%% Using my convention
R.vpsi_manu=mod(vpsi_hawc-90,360);

MI=getIntervals(R.vpsi_manu>=0 & R.vpsi_manu<=90);
R.nPeriods=size(MI,1)-1; % we don't know if the last perdio is complete

R.IT_manu=MI(end-1,2):(MI(end,2)-1); % Selecting the last period
R.psiT_manu=mod(R.vpsi_manu(R.IT_manu),360); % azimuth vector for the last period

