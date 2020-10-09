function R=fLoadR_AL(AL_DIR,RootName,IB,nB,Blade)
% Loads Actuator Line Data for blade index(es) IB


% AL_DIR: directory containing results, with or without final slash
% nB: Number of blades
% IB: vector of blade index. If only one index provided, results are not structured
% R: rotor radius

if(AL_DIR(end)=='/'), AL_DIR=AL_DIR(1:end-1); end


% Just loading the data
disp(['Loading Actuator Line Data in ' AL_DIR])
MCT=load(sprintf('%s/%s.CT01',AL_DIR,RootName));
n1=size(MCT,1);

for iB=IB
    MFp{iB}=load(sprintf('%s/%s.Fp01.bl%02d',AL_DIR,RootName,iB));
    n1=min(size(MFp{iB},1),n1);
    try
        MTp{iB}=load(sprintf('%s/%s.Tp01.bl%02d',AL_DIR,RootName,iB));
        bTp=true;
        n1=min(size(MTp{iB},1),n1);
    catch
        disp(sprintf('Impossible to load %s/%s.Tp01.bl%02d',AL_DIR,RootName,iB));
        bTp=false;
    end
end
MCT=MCT(1:n1,:);


% HACk
if(size(MCT,1)==2044)
    MCT2=zeros(3035,size(MCT,2));
    MCT2(1:2044,:)=MCT;
    MCT=MCT2;
    MCT(2045:end,3)=MCT(2044,3);
    MCT(2045:end,4)=MCT(2044,4);
    dt= MCT(2,1)-MCT(1,1) 
    omega=424.5*2*pi/60;
    dpsi= omega*dt*180/pi;
    MCT(2045:end,2)=mod(MCT(2044,2)+ dpsi*(1:991),360);
    MCT(2045:end,1)=MCT(2044,1)+ dt*(1:991);
end


% Simple data
R.vt   = MCT(:,1) ; 
R.dt   = MCT(2,1)-MCT(1,1) ; 
R.vpsi = MCT(:,2) ; 
R.CT   = MCT(:,3) ; 
R.CP   = MCT(:,4) ; 
R.nB=nB;
R.nr=length(MFp{iB}(1,4:end))/nB;
R.nt=length(R.vt);
R.vr=linspace(0,Blade.R,R.nr)';

% Interpolating blade data to CFD vr data 
twist=interp1(Blade.r,Blade.twist,R.vr); % deg
twist=twist+Blade.pitch;

% Blade Loads and velocities
for iB=IB
    Loads=reshape(MFp{iB}(1:n1,4:end),R.nt,R.nr,3);

    R.vpsi_b_t0{iB} = MFp{iB}(1,3)          ;  % Blade azimuthal position at t_0
    R.Fs{iB}        = squeeze(Loads(:,:,1)) ; 
    R.Ft{iB}        = squeeze(Loads(:,:,2)) ; 
    R.Fn{iB}        = squeeze(Loads(:,:,3)) ; 
    
    R.Ftc{iB}=zeros(size(R.Ft{iB}));
    R.Fnc{iB}=zeros(size(R.Fn{iB}));

    for ir=1:R.nr
        R.Ftc{iB}(:,ir)  = R.Ft{iB}(:,ir)*cosd(twist(ir))-R.Fn{iB}(:,ir)*sind(twist(ir));
        R.Fnc{iB}(:,ir)  = R.Ft{iB}(:,ir)*sind(twist(ir))+R.Fn{iB}(:,ir)*cosd(twist(ir));
    end

    if bTp
        Vel  =reshape(MTp{iB}(1:n1,4:end),R.nt,R.nr,3);
        R.Us{iB}        = squeeze(Vel(:,:,1))   ; 
        R.Ut{iB}        = squeeze(Vel(:,:,2))   ; 
        R.Un{iB}        = squeeze(Vel(:,:,3))   ; 
    end
end
% If user request only one blade, we remove the cell structure..
psiOff=270;
if length(IB)==1
    if(IB==1)
        psiOff=270;
    elseif(IB==2)
        psiOff=270+120
    else
        psiOff=270-120;
    end

    R.vpsi_b_t0=R.vpsi_b_t0{IB};
    R.Fs       =R.Fs{IB}       ;
    R.Ft       =R.Ft{IB}       ;
    R.Fn       =R.Fn{IB}       ;
    R.Ftc      =R.Ftc{IB}      ;
    R.Fnc      =R.Fnc{IB}      ;
    if bTp
        R.Us       =R.Us{IB}       ;
        R.Ut       =R.Ut{IB}       ;
        R.Un       =R.Un{IB}       ;
    else
        R.Us       =0*R.Fs       ;
        R.Ut       =0*R.Ft       ;
        R.Un       =0*R.Fn       ;
    end
end

% Finding last psi period
vpsiR=mod(R.vpsi+psiOff,360);
MI=getIntervals(vpsiR>=0 & vpsiR<=90); % detecting 0 crossing towards 90 deg
R.nPeriods=size(MI,1)-1; % we don't know if the last period is complete, so we'll take the one before
if(R.nPeriods>0)
    AL_DIR
    if isequal(AL_DIR, 'data/ActuatorLine/Runsdt1e-4Shen')
        disp('HACK AL perdiods')
        nn=R.nPeriods-20;
        It_R=MI(nn-1,2):(MI(nn,2)-1); % Selecting the last period
    else
        It_R=MI(end-1,2):(MI(end,2)-1); % Selecting the last period
    end
    R.psiT=mod(vpsiR(It_R),360); % azimuth vector for the last period
    R.IT=It_R;
else
    R.psiT=R.vpsi;
    R.IT=1:nt;
end

% Finding positions where the blade is at NSEW
R.ipsi0   = whichvalue(R.psiT,0) ; 
R.ipsi90  = whichvalue(R.psiT,90) ; 
R.ipsi180 = whichvalue(R.psiT,180) ; 
R.ipsi270 = whichvalue(R.psiT,270) ; 
R.ipsiT=[R.ipsi0 R.ipsi90  R.ipsi180 R.ipsi270 ];
