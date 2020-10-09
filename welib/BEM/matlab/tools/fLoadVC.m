function R=fLoadVC(SIM_FILE,Blade,U0,yaw)
% yaw her is given as "yaw_hawc" or "yaw_manu" in degrees (it is positive around yaw axis z)
% (ie the mexico experiment 30deg will be -30)
% NEED FOR SIM INPUT



disp('Loading VortexCode output..');
run(SIM_FILE);
[nt,R.nr]=size(Wings(1).alpha);
I=1:nt;
I=I(Wings(1).Cl(:,floor(R.nr/2))~=0); % looking mid span, if Cl==0, then probably no export was done at this time 

% OK, that's the current stupid format
R.nt=length(I); % that's the real length where stuff have been output
R.alpha     = zeros ( R.nt , R.nr )   ; 
R.alpha0    = zeros ( R.nt , R.nr )   ; 
R.Cl        = zeros ( R.nt , R.nr )   ; 
R.Cd        = zeros ( R.nt , R.nr )   ; 
R.Cm        = zeros ( R.nt , R.nr )   ; 
R.L         = zeros ( R.nt , R.nr )   ; 
R.D         = zeros ( R.nt , R.nr )   ; 
R.M         = zeros ( R.nt , R.nr )   ; 
R.Gamma     = zeros ( R.nt , R.nr )   ; 
R.Vrel_norm      = zeros ( R.nt , R.nr)  ; 
R.Vrel_orth_norm = zeros ( R.nt , R.nr)  ; 
R.Vrel      = zeros ( R.nt , R.nr , 3 )  ; 
R.Vrel_orth = zeros ( R.nt , R.nr , 3 )  ; 
R.V0        = zeros ( R.nt , R.nr , 3 )  ; 
R.Loads     = zeros ( R.nt , R.nr , 3 )  ; 
R.Moments   = zeros ( R.nt , R.nr , 3 )  ; 
R.Loads_s   = zeros ( R.nt , R.nr , 3 )  ; 
R.Moments_s = zeros ( R.nt , R.nr , 3 )  ; 
R.vr=Wings(1).scp'; % because all the other ones are column, but I like it better inline
R.s=Wings(1).scoord;
R.nr=length(R.vr); 

R.dt=time/R.nt; % should be ok
R.vt=(R.dt:R.dt:time)'; % note the prime for harmony with other codes...


% Hawc convention : psi=0 blade downward
% Measurement convention : psi=0 blade upward
R.vpsi_hawc=Blade.omega*R.vt*180/pi; % [deg!]
R.vpsi=mod(R.vpsi_hawc+180,360);

%% Using my convention (yaw articles)
R.vpsi_manu=mod(R.vpsi_hawc-90,360);
MI=getIntervals(R.vpsi_manu>=0 & R.vpsi_manu<=90);
R.nPeriods=size(MI,1)-1; % we don't know if the last perdio is complete
R.IT_manu=MI(end-1,2):(MI(end,2)-1); % Selecting the last period
R.psiT_manu=mod(R.vpsi_manu(R.IT_manu),360); % azimuth vector for the last period


% R.alpha_BEM = zeros ( nt , 40 )  ; 
% R.Cl_BEM    = zeros ( nt , 40 )  ; 
% R.Cd_BEM    = zeros ( nt , 40 )  ; 
% R.Cm_BEM    = zeros ( nt , 40 )  ; 
for i=1:R.nt
    R.alpha     (i,:)   = Wings(1).alpha     (I(i),:);
    R.alpha0    (i,:)   = Wings(1).alpha0    (I(i),:);
    R.Cl        (i,:)   = Wings(1).Cl        (I(i),:);
    R.Cd        (i,:)   = Wings(1).Cd        (I(i),:);
    R.Cm        (i,:)   = Wings(1).Cm        (I(i),:);
    R.L         (i,:)   = Wings(1).L         (I(i),:);
    R.D         (i,:)   = Wings(1).D         (I(i),:);
    R.M         (i,:)   = Wings(1).M         (I(i),:);
    R.Gamma     (i,:)   = Wings(1).Gamma     (I(i),:);
    R.Vrel      (i,:,:) = Wings(1).Vrel      (I(i),:,:);
    R.Vrel_orth (i,:,:) = Wings(1).Vrel_orth (I(i),:,:);
    R.V0        (i,:,:) = Wings(1).V0        (I(i),:,:);
    R.Loads     (i,:,:) = Wings(1).Loads     (I(i),:,:);
    R.Moments   (i,:,:) = Wings(1).Moments   (I(i),:,:);
    Loads_s   (i,:,:) = Wings(1).Loads_s   (I(i),:,:); % ClledFnc and Ftc
    Moments_s (i,:,:) = Wings(1).Moments_s (I(i),:,:);
end
for ir=1:R.nr
    R.Vrel_orth_norm(:,ir) = norm3d(squeeze(R.Vrel_orth(:,ir,:)));
    R.Vrel_norm(:,ir)      = norm3d(squeeze(R.Vrel(:,ir,:)));
end
R.Ftc=Loads_s(:,:,1);
R.Fnc=Loads_s(:,:,2); % Hawc/Omnivor Coordinate y is normal to the profile section
R.Mc= Moments_s(:,:,3);
R.Moments   =NaN; % TODO in Omnivor
% R=Wings(1);

if yaw~=0 
    warning('Omnivor postpro with yaw/tilt is clearly incomplete...')
end

%% Induced velocities

%
R.Un     =squeeze(R.Vrel(:,:,2))*cosd(yaw); % TODO TODO TODO tilt?
R.Uin    =R.Un-U0*cosd(yaw);
R.Ut_rel = zeros(R.nt,R.nr) ; 
R.Utw    = zeros(R.nt,R.nr) ; 
R.Uit    = zeros(R.nt,R.nr) ; 
% THIS BELOW IS NOT GOOD FOR YAW
for ir=1:R.nr
    % Using psi manu, because that's what I used for a scheme on paper, might as well use psi_racc (but different formula)
    R.Ut_rel(:,ir)  = -R.Vrel(:,ir,1).*sind(R.vpsi_manu)-R.Vrel(:,ir,3).*cosd(R.vpsi_manu) ;  % Should be negative!!!!!!
    R.Ut(:,ir) = R.Ut_rel(:,ir)-(-Blade.omega*R.vr(ir))                                           ;  % The induction due to wind and induction % 
    R.Uit(:,ir) = R.Ut(:,ir)- (-U0*sind(yaw)*sind(R.vpsi_manu))                   ; 
end


% TEMPORARY HACK DUE TO LOAD EXPORT BUG
if isequal(SIM_FILE,'data/VortexCode/results_15_00_4_40sec_3.0s_dt3_p2.3/sim_u15.00_y0.0_4/mexico_out.m');
    disp('Using temporary hack for this simulation folder')
    R.Ftc=-R.D.*cosd(R.alpha)+R.L.*sind(R.alpha);
    R.Fnc= R.L.*cosd(R.alpha)+R.D.*sind(R.alpha);
    R.Loads     =NaN;
end


% Interpolating blade data to hawc2 data
twist=interp1(Blade.r,Blade.twist,R.vr); % deg
% 
R.Ft(:,:)  = R.Ftc(:,:)*cosd(twist(:))+R.Fnc(:,:)*sind(twist(:));
R.Fn(:,:)  =-R.Ftc(:,:)*sind(twist(:))+R.Fnc(:,:)*cosd(twist(:));





% Finding last psi period % TODO TODO
vpsiR=mod(R.vpsi,360);
MI=getIntervals(vpsiR>=0 & vpsiR<=90); % detecting 0 crossing towards 90 deg
R.nPeriods=size(MI,1)-1; % we don't know if the last period is complete, so we'll take the one before
if(R.nPeriods>0)
    It_R=MI(end-1,2):(MI(end,2)-1); % Selecting the last period
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
