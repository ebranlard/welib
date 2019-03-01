function [RES] = fBEMsteadyTipLossStandAlone(WT,Sim,Wind,Algo)
cone=WT.Rotor.cone;

V0=Wind.V0(3)*cosd(cone);
VHubHeight=Wind.V0(3);

nB=WT.Rotor.nB;
ne=WT.Rotor.ne;
r=WT.Rotor.r(:)';
dr=WT.Rotor.dr(:)';
rhub=WT.Rotor.rhub;
R=WT.Rotor.R;
chord=WT.Rotor.chord;
twist=WT.Rotor.twist;
omega=Sim.Run.RPM*2*pi/60;
Rotor=WT.Rotor;

sigma=chord*nB./(2*pi*r*cosd(cone));  %!!! CONE

% Environment
rho=Sim.rho;
KinVisc=Sim.KinVisc;

pitch=Sim.Run.PITCH;
omega=Sim.Run.RPM*2*pi/60;
lambda_r=omega*r*cosd(cone)/V0;       %!!! cone
lambda=omega*R*cosd(cone)/V0;         %!!! cone


%algorithm internal paramter
nbIt=Algo.nbIt;
aTol=Algo.aTol;

% initialize vectors
a=zeros(1,ne)+0.3;       %<-------------- Using Converged solution from previous RES -> No
aprime=zeros(1,ne)+0.01; %<-------------- Using Converged solution from previous RES -> No

Pn=zeros(1,ne);
Pt=zeros(1,ne);
Vrel_norm=zeros(1,ne);
Re=zeros(1,ne);
F=zeros(1,ne);
phi=zeros(1,ne);
alpha=zeros(1,ne);
Cl=zeros(1,ne);
Cd=zeros(1,ne);
Cn=zeros(1,ne);
Ct=zeros(1,ne);
CT=zeros(1,ne);

Ftip_aprime=zeros(1,ne);
Fperf=1;
%% preparing tip-loss correction
x=(r-rhub)/(R-rhub);
if(isequal(Algo.BEM.TipLossMethod,'TipLossDB'))
    global PATH;
    load([PATH.TIPLOSSDB 'TipLDBGamma']);
    load([PATH.TIPLOSSDB 'TipLDBF']);
    nGamma=size(TipLDBGamma,2)-4;
    rGamma=single(cos(linspace(1,0,nGamma)*pi/2));
    nF=size(TipLDBF,2)-6;
    rF=cos(linspace(1,0,nF)*pi/2);
%     Gamma_fit=x*0;
%     vSSE2=zeros(1,size(TipLDBGamma,1));
%     jselect=1;
else
    % Position of helical vortex
    rhelix=[WT.Rotor.rhub r(1:end-1)+diff(r)/2 WT.Rotor.R];
    nhelix=length(rhelix);
    lambda_r_helix=omega*rhelix*cosd(cone)/V0;       %!!! cone
    Algo.bGridIn=1;
    Algo.bPolarIn=1;
    Algo.bPolarOut=1;
end



for i=1:nbIt % !!!!!!!!!!!!!!!!!!! in this order for my tip loss correction
    % --------------------------------------------------------------------------------
    % --- Step 0: Relative wind
    % --------------------------------------------------------------------------------

    % --------------------------------------------------------------------------------
    % --- Step 1: Wind Components
    % --------------------------------------------------------------------------------
    Ut = omega*r.*(1+aprime);
    Un = V0*(1-a);
    Vrel_norm = sqrt(Un.^2+Ut.^2); %more stable than  Vrel=V0*(1-a)/sind(phi);
    Re=Vrel_norm.*chord/KinVisc/10^6; % Reynolds number in Millions

    % --------------------------------------------------------------------------------
    % --- Step 2: Flow Angle
    % --------------------------------------------------------------------------------
    phi = atan2(Un,Ut)*180/pi; %(this returns a positive flow angle) [deg]
    IbadCV=find(imag(phi)~=0);% index of bad convergence
    if(~isempty(IbadCV))
        fprintf('Setting %d bad phi zero',length(IbadCV));
        phi(IbadCV)=0;
    end

    % --------------------------------------------------------------------------------
    % --- Tip loss
    % --------------------------------------------------------------------------------
    Ftip=a*0+1;
    Fperf=a*0+1; 
    if(i>1) % first time pass, no tip-losses
        switch(Algo.BEM.TipLossMethod)
            case 'TipLossDB'
                if(i>2 && vSSE2(jselect)*(1.01)>sum(x(ii:end).*(Gamma_fit(ii:end) - GammaNorm(ii:end)).^2)) % there is a trick here because the second condition is not evaluated, and hence the variables are initialized afterwards...
%                              fprintf('%d keep same Ftip\n',i);
                else
                    vSSE2=zeros(1,size(TipLDBGamma,1));
                    ii=whichvalue(x,min(x(GammaNorm==max(GammaNorm)),0.25));
                    for j=1:size(TipLDBGamma,1);
                        Gamma_fit=interp1(rGamma,TipLDBGamma(j,5:nGamma+4),x);
                        x0=TipLDBGamma(j,1);
                        ii=whichvalue(x,min(x0,0.25));
                        vSSE2(j)=sum(x(ii:end).*(Gamma_fit(ii:end) - GammaNorm(ii:end)).^2);
                    end
                    jselect=whichmin(vSSE2);
                    Gamma_fit=interp1(rGamma,TipLDBGamma(jselect,5:nGamma+4),x);
                    params=TipLDBGamma(jselect,1:4);
                    vSSE2(jselect);
                    Ifound=mfind(int16(round(TipLDBF(:,1:4)*1000)),int16(round(params*1000)));
                    LAMBDAS=unique(TipLDBF(Ifound,5)); 
                    CTS=unique(TipLDBF(Ifound,6));
                    lambdaDB=LAMBDAS(whichvalue(LAMBDAS,lambda));
                    CTDB=CTS(whichvalue(CTS,CT));
                    Ftip=TipLDBF(mfind(int16(round(TipLDBF(:,1:6)*1000)), int16(round([params lambdaDB CTDB]*1000))),7:end);
                    if(isempty(Ftip))
                        fprintf('Holes in the DB - this should not really happen...\n');
                        IgoodCT=(TipLDBF(Ifound,6)== CTDB);
                        Ifound=Ifound(IgoodCT);
                        IgoodLambda=whichvalue(TipLDBF(Ifound,5),lambda);
                        Ftip=TipLDBF(Ifound(IgoodLambda),7:end);
                    end
                        %well CT is the most important so we selct for CT, but this hole in the database is not normal
%                     [params lambdaDB CTDB lambda CT]
                    fprintf('%.3f %.3f %.3f %.3f-lDB %.3f %.3f - CTDB %.3f %.3f\n',params,lambdaDB, lambda, CTDB, CT );  % <<<<<<<<<<<<<<<<<<<<<<<< guilty of printing
                    Ftip=interp1(rF,Ftip,x);    
                  figure(), hold all, plot(rGamma,TipLDBGamma(jselect,5:nGamma+4),'r-');, plot(x,GammaNorm,'k'); figure(4), plot(r,Ftip),   set(gcf, 'Position',[  1667 432   560   420])
                                                       
                end
            case 'PrescribedWake'
                GammaBound=Gamma;
                %GammaBound=Gamma*Algo.relaxation+(1-Algo.relaxation)*GammaBoundLast;
                %figure (1200),clf,   plot(r,Gamma,r,GammaBound,r,GammaBoundLast)
                GammaBoundLast=GammaBound;

                Algo.bIt=i; %% nasty for theodorsen

                a_b=a*0; a_avg=a*0; aprime_b=a*0; aprime_avg=a*0; a_inf=a*0; aprime_inf=a*0;
                %% Computation for nb blades on lifting line
                % just on the ll
                vr_CP=r; 
                vpsi_CP=0;
                [ Wake ] = fGeometryPrescribedWake(lambda,nB,V0,R,rhelix, r, GammaBound, CTloc, a, aprime , [] ,[],[],[],[],[],[], Algo );
                [ ui ] = fUi_PrescribedWake( Wake, vr_CP,vpsi_CP,0,1,1,1,Algo); %grid,polin,polout
                a_ll=squeeze(-ui(:,1,:,3))'/V0;
                aprime_ll=squeeze(-ui(:,1,:,2))'./(lambda_r*V0);
                
                %%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%% Tip-loss factor one: Finite vs infinite number of blades 
                %%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% Computation for Infinite blades on lifting line
                if Algo.bTipLossInfiniteB
                    FactInf=60; %infinity factor
                    GammaInf=GammaBound/FactInf;
                    nBInf=nB*FactInf; % trick because this code is nasty and
                    [ WakeInf ] = fGeometryPrescribedWake(lambda,nBInf,V0,R,rhelix, r, GammaInf, CTloc, a, aprime , [] ,[],[],[],[],[],[], Algo );

                    [ uiinf ] = fUi_PrescribedWake( WakeInf, vr_CP,vpsi_CP,0,1,1,1,Algo); %grid,polin,polout
                    a_inf=squeeze(-uiinf(:,1,:,3))'/V0;
                    aprime_inf=squeeze(-uiinf(:,1,:,2))'./(lambda_r*V0);

                    % This seems to make results sensitive
                    nr=floor(length(r)/2)*2;
                    istart=find(a_inf(1:end-3)./a_ll(1:end-3)>1);
                    if(isempty(istart))
                        istart=1;
                    end
                    istart=istart(end);
                end
                %                         keyboard
                % %    
                % weirdGamma=[Wake.Gamma Wake.Gamma(:,end)];
                % figure(1227),clf, set(1227,'Renderer','OpenGL'), hold on
                % surf(Wake.X/R, -Wake.Z/R, Wake.Y/R,weirdGamma,'FaceAlpha',0.5);colorbar
                % grid on,box on; xlabel('x'); ylabel('z'); zlabel('y'); rotate3d on;
                % xlim(1.5*[-1 1])
                % zlim(1.5*[-1 1])
                % view([115 30]);
                %
                % weirdGamma=[WakeInf.Gamma WakeInf.Gamma(:,end)];
                % figure(1229),clf, set(1229,'Renderer','OpenGL'), hold on
                % surf(WakeInf.X/R, -WakeInf.Z/R, WakeInf.Y/R,weirdGamma,'FaceAlpha',0.5);colorbar
                % grid on,box on; xlabel('x'); ylabel('z'); zlabel('y'); rotate3d on;
                % xlim(1.5*[-1 1])
                % zlim(1.5*[-1 1])
                % view([115 30]);

                %%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%% Tip-loss factor two: Comparison with azimuthally average
                %%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% simulation on the rotor plane
                if Algo.bTipLossPlaneAvg
                    vr_CP=r; 
                    vpsi_CP=linspace(0,2*pi,60);
                    [ Wake ] = fGeometryPrescribedWake(lambda,nB,V0,R,rhelix, r, Gamma, CTloc , a, aprime , [] ,[],[],[],[],[],[], Algo );
                    [ uiplane ] = fUi_PrescribedWake( Wake, vr_CP,vpsi_CP,0,1,1,1,Algo); %grid,polin,polout
                    WI=squeeze(uiplane(:,:,1,3));
                    TI=squeeze(uiplane(:,:,1,2));
                    a_b=-WI(:,1)'/V0;
                    a_avg=-mean(WI')/V0;
                    aprime_b=-TI(:,1)'./(lambda_r*V0);
                    aprime_avg=-mean(TI')./(lambda_r*V0);
                end

                %%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%% Tip-loss factor 
                %%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% tip-loss by comparison with infinite number of blades
                nr=floor(length(r)/2)*2;
                istart=nr/2;
                %istart=find(abs((a_inf(1:end-3)-a_ll(1:end-3))./(a_inf(1:end-3)))<0.01);
                % istart=max(istart)-1;

                % the tip losses by comparison with infinite number of blades
                Ftip=0*a+1;
                Ftip(istart:end)=a_inf(istart:end)./a_ll(istart:end); % I'm suspicious about this 
%                 Ftip(1:istart)=Ftip(istart);
                %Ftip=Ftip/max(Ftip); 

                Ftip_aprime=0*a+1;
                Ftip_aprime(istart:end)=aprime_inf(istart:end)./aprime_ll(istart:end); % I'm suspicious about this 
%                 Ftip_aprime(1:istart)=Ftip_aprime(istart);
                %                         Ftip_aprime=Ftip_aprime/max(Ftip_aprime); 
                %% tip loss by comparison with avg
                F_bar_a=a_avg./a_b;
                F_bar_aprime=aprime_avg./aprime_b;
                F_bar_a(1:istart)=1;
                F_bar_aprime(1:istart)=1;
%                 F_bar_a(1:istart)=F_bar_aprime(istart);
%                 F_bar_aprime(1:istart)=F_bar_aprime(istart);
                % F_bar_a=F_bar_a/max(F_bar_a); 
                % F_bar_aprime=F_bar_aprime/max(F_bar_aprime); 
                if Algo.bTipLossPlaneAvg && ~Algo.bTipLossInfiniteB
                    Ftip=F_bar_a;
                    Ftip_aprime=F_bar_aprime;
                end


                istart=max(whichvalues(Ftip>=1,1));
                if(~isempty(istart))
                    Ftip(1:istart)=1;
                    Ftip_aprime(1:istart)=1;
                end


                figure(1),clf
                hold all
                plot(r/R,a_ll,'b')
                plot(r/R,a_inf,'r')
                plot(r/R,a,'g')
                plot(r/R,a_b,'b+')
                plot(r/R,a_avg,'r+')
                plot(r(istart)/R,a_ll(istart),'ko')
                ylabel('a [.]')
                legend('a ll','a infB','a RES','a 2','a avg',0)

                figure(2),clf, hold all
                plot(r/R,Ftip)
                plot(r/R,F_bar_a)
                % Ftip=smooth(Ftip,nhelix*4/100,'rloess')';
                % plot(r/R,Ftip)
                ylim([0 1])
                set(gcf, 'Position',[  1667 432   560   420])
                legend('a Ftip','a Fbara',0)

                figure(3),clf,hold all
                plot(r/R,Ftip_aprime)
                plot(r/R,F_bar_aprime)
                ylim([0 1])
                set(gcf, 'Position',[  467 432   560   420])
                legend('aprime Ftip','aprime Fbara',0)

                figure(4)
                clf
                hold all
                plot(r/R,aprime_ll,'b')
                plot(r/R,aprime_inf,'r')
                plot(r/R,aprime,'g')
                plot(r/R,aprime_b,'b+')
                plot(r/R,aprime_avg,'r+')
                ylim([-0.05 0.05])
                set(gcf, 'Position',[  967 432   560   420])
                legend('ap ll','ap infB','ap RES','ap 2','ap avg',0)

                %warning('good time tyo stop')
                %keyboard


            otherwise
                error('No other tiploss method')
        end
    end

    % --------------------------------------------------------------------------------
    % --- Hub loss
    % --------------------------------------------------------------------------------
    Fhub=a*0+1;
    if(Algo.BEM.bHubLoss)
        %prandtl hub loss correction
        Fhub=2/pi*acos(exp(-nB/2*(r-rhub)./(rhub*sind(phi))));
    end
    F=Ftip.*Fhub;

    % --------------------------------------------------------------------------------
    % --- Step 3: Angle of attack
    % --------------------------------------------------------------------------------
    if(Algo.bAlphaHacked)
        alpha=interp1(Rotor.AlphaHacked(:,1),Rotor.AlphaHacked(:,2),r);
        phi=alpha+(twist+pitch);
    else
        alpha=phi-(twist+pitch); %[deg]
    end

    % --------------------------------------------------------------------------------
    % --- Step 4: Profile Data
    % --------------------------------------------------------------------------------
    [Cl Cd Cn Ct CnForAI CtForTI ] = fAeroCoeffWrap(0,alpha,phi,chord,Vrel_norm,Re,Fperf,WT,Algo); %Includes prescibed circulation


    % --------------------------------------------------------------------------------
    % --- Induction Factor in coordinate system
    % --------------------------------------------------------------------------------

    % --------------------------------------------------------------------------------
    % --- Step 5: Induction Coefficients
    % --------------------------------------------------------------------------------
    %a=(V0-Un)/V0;   %induction factor %abs for compatibility with unsteady code
    % Storing last values
    a_last=a;
    VrelVec=zeros(3,length(r));
    VrelVec(3,:)=Vrel_norm;
    nWVec=zeros(3,length(r));
    nWVec(3,:)=-a_last*V0;
    [ a aprime CTloc ] = fInductionCoefficients(a_last,VrelVec,Un,Ut,[0;0;V0],[0;0;V0],nWVec,omega,chord,F,Ftip,CnForAI,CtForTI,lambda_r,sigma,phi,Algo) ;

    % --------------------------------------------------------------------------------
    % --- Convergence Criteria
    % --------------------------------------------------------------------------------
    if (i>3 && max(abs(a-a_last))<aTol)
        break;
    end
    % tip Loss requirements
    Ftip_previous=Ftip;
    Gamma=0.5*Vrel_norm.*Cl.*chord;
    GammaNorm=Gamma/max(Gamma);
    Thrust = nB*sum(dr.*(Pn*cosd(cone)));    %Rotor shaft thrust at t in Newton
    CT=Thrust/(0.5*rho*V0^2*pi*R^2);
    if(i==1) 
        GammaBoundLast=Gamma;
    end
%
end %end iterative loop for one element
if(i==nbIt)
    fprintf('Maximum iterations reached\n');
else
     fprintf('Converged after %d iterations\n',i);
end


    % --------------------------------------------------------------------------------
    % --- Step 6: Aerodynamic Forces PER LENGTH
    % --------------------------------------------------------------------------------
    %     L=0.5*rho*Vrel_norm.^2*chord(e)*Cl;
    %     D=0.5*rho*Vrel_norm.^2*chord(e)*Cd;
    %     Pn(e) = L*cosd(phi) + D*sind(phi);   %load normal to the rotor plane
    %     Pt(e) = L*sind(phi) - D*cosd(phi);   %load tangential to the rotor plane
    Pn =   0.5*rho*Vrel_norm.^2*chord.*Cn;
    Pt =   0.5*rho*Vrel_norm.^2*chord.*Ct;
    



RES.Vrel=Vrel_norm;
RES.Re=Re;
RES.F=F;
RES.a=a;
RES.aprime=aprime;
RES.phi=phi;
RES.alpha=alpha;
RES.Cl=Cl;
RES.Cd=Cd;
RES.Cn=Cd;
RES.Ct=Ct;

RES.Pn=Pn;
RES.Pt=Pt;

RES.ThrLoc=dr.*Pn*cosd(cone);
RES.ThrLocLn=Pn*cosd(cone);
RES.CTloc=nB*RES.ThrLoc./(0.5*rho*VHubHeight^2*(2*pi*r.*cosd(cone).*dr)) ;

RES.TqLoc=dr.*r.*Pt*cosd(cone);
RES.TqLocLn=r.*Pt*cosd(cone);
RES.CQloc=nB.*RES.TqLoc./(0.5*rho*VHubHeight^2*(2*pi*r.*cosd(cone).*dr.*r*cosd(cone))) ;

%%%% Returning Aerodynamic Forces
if(isequal(WT.Sources.Format,'wtperf'))
    RES.Thrust = nB*sum(Rotor.dr.*(Pn*cosd(cone)));    %Rotor shaft thrust at t in Newton
    RES.Torque = nB*sum(Rotor.dr.*Pt.*(r*cosd(cone))); %Rotor shaft torque at t in Newton
else
    RES.Torque = nB*getTorqueFromBlade(r,Pt*cosd(cone),R);   %Rotor shaft torque at t in Newton
    RES.Thrust = nB*getThrustFromBlade(r,Pn*cosd(cone),R);   %Rotor shaft thrust at t in Newton
end
RES.Flap   = sum(dr.*(Pn*cosd(cone)).*(r-rhub));% approximate
RES.Edge   = sum(dr.*Pt.*(r*cosd(cone)).*(r-rhub));% approximate

RES.Power=omega*RES.Torque;
RES.CP=RES.Power/(0.5*rho*V0^3*WT.Rotor.SweptArea);
RES.CT=RES.Thrust/(0.5*rho*V0^2*pi*R^2);
RES.CQ=RES.Torque/(0.5*rho*V0^2*pi*R^3);
RES.Gamma=0.5*RES.Re.*RES.Cl*KinVisc*10^6;
if(isequal(Algo.BEM.TipLossMethod,'TipLossDB'))
    RES.Gamma_fit=Gamma_fit;
    RES.FParams=params;
else
RES.Ftip_aprime=Ftip_aprime;
RES.F_bar_a=F_bar_a;
RES.F_bar_aprime=F_bar_aprime;
end
RES.r=r;
RES.uia=V0*RES.a;
RES.uit=omega*r.*RES.aprime;
end

