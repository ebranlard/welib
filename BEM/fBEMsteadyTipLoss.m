function [RES] = fBEMsteadyTipLoss(WT,Sim,Wind,Algo.BEM)
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
cone=WT.Rotor.cone;
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
a=RES.a;           %<-------------- Using Converged solution from previous RES
aprime=RES.aprime; %<-------------- Using Converged solution from previous RES

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

%% preparing tip-loss correction
x=(r-rhub)/(R-rhub);
Gamma=0.5*RES.Vrel_norm.*RES.Cl.*chord;
Gamma=Gamma./max(Gamma);
global PATH;
load([PATH.TIPLOSSDB 'TipLDBGamma']);
load([PATH.TIPLOSSDB 'TipLDBF']);
nGamma=size(TipLDBGamma,2)-4;
rGamma=single(cos(linspace(1,0,nGamma)*pi/2));
nF=size(TipLDBF,2)-6;
rF=cos(linspace(1,0,nF)*pi/2);
% load([MAINPATH 'data/TipLossDB/TipLDBParams']);
%         LAMBDAS=unique(TipLDBF(:,5));
%         CTS=unique(TipLDBF(:,6));
%         lambdaDB=LAMBDAS(whichvalue(LAMBDAS,lambda));
%         CTDB=CTS(whichvalue(CTS,RES.CT));
%         TipLDBF=TipLDBF(mfind(int16(round(TipLDBF(:,5:6)*1000)), int16(round([lambdaDB CTDB]*1000))),:);







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
        if(i>2 && vSSE2(jselect)*(1.01)>sum(x(ii:end).*(Gamma_fit(ii:end) - GammaNorm(ii:end)).^2)) % there is a trick here because the second condition is not evaluated, and hence the variables are initialized afterwards...
            %         fprintf('%d keep same Ftip\n',i);
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
            fprintf('%f %f %f %f\n',params);
            
            
            vSSE2(jselect);
            
            Ifound=mfind(int16(round(TipLDBF(:,1:4)*1000)),int16(round(params*1000)));
            
            LAMBDAS=unique(TipLDBF(Ifound,5));
            CTS=unique(TipLDBF(Ifound,6));
            lambdaDB=LAMBDAS(whichvalue(LAMBDAS,lambda));
            CTDB=CTS(whichvalue(CTS,RES.CT));
            Ftip=TipLDBF(mfind(int16(round(TipLDBF(:,1:6)*1000)), int16(round([params lambdaDB CTDB]*1000))),7:end);
            TipLDBF(Ifound,1:6);
            [lambda RES.CT];
            [params lambdaDB CTDB];
            
            %         if(length(Ifound)~=1)
            %             fprintf('not found')
            %             params
            %             Ifound
            %             error();
            %         end
            
            Ftip=interp1(rF,Ftip,x);    
            %         figure
            %         hold all
            %         plot(rGamma,TipLDBGamma(jselect,5:nGamma+4),'k+');
            %         plot(x,Gamma,'k');
            %     plot(x,Fitted_Curveb,'bo');
            % 
            %          figure
            %          plot(r,Ftip)
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
RES.Gamma_fit=Gamma_fit;
RES.FParams=params;
RES.r=r;
end

