function [BEM] = fBEMsteadyTipLoss(WT,Sim,Wind,Algo,BEM)
global MAINPATH
addpath([MAINPATH 'code/CirculationFamilyCurves/']);

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
omega=Sim.Omega;

Format=WT.Sources.Format;

rho=Sim.rho;
KinVisc=Sim.KinVisc;

lambda_r=omega*r*cosd(cone)/V0;       %!!! cone
lambda=omega*R*cosd(cone)/V0;         %!!! cone
sigma=chord*nB./(2*pi*r*cosd(cone));  %!!! CONE

pitch=Sim.PITCH;
%algorithm internal paramter
nbIt=Algo.nbIt;
aTol=Algo.aTol;
bReInterp=Algo.ReInterp;
bRough=Algo.RoughProfiles;
bCl2piAlpha=Algo.Cl2piAlpha;


% initialize vectors
a=BEM.a;           %<-------------- Using Converged solution from previous BEM
aprime=BEM.aprime; %<-------------- Using Converged solution from previous BEM

Pn=zeros(1,ne);
Pt=zeros(1,ne);
Vrel=zeros(1,ne);
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
Gamma=0.5*BEM.Vrel.*BEM.Cl.*chord;
Gamma=Gamma./max(Gamma);
load([MAINPATH 'data/TipLossDB/TipLDBGamma']);
load([MAINPATH 'data/TipLossDB/TipLDBF']);
% load([MAINPATH 'data/TipLossDB/TipLDBParams']);
%         LAMBDAS=unique(TipLDBF(:,5));
%         CTS=unique(TipLDBF(:,6));
%         lambdaDB=LAMBDAS(whichvalue(LAMBDAS,lambda));
%         CTDB=CTS(whichvalue(CTS,BEM.CT));
%         TipLDBF=TipLDBF(mfind(int16(round(TipLDBF(:,5:6)*1000)), int16(round([lambdaDB CTDB]*1000))),:);



nGamma=size(TipLDBGamma,2)-4;
rGamma=single(cos(linspace(1,0,nGamma)*pi/2));
nF=size(TipLDBF,2)-6;
rF=cos(linspace(1,0,nF)*pi/2);
x=(r-rhub)/(R-rhub);



for i=1:nbIt % !!!!!!!!!!!!!!!!!!! in this order for my tip loss correction
    %%% Step 2 : flow angle
%     phi=atan( (1-a) ./((1+aprime).*lambda_r) )*180/pi;  %[deg]
    phi=atan2( (1-a) ,((1+aprime).*lambda_r) )*180/pi;   %[deg]
    IbadCV=find(imag(phi)~=0);% index of bad convergence
    %         if(imag(phi)~=0)
    %             disp(['Algorithm did not converge :',sprintf('lambda=%.2f beta=%.2f V0=%.2f r=%.2f - it %d',lambda,twist(e),V0,r(e),i)])
    %             phi=0;
    %             break
    %         end
    if(~isempty(IbadCV))
        fprintf('Setting %d bad phi zero',length(IbadCV));
        phi(IbadCV)=0;
    end
    %%% Tip Loss

    if(i>1 && vSSE2(jselect)*(1.01)>sum(x(ii:end).*(Gamma_fit(ii:end) - Gamma(ii:end)).^2))
%         fprintf('%d keep same Ftip\n',i);
    else
        Ftip=a*0+1;
        Fhub=a*0+1;
        vSSE2=zeros(1,size(TipLDBGamma,1));
        ii=whichvalue(x,min(x(Gamma==max(Gamma)),0.25));
        for j=1:size(TipLDBGamma,1);
            Gamma_fit=interp1(rGamma,TipLDBGamma(j,5:nGamma+4),x);
            x0=TipLDBGamma(j,1);
            ii=whichvalue(x,min(x0,0.25));
            vSSE2(j)=sum(x(ii:end).*(Gamma_fit(ii:end) - Gamma(ii:end)).^2);
        end
        jselect=whichmin(vSSE2);
        Gamma_fit=interp1(rGamma,TipLDBGamma(jselect,5:nGamma+4),x);
        params=TipLDBGamma(jselect,1:4);
        
        %hack
%         params=[0.4 0.95 1 0.3];
%         [sse Gamma_fit]=fFitGamma(params,x,Gamma);

        BEM.FParams=params;
        vSSE2(jselect);

        Ifound=mfind(int16(round(TipLDBF(:,1:4)*1000)),int16(round(params*1000)));
        
        LAMBDAS=unique(TipLDBF(Ifound,5));
        CTS=unique(TipLDBF(Ifound,6));
        lambdaDB=LAMBDAS(whichvalue(LAMBDAS,lambda));
        CTDB=CTS(whichvalue(CTS,BEM.CT));
        Ftip=TipLDBF(mfind(int16(round(TipLDBF(:,1:6)*1000)), int16(round([params lambdaDB CTDB]*1000))),7:end);
        TipLDBF(Ifound,1:6);
        [lambda BEM.CT];
        [params lambdaDB CTDB];
        
%         if(length(Ifound)~=1)
%             fprintf('not found')
%             params
%             Ifound
%             error();
%         end

%         Ftip=TipLDBF(Ifound,7:end);
        Ftip=interp1(rF,Ftip,x);    
%     figure
%     hold all
%     plot(rGamma,TipLDBGamma(jselect,5:nGamma+4),'k+');
%     plot(x,Gamma,'k');
%     plot(x,Fitted_Curveb,'bo');
% 
%     figure
%     plot(rF,Ftip)
    end
    
    if(Algo.HubLoss)
        %prandtl hub loss correction
        Fhub=2/pi*acos(exp(-nB/2*(r-rhub)./(rhub*sind(phi))));
    end
    F=Ftip.*Fhub;
    %%% Step 3 : angle of attack
    alpha=phi-(twist+pitch); %[deg]
    %%% Step 4 : profiles data
    % Relative wind estimate
    Vrel=V0*(1-a)./sind(phi);
    Re=Vrel.*chord/KinVisc/10^6;  % [Milllions]
    if(bCl2piAlpha)
        Cl=2*pi*sind(alpha).*(1-exp(-((Rotor.r)-Rotor.rhub)/(Rotor.R-Rotor.rhub)*1/0.1 ));
        %;*acos(exp(-((Rotor.r(e)-Rotor.rhub)/(Rotor.R-Rotor.rhub))));
        Cd=Cl/(Algo.ClOverCd);%/((Rotor.r(e)-Rotor.rhub)/(Rotor.R-Rotor.rhub))^(1/20);
    else
        if(~isequal(Format,'flex')) % bHawc, Hawc2 format with Pc files
            ClCdCm= fAeroCoeff(alpha,WT.Profiles,WT.Rotor.ProfileSet(:,:),WT.Rotor.thickness_rel(:),Re,bReInterp,bRough);
            Cl=ClCdCm(:,1)';
            Cd=ClCdCm(:,2)';
        else
            % to be done
            ee=WT.Rotor.ProfileSet(2,e);
            % Badly programmed, what if all the alphas are not the same,
            % then the use of a table is bad
            % static aerodynamic coefficients
            Cd= interp1(WT.Profiles.alpha(:,ee) , WT.Profiles.Cd(:,ee)  , alpha);
            Cl= interp1(WT.Profiles.alpha(:,ee) , WT.Profiles.Cl(:,ee)  , alpha);
        end
    end
    
    %%% Step 5 : Aerodynamic coefficients
    if(Algo.AIDrag)
        Cn=Cl.*cosd(phi)+Cd.*sind(phi);
    else
        Cn=Cl.*cosd(phi);
    end
    Ct=Cl.*sind(phi)-Cd.*cosd(phi);
    
    %%% Induction factor
    %by default the next a is :
    a_last=a;
    %normal expression, the default one
    a=1./( (4*F.*sind(phi).^2) ./(sigma.*Cn)+1 );
    % Thrust coefficient from the momentum theory => alast
    CT=(1-a_last).^2.*sigma.*Cn./((sind(phi)).^2);
    
    if(isequal(Algo.correction,'Shen'))
        ac=1/3;
        Iac=a>ac;
        a(Iac)=(CT(Iac)/4*Ftip(Iac)-Ftip(Iac)*ac^2)./(1-2*ac*Ftip(Iac)) ;
        %%% Glauert Eaxct correction
    elseif(isequal(Algo.correction,'Glauert'))
        ac=0.3;
        Iac=a>ac;
        error();
        A=sigma.*Cn./sind(phi)^2;
        a=fzero(@(aa) -A+aa*(4*F+2*A)+aa.^2*(-5*F-A)+3*F*aa.^3    ,[0 1]);
        %%% Glauert Exact correction
    elseif(isequal(Algo.correction,'GlauertExact'))
        ac=0.3;
        error();
        if a>ac
            A=sigma(e)*Cn/sind(phi)^2;
            asolutions=GlauertSolutions(F,A);
            a=asolutions(whichmin(abs(asolutions-a_last)));
        end
        %%% Glauert correction REQUIRES RELAXATION
    elseif(isequal(Algo.correction,'GlauertCT'))
        ac=0.3;
        Iac=a>ac;
        fg=0.25*(5-3*a(Iac));
        a(Iac)=CT(Iac)./(4*F(Iac).*(1-fg.*a(Iac)));
        %SperaExact correction
    elseif(isequal(Algo.correction,'Spera'))
        ac=0.34;
        Iac=a>ac;
        K=4*F(Iac).*(sind(phi(Iac))).^2./(sigma(Iac).*Cn(Iac));
        a(Iac)=0.5*(2+K*(1-2*ac)-sqrt((K*(1-2*ac)+2 ).^2 + 4*(K*ac^2-1)    )  );
        %Spera correction REQUIRES RELAXATION
    elseif(isequal(Algo.correction,'SperaCT'))
        ac=0.34;
        Iac=a>ac;
        fgs=ac./a(Iac)*(2-ac./a(Iac));
        a=CT(Iac)./(4*F(Iac).*(1-fgs.*a(Iac)));
        %WE handbook correction
    elseif(isequal(Algo.correction,'HandbookCT'))
        Ict=CT>0.96;
        a(Ic)=1./F(Ic).*(0.143 + sqrt( 0.0203-0.6427 *(0.889-CT(Ic) ) ));
        %Aerodyn correction, REQUIRES RELAXATION
    elseif(isequal(Algo.correction,'AeroDyn'))
        error();
        CT=min(max(-2.0,CT),2);
        if CT>0.96*F    %
            a=0.1432+sqrt(-0.55106+0.6427*CT/F);
        else
            a=0.5*(1-sqrt(1-CT/F));
        end
    end
    
    %relaxation
    a=a*Algo.relaxation+(1-Algo.relaxation)*a_last;
    
    %%% Swirl
    if(Algo.Swirl)
        aprime=1./( (4*F.*sind(phi).*cosd(phi)) ./(sigma.*Ct)  -1 )   ;
        if(isequal(Algo.correction,'AeroDyn'))
            SwlAng=1+4*a.*F.*(1-a)./lambda_r.^2;
            aprime=0.5*(sqrt(SwlAng)-1);
        end
    else
        aprime=a*0;
    end
    if(sum(isnan(a))>0 )
        disp('BEM is crashing')
        error('Beak')
    end

    %%% Step 5 : Aerodynamic forces PER LENGTH and coefficients
    Cn=Cl.*cosd(phi)+Cd.*sind(phi);
    Ct=Cl.*sind(phi)-Cd.*cosd(phi); %redundant   
%     L=0.5*rho*Vrel.^2*chord.*Cl;
%     D=0.5*rho*Vrel.^2*chord.*Cd;
%     
%     Pn = L.*cosd(phi) + D.*sind(phi);   %load normal to the rotor plane
%     Pt = L.*sind(phi) - D.*cosd(phi);   %load tangential to the rotor plane
    Pn =   0.5*rho*Vrel.^2.*chord.*Cn;  %load normal to the rotor plane
    Pt =   0.5*rho*Vrel.^2.*chord.*Ct;  %load tangential to the rotor plane
    
     %%% convergence criteria and preparation for next step
    if (i>3 && max(abs(a-a_last))<aTol)
        break;
    end    
    Ftip_previous=Ftip;
    Gamma=0.5*Vrel.*Cl.*chord;
    Gamma=Gamma/max(Gamma);
    %

end %end iterative loop for one element
if(i==nbIt)
    fprintf('Maximum iterations reached\n');
else
     fprintf('Converged after %d iterations reached\n',i);
end





BEM.Vrel=Vrel;
BEM.Re=Re;
BEM.F=F;
BEM.a=a;
BEM.aprime=aprime;
BEM.phi=phi;
BEM.alpha=alpha;
BEM.Cl=Cl;
BEM.Cd=Cd;
BEM.Cn=Cd;
BEM.Ct=Ct;

BEM.Pn=Pn;
BEM.Pt=Pt;

BEM.ThrLoc=dr.*Pn*cosd(cone);
BEM.ThrLocLn=Pn*cosd(cone);
BEM.CTLoc=nB*BEM.ThrLoc./(0.5*rho*VHubHeight^2*(2*pi*r.*cosd(cone).*dr)) ;

BEM.TqLoc=dr.*r.*Pt*cosd(cone);
BEM.TqLocLn=r.*Pt*cosd(cone);
BEM.CQLoc=nB.*BEM.TqLoc./(0.5*rho*VHubHeight^2*(2*pi*r.*cosd(cone).*dr.*r*cosd(cone))) ;

%%%% Returning Aerodynamic Forces
if(isequal(Format,'wtperf'))
    BEM.Thrust = nB*sum(Rotor.dr.*(Pn*cosd(cone)));    %Rotor shaft thrust at t in Newton
    BEM.Torque = nB*sum(Rotor.dr.*Pt.*(r*cosd(cone))); %Rotor shaft torque at t in Newton
else
    BEM.Torque = nB*getTorqueFromBlade(r,Pt,R);   %Rotor shaft torque at t in Newton
    BEM.Thrust = nB*getThrustFromBlade(r,Pn,R);   %Rotor shaft thrust at t in Newton
end
BEM.Thrust = nB*sum(dr.*(Pn*cosd(cone)));    %Rotor shaft thrust at t in Newton
BEM.Torque = nB*sum(dr.*Pt.*(r*cosd(cone))); %Rotor shaft torque at t in Newton
BEM.Flap   = sum(dr.*(Pn*cosd(cone)).*(r-rhub));% approximate
BEM.Edge   = sum(dr.*Pt.*(r*cosd(cone)).*(r-rhub));% approximate

BEM.Power=omega*BEM.Torque;
BEM.CP=BEM.Power/(0.5*rho*V0^3*WT.Rotor.SweptArea);
BEM.CT=BEM.Thrust/(0.5*rho*V0^2*pi*R^2);
BEM.CQ=BEM.Torque/(0.5*rho*V0^2*pi*R^3);
BEM.Gamma=0.5*BEM.Re.*BEM.Cl*KinVisc*10^6;
BEM.Gamma_fit=Gamma_fit;
BEM.r=r;
end

