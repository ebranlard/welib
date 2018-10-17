function [BEM] = fBEMsteady(WT,Sim,Wind,Algo)
cone=WT.Rotor.cone;

V0=Wind.V0(3)*cosd(cone);
VHubHeight=Wind.V0(3);

nB=WT.Rotor.nB;
ne=WT.Rotor.ne;
r=WT.Rotor.r;
dr=WT.Rotor.dr;
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
Pn=zeros(1,ne);
Pt=zeros(1,ne);
BEM.Vrel=zeros(1,ne);
BEM.Re=zeros(1,ne);
BEM.F=zeros(1,ne);
BEM.a=zeros(1,ne);
BEM.aprime=zeros(1,ne);
BEM.phi=zeros(1,ne);
BEM.alpha=zeros(1,ne);
BEM.Cl=zeros(1,ne);
BEM.Cd=zeros(1,ne);
BEM.Cn=zeros(1,ne);
BEM.Ct=zeros(1,ne);
BEM.CT=zeros(1,ne);


if(isequal(Algo.TipLossMethod,'GoldsteinSimple'))
    l_bar=1/lambda;
    w=1;
    Vr=linspace(0,R,length(r));
    FGoldstein=fTipLossGoldsteinOkulov( l_bar,w,R,nB,Vr,1 );
    figure
    plot(Vr,FGoldstein)
end

for e=ne:-1:1 % !!!!!!!!!!!!!!!!!! backward for Xu and Sankar correction
    %initialization
    a=0.3;
    aprime=0.01;
    
    Ftip_previous=1;
    %BEM loop
    for i=1:nbIt
        %%% Step 2 : flow angle
        phi=atan( (1-a) /((1+aprime)*lambda_r(e)) )*180/pi;  %[deg]
        if(imag(phi)~=0)
            disp(['Algorithm did not converge :',sprintf('lambda=%.2f beta=%.2f V0=%.2f r=%.2f - it %d',lambda,twist(e),V0,r(e),i)])
            phi=0;
            break
        end
        
        Ftip=1;
        Fhub=1;
        Fshen=1;
        
        if(Algo.TipLoss)
            % tip loss correction
            if(isequal(Algo.TipLossMethod,'Glauert'))
                if(abs(sind(phi))>0.01) %singularity if sin phi close to zero
                    Ftip=2/pi*acos(exp(-nB/2*(R-r(e))/(r(e)*sind(phi))));
                end
            elseif(isequal(Algo.TipLossMethod,'Prandtl'))
                Ftip=2/pi*acos(exp(-nB/2*(1-lambda_r(e)/lambda) *sqrt(1+lambda^2 )));
            elseif(isequal(Algo.TipLossMethod,'GoldsteinSimple'))
                Ftip=interp1(Vr,FGoldstein,r(e));
                Ftip=min(Ftip,1);
            elseif(isequal(Algo.TipLossMethod,'XuSankar'))
                if(abs(sind(phi))>0.01) %singularity if sin phi close to zero
                    if(r(e)/R>=0.7)
                        FGl=2/pi*acos(exp(-nB/2*(R-r(e))/(r(e)*sind(phi))));
                        Ftip=0.5*(FGl^0.86+0.5);
                    else
                        i07=whichvalue(r,0.7*R);
                        FGl=2/pi*acos(exp(-nB/2*(R-0.7*R)/(0.7*R*sind(BEM.phi(i07)))));
                        Ftip=1-r(e)/R*(1-FGl)/0.7 ;
                    end
                end
            elseif(isequal(Algo.TipLossMethod,'Lindenburg'))
                Ftip=2/pi*acos(exp(-nB/2*(1-lambda_r(e)/lambda)*sqrt(1+lambda_r(e)^2*   ( (1+2*sqrt(Ftip_previous)*aprime/2)/(1-sqrt(Ftip_previous)*a/2) )^2    )));
            elseif(isequal(Algo.TipLossMethod,'Shen'))
                if(abs(sind(phi))>0.01) %singularity if sin phi close to zero
                    Ftip=2/pi*acos(exp(-nB/2*(R-r(e))/(r(e)*sind(phi))));
                    g=exp(-0.125*(nB*lambda-21))+0.1;
                    Fshen  =2/pi*acos(exp(-nB/2*(R-r(e))/(r(e)*sind(phi)) *g));
                end
            elseif(isequal(Algo.TipLossMethod,'Manu'))
                
            end
        end
        if(Algo.HubLoss)
            %prandtl hub loss correction
            Fhub=2/pi*acos(exp(-nB/2*(r(e)-rhub)/(rhub*sind(phi))));
        end
        F=Ftip*Fhub;
        %%% Step 3 : angle of attack
        if(Algo.AlphaHacked)
            alpha=interp1(Rotor.AlphaHacked(:,1),Rotor.AlphaHacked(:,2),r(e));
            phi=alpha+(twist(e)+pitch);
        else
            alpha=phi-(twist(e)+pitch); %[deg]
        end
        %%% Step 4 : profiles data
        % Relative wind estimate
        Vrel=V0*(1-a)/sind(phi);
        Re=Vrel*chord(e)/KinVisc/10^6;
        if(bCl2piAlpha)
            Cl=2*pi*sind(alpha)*(1-exp(-((Rotor.r(e))-Rotor.rhub)/(Rotor.R-Rotor.rhub)*1/0.1 ));
            %;*acos(exp(-((Rotor.r(e)-Rotor.rhub)/(Rotor.R-Rotor.rhub))));
            Cd=Cl/(Algo.ClOverCd);%/((Rotor.r(e)-Rotor.rhub)/(Rotor.R-Rotor.rhub))^(1/20);
        else
            if(~isequal(Format,'flex'))
                ClCdCm= fAeroCoeff(alpha,WT.Profiles,WT.Rotor.ProfileSet(:,e),WT.Rotor.thickness_rel(e),Re,bReInterp,bRough);
                Cl=ClCdCm(1);
                Cd=ClCdCm(2);
                
            else
                ee=WT.Rotor.ProfileSet(2,e);
                % Badly programmed, what if all the alphas are not the same,
                % then the use of a table is bad
                % static aerodynamic coefficients
                Cd= interp1(WT.Profiles.alpha(:,ee) , WT.Profiles.Cd(:,ee)  , alpha);
                Cl= interp1(WT.Profiles.alpha(:,ee) , WT.Profiles.Cl(:,ee)  , alpha);
            end
            
            
        end
        
        %%% Step 5 : Aerodynamic  and coefficients
        if(Algo.AIDrag)
            Cn=Cl*cosd(phi)+Cd*sind(phi);
        else
            Cn=Cl*cosd(phi);
        end
        Ct=Cl*sind(phi)-Cd*cosd(phi);
        
        Cn=Fshen*Cn;
        Ct=Fshen*Ct;
        
        %%% Induction factor
        %by default the next a is :
        a_last=a;
        %normal expression, the default one
        a=1/( (4*F*sind(phi)^2) /(sigma(e)*Cn)+1 );
        % Thrust coefficient from the momentum theory => alast
        CT=(1-a_last)^2*sigma(e)*Cn/((sind(phi))^2);
        
        if(Algo.TipLoss && isequal(Algo.TipLossMethod,'Shen'))
            Y1=4*Ftip*sind(phi)^2/(sigma(e)*Cn*Fshen);
            Y2=4*Ftip*sind(phi)*cosd(phi)/(sigma(e)*Ct*Fshen);
            a=(2+Y1-sqrt(4*Y1*(1-Ftip)+Y1^2) ) /(2*(1+Ftip*Y1));
        end
        
        
        if(isequal(Algo.correction,'Shen'))
            ac=1/3;
            if a>ac
                a=(CT/4*Ftip-Ftip*ac^2)/(1-2*ac*Ftip) ;
            end
            
            %%% Glauert Eaxct correction
        elseif(isequal(Algo.correction,'Glauert'))
            ac=0.3;
            if a>ac
                A=sigma(e)*Cn/sind(phi)^2;
                a=fzero(@(aa) -A+aa*(4*F+2*A)+aa.^2*(-5*F-A)+3*F*aa.^3    ,[0 1]);
            end
            %%% Glauert Exact correction
        elseif(isequal(Algo.correction,'GlauertExact'))
            ac=0.3;
            if a>ac
                A=sigma(e)*Cn/sind(phi)^2;
                asolutions=GlauertSolutions(F,A);
                a=asolutions(whichmin(abs(asolutions-a_last)));
            end
            %%% Glauert correction REQUIRES RELAXATION
        elseif(isequal(Algo.correction,'GlauertCT'))
            ac=0.3;
            if a>ac
                fg=0.25*(5-3*a);
                a=CT/(4*F*(1-fg*a));
            end
            %SperaExact correction
        elseif(isequal(Algo.correction,'Spera'))
            ac=0.34;
            if a>ac
                K=4*F*(sind(phi))^2/(sigma(e)*Cn);
                a=0.5*(2+K*(1-2*ac)-sqrt( (K*(1-2*ac)+2 )^2 + 4*(K*ac^2-1)    )  );
            end
            %Spera correction REQUIRES RELAXATION
        elseif(isequal(Algo.correction,'SperaCT'))
            ac=0.34;
            if a>ac
                fgs=ac/a*(2-ac/a);
                a=CT/(4*F*(1-fgs*a));
            end
            %WE handbook correction
        elseif(isequal(Algo.correction,'HandbookCT'))
            if CT>0.96    %
                a=1/F*(0.143 + sqrt( 0.0203-0.6427 *(0.889-CT ) ));
            end
            %Aerodyn correction, REQUIRES RELAXATION
        elseif(isequal(Algo.correction,'AeroDyn'))
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
            if(Algo.TipLoss && isequal(Algo.TipLossMethod,'Shen'))
                aprime=1/((1-a*Ftip)*Y2/ ((1-a)-1));  % !!!!!!!!!!!!!!!!!!!!!! doubts on this one
            else
                
                aprime=1/( (4*F*sind(phi)*cosd(phi)) /(sigma(e)*Ct)  -1 )   ;
                if(isequal(Algo.correction,'AeroDyn'))
                    SwlAng=1+4*a*F*(1-a)/lambda_r(e)^2;
                    aprime=0.5*(sqrt(SwlAng)-1);
                end
            end
        else
            aprime=0;
        end
        if(isnan(a) )
            disp('BEM is crashing')
            error('Beak')
        end
        
        %%% convergence criteria
        if (i>3 && abs(a-a_last)<aTol)
            break;
        end
        %if (i>3 && abs(alpha-alpha_last)<alphaCrit  && abs(alpha-alpha_last_last) < alphaCrit)
        %    break;
        %end
        %alpha_last_last=alpha_last;
        %alpha_last=alpha;
        
        Ftip_previous=Ftip;
        
    end %end iterative loop for one element
    if(i==nbIt)
        fprintf('Maximum iterations reached : lambda=%.2f beta=%.2f V0=%.2f r=%.2f\n',lambda,twist(e),V0,r(e));
    end
    
    
    %%% Step 5 : Aerodynamic forces PER LENGTH and coefficients
    Cn=Cl*cosd(phi)+Cd*sind(phi);
    Ct=Cl*sind(phi)-Cd*cosd(phi);
    
    Cn=Fshen*Cn;
    Ct=Fshen*Ct;
    
    L=0.5*rho*norm(Vrel).^2*chord(e)*Cl;
    D=0.5*rho*norm(Vrel).^2*chord(e)*Cd;
    
    Pn(e) = L*cosd(phi) + D*sind(phi);   %load normal to the rotor plane
    Pt(e) = L*sind(phi) - D*cosd(phi);   %load tangential to the rotor plane
    
    Pn(e) =   0.5*rho*norm(Vrel).^2*chord(e)*Cn;
    Pt(e) =   0.5*rho*norm(Vrel).^2*chord(e)*Ct;
    
    BEM.Vrel(e)=Vrel;
    BEM.Re(e)=Re;
    BEM.F(e)=F;
    BEM.a(e)=a_last;
    BEM.aprime(e)=aprime;
    BEM.phi(e)=phi;
    BEM.alpha(e)=alpha;
    BEM.Cl(e)=Cl;
    BEM.Cd(e)=Cd;
    BEM.Cn(e)=Cn;
    BEM.Ct(e)=Ct;
    BEM.CT(e)=CT;
end %loop on elements

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
BEM.r=r;
end

