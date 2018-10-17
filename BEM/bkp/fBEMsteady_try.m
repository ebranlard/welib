function [BEM] = fBEMsteady()
global Aero Rotor Algo Environment Profiles Controller Nacelle Shaft Tower
nB=Rotor.nB;
ne=Rotor.ne;
nS=Algo.NumSect;

r=Rotor.r;
rhub=Rotor.rhub;
R=Rotor.R;

cone=Rotor.cone;

chord=Rotor.chord;
beta=Rotor.beta;
omega=Rotor.Omega;

rho=Environment.rho;
KinVisc=Environment.KinVisc;

%lambda_r=omega*r*cosd(cone)/V0;       %!!! cone
%lambda=omega*R*cosd(cone)/V0;         %!!! cone
sigma=chord*nB./(2*pi*r*cosd(cone));  %!!! CONE

pitch=Controller.pitch;

BEM.Thrust = 0;
BEM.Torque = 0;
BEM.Flap   = 0;
BEM.Edge   = 0;


%algorithm internal paramter
nbIt=Algo.nbIt;
aTol=Algo.aTol;

% normal vectors
n_thrust_in3=[0; 0; -1];       %vector along the thrust in system 3
n_rotor_in2=[0; 0; 1];         %normal vector to the rotor plan in 2
% induced velocities
W=zeros(3,ne,nS);            %induced velocity
W_qs=zeros(3,ne,nS);         %quasistatic induced velocities


%aero force on element
Pn=zeros(ne,1);
Pt=zeros(ne,1);


% yaw model initialization
khi=max([Controller.yaw; Nacelle.tilt]);
psi0=0;

% Angles between blades updated depending on Dynamic degree of freedom 2
Vpsi=mod(0:(360/nS):(360/nS)*(nS-1),360);


% Transformation matrices
[a12 a23 a34]=getTransfoMatrices(Controller.yaw, Nacelle.tilt, 0, Rotor.cone);
% Shaft vector coordinates
rs_in1=a12'*Shaft.rs_in2;

%%% Loop on Sectors
for idS=1:nS
    psi=Vpsi(idS);  %azimutal position of the blade
    % Transformation matrix
    a23=[cosd(psi) sind(psi) 0;  -sind(psi) cosd(psi) 0;  0 0 1];
    % transformation of coordinates
    a41=a12'*a23'*a34';
    a31=a12'*a23';
    a14=a41';
    a13=a31';
    
    for e=1:ne
        %initialization
       
        %%% Step1 : Relative wind estimate
        rb_in4=[r(e); 0; 0];
        rb_in1=a41*rb_in4;
        r_position=Tower.rt_in1+rs_in1+rb_in1;
        % Incoming wind
        V0_in1=getPointIncomingWind(r_position,psi);
        V0_in4=a14*V0_in1;
        V0_in3=a13*V0_in1;
        % Neglecting component on x ?????
        V0_in4=[0 ; V0_in4(2) ; V0_in4(3) ];
        V0_in3=[0 ; V0_in3(2) ; V0_in3(3) ];
        % Velocity seen by the blade
        Vb_in4=[0;  -omega*r(e)*cosd(Rotor.cone); 0]; %blade speed !!!CONE
        
        %BEM loop
        for i=1:nbIt
            % Relative speed speed
            Vrel_in4=V0_in4+Aero.last.W(:,e,idS)+Vb_in4 + 0;    %relative speed
            
            
            
            %%% Step 2 : flow angle
            phi=atan2(Vrel_in4(3),-Vrel_in4(2))*180/pi;
            
            F=1;
            if(imag(phi)~=0)
                disp(['Algorithm did not converge :',sprintf('lambda=%.2f beta=%.2f V0=%.2f r=%.2f\n',omega,beta(e),V0,r(e))])
                break
            end
            Ftip=1;
            Fhub=1;
            if(Algo.TipLoss)
                %prandtl tip loss correction
                Ftip=2/pi*acos(exp(-nB/2*(R-r(e))/(r(e)*sind(phi))));
            end
            if(Algo.HubLoss)
                %prandtl hub loss correction
                Fhub=2/pi*acos(exp(-nB/2*(r(e)-rhub)/(rhub*sind(phi))));
            end
            F=Ftip*Fhub;
            %%% Step 3 : angle of attack
            alpha=phi-(beta(e)+pitch); %[deg]
            %%% Step 4 : profiles data
            % Relative wind estimate
            Re=norm(Vrel_in4)*chord(e)/KinVisc/10^6;
            if(~isequal(Algo.Format,'flex'))
                ClCdCm= fAeroCoeff(alpha,Profiles,Rotor.ProfileSet(:,e),Rotor.thickness_rel(e),Re);
                Cl=ClCdCm(1);
                Cd=ClCdCm(2);
            else
                ee=Rotor.ProfileSet(2,e);
                % Badly programmed, what if all the alphas are not the same,
                % then the use of a table is bad
                % static aerodynamic coefficients
                Cd= interp1(Profiles.alpha(:,ee) , Profiles.Cd(:,ee)  , alpha);
                Cl= interp1(Profiles.alpha(:,ee) , Profiles.Cl(:,ee)  , alpha);
            end
            
            %%% Step 5 : Aerodynamic forces PER LENGTH and coefficients
            Cn=Cl*cosd(phi)+Cd*sind(phi);
            Ct=Cl*sind(phi)-Cd*cosd(phi);
            
            if(Algo.AIDrag)
                Cn=Cl*cosd(phi)+Cd*sind(phi);
            else
                Cn=Cl*cosd(phi);
            end
            
        end
        
        L=0.5*rho*norm(Vrel_in4).^2*chord(e)*Cl;
        %%% Induction factor
        Wn_in4=[0 ; 0 ; Aero.last.W(3,e,idS) ];
        Wn_in3=a34'*Wn_in4;
        nnW_in3=n_thrust_in3.*(n_thrust_in3.*Wn_in3);
        nnW_in4=Wn_in4;
        V_prime_induction_in3=V0_in3+nnW_in3;
        a=(norm(V0_in3)-norm(V_prime_induction_in3))/norm(V0_in3);   %induction factor
        if(isnan(a) || isnan(phi))
            disp('BEM is crashing')
            Algo.breakdown=1;
            break;
        end
        
        
        %%% Induction factor
        %by default the next a is :
        a_last=a;
        %normal expression, the default one
        a=1/( (4*F*sind(phi)^2) /(sigma(e)*Cn)+1 );
        % Thrust coefficient from the momentum theory => alast
        CT=(1-a_last)^2*sigma(e)*Cn/((sind(phi))^2);
        
        if(isequal(Algo.correction,'Spera'))
            ac=0.34;
            if a>ac
                K=4*F*(sind(phi))^2/(sigma(e)*Cn);
                a=0.5*(2+K*(1-2*ac)-sqrt( (K*(1-2*ac)+2 )^2 + 4*(K*ac^2-1)    )  );
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
        
        
        
        %%% Dynamic wake model
        % quastistatic induced velocity
        W_z_qs=-nB*L*cosd(phi)/(4*pi*rho*r(e)*F*norm(V0_in4+nnW_in4 ));
        W_y_qs=-nB*L*sind(phi)/(4*pi*rho*r(e)*F*norm(V0_in4+nnW_in4 ));
        W_qs(:,e,idS)=[0 ; W_y_qs; W_z_qs];
        
        
        
        %%% Yaw model, repartition of the induced velocity
        if(Algo.YawModel)
            %%% Yaw model, Skew angle and psi0
            if(e==Rotor.e_ref_for_khi)
                %%% Determination of psi0
                r_hub=Tower.rt_in1+rs_in1;
                % Incoming wind at hub
                V0_in1=getPointIncomingWind(r_hub,psi);
                V0_in2=a12*V0_in1;
                % psi0
                psi0=atan2(V0_in2(2),V0_in2(1))*180/pi;
                %%% Determination of skew angle
                % Averaging Wn on each blade
                meanWn_in4=[0;0;mean(W_qs(3,Rotor.e_ref_for_khi,:)) ];
                meanWn_in2=a34'*meanWn_in4;
                V_prime_for_khi_in2=V0_in2+meanWn_in2;
                khi=acosd(dot(n_rotor_in2,V_prime_for_khi_in2)/norm(V_prime_for_khi_in2));
            end
            W(:,e,idS)=W_qs(:,e,idS)*(1+r(e)/R *tand(khi/2)*cosd(Vpsi(idS)-psi0));
        else
            W(:,e,idS)=W_qs(:,e,idS);
        end

        %%% Swirl
        if(Algo.Swirl)
            aprime=1/( (4*F*sind(phi)*cosd(phi)) /(sigma(e)*Ct)  -1 )   ;
            if(isequal(Algo.correction,'AeroDyn'))
                lambda_r=omega*r(e)/norm(V0_in4);
                SwlAng=1+4*a*F*(1-a)/lambda_r^2;
                aprime=0.5*(sqrt(SwlAng)-1);
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
        if(i==nbIt)
            disp(['Maximum iterations reached :',sprintf('lambda=%.2f beta=%.2f V0=%.2f r=%.2f\n',omega,beta(e),V0_in4,r(e))])
        end
        
        
        %%% Step 5 : Aerodynamic forces PER LENGTH and coefficients
        Cn=Cl*cosd(phi)+Cd*sind(phi);
        Ct=Cl*sind(phi)-Cd*cosd(phi);
        
        L=0.5*rho*norm(Vrel_in4).^2*chord(e)*Cl;
        D=0.5*rho*norm(Vrel_in4).^2*chord(e)*Cd;
        Pn(e) = L*cosd(phi) + D*sind(phi);   %load normal to the rotor plane
        Pt(e) = L*sind(phi) - D*cosd(phi);   %load tangential to the rotor plane
        
        BEM.Vrel(e)=norm(Vrel_in4);
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
    %%%% Returning Aerodynamic Forces
    if(isequal(Algo.Format,'wtperf'))
        BEM.Thrust = BEM.Thrust+nB*sum(Rotor.dr.*(Pn*cosd(cone)))/nS;    %Rotor shaft thrust at t in Newton
        BEM.Torque = BEM.Torque+nB*sum(Rotor.dr.*Pt.*(r*cosd(cone)))/nS; %Rotor shaft torque at t in Newton
        BEM.Flap   = BEM.Flap+sum(Rotor.dr.*(Pn*cosd(cone)).*(r-rhub))/nS;% approximate
        BEM.Edge   = BEM.Edge+sum(Rotor.dr.*Pt.*(r*cosd(cone)).*(r-rhub))/nS;% approximate
    else
        BEM.Torque = nB*getTorqueFromBlade(r,Pt,R);   %Rotor shaft torque at t in Newton
        BEM.Thrust = nB*getThrustFromBlade(r,Pn,R);   %Rotor shaft thrust at t in Newton
    end
end % loop on sectors




BEM.Power=omega*BEM.Torque;
% BEM.CP=BEM.Power/(0.5*rho*V0^3*Rotor.SweptArea);
% BEM.CT=BEM.Thrust/(0.5*rho*V0^2*pi*R^2);
% BEM.CQ=BEM.Torque/(0.5*rho*V0^2*pi*R^3);
end

