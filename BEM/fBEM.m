function [RES WT]=fBEM(x,v,update, WT,Sim,Wind,Algo)
% x and v have minimum length of two 
% the first coordinate is the nacelle displacement 
% the second coordinate is the rotor rotational speed


nB=WT.Rotor.nB;
ne=WT.Rotor.ne;
r=WT.Rotor.r;
dr=WT.Rotor.dr;
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

pitch=WT.Controller.pitch;

%% Initializations
% normal vectors
n_thrust_in3=[0; 0; -1];       %vector along the thrust in system 3
%n_rotor_in4=[0; 0; 1];         %normal vector to the rotor plan in 4
n_rotor_in2=[0; 0; 1];         %normal vector to the rotor plan in 2
% induced velocities
W0=ones(3,ne,nB);            %induced velocity without yaw/tilt
W=zeros(3,ne,nB);            %induced velocity
W_qs=zeros(3,ne,nB);         %quasistatic induced velocities
W_int=zeros(3,ne,nB);        %intermediate induced velocity (used for filtering)

% Forces and Power
RES.BladeTorque=zeros(1,nB);
RES.BladeThrust=zeros(1,nB);
RES.BladeEdge=zeros(1,nB);
RES.BladeFlap=zeros(1,nB);


%aero force on element
Pn=zeros(ne,1);
Pt=zeros(ne,1);
RES.Pn=zeros(nB,ne);
RES.Pt=zeros(nB,ne);
RES.Un=zeros(nB,ne);
RES.Ut=zeros(nB,ne);
RES.Cl=zeros(nB,ne);
RES.Cd=zeros(nB,ne);
RES.Gamma=zeros(nB,ne);
RES.alpha=zeros(nB,ne);

%weight force on elements
if(Algo.bWeight==0)
     Pz=0;
     Py=0;
else
    Py=zeros(ne,1);
    Pz=zeros(ne,1);
end
%Generalized force on each blade
if(Algo.DOF~=0)
    RES.GF=zeros(nB,3);
end

%% Preliminary calculations


% yaw model initialization
khi=WT.Aero.last.chi;
psi0=0;

%% Elasticity dependence
% Angles between blades updated depending on Dynamic degree of freedom 2
Vpsi0=mod(0:(360/nB):(360/nB)*(nB-1),360);
Vpsi= mod(Vpsi0 + x(2)*180/pi,360); % [deg]
omega=v(2);





% Transformation matrices
[a12 a23 a34]=getTransfoMatrices(WT.Controller.yaw, WT.Nacelle.tilt, 0, Rotor.cone);
% Shaft vector coordinates
rs_in1=a12'*WT.Shaft.rs_in2;

%%% Loop on blades
for idB=1:nB
    psi=Vpsi(idB);  %azimutal position of the blade
    % Transformation matrix
    a23=[cosd(psi) sind(psi) 0;  -sind(psi) cosd(psi) 0;  0 0 1];
    % transformation of coordinates
    a41=a12'*a23'*a34';
    a31=a12'*a23';
    a14=a41';
    a13=a31';
    % loop on elements
    for e=1:ne
        % --------------------------------------------------------------------------------
        % --- Step 0: Relative wind
        % --------------------------------------------------------------------------------
        rb_in4=[r(e); 0; 0];
        rb_in1=a41*rb_in4;
        r_position=WT.Tower.rt_in1+rs_in1+rb_in1;
        % Incoming wind
        V0_in1=getPointIncomingWindLegacy02(r_position,psi,WT,Wind,Algo);
        V0_in4=a14*V0_in1;
        V0_in3=a13*V0_in1;
        % Neglecting component on x ?????
        V0_in4=[0 ; V0_in4(2) ; V0_in4(3) ];
        V0_in3=[0 ; V0_in3(2) ; V0_in3(3) ];
        % Velocity seen by the blade
        Vb_in4=[0;  -omega*r(e)*cosd(Rotor.cone); 0]; %blade speed
        % Relative speed speed
        if(Algo.DOF>3)
            % Velocity of the blade from elasticity
            ee=Rotor.ee(e);
            Velast_in4(2)=v(idB*3+1)*Rotor.Blade.eigen1f(ee,1)+v(idB*3+2)*Rotor.Blade.eigen1e(ee,2)+v(idB*3+3)*Rotor.Blade.eigen2f(ee,2);
            Velast_in4(3)=v(idB*3+1)*Rotor.Blade.eigen1f(ee,2)+v(idB*3+2)*Rotor.Blade.eigen1e(ee,3)+v(idB*3+3)*Rotor.Blade.eigen2f(ee,3);
            Velasicity_in4=a14*[0;0; v(1)];
            Vrel_in4=V0_in4+WT.Aero.last.W(:,idB,e)+Vb_in4 - Velasicity_in4-Velast_in4;   
        else
            Velasicity_in4=a14*[0;0; 0*-v(1)];
            Vrel_in4=V0_in4+WT.Aero.last.W(:,e,idB)+Vb_in4 + Velasicity_in4;    %relative speed
        end
        lambda_r=omega*r(e)*cosd(cone)/norm(V0_in3);       %!!!!!!!!!!!!!!!!!!!! cone
        % --------------------------------------------------------------------------------
        % --- Step 1: Wind Components
        % --------------------------------------------------------------------------------
        % pas sur, it might be in 3...
        Ut = -Vrel_in4(2) ; 
        Un = Vrel_in4(3)  ; 
        Vrel_norm = sqrt(Un.^2+Ut.^2); %more stable than  Vrel=V0*(1-a)/sind(phi);
        Re=Vrel_norm*chord(e)/KinVisc/10^6; % Reynolds number in Millions
        RES.Un(idB,e)=Un;
        RES.Ut(idB,e)=Ut+Vb_in4(2) ;

        % --------------------------------------------------------------------------------
        % --- Step 2: Flow Angle
        % --------------------------------------------------------------------------------
        phi = atan2(Un,Ut)*180/pi; %(this returns a positive flow angle) [deg]
        if(imag(phi)~=0)
            disp(['Algorithm did not converge :',sprintf('lambda=%.2f beta=%.2f V0=%.2f r=%.2f - it %d',lambda,twist(e),V0,r(e),i)])
            phi=0;
            break
        end

        % --------------------------------------------------------------------------------
        % --- Tip loss
        % --------------------------------------------------------------------------------
        Ftip=1;
        Fperf=1;
        if(Algo.BEM.bTipLoss)
            if(sind(phi)>0.01)
                %prandtl tip correction
                Ftip = (2/pi)*(acos(exp(-(nB*(R-r(e)))/(2*r(e)*sind(phi)))));
            end
        end
        % --------------------------------------------------------------------------------
        % --- Hub loss
        % --------------------------------------------------------------------------------
        Fhub=1;
        if(Algo.BEM.bHubLoss)
            %prandtl hub loss correction
            Fhub=2/pi*acos(exp(-nB/2*(r(e)-rhub)/(rhub*sind(phi))));
        end
        F=Ftip.*Fhub;

        % --------------------------------------------------------------------------------
        % --- Step 3: Angle of attack
        % --------------------------------------------------------------------------------
        if(Algo.bAlphaHacked)
            alpha=interp1(Rotor.AlphaHacked(:,1),Rotor.AlphaHacked(:,2),r(e));
            phi=alpha+(twist(e)+pitch);
        else
            alpha=phi-(twist(e)+pitch); %[deg]
        end
        
        % --------------------------------------------------------------------------------
        % --- Step 4: Profile Data
        % --------------------------------------------------------------------------------
        [Cl Cd Cn Ct CnForAI CtForTI fs] = fAeroCoeffWrap(e,alpha,phi,Rotor.chord,norm(Vrel_in4),Re,1,WT,Algo); %Includes prescibed circulation


        L=0.5*rho*norm(Vrel_in4).^2*Rotor.chord(e)*Cl;
        
        % --------------------------------------------------------------------------------
        % --- Weight
        % --------------------------------------------------------------------------------
        if(Algo.bWeight)
            P_in1=[-Environment.g*Rotor.Blade.Mass(e);0;0];
            P_in4=a41'*P_in1;
            Py(e)=P_in4(2);
            Pz(e)=P_in4(3);
        end


        % --------------------------------------------------------------------------------
        % --- Induction Factor in coordinate system
        % --------------------------------------------------------------------------------
        Wn_in4=[0 ; 0 ; WT.Aero.last.W(3,e,idB) ];
        Wn_in3=a34'*Wn_in4;
        nnW_in3=n_thrust_in3.*(n_thrust_in3.*Wn_in3);
        nnW_in4=Wn_in4;
        V_prime_induction_in3=V0_in3+nnW_in3;
        sign=1;
        if(V_prime_induction_in3(3)<0)
            sign=-1;
        end
        a=(norm(V0_in3)-sign*norm(V_prime_induction_in3))/norm(V0_in3);   %induction factor
        if(isnan(a) || isnan(phi))            
            disp('RES is crashing')
            keyboard
            Algo.breakdown=1;
            break;
        end

        % --------------------------------------------------------------------------------
        % --- Step 5: Induction Coefficients
        % --------------------------------------------------------------------------------
        %a=(V0-Un)/V0;   %induction factor %abs for compatibility with unsteady code
        % Storing last values
        a_last=a;

        [ a aprime CT] = fInductionCoefficients(a_last,Vrel_in4,Un,Ut,V0_in3,V0_in4,nnW_in4,omega,chord(e),F,Ftip,CnForAI,CtForTI,lambda_r,sigma(e),phi,Algo) ;
        
        % --------------------------------------------------------------------------------
        % --- Dynamic wake model
        % --------------------------------------------------------------------------------
        W_y_qs = - omega*r(e)*aprime;  %!!!!!!!!1 and here the cone disappered...
        W_z_qs = -norm(V0_in3)*a;
        W_qs(:,e,idB)=[0 ; W_y_qs; W_z_qs];
        if(sum(isnan(squeeze(W_qs(:,e,idB))))>0)
            warning('here crash');
            keyboard
        end
        if(Algo.BEM.bDynaWake)
%             %%% Glauert correction
%             if RES.a<=0.2
%                 fg=1;
%             elseif RES.a>0.2
%                 fg=0.2/RES.a*(2-0.2/RES.a);
%             end
%             W_z_qs=-nB*RES.L*cosd(phi)/(4*pi*rho*r(e)*RES.F*norm(V0_in4+fg*nnW_in4 ));
%             W_y_qs=-nB*RES.L*sind(phi)/(4*pi*rho*r(e)*RES.F*norm(V0_in4+fg*nnW_in4 ));
            % dynamic wake
            tau1=1.1/(1-1.3*min(a,0.5))*R/norm(V0_in4);     %time constant 1
            tau2=(0.39-0.26*(r(e)/R)^2)*tau1;       %time constant 2
            H=W_qs(:,e,idB)+0.6*tau1*(W_qs(:,e,idB)-WT.Aero.last.W_qs(:,e,idB))/Algo.dt;

            W_int(:,e,idB)=H+(WT.Aero.last.W_int(:,e,idB)-H)*exp(-Algo.dt/tau1);  %intermediate W
            W0(:,e,idB)=W_int(:,e,idB)+(WT.Aero.last.W0(:,e,idB)-W_int(:,e,idB))*exp(-Algo.dt/tau2);  %W without yaw/tilt
        else
             W0(:,e,idB)=W_qs(:,e,idB);
        end
        
        % --------------------------------------------------------------------------------
        % ---  Yaw model, repartition of the induced velocity
        % --------------------------------------------------------------------------------
        if(Algo.BEM.bYawModel)
            %%% Yaw model, Skew angle and psi0
            if(e==Rotor.e_ref_for_khi && idB==1)
                %%% Determination of psi0
                r_hub=WT.Tower.rt_in1+rs_in1;
                % Incoming wind at hub
                V0_in1=getPointIncomingWindLegacy02(r_hub,psi,WT,Wind,Algo);
                V0_in2=a12*V0_in1;
                % psi0
                psi0=atan2(V0_in2(2),V0_in2(1))*180/pi;
                %%% Determination of skew angle
                % Averaging Wn on each blade
                meanWn_in4=[0;0;mean(W0(3,Rotor.e_ref_for_khi,:)) ];
                meanWn_in2=a34'*meanWn_in4;
                V_prime_for_khi_in2=V0_in2+meanWn_in2;
                khi=acosd(dot(n_rotor_in2,V_prime_for_khi_in2)/norm(V_prime_for_khi_in2)); %deg
            end
            W(:,e,idB)=W0(:,e,idB)*(1+r(e)/R *tand(khi/2)*cosd(Vpsi(idB)-psi0));
        else
            W(:,e,idB)=W0(:,e,idB);
        end

        % --------------------------------------------------------------------------------
        % --- Convergence Criteria
        % --------------------------------------------------------------------------------

        
        % --------------------------------------------------------------------------------
        % --- Step 6: Aerodynamic Forces PER LENGTH
        % --------------------------------------------------------------------------------
        L=0.5*rho*Vrel_norm.^2*chord(e)*Cl;
        D=0.5*rho*Vrel_norm.^2*chord(e)*Cd;
        Pn(e) = L*cosd(phi) + D*sind(phi);   %load normal to the rotor plane
        Pt(e) = L*sind(phi) - D*cosd(phi);   %load tangential to the rotor plane    
        Pn(e) =   0.5*rho*Vrel_norm.^2*chord(e)*Cn;
        Pt(e) =   0.5*rho*Vrel_norm.^2*chord(e)*Ct;
    %nasty
%     valpha(e)=alpha;
%     vphi(e)=phi;
%     va(e)=a;
%     vaprime(e)=aprime;
%     vCl(e)=RES.Cl;
%     vCd(e)=RES.Cd;
% %     va(e)=a;
% %     vaprime(e)=aprime;
%     vCn(e)=Cn;
%     vCnAI(e)=CnForAI;
%     vUn(e)=Un;
%     vUt(e)=Ut;
%     vF(e)=RES.F;
% 

        RES.Cl(idB,e)=Cl;
        RES.Cd(idB,e)=Cd;
        RES.alpha(idB,e)=alpha;
        RES.Gamma(idB,e)=0.5*Re.*Cl*KinVisc*10^6;
    end %loop on elements
	%%%
    if(Algo.DOF>3)
        RES.GF(idB,1)=trapz([Rotor.r;R],[(Pz+Pn).*Rotor.Blade.eigen1f(Rotor.ee,3)+(Py+Pt).*Rotor.Blade.eigen1f(Rotor.ee,2);0]);
        RES.GF(idB,2)=trapz([Rotor.r;R],[(Pz+Pn).*Rotor.Blade.eigen1e(Rotor.ee,3)+(Py+Pt).*Rotor.Blade.eigen1e(Rotor.ee,2);0]);
        RES.GF(idB,3)=trapz([Rotor.r;R],[(Pz+Pn).*Rotor.Blade.eigen2f(Rotor.ee,3)+(Py+Pt).*Rotor.Blade.eigen2f(Rotor.ee,2);0]);
    end

    %%%% Returning Aerodynamic Forces
    if(isequal(WT.Sources.Format,'wtperf'))
        RES.BladeThrust(idB) = sum(Rotor.dr.*(Pn*cosd(cone)));    %Rotor shaft thrust at t in Newton
        RES.BladeTorque(idB) = sum(Rotor.dr.*Pt.*(r*cosd(cone))); %Rotor shaft torque at t in Newton
    else
        RES.BladeTorque(idB) = getTorqueFromBlade(r,Pt*cosd(cone),R);   %Rotor shaft torque at t in Newton
        RES.BladeThrust(idB) = getThrustFromBlade(r,Pn*cosd(cone),R);   %Rotor shaft thrust at t in Newton
    end
    RES.BladeFlap(idB)   = sum(Rotor.dr.*(Pn'*cosd(cone)).*(r-rhub));% approximate
    RES.BladeEdge(idB)   = sum(Rotor.dr.*Pt'.*(r*cosd(cone)).*(r-rhub));% approximate
    
    RES.Pn(idB,:)=Pn;
    RES.Pt(idB,:)=Pt;


    %%% Torque momentum at hub
    %RES.BladeTorque(idB)=getTorqueFromBlade(r,(Py*0+Pt),R);
    %RES.BladeThrust(idB)=getThrustFromBlade(r,(Pz*0+Pn),R);
    %RES.BladeFlap(idB)=trapz([Rotor.r;R],[Rotor.r.*(Pz*0+Pn);0]);
    %RES.BladeEdge(idB)=trapz([Rotor.r;R],[Rotor.r.*(Py*0+Pt);0]);
end %loop on blades


%%%% Returning Aerodynamic Forces
RES.Torque = sum(RES.BladeTorque);   %Rotor shaft torque at t [N]
RES.Thrust = sum(RES.BladeThrust);   %Rotor shaft thrust at t [N]
RES.Flap   = sum(RES.BladeFlap);       %  [N]
RES.Edge   = sum(RES.BladeEdge);       %  [N]
RES.Power  = omega*RES.Torque;          %  [W]   

%   keyboard

%%% updating the aero values
if(update)
    WT.Aero.last.W_qs=W_qs;      % quasistatic induced velocity
    WT.Aero.last.W=W;            % induced velocity
    WT.Aero.last.W0=W0;          % induced velocity
    WT.Aero.last.W_int=W_int;    % intermediate induced velocity
    WT.Aero.last.fs=fs;
    WT.Aero.last.chi=khi;
    WT.Aero.Torque=RES.Torque;   %  [Nm]
    WT.Aero.Thrust=RES.Thrust;   %  [N]
    WT.Aero.Flap=RES.BladeFlap;  %  [N]
    WT.Aero.Edge=RES.BladeTorque;%  [Nm]
    WT.Aero.Power=RES.Power;     %  [W] 
end

