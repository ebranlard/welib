function [BEM]=fBEM(x,v,update)
global Aero Rotor Angles Controller Nacelle Algo Shaft Environment Tower
nB=Rotor.nB;
ne=Rotor.ne;
r=Rotor.r;
R=Rotor.R;
rho=Environment.rho;
%% Initializations
% normal vectors
n_thrust_in3=[0; 0; -1];       %vector along the thrust in system 3
n_rotor_in4=[0; 0; 1];         %normal vector to the rotor plan in 4
n_rotor_in2=[0; 0; 1];         %normal vector to the rotor plan in 2
% induced velocities
W0=ones(3,nB,ne);            %induced velocity without yaw/tilt
W=zeros(3,nB,ne);            %induced velocity
W_qs=zeros(3,nB,ne);        %quasistatic induced velocities
W_int=zeros(3,nB,ne);        %intermediate induced velocity (used for filtering)

% Forces and Power
BEM.BladeTorque=zeros(1,nB);
BEM.BladeThurst=zeros(1,nB);
BEM.BladeEdge=zeros(1,nB);
BEM.BladeFlap=zeros(1,nB);


%aero force on element
Pn=zeros(ne,1);
Pt=zeros(ne,1);
%weight force on elements
Py=zeros(ne,1);
Pz=zeros(ne,1);
%Generalized force on each blade
BEM.GF=zeros(nB,3);



%% Preliminary calculations


% yaw model initialization
khi=max([Controller.yaw; Nacelle.tilt]);
psi0=0;

%% Elasticity dependence
% Angles between blades updated depending on Dynamic degree of freedom 2
Vpsi0=mod([0:(360/nB):(360/nB)*(nB-1)],360);
Vpsi= mod(Vpsi0 + x(2)*180/pi,360); % [deg]
omega=v(2);

% Transformation matrices
[a12 a23 a34]=getTransfoMatrices(Controller.yaw, Nacelle.tilt, 0, Rotor.cone);
% Shaft vector coordinates
rs_in1=a12'*Shaft.rs_in2;

%%% Loop on blades
for idB=1:nB
    psi=Vpsi(idB);  %azimutal position of the blade
    % Transformation matrix
    a23=[cosd(psi) sind(psi) 0;  -sind(psi) cosd(psi) 0;  0 0 1];
    % transformation of coordinates
    a41=a12'*a23'*a34';
    a31=a12'*a23';
    
    % loop on elements
    for e=1:ne
        %%% Step1 : Relative wind estimate
        rb_in4=[r(e); 0; 0];
        rb_in1=a41*rb_in4;
        r_position(e,:)=Tower.rt_in1+rs_in1+rb_in1;
        % Incoming wind
        V0_in1=getPointIncomingWind(r_position(e,:),psi);
        V0_in4=a41'*V0_in1;
        V0_in3=a31'*V0_in1;
        % Neglecting component on x ?????
        V0_in4=[0 ; V0_in4(2) ; V0_in4(3) ];
        V0_in3=[0 ; V0_in3(2) ; V0_in3(3) ];
        % Velocity seen by the blade
        Vb_in4=[0;  -omega*r(e)*cosd(Rotor.cone); 0]; %blade speed
        % Velocity of the blade from elasticity
        Velast_in4=[0;0;0];
        ee=Rotor.ee(e);
        Velast_in4(2)=v(idB*3+1)*Rotor.Blade.eigen1f(ee,1)+v(idB*3+2)*Rotor.Blade.eigen1e(ee,2)+v(idB*3+3)*Rotor.Blade.eigen2f(ee,2);
        Velast_in4(3)=v(idB*3+1)*Rotor.Blade.eigen1f(ee,2)+v(idB*3+2)*Rotor.Blade.eigen1e(ee,3)+v(idB*3+3)*Rotor.Blade.eigen2f(ee,3);
        Velasicity_in4=a41'*[0;0; v(1)];
        Vrel_in4=V0_in4+Aero.last.W(:,idB,e)+Vb_in4 - Velasicity_in4-Velast_in4;    %relative speed
        
        %%% Step 2 : flow angle
        phi=atan2(Vrel_in4(3),-Vrel_in4(2))*180/pi;
        if(sind(phi)<0.01)
            % To avoid complex values
            BEM.F=1;
        else
            %prandtl tip correction
            f = (3*(R-r(e)))/(2*r(e)*sind(phi));
            BEM.F = (2/pi)*(acos(exp(-f)));
        end
        %%% Step 3 : angle of attack
        alpha=phi-(Rotor.beta(e)+Controller.pitch);  %[deg]
        %%% Step 4 : profiles data
        BEM.Cd= interp1(Rotor.Profiles.alpha(e,:) , Rotor.Profiles.Cd(e,:)  , alpha);
        fs=0;
        if(Algo.dynastall)
            % dynamic stall
            % interpolation from data
            f_st=interp1(Rotor.Profiles.alpha(e,:) , Rotor.Profiles.f_st(e,:) , alpha);
            Clinv=interp1(Rotor.Profiles.alpha(e,:) , Rotor.Profiles.Cl_inv(e,:) , alpha);
            Clfs= interp1(Rotor.Profiles.alpha(e,:) , Rotor.Profiles.Cl_fs(e,:)  , alpha);
            % dynamic stall model
            tau=4 * Rotor.chord(e) / norm(Vrel_in4);            
            fs=f_st + ( Aero.last.fs-f_st )*exp(- Algo.dt / tau);
            BEM.Cl=fs*Clinv+(1-f_st)*Clfs;
            
        else
            % static aerodynamic coefficients
            BEM.Cl= interp1(Rotor.Profiles.alpha(e,:) , Rotor.Profiles.Cl(e,:)  , alpha);
        end        
        
        
        
        %%% Step 5 : Aerodynamic forces
        BEM.L=0.5*rho*norm(Vrel_in4).^2*Rotor.chord(e)*BEM.Cl;
        BEM.D=0.5*rho*norm(Vrel_in4).^2*Rotor.chord(e)*BEM.Cd;
        Pn(e) = BEM.L*cosd(phi) + BEM.D*sind(phi);   %load normal to the rotor plane
        Pt(e) = BEM.L*sind(phi) - BEM.D*cosd(phi);   %load tangential to the rotor plane
        
        %%% Weight
        P_in1=[-9.81*Rotor.Blade.Mass(e);0;0];
        P_in4=a41'*P_in1;
        Py(e)=P_in4(2);
        Pz(e)=P_in4(3);
        
        %%% Induction factor
        Wn_in4=[0 ; 0 ; Aero.last.W(3,idB,e) ];
        Wn_in3=a34'*Wn_in4;
        nnW_in3=n_thrust_in3.*(n_thrust_in3.*Wn_in3);
        nnW_in4=Wn_in4;
        V_prime_induction_in3=V0_in3+nnW_in3;
        BEM.a=(norm(V0_in3)-norm(V_prime_induction_in3))/norm(V0_in3);   %induction factor
        %%% Glauert correction
        if BEM.a<=0.2
            fg=1;
        elseif BEM.a>0.2
            fg=0.2/BEM.a*(2-0.2/BEM.a);
        end
        %%% Dynamic wake model
        % quastistatic induced velocity
        W_z_qs=-nB*BEM.L*cosd(phi)/(4*pi*rho*r(e)*BEM.F*norm(V0_in4+fg*nnW_in4 ));
        W_y_qs=-nB*BEM.L*sind(phi)/(4*pi*rho*r(e)*BEM.F*norm(V0_in4+fg*nnW_in4 ));
        W_qs(:,idB,e)=[0 ; W_y_qs; W_z_qs];
        % dynamic wake
        tau1=1.1/(1-1.3*min(BEM.a,0.5))*R/norm(V0_in4);      %time constant 1
        tau2=(0.39-0.26*(r(e)/R)^2)*tau1;           %time constant 2
        H=W_qs(:,idB,e)+0.6*tau1*(W_qs(:,idB,e)-Aero.last.W_qs(:,idB,e))/Algo.dt;
        W_int(:,idB,e)=H+(Aero.last.W_int(:,idB,e)-H)*exp(-Algo.dt/tau1);  %intermediate W
        W0(:,idB,e)=W_int(:,idB,e)+(Aero.last.W0(:,idB,e)-W_int(:,idB,e))*exp(-Algo.dt/tau2);  %W without yaw/tilt
        
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
            meanWn_in4=[0;0;mean(W0(3,:,Rotor.e_ref_for_khi)) ];
            meanWn_in2=a34'*meanWn_in4;
            V_prime_for_khi_in2=V0_in2+meanWn_in2;
            khi=acosd(dot(n_rotor_in2,V_prime_for_khi_in2)/norm(V_prime_for_khi_in2));
        end
        %%% Yaw model, repartition of the induced velocity
        if(Algo.YawModel)
            W(:,idB,e)=W0(:,idB,e)*(1+r(e)/R *tand(khi/2)*cosd(Vpsi(idB)-psi0));
        else
            W(:,idB,e)=W0(:,idB,e);
        end
    end %loop on elements
    if(Algo.Weight==0)
      Pz=0;
      Py=0;
    end
    BEM.GF(idB,1)=trapz([Rotor.r;R],[(Pz+Pn).*Rotor.Blade.eigen1f(Rotor.ee,3)+(Py+Pt).*Rotor.Blade.eigen1f(Rotor.ee,2);0]);
    BEM.GF(idB,2)=trapz([Rotor.r;R],[(Pz+Pn).*Rotor.Blade.eigen1e(Rotor.ee,3)+(Py+Pt).*Rotor.Blade.eigen1e(Rotor.ee,2);0]);
    BEM.GF(idB,3)=trapz([Rotor.r;R],[(Pz+Pn).*Rotor.Blade.eigen2f(Rotor.ee,3)+(Py+Pt).*Rotor.Blade.eigen2f(Rotor.ee,2);0]);
    
    %%% Torque momentum at hub
    BEM.BladeTorque(idB)=getTorqueFromBlade(r,(Py+Pt),R);
    BEM.BladeThrust(idB)=getThrustFromBlade(r,(Pz+Pn),R);
    %%%
    BEM.BladeFlap(idB)=trapz([Rotor.r;R],[Rotor.r.*(Pz+Pn);0]);
end %loop on blades


%%%% Returning Aerodynamic Forces
BEM.Torque = sum(BEM.BladeTorque);  %Rotor shaft torque at t in N
BEM.Thrust = sum(BEM.BladeThrust);   %Rotor shaft thrust at t in N
BEM.Flap = sum(BEM.BladeFlap);
BEM.Power=omega*BEM.Torque;

%%% updating the aero values
if(update)
    Aero.last.W_qs=W_qs;      % quasistatic induced velocity
    Aero.last.W=W;            % induced velocity
    Aero.last.W0=W0;          % induced velocity
    Aero.last.W_int=W_int;    % intermediate induced velocity
    Aero.last.fs=fs;
    Aero.Torque=BEM.Torque;
    Aero.Thrust=BEM.Thrust;
    Aero.Power=BEM.Power;
    Aero.Flap=BEM.BladeFlap;
    Aero.Edge=BEM.BladeTorque;
end

