%% Initializations
nt=length(Vtime);
% normal vectors
n_thrust_in3=[0; 0; -1];       %vector along the thrust in system 3
n_rotor_in4=[0; 0; 1];         %normal vector to the rotor plan in 4
n_rotor_in2=[0; 0; 1];         %normal vector to the rotor plan in 2
% induced velocities
W0=ones(3,nB,ne);            %induced velocity without yaw/tilt
W=zeros(3,nB,ne);            %induced velocity
W_qs=zeros(3,nB,ne);        %quasistatic induced velocities
W_int=zeros(3,nB,ne);        %intermediate induced velocity (used for filtering)


W_previous=ones(3,nB,ne).*w_guess;       %temporary induced velocity
W0_previous=ones(3,nB,ne).*w_guess;       %temporary induced velocity
W_qs_previous=ones(3,nB,ne).*w_guess; %temporary quasistatic induced velocity
W_int_previous=ones(3,nB,ne).*w_guess; %temporary intermediate induced velocity

% Forces and Power 
BladeTorque=zeros(nt,nB);
BladeThurst=zeros(nt,nB);
Torque=zeros(1,nt);
Thrust=zeros(1,nt);
Power=zeros(1,nt);

if(exist('Wguess'))
    disp('Using equilibrium')
    W_previous=Wguess;
    W_qs_previous=Wguess;
    W_int_previous=Wguess;    
end
if(exist('W0guess'))
    disp('Using equilibrium')
    W0_previous=W0guess;
end


Pn=zeros(ne,1);
Pt=zeros(ne,1);

if(BigStorage)
    Mpsi=zeros(nt,nB);
    MW=zeros(nt,3,nB,ne);
    MW0=zeros(nt,3,nB,ne);
    MPn=zeros(nt,nB,ne);
    MKhi=zeros(nt,nB);
end
%% Preliminary calculations
% Angles
Angles.pitch=Vpitch_of_t(1);
Angles.yaw=Vyaw_of_t(1);
Angles.tilt=tilt;
Angles.cone=cone;
Angles.psi=0;

% yaw model initialization
khi=max([Angles.yaw; Angles.tilt]);
psi0=0;
%Angles between blades
Vpsi=mod([0:(360/nB):(360/nB)*(nB-1)],360);


%% Time loop
for idT=1:length(Vtime)
    t=Vtime(idT);
    if(abs(mod(t,2*pi/omega))<=dt/2)
        disp(['Turn ' num2str(round(t/(2*pi/omega))) ' - t=' num2str(t)])
    end
    
    
    %%% Time dependent parameters
    Angles.pitch=Vpitch_of_t(mod(idT,length(Vpitch_of_t))+1);
    Angles.yaw=Vyaw_of_t(mod(idT,length(Vyaw_of_t))+1);
    VelocityParams.V0=VV0_of_t(mod(idT,length(VV0_of_t(:,1)))+1,:);
    
    % Transformation matrices
    [a12 a23 a34]=getTransfoMatrices(Angles.yaw,Angles.tilt,Angles.psi,Angles.cone);
    % Shaft vector coordinates
    rs_in1=a12'*rs_in2;
    
    %%% Loop on blades
    for idB=1:nB
        Angles.psi=Vpsi(idB);  %azimutal position of the blade
        % Transformation matrix
        a23=[cosd(Angles.psi) sind(Angles.psi) 0;  -sind(Angles.psi) cosd(Angles.psi) 0;  0 0 1];
        % transformation of coordinates
        a41=a12'*a23'*a34';
        a31=a12'*a23';
        
        % loop on elements
        for e=Velements
            %%% Step1 : Relative wind estimate
            rb_in4=[r(e); 0; 0];
            rb_in1=a41*rb_in4;
            r_position(e,:)=rt_in1+rs_in1+rb_in1;
            % Incoming wind
            V0_in1=getPointIncomingWindOld(r_position(e,:), Tower, Angles, Model,VelocityParams);
            V0_in4=a41'*V0_in1';
            V0_in3=a31'*V0_in1';
            % Neglecting component on x ?????
             V0_in4=[0 ;V0_in4(2) ;V0_in4(3) ];
             V0_in3=[0 ;V0_in3(2) ;V0_in3(3) ];
            % Velocity seen by the blade
            Vb_in4=[0;  -omega*r(e)*cosd(Angles.cone); 0]; %blade speed
            Vrel_in4=V0_in4+W_previous(:,idB,e)+Vb_in4;    %relative speed
            
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
            alpha=phi-(beta(e)+Angles.pitch);  %[deg]
            %%% Step 4 : profiles data
            Cl_data = Profiles.Cl(e,:);
            Cd_data = Profiles.Cd(e,:);
            alpha_data= Profiles.alpha(e,:);
            BEM.Cl = interp1(alpha_data,Cl_data,alpha);
            BEM.Cd = interp1(alpha_data,Cd_data,alpha);
            %%% Step 5 : Aerodynamic forces
            BEM.L=0.5*rho*norm(Vrel_in4).^2*chord(e)*BEM.Cl;
            BEM.D=0.5*rho*norm(Vrel_in4).^2*chord(e)*BEM.Cd;
            Pn(e) = BEM.L*cosd(phi) + BEM.D*sind(phi);   %load normal to the rotor plane
            Pt(e) = BEM.L*sind(phi) - BEM.D*cosd(phi);   %load tangential to the rotor plane
            
            %%% Induction factor
            Wn_in4=[0 ; 0 ; W_previous(3,idB,e) ];
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
            tau2=(0.39-0.26*(r(e)/R)^2)*tau1;            %time constant 2
            H=W_qs(:,idB,e)+0.6*tau1*(W_qs(:,idB,e)-W_qs_previous(:,idB,e))/dt;
            W_int(:,idB,e)=H+(W_int_previous(:,idB,e)-H)*exp(-dt/tau1);  %intermediate W
            W0(:,idB,e)=W_int(:,idB,e)+(W0_previous(:,idB,e)-W_int(:,idB,e))*exp(-dt/tau2);  %W without yaw/tilt
            %%% Yaw model, Skew angle and psi0
            if(e==e_ref_for_khi)
                %%% Determination of psi0
                r_hub=rt_in1+rs_in1;
                % Incoming wind at hub
                V0_in1=getPointIncomingWindOld(r_hub, Tower, Angles, Model,VelocityParams);
                V0_in2=a12*V0_in1';
                % psi0
                psi0=atan2(V0_in2(2),V0_in2(1))*180/pi;
                %%% Determination of skew angle
                % Averaging Wn on each blade
                meanWn_in4=[0;0;mean(W0(3,:,e_ref_for_khi)) ];
                meanWn_in2=a34'*meanWn_in4;
                V_prime_for_khi_in2=V0_in2+meanWn_in2;
                khi=acosd(dot(n_rotor_in2,V_prime_for_khi_in2)/norm(V_prime_for_khi_in2));
            end
            %%% Yaw model, repartition of the induced velocity
            if(YawModel)
                W(:,idB,e)=W0(:,idB,e)*(1+r(e)/R *tand(khi/2)*cosd(Vpsi(idB)-psi0));
            else
                W(:,idB,e)=W0(:,idB,e);
            end
        end %loop on elements
        
        %%% Torque momentum at hub
        BladeTorque(idT,idB)=getTorqueFromBlade(r,Pt);
        BladeThrust(idT,idB)=getThrustFromBlade(r,Pn);        
        
        if(BigStorage) 
            MPn(idT,idB,:)=Pn;
            Mpsi(idT,idB)=Vpsi(idB);
            MW(idT,:,idB,:)=W(:,idB,:);
            MW0(idT,:,idB,:)=W0(:,idB,:);
        end
    end %loop on blades
    
    
    %%% Rotation of the blades
    Vpsi=mod(Vpsi+omega*dt*180/pi,360) ;
    
    %%%% temporary induced velocity matrix storage
    W_qs_previous=W_qs;         %temporary quasistatic induced velocity
    W_previous=W;                %temporary induced velocity
    W0_previous=W0;                %temporary induced velocity
    W_int_previous=W_int;        %temporary intermediate induced velocity
    
    %%%% Storing Total Aerodynamic Forces
    Torque(idT) = sum(BladeTorque(idT,:))/1000;   %Rotor shaft torque at t in kN
    Thrust(idT) = sum(BladeThrust(idT,:))/1000;   %Rotor shaft thrust at t in kN    
end % time loop
Power=omega*Torque;