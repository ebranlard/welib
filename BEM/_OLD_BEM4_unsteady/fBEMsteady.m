function [BEM CP CT CQ] = fBEMsteady()
%Param,Profiles,r,x,chord,beta,a0,aprime0,nB,R,lambda,V0,rho)
    global Aero Rotor BEMParam Environment Iblade Profiles opt
    V0=Aero.Wind.V0(3);
    nB=Rotor.nB;
    ne=Rotor.ne;
    r=Rotor.r;
    R=Rotor.R;
    chord=Rotor.chord;
    beta=Rotor.beta;
    omega=Rotor.Omega;
    rho=Environment.rho;
    lambda_r=omega*r/V0;
    lambda=omega*R/V0;
    sigma=chord*nB./(2*pi*r);
    pitch=0;
    %algorithm internal paramter
    nbIt=BEMParam.nbIt;
    alphaCrit=BEMParam.alphaCrit;
    
    %aero force on element
    Pn=zeros(ne,1);
    Pt=zeros(ne,1);

    for e=1:ne
        %initialization
        a=0.2;
        alpha_last=0;
        alpha_last_last=0;
        aprime=0.01;
        
        %BEM loop
        for i=1:nbIt
          
            %%% Step 2 : flow angle
            phi=atan( (1-a) /((1+aprime)*lambda_r(e)) )*180/pi;  %[deg]
            F=1;% To avoid complex values
            if(imag(phi)~=0)
                    disp(['Algorithm did not converge :',sprintf('lambda=%.2f beta=%.2f V0=%.2f r=%.2f\n',lambda,beta(e),V0,r(e))])
                    break
            end
            if(sind(phi)>0.01)
                if(BEMParam.Tiploss)
                    %prandtl tip correction
                    f=nB/2*(R-r(e))/(r(e)*sind(phi));
                    F=2/pi*acos(exp(-f));
                end
            end

            %%% Step 3 : angle of attack
            alpha_=phi-(beta(e)+pitch); %[deg]
            %%% Step 4 : profiles data
%             alpha_data= Rotor.Profiles.alpha(Rotor.pe(e),:);
%             Cl = interp1(alpha_data,Rotor.Profiles.Cl(Rotor.pe(e),:),alpha);
%             Cd = interp1(alpha_data,Rotor.Profiles.Cd(Rotor.pe(e),:),alpha);
%             Rotor.ProfileSet(e+(length(Iblade(Iblade==0))),:)
%             opt.t_c
 
            % method interp Risoe
            ClCdCm_= fAeroCoeff(alpha_,Profiles,Rotor.ProfileSet(e+(length(Iblade(Iblade==0))),:),Rotor.opt_rel_thickness(e+(length(Iblade(Iblade==0)))));
            Cl=ClCdCm_(:,1);
            Cd=ClCdCm_(:,2);
            
            %%% Step 5 : Aerodynamic forces PER LENGTH and coefficients
            Cn=Cl*cosd(phi)+Cd*sind(phi);
            Ct=Cl*sind(phi)-Cd*cosd(phi);
            % Relative wind estimate
            Vrel=V0*(1-a)/sind(phi);
            L=0.5*rho*norm(Vrel).^2*chord(e)*Cl;
            D=0.5*rho*norm(Vrel).^2*chord(e)*Cd;
            Pn(e) = L*cosd(phi) + D*sind(phi);   %load normal to the rotor plane
            Pt(e) = L*sind(phi) - D*cosd(phi);   %load tangential to the rotor plane
                  
            %%% Induction factor
            %by default the next a is :
            a_last=a;
            %normal expression, the default one
            a=1/( (4*F*sind(phi)^2) /(sigma(e)*Cn)+1 ); 
            % Thrust coefficient from the momentum theory => alast
            CT=(1-a_last)^2*sigma(e)*Cn/((sind(phi))^2); 
            
            %%% Glauert Eaxct correction
            if(isequal(BEMParam.correction,'Glauert'))  
                ac=0.3;
                if a>ac
                   A=sigma(e)*Cn/sind(phi)^2;
                   a=fzero(@(aa) -A+aa*(4*F+2*A)+aa.^2*(-5*F-A)+3*F*aa.^3    ,[0 1]);
                end
            end
            
            %%% Glauert Eaxct correction
            if(isequal(BEMParam.correction,'GlauertExact'))  
                ac=0.3;
                if a>ac
                    A=sigma(e)*Cn/sind(phi)^2;
                    asolutions=GlauertSolutions(F,A);
                   a=asolutions(whichmin(abs(asolutions-a_last)));
                end
            end
           
            %%% Glauert correction REQUIRES RELAXATION
            if(isequal(BEMParam.correction,'GlauertCT'))  
                ac=0.3;
                if a>ac
                    fg=0.25*(5-3*a); 
                    a=CT/(4*F*(1-fg*a)); 
                end
      
            end
            %SperaExact correction
            if(isequal(BEMParam.correction,'Spera'))
                ac=0.34;
                if a>ac
                    K=4*F*(sind(phi))^2/(sigma(e)*Cn);
                    a=0.5*(2+K*(1-2*ac)-sqrt( (K*(1-2*ac)+2 )^2 + 4*(K*ac^2-1)    )  );
                end
            end
            %Spera correction REQUIRES RELAXATION 
            if(isequal(BEMParam.correction,'SperaCT'))
                ac=0.34;
                if a>ac
                    fgs=ac/a*(2-ac/a);
                    a=CT/(4*F*(1-fgs*a)); 
                end
            end
            %WE handbook correction
            if(isequal(BEMParam.correction,'HandbookCT'))
                if CT>0.96    %
                    a=1/F*(0.143 + sqrt( 0.0203-0.6427 *(0.889-CT ) ));
                end
            end             
                       
            %relaxation
            a=a*BEMParam.relaxation+(1-BEMParam.relaxation)*a_last;
            aprime=1/( (4*F*sind(phi)*cosd(phi)) /(sigma(e)*Ct)  -1 )   ;     
            

            %%% Storage
            if BEMParam.BigStorage
                BEM.phi(e,i)=phi;
                BEM.alpha_(e,i)=alpha_;
                BEM.a(e,i)=a_last;
                BEM.aprime(e,i)=aprime;
                BEM.Cd(e,i)=Cd;
                BEM.Cl(e,i)=Cl;
                BEM.Cn(e,i)=Cn;
                BEM.Ct(e,i)=Ct;
                BEM.CT(e,i)=CT;
                BEM.Vrel(e,i)=Vrel;
                BEM.Pn(e,i)=Pn(e);
                BEM.Pt(e,i)=Pt(e);
            else
                BEM.phi(e)=phi;
                BEM.alpha_(e)=alpha_;
                BEM.a(e)=a_last;
                BEM.aprime(e)=aprime;
                BEM.Cd(e)=Cd;
                BEM.Cl(e)=Cl;
                BEM.Cn(e)=Cn;
                BEM.Ct(e)=Ct;
                BEM.CT(e)=CT;
                BEM.Vrel(e)=Vrel;
                BEM.Pn(e)=Pn(e);
                BEM.Pt(e)=Pt(e);
            end
            %%% convergence criteria
            if (i>3 && abs(alpha_-alpha_last)<alphaCrit  && abs(alpha_-alpha_last_last) < alphaCrit)
                break;
            end
            
            alpha_last_last=alpha_last;
            alpha_last=alpha_;
        end %end iterative loop for one element
        if(i==nbIt)
             disp(['Algorithm did not converge :',sprintf('lambda=%.2f beta=%.2f V0=%.2f r=%.2f\n',lambda,beta(e),V0,r(e))])
        end
    end %loop on elements

%%%% Returning Aerodynamic Forces
BEM.Torque = nB*getTorqueFromBlade(r,Pt,R);   %Rotor shaft torque at t in Newton
BEM.Thrust = nB*getThrustFromBlade(r,Pn,R);   %Rotor shaft thrust at t in Newton
% BEM.Flap = nB;
BEM.Power=omega*BEM.Torque;
BEM.CP=BEM.Power/(0.5*rho*V0^3*pi*R^2);
BEM.CT=BEM.Thrust/(0.5*rho*V0^2*pi*R^2);
BEM.CQ=BEM.Torque/(0.5*rho*V0^2*pi*R^3);
end

