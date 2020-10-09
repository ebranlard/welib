function [BEM CP CT CQ] = BEMfunction(Param,Profiles,r,dr,x,c,beta,a,aprime,B,R,lambda,V0,rho)
    %algorithm internal paramter
    nbIt=Param.nbIt;
    alphaCrit=Param.alphaCrit;
    
    %variable needed  
    ne=length(r);
    omega=lambda*V0/R;
    sigma=c*B./(2*pi*r);
    %% matrices that will store the paramters value at each iterational step
    BEM.a=zeros(ne,nbIt);
    BEM.aprime=zeros(ne,nbIt);
    BEM.alpha=zeros(ne,nbIt);
    BEM.phi=zeros(ne,nbIt);
    BEM.Cl=zeros(ne,nbIt);
    BEM.Cd=zeros(ne,nbIt);
    BEM.Cn=zeros(ne,nbIt);
    BEM.Ct=zeros(ne,nbIt);
    BEM.F=zeros(ne,nbIt);
    %% vectors that will store the loads paramters for each elements
    BEM.Vr=zeros(ne,1);
    BEM.Ct=zeros(ne,1);
    BEM.Cp=zeros(ne,1);
    BEM.L=zeros(ne,1);
    BEM.D=zeros(ne,1);
    BEM.T=zeros(ne,1);
    BEM.Q=zeros(ne,1);
    BEM.P=zeros(ne,1);
    Pn=zeros(ne,1);
    Pt=zeros(ne,1);
    %% initialization
   
    for e=1:ne
        %initialization
        BEM.a(e,1)=0.2;
        BEM.aprime(e,1)=0.01;
        %BEM loop
        for i=1:nbIt
            %step 2
            BEM.phi(e,i)=atan( (1-BEM.a(e,i))/((1+BEM.aprime(e,i))*x(e)) );  %rad
            if(sin(BEM.phi)<0.01)
                % To avoid complex values
                F=1;
            else
                %prandtl tip correction
                 f=B/2*(R-r(e))/(r(e)*sin(BEM.phi(e,i)));
                 BEM.F(e,i)=2/pi*acos(exp(-f));
            end  
            %step 3
            BEM.alpha(e,i)=(BEM.phi(e,i)*180/pi)-(beta(e)); %deg
            %step 4
            alpha_data= Profiles.alpha;
            BEM.Cl(e,i) = interp1(alpha_data,Profiles.Cl,BEM.alpha(e,i));
            BEM.Cd(e,i) = interp1(alpha_data,Profiles.Cd,BEM.alpha(e,i));
            %step 5

            BEM.Cn(e,i)=BEM.Cl(e,i)*cos(BEM.phi(e,i))+BEM.Cd(e,i)*sin(BEM.phi(e,i));
            BEM.Ct(e,i)=BEM.Cl(e,i)*sin(BEM.phi(e,i))-BEM.Cd(e,i)*cos(BEM.phi(e,i));
            BEM.Vr(e)=V0*(1-BEM.a(e,i))/sin(BEM.phi(e,i));
            L=0.5*rho*norm(BEM.Vr(e)).^2*c(e)*BEM.Cl(e,i);
            D=0.5*rho*norm(BEM.Vr(e)).^2*c(e)*BEM.Cd(e,i);
            BEM.Pn(e) = L*cos(BEM.phi(e,i)) + D*sin(BEM.phi(e,i));   %load normal to the rotor plane
            BEM.Pt(e) = L*sin(BEM.phi(e,i)) - D*cos(BEM.phi(e,i));   %load tangential to the rotor plane
            

            %preparing for next iteration
           
            %by default the next a is :
            BEM.a(e,i+1)=1/( (4*BEM.F(e,i)*sin(BEM.phi(e,i))^2) /(sigma(e)*BEM.Cn(e,i))+1 );
            %Spera correction
            if(isequal(Param.correction,'Spera'))
                ac=0.34;
                if BEM.a(e,i)>ac
                    K=4*BEM.F(e,i)*(sin(BEM.phi(e,i)))^2/(sigma(e)*BEM.Cn(e,i));
                    BEM.a(e,i+1)=0.5*(2+K*(1-2*ac)-sqrt( (K*(1-2*ac)+2 )^2 + 4*(K*ac^2-1)    )  );
                end
            end             
     
            %relaxation
            BEM.a(e,i+1)=BEM.a(e,i+1)*Param.relaxation+(1-Param.relaxation)*BEM.a(e,i);
            
            BEM.aprime(e,i+1)=1/( (4*BEM.F(e,i)*sin(BEM.phi(e,i))*cos(BEM.phi(e,i))) /(sigma(e)*BEM.Ct(e,i))-1 )   ;     
                        %convergence criteria
            if (i>3 && abs(BEM.alpha(e,i)-BEM.alpha(e,i-1))<alphaCrit  && abs(BEM.alpha(e,i)-BEM.alpha(e,i-2)) < alphaCrit)
                break;
            end
        end %end while BEM for one element
        if(i==nbIt)
             BEM.Cp(e)=NaN;
             disp(['Algorithm did not converge :',sprintf('lambda=%.2f beta=%.2f V0=%.2f r=%.2f\n',lambda,beta(e),V0,r(e))])
        end
    end %end for loop on elements
    %%%% Returning Aerodynamic Forces
    BEM.Torque = B*getTorqueFromBlade(r,BEM.Pt,R);   %Rotor shaft torque at t in Newton
    BEM.Thrust = B*getThrustFromBlade(r,BEM.Pn,R);   %Rotor shaft thrust at t in Newton
    BEM.Power=omega*BEM.Torque;
    CP=BEM.Power/(0.5*rho*V0^3*pi*R^2);
    CT=BEM.Thrust/(0.5*rho*V0^2*pi*R^2);
    CQ=BEM.Torque/(0.5*rho*V0^2*pi*R^3);

end

