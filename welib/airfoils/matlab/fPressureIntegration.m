function [Force dF dFn]=fPressureIntegration(CP,P,N,Ncp,Pressure)
    dF=N*0; %initialization
    dFn=N*0;

    %% Using rectangle centered on P's between CPs where pressure is known
%     for i = 1:(length(Pressure)-1)
%         dx=norm(P(:,i+1) - P(:,i));
%         dp=Pressure(i+1);
%         dF(:,i+1)=dx*dp*N(:,i+1);
%         dFn(:,i+1)=dF(:,i+1)/dx;
%     end
%     %First Point
%     dx=norm(P(:,1) - P(:,end));
%     dp=Pressure(1);
%     dF(:,1)=-dx*dp*N(:,1);
%     dFn(:,1)=dF(:,1)/dx;
%     
%     Force=sum(dF,2);
    
    
    %% Using Trapezoidal method 
    % CP ; points where pressure is known
    % N normal in between CPs
    for i = 1:(length(Pressure)-1)
        dx=norm(CP(:,i+1) - CP(:,i));
        ptrapz=(Pressure(i)+Pressure(i+1))/2;
        dF(:,i)=- dx*ptrapz* N(:,i);
        dFn(:,i)=dF(:,i)/dx;
    end
    %Last Point
    dx=norm(CP(:,1) - CP(:,end));
    ptrapz=(Pressure(1)+Pressure(end))/2;
    dF(:,end)=dx*ptrapz*N(:,end);
    dFn(:,end)=dF(:,end)/dx;
    
    
    
    
    Force=sum(dF,2);
    %trapz(CP,dFn)
    
%     Lift=dot(Force,[-Vrel(2);Vrel(1)])/norm(Vrel);
%     Drag=dot(Force,Vrel)/norm(Vrel);
%     Cl=norm(Lift)/(0.5*rho*chord*norm(Vrel)^2);
%     Cd=norm(Drag)/(0.5*rho*chord*norm(Vrel)^2);    
% %    Lift=trapz(P(1,:),pinf-psup)
%     
%     Lift=[-Vrel(2);Vrel(1)]/norm(Vrel)*Lift;   
%     Drag=Vrel/norm(Vrel)*Drag;

end