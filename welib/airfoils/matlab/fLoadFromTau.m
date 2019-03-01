function [Force dF dFn CP]=fLoadFromTau(x,y,Tau)

CP(1,:)=x;
CP(2,:)=y;

    dF=Tau*0; %initialization
    dFn=Tau*0;

    %% Using rectangle centered on P's between CPs where pressure is known
    for i = 1:(length(x)-1)
        dx=norm(CP(:,i+1) - CP(:,i));
        dF(:,i)=dx*Tau(:,i);
        dFn(:,i)=dF(:,i)/dx;
    end
    %Last Point
    dx=norm(CP(:,1) - CP(:,end));
    dF(:,end)=dx*Tau(:,end);
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