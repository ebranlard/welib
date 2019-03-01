function [U Grad]=getPointIncomingWindLegacy03(P,psi,WT,Wind,t,bComputeGrad)
% Now I go to line vectors....
% input P should be a matrix of size np*3, checks are done at the begining though
Model=Wind.Model;
Vhub=Wind.V0(:)';
nu=Wind.nu;
V0=norm(Vhub); 

np=length(P)/3;
if floor(np)~=np
    P=[0 0 P(:)];
    np=length(P);
end
if np~=3 
    if size(P,1)==3
        P=P';
    end
end


U=zeros(np,3);
if bComputeGrad
    Grad=zeros(np, 9);
else
    Grad=[];
end
if length(findstr('Constant',Model))>0
    U(:,1:3)=repmat(Vhub,np,1);
    %grad will be zero
else
    warning('To redo')
    keyboard

% if length(findstr('PowerLaw',Model))>0
% %     U(:,1:3)=rempat(Vhub,*(P(:,1)/60).^nu;
% %     V0=sqrt(sum((U(:,1:3).^2)'))' ;
% %     if bComputeGrad
% % 
% %     end
% end
% if length(findstr('TowerEffect',Model))>0
%     warning('Toredo')
%     keyboard
% %     if psi>80 & psi<280
% %           a=WT.Tower.r1-(P(1)-49)*(WT.Tower.r1/WT.Tower.H1);
% %     elseif psi>120 & psi<240
% %           a=2.375-(P(1)-28)*7.8125e-3;
% %     else
% %           a=0;
% %     end
% %     % projection of 0P on the ground
% %     r=sqrt(P(3)^2+P(2)^2);
% %     cosT=P(3)/r;
% %     sinT=-P(2)/r;
% %     %Polar velocity
% %     Vr=V0*(1-(a/r)^2)*cosT;
% %     V_theta=-V0*(1+(a/r)^2)*sinT;
% %     %Cartesian Velocity
% %     % Wrong in case of slanted flow!!!!!!!!!1
% %     Vpoint(3)=Vr*cosT-V_theta*sinT;
% %     Vpoint(2)=-Vr*sinT-V_theta*cosT;
% end
% 
% if length(findstr('Stochastic',Model))>0
% %     Vpoint=Vpoint+Wind.fTurb(P(3),-P(2),P(1),t);
end


end


