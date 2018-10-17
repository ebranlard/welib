function [Mass,step]=getMassMatrix(Data2,Eigen_M,gamma,theta_cone,mtot)
step=Data2(2:end,1)-Data2(1:end-1,1);
% Integration of Moment of inertia that goes into 2,2 position in the mass
% matrix
IM=3*sum(step.*(Data2(1:end-1,end).*Data2(1:end-1,1).^2.*cos(gamma).^2+...
     Data2(2:end,end).*Data2(2:end,1).^2.*cos(gamma).^2)/2);
% The coordinate system of the rotor rotates with blades therefore it must be transformed to the nacell coordinate system. It is done by the transformation function.
% The load force has two components "y" (paralel to the rotor plane) and "z"
% (tangential to the rotor plane). The integral is calcculated for both
% directions of each eigenmode in order to extract the generalized force for
% a unit displacement for the whole blade.
P=[zeros(length(Data2),2) Data2(:,end)];
for i=1:length(Data2)
     P_trans(i,:)=transformation(theta_cone,P(i,:));
end
for m=1:3
     m_1(m)=sum(step.*(P_trans(1:end-1,2).*Eigen_M(1:end-1,3,m)+...
         P_trans(1:end-1,3).*Eigen_M(1:end-1,2,m)+...
         P_trans(2:end,2).*Eigen_M(2:end,3,m)+...
         P_trans(2:end,3).*Eigen_M(2:end,2,m))/2);
end
% This is for the second column in the Mass matrix starting from the
3rd
% row and down
P=[zeros(length(Data2),1) Data2(:,end).*Data2(:,1)
zeros(length(Data2),1)];
for m=1:3
     m_2(m)=sum(step.*(P(1:end-1,2).*Eigen_M(1:end-1,3,m)+...
         P(1:end-1,3).*Eigen_M(1:end-1,2,m)+...
         P(2:end,2).*Eigen_M(2:end,3,m)+...
         P(2:end,3).*Eigen_M(2:end,2,m))/2);
end
% The diagonal positions 3,3 4,4 and 5,5 of the Mass matrix are
filled.
for m=1:3
     P=[zeros(length(Data2),1) Data2(:,end).*Eigen_M(:,3,m) Data2(:,end).*Eigen_M(:,2,m)];
     m_3(m)=sum(step.*(P(1:end-1,2).*Eigen_M(1:end-1,3,m)+...
         P(1:end-1,3).*Eigen_M(1:end-1,2,m)+...
         P(2:end,2).*Eigen_M(2:end,3,m)+...
         P(2:end,3).*Eigen_M(2:end,2,m))/2);
end
% First raw is the same as the first column in the Mass matrix
% It is repeated for other diagonal positions from 6 to 11 because
blades
% are considered to be the same.
Mass=zeros(11,11);          % Mass matrix for 11 DOF system
Mass(1,1)=mtot;
Mass(2,2)=IM;
for i=0:2
    Mass(3+i:3:9+i,1)=m_1(1+i);
    Mass(1,3+i:3:9+i)=m_1(1+i);
    Mass(3+i:3:9+i,2)=m_2(1+i);
    Mass(2,3+i:3:9+i)=m_2(1+i);
    Mass(3+i,3+i)=m_3(i+1);
    Mass(6+i,6+i)=m_3(i+1);
    Mass(9+i,9+i)=m_3(i+1);
end
end
 

