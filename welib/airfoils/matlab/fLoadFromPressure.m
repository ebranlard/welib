function [Force dF dFn AeDataCalc CP N T P]=fLoadFromPressure(x,y,Pressure,Vrel)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Load from Cp
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Rot=[0 -1 ;1 0];

n=length(x);

CP=zeros(2,n); % CP control points!!
P=zeros(2,n);
T=zeros(2,n);
N=zeros(2,n);


CP(1,:)=x;
CP(2,:)=y;
Ncp=zeros(2,n);
Tcp=zeros(2,n);
for i=1:(n-1)
    Tcp(:,i)= (CP(:,i+1)-CP(:,i) ) /norm(CP(:,i+1)-CP(:,i));
    Ncp(:,i)=Rot*Tcp(:,i);
end
%last point
Tcp(:,n)= (CP(:,1)-CP(:,end) ) /norm(CP(:,1)-CP(:,end));
Ncp(:,n)=Rot*Tcp(:,n);


Vdelta_thickness=0*norm(CP(:,1)-CP(:,2))/3;

%computation of singularities Point positions
P(:,1:n-1)=(CP(:,2:n)+CP(:,1:n-1))/2-Ncp(:,1:n-1).*Vdelta_thickness;
%last point is extrapolated
P(:,n)=(CP(:,1)+CP(:,end))/2-Ncp(:,n).*Vdelta_thickness;


%projection
for i=2:n
    T(:,i)= (P(:,i)-P(:,i-1) ) /norm(P(:,i)-P(:,i-1));
    N(:,i)=Rot*T(:,i);
end
%First point
T(:,1)= (P(:,1)-P(:,end) ) /norm(P(:,1)-P(:,end));
N(:,1)=Rot*T(:,1);
N(:,floor(n/2)+1)=Ncp(:,floor(n/2)+1);
N(:,floor(n/2)+2)=Ncp(:,floor(n/2)+2);
T(:,floor(n/2)+1)=Tcp(:,floor(n/2)+1);
T(:,floor(n/2)+2)=Tcp(:,floor(n/2)+2);
T(:,1)=Tcp(:,1);
N(:,1)=Ncp(:,1);

[Force dF dFn]=fPressureIntegration(CP,P,N,Ncp,Pressure);

% calculation for verifications
e_chord=P(:,floor(n/2))-P(:,1);

AeDataCalc.chord=norm(CP(:,1)-CP(:,floor(n/2)));
AeDataCalc.thickness=max(CP(2,:))-min(CP(2,:));
AeDataCalc.thickness_rel=AeDataCalc.thickness/AeDataCalc.chord;
AeDataCalc.phi=atan2(Vrel(2),Vrel(1))*180/pi;
AeDataCalc.twist=atan2(e_chord(2),e_chord(1))*180/pi;
AeDataCalc.alpha=AeDataCalc.phi-AeDataCalc.twist;
