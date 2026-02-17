function [pt,co,cp,Cl,Gamma_tot,Mu,sig,AI ] = constantDoubletSourcePanelsPot( AirfoilPoints,U0,alpha ,inverseMethod)
% inspired from:
% PROGRAM N0. 8= CONsTANT STRENGTH SOURCEIDOUBLET POTENTIAL- KAtz

%--- Inputs
% AirfoilPoints: nx2 matrix with first and last point equal to [1 0]
% alpha: in degrees

if nargin==0
    addpath('../../../airfoils/matlab/')
    %InitClear
    %require('PROFILES','v00');
    param.tau   = 17   ;  % Tail Angle
    param.t_rel = 0.15 ;  % approx t/c (approx ymax_SS-ymin_PS)
    U0=10;
    % Inputs and param
    param.n=91;
    param.alpha=05;
    chord=2;
    off=0.5;

    % ---- Reference
    [AP_ref,~,~,Cp_ref] = fProfileVanDeVooren( param.t_rel,param.tau,1,param.n,2,U0,param.alpha); % to get Cp
    x_ref = AP_ref(:,1);

    % --- Points Used below
    % THIS
    % AirfoilPoints PS SS ~] = fProfileVanDeVooren( param.t_rel,param.tau,1,param.n,2,U0,param.alpha); % to get Cp
    % OR 
    AirfoilPoints=readmatrix('../data/VonDeVooren_esp0.075_k1.906_AFOIL2.csv'); % HACK important, they have chord 0f 2 => different everything, AI, Cl, g etc
    %x_ref =x_ref-0.5;

    alpha=param.alpha;
    inverseMethod=2; % 1= pivot fortran, 2=matlab
end



%% checks
ept=AirfoilPoints;
% if ~(norm(ept(1,:)+ept(end,:)-[2 0])<10^-6)
%     norm(ept(1,:)+ept(end,:)-[2 0])
%     error('Airfoil Point Format is nx2 with first and last point equal to [1 0]')
% end

if nargin==2
    inverseMethod=2; % 1= pivot fortran, 2=matlab
end


%% init 
m= size(ept,1)-1; % 'NUMBER OF PANELS');
n=m+1;

% warning('hack')
% ept=load(sprintf('AFOIL2_%d_%d.DAT',n-1,alpha));

ep=zeros(n,2); 
pt1=zeros(m,2);
pt2=zeros(m,2);
co=zeros(m,2); 
th=zeros(1,m);
a=zeros(n,n);
rhs=zeros(n,1);
b=zeros(n,n);
g=zeros(1,n); 
cp=zeros(1,m); 
sig=zeros(1,m); 


al=alpha*pi/180;

% Make sure PANELING is CLOCKWISE
if ept(floor(m/4),2)>0
    for i=1:m+1;
        ep(i,1)=ept(n-i+1,1);
        ep(i,2)=ept(n-i+1,2);
    end;
else
    ep=ept;
end

% ESTABLISH COORDINATES OF PANEL END POINTS
for i=1:m;
    pt1(i,1)=ep(i,1);
    pt2(i,1)=ep(i+1,1);
    pt1(i,2)=ep(i,2);
    pt2(i,2)=ep(i+1,2);
end;

% FIND PANEL ANGLES TH(J)
for i=1:m;
    dz=pt2(i,2)-pt1(i,2);
    dx=pt2(i,1)-pt1(i,1);
    th(i)=atan2(dz,dx);
end; 

% ESTABLISH COLOCATION POINTS
co=(pt2+pt1)./2;

% ESTABLISH SOURCE STRENGTHS (SIGMA=V DOT N)
for i=1:m;
    sig(i)=U0*(cos(al).*sin(th(i))-sin(al).*cos(th(i)));
end; i=m+1;

accum=0;
% ESTABLISH INFLUENCE COEFFICIENTS
for i=1:m;
    temp=0;
    for j=1:m;
        accum=accum+1;
        %         CONVERT THE COLOCATION POINT TO LOCAL PANEL COORDS
        xt=co(i,1)-pt1(j,1);
        zt=co(i,2)-pt1(j,2);
        x2t=pt2(j,1)-pt1(j,1);
        z2t=pt2(j,2)-pt1(j,2);
        x=xt.*cos(th(j))+zt.*sin(th(j));
        z=-xt.*sin(th(j))+zt.*cos(th(j));
        x2=x2t.*cos(th(j))+z2t.*sin(th(j));
        z2=0;
        %         SAVE PANEL LENGTHS
        if(i==1)
            dl(j)=x2;
        end;

        %         COMPUTE R AND THETA VALUES FOR THE COLOC. PoINT
        r1=sqrt(x.^2+z.^2);
        r2=sqrt((x-x2).^2+z.^2);
        th1=atan2(z,x);
        th2=atan2(z,x-x2);
        % DOUBLET POTENTIAL % Katz 11.64 (like 11.40 for source velocity)
        if(i==j)
            a(i,j)=0.5;
        else;
            a(i,j)=-1/(2*pi).*(th2-th1); 
        end;

        % Source Potential ( Katz 11.63)
        % ADD THEM UP TO GIVE THE RHS

%         DEB(accum,:)=[i,j,sig(j),temp,x,x2,th1,th2,z,r1,r2,log(r1),log(r2),sig(j)/3.14159265*(x*log(r1)),th2,th1,(x-x2)*log(r2)+z*(th2-th1),sig(j)/6.28319*(x*log(r1)-(x-x2)*log(r2)+z*(th2-th1)) ];
        if(i==j)
            temp=temp+sig(j)./pi.*(x.*log(r1));
        else;
            temp=temp+sig(j)./(2*pi).*(x.*log(r1)-(x-x2).*log(r2)+z.*(th2-th1));
        end;
%             IF(I.EQ.J) THEN
%                 TEMP=TEMP+SIG(J)/3.14159265*(X*LOG(R1))
%             ELSE
%                 TEMP=TEMP+SIG(J)/6.28319*(X*LOG(R1)-(X-X2)*LOG(R2)+Z*(TH2-TH1))
%             END IF
    end;
    %     ADD WAKE INFLUENCE COEFF.
    xw=co(i,1)-pt2(m,1);
    zw=co(i,2)-pt2(m,2);
    dthw=-atan(zw./xw);

    a(i,n)=-1/(2*pi).*(dthw);
    rhs(i,1)=temp; % MANU modif _>since we put it on the other side I thought I should have put a minus sign..
    if inverseMethod==1
        a(i,n+1)=temp;
    end
end 

% ADD AN EXPLICIT KUTTA CONDITION
if inverseMethod==2
    upperBound=n; % I would put n, but n+1 works...
else
    upperBound=n+1;
end
for i=1:upperBound; % !!! MANU modif
    a(n,i)=0;
end; i=n+1+1;
a(n,1)=-1;
a(n,m)=1;
a(n,n)=-1;

% SOLVE FOR THE SOLUTION VECTOR OF DOUBLET STRENGTHS
if inverseMethod==2
    g=a\rhs;
else
    [a2,~,g]=matrx(a,n+1,g);
end
% kbd
% 
% amat=a;
% amat(:,n+1)=rhs;
% g=a\rhs;
% [a2,~,g2]=matrx(amat,n+1,g*0);
% amat2=amat;
% amat2(n+1,:)=0;
% [a2,~,g3]=matrx(amat2,n+1,g*0);
% G=load('GDS.DAT');
% AIDS=load('AIDS.DAT');
% [a2,~,G2]=matrx(AIDS,n+1,g*0);
% figure, plot(1:n,g,'+',1:n,g2,'.',1:n,g3,'o',1:n,G,'k',1:n,G,'k+')




% CONVERT DOUBLET STRENGTHS INTO TANGENTIAL
% VELOCITIES ALONG THE AIRFOIL SURFACE AND CP'S
% ON EACH PANEL.
for i=1:m;
    phi(i)=U0*(co(i,1).*cos(al)+co(i,2).*sin(al))+g(i);
    phiU0(i)=U0*(co(i,1).*cos(al)+co(i,2).*sin(al));
end; i=m+1;

for i=1:m-1;
    r(i)=(dl(i+1)+dl(i))./2;
    vel(i)=(phi(i)-phi(i+1))./r(i);
    cp(i)=1-vel(i)^2/U0^2;
end;

% ----- Outputs
chord=max(AirfoilPoints(:,1))-min(AirfoilPoints(:,1));
Gamma_tot=g(m+1);
Cl=2*Gamma_tot/(U0*chord);
Mu=g; % used to be 1:end-2..
AI=a;
% due to the finite difference, the solution is well approximated at panel points
cp=cp(1:end-1);
pt=pt2(1:end-1,1);

end

