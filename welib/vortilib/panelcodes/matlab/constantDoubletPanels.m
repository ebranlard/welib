function [pt, co, cp,Cl,Gamma_tot,Mu,AI,rhs ] = constantDoubletPanels( AirfoilPoints,U0,alpha,inverseMethod )
% inspired from:
%    PROGRAM N0. 3: CONSTANT STRENGTH DOUBLET
% Should reproduce Figure 11.22 from Katz and Plotkin


%--- Inputs
% AirfoilPoints: nx2 matrix with first and last point equal to [1 0]
% alpha: in degrees

if nargin==0
    %addpath(genpath('C:/Config/path/MatlabPath/'))
    %InitClear
    %require('PROFILES','v00');

    addpath('../../../airfoils/matlab/')

    param.tau   = 17   ;  % Tail Angle
    param.t_rel = 0.15 ;  % approx t/c (approx ymax_SS-ymin_PS)
    U0=1;
    % Inputs and param
    param.n=361;
    param.alpha=5;
    chord=2;
    off=0.5;

    % ---- Reference
    [AP_ref,~,~,Cp_ref] = fProfileVanDeVooren( param.t_rel,param.tau,1,param.n,2,U0,param.alpha); % to get Cp
    x_ref = AP_ref(:,1);

    % --- Points Used below
    % THIS
    [ AirfoilPoints PS SS ~] = fProfileVanDeVooren( param.t_rel,param.tau,1,param.n,2,U0,param.alpha); % to get Cp
    % OR 
    AirfoilPoints=readmatrix('../data/VonDeVooren_esp0.075_k1.906_AFOIL2.csv'); % HACK important, they have chord 0f 2 => different everything, AI, Cl, g etc

    chord=1;off=0.5;
    alpha=param.alpha;
    inverseMethod=2;
end


%% checks
ept=AirfoilPoints;
% if ~(norm(ept(1,:)+ept(end,:)-[2 0])<10^-6)
%     norm(ept(1,:)+ept(end,:)-[2 0])
%     error('Airfoil Point Format is nx2 with first and last point equal to [1 0]')
% end
if nargin==3
    inverseMethod=2;
end



%% init 
m= size(ept,1)-1; % 'NUMBER OF PANELS');
n=m+1; % number of points
ep=zeros(n,2); 
pt1=zeros(m,2);
pt2=zeros(m,2);
co=zeros(m,2); 
th=zeros(m,1);
aDvel=zeros(n,n);
bDvel=zeros(n,n);
rhs=zeros(n,1);
g=zeros(1,n); 
cp=zeros(1,m); 


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
Normal(:,1:2)=[-sin(th)  cos(th)];
Tangent(:,1:2)=[cos(th)  sin(th)];

% ESTABLISH COLOCATION POINTS
co=(pt2+pt1)./2;


% ESTABLISH INFLUENCE COEFFICIENTS
for i=1:m;
    for j=1:m;
        %         CONVERT THE COLOCATION POINT
        %         TO LOCAL PANEL COORDS.
        xt=co(i,1)-pt1(j,1);
        zt=co(i,2)-pt1(j,2);
        x2t=pt2(j,1)-pt1(j,1);
        z2t=pt2(j,2)-pt1(j,2);
        x=xt.*cos(th(j))+zt.*sin(th(j));
        z=-xt.*sin(th(j))+zt.*cos(th(j));
        x2=x2t.*cos(th(j))+z2t.*sin(th(j));
        z2=0;

        % SAVE PANEL LENGTHS
        if(i==1)
            dl(j)=x2;
        end;
        r1=sqrt(x.^2+z.^2);
        r2=sqrt((x-x2).^2+z.^2);

        %         COMPUTE THE VELOCITY INDUCED AT THE ITH
        %         COLOCATION POINT BY THE JTH PANEL
        if(i==j)
            ul=0;
%             wl=-1/(pi*x);  % this is only true is the control point is at the middle
            wl=-1/(2*pi)*(1/x-1/(x-x2)); % Katz 10.33
%             x
%             wl=-1/(3.14159265*x);
        else;
            ul=1/(2*pi)*(z./(r1.^2)-z./(r2.^2));
            wl=-1/(2*pi)*(x./(r1.^2)-(x-x2)./(r2.^2));
        end;
        u=ul.*cos(-th(j))+wl.*sin(-th(j));
        w=-ul.*sin(-th(j))+wl.*cos(-th(j));

        %         A(I,J) IS THE COMPONENT OF VELOCITY INDUCED IN THE
        %         DIRECTION NORMAL TO PANEL I BY PANEL J AT THE ITH
        %         COLOCATION POINT
        aDvel(i,j)=-u.*sin(th(i))+w.*cos(th(i));
        bDvel(i,j)=u.*cos(th(i))+w.*sin(th(i));
    end; j=m+1;

    %     INCLUDE THE INFLUENCE OF THE WAKE PANEL
    r=sqrt((co(i,1)-pt2(m,1)).^2 +(co(i,2)-pt2(m,2)).^2);

    u=1/(2*pi).*(co(i,2)./(r.^2));
    w=-1/(2*pi).*(co(i,1)-pt2(m,1))./(r.^2);
    aDvel(i,n)=-u.*sin(th(i))+w.*cos(th(i));
    bDvel(i,n)=u.*cos(th(i))+w.*sin(th(i));
    rhs(i,1)=U0*(cos(al).*sin(th(i))-sin(al).*cos(th(i)));   % MANU modif, this is:  -U0.n 
end;

% PREPARE THE MATRIX FOR SOLUTION BY PROVIDING
% A KUTTA CONDITION
for i=1:n; % !!! MANU modif
    aDvel(n,i)=0;
end; i=n+1+1;
aDvel(n,1)=-1;
aDvel(n,m)=1;
aDvel(n,n)=-1;
%     OPEN(10,FILE='A.DAT',STATUS='NEW')
%     DO I=1,N+1
%         DO j=1,N+1
%             write(10,*)
%     enddo
%     close(10)

% SOLVE FOR THE SOLUTION VECTOR OF DOUBLET STRENGTHS
if inverseMethod==2
    rhs1=rhs;
    a1=aDvel;
    g1=a1\rhs1;
     a2=aDvel;
%     a3=a;
%     a4=a;
%      a5=a;
% 
%     amat2=a; 
%     amat2(:,n+1)=rhs;
%     amat2(n+1,:)=0;

%     a1(n+1,1:n-1)=0; % adding an extra line of zero to make the matrix square
     a2(n+1,1:n-1)=1; % adding an extra line of zero to make the matrix square
%     a3(n+1,:)=0; % Not good
%     a4(:,n+1)=0; % Not good
%      a5(1:n,n+1)=1; % OK
%     rhs4=rhs;
    rhs1(n+1)=0;
%     g=a1\rhs1;
     g=a2\rhs1;
%     g3=a3\rhs1;
%     g4=a4\rhs4;
%      g5=a5\rhs4;
%     G=load('GD.DAT');
%     A=load('AID.DAT');
%     [A2,~,G2]=matrx(A,n+1,g);
%     figure,plot(1:n,g-pi,1:n,g2,1:n+1,g5-pi,1:n,G(1:n),'k')
%     kbd
%     [~,~,gmat]=matrx(amat2,n+1,g*0);
%     figure,plot(1:n,g,1:n,g2,1:n,g3,1:n+1,g4)
    aDvel=a2;
    rhs=rhs1;
else
    aDvel(:,n+1)=rhs;
    aDvel(n+1,:)=0; % adding an extra line of zero to make the matrix square
%     A=load('AID.DAT');
%     G=load('GD.DAT');
%     B=a-A; 
    warning('hack')
    [A2,~,G2]=matrx(A,n+1,g);
    kbd
    [a2,~,g2]=matrx(a,n+1,g);
end

mu=g;


% CONVERT DOUBLET STRENGTHS INTO TANGENTIAL
% VELOCITIES ALONG THE AIRFOIL SURFACE AND CP'S
% ON EACH OF THE PANELS
%% NEW the vectorial way
vtU0= U0*( cos(al)*cos(th) + sin(al)*sin(th));
vtDothers=bDvel(:,(1:(m+1)))*g(1:(m+1)); % !!!!!!!!!! here the m+1 should be important, you need the influence of the wake as well in the total velocity
vtDothers=vtDothers(1:m)';

for i=1:m;
    if(i~=1&&i~=m)
        % normal case
        r=sqrt((co(i+1,1)-co(i-1,1)).^2  +(co(i+1,2)-co(i-1,2)).^2);
        vtD(i)=(g(i+1)-g(i-1))./r;
    elseif(i==1) ;
        r=sqrt((co(2,1)-co(1,1)).^2 +(co(2,2)-co(1,2)).^2);
        vtD(i)=(g(2)-g(1))./r;
    elseif(i==m) ;
        r=sqrt((co(m,1)-co(m-1,1)).^2 +(co(m,2)-co(m-1,2)).^2);
        vtD(i)=(g(m)-g(m-1))./r;
    end;
end; i=m+1;
vtall=vtU0'+vtDothers+vtD/2;
cp=1-vtall.^2/U0^2;

%% Another try of differentiation
for i=1:m;
    phi(i)=U0*(co(i,1).*cos(al)+co(i,2).*sin(al))+g(i);
    phiU0(i)=U0*(co(i,1).*cos(al)+co(i,2).*sin(al));
end; i=m+1;
for i=1:m-1;
    r(i)=(dl(i+1)+dl(i))./2;
    velD(i)=-(g(i)-g(i+1))./r(i);
    velU0(i)=-(phiU0(i)-phiU0(i+1))./r(i);
end;
vtDothers2=(vtDothers(2:end)+vtDothers(1:end-1))/2;
vtDothers3=vtDothers(1:end-1);
vel=velU0 + 1*vtDothers3 + velD./2;
cp2=1-(vel).^2/U0^2;

%% For fun Normal velocity
% Doublet 
vnD=aDvel(1:m,(1:(m+1)))*mu(1:(m+1)); % !!!!!!!!!! here the m+1 should be important, you need the influence of the wake as well in the total velocity
% Free stream 
vnU0= U0*( -cos(al)*sin(th) + sin(al)*cos(th));  % cf -rhs(1:m)
% total normal velocity
vnall=vnU0+vnD;

%% Vectorial velocities
Vn=Normal.*[vnall(:) vnall(:)];
Vt=Tangent.*[vtall(:) vtall(:)];
Vall=Vn+Vt;

VtD=Tangent.*[vtD(:) vtD(:)]/2;
VtDo=Tangent.*[vtDothers(:) vtDothers(:)];

Vn0=Normal.*[vnU0(:) vnU0(:)];
Vt0=Tangent.*[vtU0(:) vtU0(:)];
V0=Vn0+Vt0;


% ----- Outputs
chord=max(AirfoilPoints(:,1))-min(AirfoilPoints(:,1));
Gamma_tot=g(m+1);
Cl=2*Gamma_tot/(U0*chord);
Mu=g; % gammas actually
AI=aDvel;

% due to the finite difference, the solution might be better approximated at panel points
pt=pt2(:,1);



if nargin==0
    figure,plot(g)
    figure,plot(aDvel*g)
    figure,plot(rhs)
    figure,plot(aDvel*g-rhs)
   
 figure, quiver(co(:,1),co(:,2),Vall(:,1),Vall(:,2)); hold all, plot(pt1(:,1),pt1(:,2),'k')



    figure,hold all
    plot(x_ref, Cp_ref,'k'),axis ij
    plot(pt2(1:end-1,1)/chord+off,cp(1:end-1),'+'),axis ij
    plot(pt2(1:end-1,1)/chord+off,cp2,'.'),axis ij
    legend('Theory', 'Cp 1', 'Cp2')
    % plot(pt2(1:end-1,1),cppot,'d'),axis ij
%     plot(pt2(1:end-1,1),cppot2,'d'),axis ij
    xlim([0 1])
    ylim([-1.8 1])
    
    keyboard










end

