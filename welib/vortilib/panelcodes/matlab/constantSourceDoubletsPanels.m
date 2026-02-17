function [pt,co,cp,Cl,Gamma_tot,Mu,Sigmas,AI ] = constantSourceDoubletsPanels( AirfoilPoints,U0,alpha)
%
% NOTE: Seems buggy for now
%
% Maybe:  See Katz section 11.3.1
%
% SETUP:
%  - m Sources
%  - 1 Gamma to scale the entire mu
%
%  - m Doublet prescribed based on curvilinear coordinate
% 
% 

if nargin==0
    addpath('../../../airfoils/matlab/')
    %InitClear
    %require('PROFILES','v00');
    param.tau   = 07   ;  % Tail Angle
    param.t_rel = 0.15 ;  % approx t/c (approx ymax_SS-ymin_PS)
    U0=10;
    % Inputs and param
    param.n=361;
    param.alpha=05;
    chord=2; off=0.5;

    % ---- Reference
    [AP_ref,~,~,Cp_ref] = fProfileVanDeVooren( param.t_rel,param.tau,1,param.n,2,U0,param.alpha); % to get Cp
    x_ref = AP_ref(:,1);

    % --- Points Used below
    % THIS
    %[ AirfoilPoints PS SS ~] = fProfileVanDeVooren( param.t_rel,param.tau,1,param.n,2,U0,param.alpha); % to get Cp
    % OR 
    AirfoilPoints=readmatrix('../data/VonDeVooren_esp0.075_k1.906_AFOIL2.csv'); % HACK important, they have chord 0f 2 => different everything, AI, Cl, g etc


    [ptDS, coDS,CpDS,Cl,GammaDS,SigDS,AI ] = fPanelCode2DDoubletSourcePot( AirfoilPoints,U0,param.alpha,2 );
    [ptD, coD,CpD,Cl,GammaD,AID,rhsD ] = fPanelCode2DDoublet( AirfoilPoints,U0,param.alpha );
%     AID(end+1,:)=0;
    alpha=param.alpha;
end



ept=AirfoilPoints;

%% init 
m= size(ept,1)-1; % 'NUMBER OF PANELS');
n=m+1;
ep=zeros(n,2); 
pt1=zeros(m,2);
pt2=zeros(m,2);
co=zeros(m,2); 
th=zeros(m,1);
aSvel=zeros(m,m);
bSvel=zeros(m,m);
aDvel=zeros(m,n);
bDvel=zeros(m,n);
rhs=zeros(m,1);
a=zeros(n,n);
b=zeros(n,n);
g=zeros(1,n); 

Svel_n_self=zeros(1,m); 
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


% Panel lengths
for i=1:m, 
    dl(i)=norm(pt2(i,:)-pt1(i,:)); 
end
L=sum(dl);
s_coord=cumsum(dl);%-dl(1)/2;

% ESTABLISH Doublet STRENGTHS (mu = s )
mu_nocirc=s_coord;
mu_nocirc(m+1)=mu_nocirc(m)-mu_nocirc(1); % Kutta condition for wake panel
% mu_nocirc(m+1)=0; % Kutta condition for wake panel


Q0t=U0*(cos(al).*cos(th)+sin(al).*sin(th));
Q0n=U0*(cos(al).*sin(th)-sin(al).*cos(th));
% ESTABLISH INFLUENCE COEFFICIENTS
for ic=1:m;  % loop on collocation point
    for j=1:m;
        % CONVERT THE COLOCATION POINT TO LOCAL PANEL COORDS
        % collocation point
        xt=co(ic,1)-pt1(j,1);
        zt=co(ic,2)-pt1(j,2);
        % second panel point
        x2t=pt2(j,1)-pt1(j,1);
        z2t=pt2(j,2)-pt1(j,2);
        % collocation point in panel coord
        x=xt.*cos(th(j))+zt.*sin(th(j));
        z=-xt.*sin(th(j))+zt.*cos(th(j));
        % point 2 in panel coordinates (point 1 is 0,0)
        x2=x2t.*cos(th(j))+z2t.*sin(th(j));
        z2=0;

        % COMPUTE R AND THETA VALUES FOR THE COLOC. PoINT
        r1=sqrt(x.^2+z.^2);
        r2=sqrt((x-x2).^2+z.^2);
        th1=atan2(z,x);
        th2=atan2(z,x-x2);

        %% Doublet velocity (prescribed - kept as matrix of unit strengtg for now for easy debugging) 
        % COMPUTE THE VELOCITY INDUCED AT THE ITH COLOCATION POINT BY THE JTH PANEL
        if(ic==j)
            ul=0;
            wl=-1/(2*pi)*(1/x-1/(x-x2)); % Katz 10.33
        else;
            ul=1/(2*pi)*(z./(r1.^2)-z./(r2.^2));
            wl=-1/(2*pi)*(x./(r1.^2)-(x-x2)./(r2.^2));
        end;
        % back to ref frame
        u=ul.*cos(-th(j))+wl.*sin(-th(j));
        w=-ul.*sin(-th(j))+wl.*cos(-th(j));
        % projection on collocation point normal (and tangent)
        aDvel(ic,j)=-u.*sin(th(ic))+w.*cos(th(ic));
        bDvel(ic,j)=u.*cos(th(ic))+w.*sin(th(ic));

        %% Source Velocity due to panel j at collocation point ic
        % compute source velocity in local ref. frame
        if(ic==j)
            uls=0;
            wls=0.5;   % Katz 10.24
        else
            uls=1/(2*pi)*log(r1/r2);  % Katz 10.17
            wls=1/(2*pi)*(th2-th1);   % Katz 10.18
        end 
        % return velocity to global ref. frame
        us=uls*cos(-th(j))+wls*sin(-th(j));
        ws=-uls*sin(-th(j))+wls*cos(-th(j));
        % project to collocation point normal and tangent
        aSvel(ic,j)=-us*sin(th(ic))+ws*cos(th(ic));
        bSvel(ic,j)=us*cos(th(ic))+ws*sin(th(ic));
    end

    %% Doublet velocity (WAKE PANEL)
    r=sqrt((co(ic,1)-pt2(m,1)).^2 +(co(ic,2)-pt2(m,2)).^2);
    uwd=1/(2*pi).*(co(ic,2)./(r.^2));
    wwd=-1/(2*pi).*(co(ic,1)-pt2(m,1))./(r.^2);
    aDvel(ic,n)=-uwd.*sin(th(ic))+wwd.*cos(th(ic));
    bDvel(ic,n)=uwd.*cos(th(ic))+wwd.*sin(th(ic));
end 

%% RHS velocity , free stream on normal surface
rhs(:,1)=- U0*( -cos(al)*sin(th) + sin(al)*cos(th));  

%% Sum Doublet contribs
viD_nocirc=aDvel*mu_nocirc';


%% Building matrix
a(1:m,1:m)=aSvel;
a(1:m,n)=viD_nocirc;
% ADD a KUTTA CONDITION for Sources, Equality of tangential velocity ADD principal value!!!!!!!!!!!!!!!!!!...
Vt1=[bSvel(1,:) bDvel(1,:)*mu_nocirc'+1/2];
VtM=[bSvel(m,:) bDvel(m,:)*mu_nocirc'+1/2];
U0t1=  U0*( cos(al)*cos(th(1)) + sin(al)*sin(th(1)));  
U0tM=  U0*( cos(al)*cos(th(m)) + sin(al)*sin(th(m)));  
a(n,1:n)=Vt1+VtM;  % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Signs
rhs(n)=-(U0t1+U0tM);

% SOLVE FOR THE SOLUTION VECTOR OF SOURCE STRENGTHS
g=a\rhs;

Sigmas=g(1:m);
Gamma=g(m+1);
mu_circ=mu_nocirc*Gamma;


%% Computing tangential velocity
% Source contributions 
vtSothers=bSvel(1:m,1:m)*Sigmas;
% Doublet (excluding principal)
vtDothers=bDvel(1:m,(1:(m+1)))*mu_circ(1:(m+1))'; % !!!!!!!!!! here the m+1 should be important, you need the influence of the wake as well in the total velocity
% Doublet principal value and free stream with centered diff on pt2
phiU0=U0*(co(:,1)*cos(al)+co(:,2)*sin(al));
vtD=zeros(m,1);
vtU0=zeros(m,1);
for i=1:m-1;
    r(i)=(dl(i+1)+dl(i))./2;
%     r(i)=dl(i);
%     r(i)=dl(i+1);
    vtD(i)=-(mu_circ(i)-mu_circ(i+1))./r(i);
    r(i)=(dl(i+1)+dl(i))./2;
    vtU0(i)=-(phiU0(i)-phiU0(i+1))./r(i);
end;
vtDreal=vtD*0+Gamma; % !!!!!!!!!!!HACK 
% vtD(1)=vtDreal(1);
vtD(m)=vtDreal(m);
% vtD=vtDreal;

% vtD=vtD*0; % !!!!!!!!!!!HACK 
vtU0real= U0*( cos(al)*cos(th) + sin(al)*sin(th));     % --------------------<< try it 
% vtU0(1)=vtU0real(1);
vtU0(m)=vtU0real(m);
% vtU0(1:m)=vtU0real(1:m);

% total tangential velocity
vtDothers2=[ (vtDothers(2:end)+vtDothers(1:end-1))/2;vtDothers(end)  ];
vtSothers2=[ (vtSothers(2:end)+vtSothers(1:end-1))/2;vtSothers(end)  ];
vtall=vtU0+1*vtDothers+1*vtD/2+1*vtSothers;
cpall=1-(vtall).^2/U0^2;


%% For fun Normal velocity
% Source contributions 
vnSothers=aSvel(1:m,1:m)*Sigmas;
% Doublet 
vnD=aDvel(1:m,(1:(m+1)))*mu_circ(1:(m+1))'; % !!!!!!!!!! here the m+1 should be important, you need the influence of the wake as well in the total velocity
% Free stream 
vnU0= U0*( -cos(al)*sin(th) + sin(al)*cos(th));  % cf -rhs(1:m)
% total normal velocity
vnall=vnU0+1*vnD+1*vnSothers;

%% Vectorial Velocities
Vn=Normal.*[vnall vnall];
Vt=Tangent.*[vtall vtall];
Vall=Vn+Vt;

Vn0=Normal.*[vnU0 vnU0];
Vt0=Tangent.*[vtU0 vtU0];
V0=Vn0+Vt0;

if nargin==0

%%
%     figure, hold all, quiver(co(:,1),co(:,2),Normal(:,1),Normal(:,2),'k'); quiver(co(:,1),co(:,2),Tangent(:,1),Tangent(:,2));
%       figure, hold all, quiver(co(:,1),co(:,2),Normal(:,1).*vnU0,Normal(:,2).*vnU0,'k'); quiver(co(:,1),co(:,2),Tangent(:,1).*vtU0,Tangent(:,2).*vtU0);
%       figure, hold all, quiver(co(:,1),co(:,2),V0(:,1),V0(:,2),'k');axis equal

     figure, quiver(co(:,1),co(:,2),Vall(:,1),Vall(:,2)); hold all, plot(pt1(:,1),pt1(:,2),'k')
     
     figure,hold all,  plot(pt1(:,1),pt1(:,2),'k')
                      quiver(co(2,1),co(2,2),Vall(2,1),Vall(2,2)); 
                      quiver(co(1,1),co(1,2),Vall(1,1),Vall(1,2),'k'); 
                      quiver(co(m,1),co(m,2),Vall(m,1),Vall(m,2),'k');
                      quiver(co(m-1,1),co(m-1,2),Vall(m-1,1),Vall(m-1,2));

%      figure, quiver(co(1,1),co(1,2),Vall(:,1),Vall(:,2)); hold all, plot(pt1(:,1),pt1(:,2),'k')


   
    figure,hold all,plot(a*g-rhs),title('Error in solving rhs')

    figure,hold all, plot(vtU0), plot(vtSothers), plot(vtDothers), plot(vtD/2)
    legend('U0','Sources','Doublets Neighbor','Doublet Principal');
    
    figure,hold all
    plot(x_ref,Cp_ref,'k'),axis ij
%     plot(pt2(1:end-1,1)/chord+0.5,cpdonly,'+'),axis ij  % this on has only doublet contrib, has to be wrong...
    % plot(pt2(1:end-1,1),cp2(1:end-1),'.'),axis ij
    % plot(co(1:end-1,1),cp2(1:end-1),'.'),axis ij
    plot(pt2(:,1)/chord+off,cpall,'b-'),axis ij
%     plot(co(:,1)/chord+off,cpall,'r-'),axis ij
    % plot(pt2(1:end-1,1),cppot,'d'),axis ij
%     plot(pt2(1:end-1,1),cppot2,'d'),axis ij
    xlim([0 1])
%     ylim([-10 1])


    keyboard
end




% ----- Outputs
chord=max(AirfoilPoints(:,1))-min(AirfoilPoints(:,1));
Gamma_tot=Gamma*L;
Cl=2*Gamma_tot/(U0*chord); % check that, chord??
Sigmas=Sigmas; % used to be 1:end-2..
Mu=mu_circ;
AI=aDvel;
% due to the finite difference, the solution is well approximated at panel points
cp=cpall;
pt=pt2(:,1);


%% END Script 2DDoubletSource
