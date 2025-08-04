function [pt,co,cpall,Cl,Gamma_tot,Mus,Sigmas,AI,U,V ,X,Y] = constantDoubletSourcePanels( AirfoilPoints,U0,alpha,Vextern,Xgrid,Ygrid)

%
% Maybe:  See Katz section 11.3.1
%
% SETUP:
%  - m Doublet to be solved for
% 
%  - m Sources prescribed based on the freestream?
% 
% 

% This script is a purged version of fPanelCode2DDoubletSource_ManyOptions
%
% Additional features were added for external flow and flow field plot

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
    %[ AirfoilPoints PS SS ~] = fProfileVanDeVooren( param.t_rel,param.tau,1,param.n,2,U0,param.alpha); % to get Cp
    % OR 
    AirfoilPoints=readmatrix('../data/VonDeVooren_esp0.075_k1.906_AFOIL2.csv'); % HACK important, they have chord 0f 2 => different everything, AI, Cl, g etc

    [ptDS, coDS,CpDS,Cl,GammaDS,SigDS,AI ] = fPanelCode2DDoubletSourcePot( AirfoilPoints,U0,param.alpha,2 );
    [ptD, coD,CpD,Cl,GammaD,AID,rhsD ] = fPanelCode2DDoublet( AirfoilPoints,U0,param.alpha );
    %     AID(end+1,:)=0;
    alpha=param.alpha;

end


ept=AirfoilPoints;

% Param that seem to work: 
% Add Extra Line Condition, Don't double source, Include U0 in the RHS, the Source velocity is kept as its negative default (NegVel=0)
% Not as good: Add Extra Column Condition, Don't double source, Include U0 in the RHS, the Source velocity is kept as its negative default (NegVel=0)
% Not as good: Add Extra line Condition, double source, No U0, NegVel=o and Double source 
bAddExtraLineCondition=1; 
bAddExtraColumnCondition=0; 
bSimplifyMatrix=0; % !!!!!!!!!!!!!!!!!! has priority


%% init 
m= size(ept,1)-1; % 'NUMBER OF PANELS');
n=m+1;
ep=zeros(n,2); 
pt1=zeros(m,2);
pt2=zeros(m,2);
co=zeros(m,2); 
th=zeros(m,1);
Sigmas=zeros(m,1);
aDvel=zeros(n,n);
bDvel=zeros(n,n);
rhsvel=zeros(m,1);
rhsU0=zeros(m,1);
b=zeros(n,n);
g=zeros(1,n); 
cp=zeros(1,m); 
al=alpha*pi/180;

if nargin<4
    Vextern=zeros(m,2);
    Xgrid=[];
    Ygrid=[];
end
U=[];
V=[];
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


vnExtern = Vextern(:,1).*Normal(:,1)  + Vextern(:,2).*Normal(:,2);
vtExtern = Vextern(:,1).*Tangent(:,1) + Vextern(:,2).*Tangent(:,2);


% ESTABLISH COLOCATION POINTS
co=(pt2+pt1)./2;

% Panel lengths
for i=1:m, 
    dl(i)=norm(pt2(i,:)-pt1(i,:)); 
end
L=sum(dl);


% ESTABLISH SOURCE STRENGTHS (SIGMA=V DOT N)
for i=1:m;
    Sigmas(i)=U0*(cos(al).*sin(th(i))-sin(al).*cos(th(i)));
end; 


% ESTABLISH INFLUENCE COEFFICIENTS
for ic=1:m;  % loop on collocation point
    temp_Spot=0;
    temp_Svel_n=0;
    temp_Svel_others_n=0;
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

        %% Doublet velocity
        % COMPUTE THE VELOCITY INDUCED AT THE ITH
        % COLOCATION POINT BY THE JTH PANEL
        if(ic==j)
            ul=0;
            %             wl=-1/(pi*x); % assumes control points at the middle
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
        aSvel(ic,j)=-us*sin(th(ic))+ws*cos(th(ic));
        bSvel(ic,j)=us*cos(th(ic))+ws*sin(th(ic));

        % project to collocation point normal  !!!!!!!! At this point we put the source intensity sig!!!
        temp_Svel_n= temp_Svel_n+ Sigmas(j)*( -us*sin(th(ic))+ws*cos(th(ic)));
        %         b(ic,j)=u*cos(th(ic))+w*sin(th(ic))


    end;
    %% Doublet velocity (WAKE PANEL)
    r=sqrt((co(ic,1)-pt2(m,1)).^2 +(co(ic,2)-pt2(m,2)).^2);
    u=1/(2*pi).*(co(ic,2)./(r.^2));
    w=-1/(2*pi).*(co(ic,1)-pt2(m,1))./(r.^2);

    aDvel(ic,n)=-u.*sin(th(ic))+w.*cos(th(ic));
    bDvel(ic,n)=u.*cos(th(ic))+w.*sin(th(ic));

    %% RHS velocity 
    % free stream on normal surface
    rhsU0(ic,1)=U0*(cos(al).*sin(th(ic))-sin(al).*cos(th(ic)));  
    rhsvel(ic,1)=-temp_Svel_n;
end 
rhs=rhsvel+rhsU0 - vnExtern;


% ADD AN EXPLICIT KUTTA CONDITION
if bSimplifyMatrix
    aDvel(:,1)=aDvel(:,1)-aDvel(:,n);
    aDvel(:,n-1)=aDvel(:,n-1)+aDvel(:,n);
    aDvel=aDvel(1:(n-1),1:(n-1));
    if bAddExtraLineCondition
        aDvel(n,1:m)=1;
        rhs(n)=0;
    end
else
    rhs(n)=0;
    if bAddExtraColumnCondition
        upperBound=n+1; % I would put n, but n+1 works...
    else
        upperBound=n;
    end
    for i=1:upperBound; % !!! MANU modif
        aDvel(n,i)=0;
    end; i=n+1+1;
    aDvel(n,1)=-1;
    aDvel(n,m)=1;
    aDvel(n,n)=-1;

    aDpot(n,1:n)=aDvel(n,1:n);
    if bAddExtraLineCondition
        aDvel(n+1,1:(n))=1;
        rhs(n+1)=0;
    end
end

g=aDvel\rhs;

mu=g(1:m+1);

%% Doublet formulation
vDothers=bDvel(:,(1:(m+1)))*g(1:(m+1)); % !!!!!!!!!! here the m+1 should be important, you need the influence of the wake as well in the total velocity
vDothers=vDothers(1:m)';
vSothers=(bSvel(:,1:m)*Sigmas(1:m));
vSothers=vSothers(1:m)';

%%  Doublet gradient and velocity
for i=1:m;
    phi(i)=U0*(co(i,1).*cos(al)+co(i,2).*sin(al))+g(i);
    phiU0(i)=U0*(co(i,1).*cos(al)+co(i,2).*sin(al));
end;
for i=1:m-1;
    r(i)=(dl(i+1)+dl(i))./2;
    vtD(i)=-(g(i)-g(i+1))./r(i);
    vtU0(i)=-(phiU0(i)-phiU0(i+1))./r(i);
end;
vtD(m)=vtD(m-1);
vtU0(m)=vtU0(m-1);

ptvel=pt2;
vDothers2=[ (vDothers(2:end)+vDothers(1:end-1))/2 vDothers(end)  ];
vSothers2=[ (vSothers(2:end)+vSothers(1:end-1))/2 vSothers(end)  ];

vtall=vtU0(:) + 1*vDothers(:) + vtD(:)./2 +1*vSothers(:) + vtExtern(:);
cpall=1-(vtall).^2/U0^2;


%% For fun Normal velocity
% Source contributions 
vnSothers=aSvel(1:m,1:m)*Sigmas;
% Doublet 
vnD=aDvel(1:m,(1:(m+1)))*mu(1:(m+1)); % !!!!!!!!!! here the m+1 should be important, you need the influence of the wake as well in the total velocity
% Free stream 
vnU0= U0*( -cos(al)*sin(th) + sin(al)*cos(th));  % cf -rhs(1:m)
% total normal velocity
vnall=vnU0+1*vnD+1*vnSothers;

%% Vectorial velocities
Vn=Normal.*[vnall(:) vnall(:)];
Vt=Tangent.*[vtall(:) vtall(:)];
Vall=Vn+Vt;

Vn0=Normal.*[vnU0(:) vnU0(:)];
Vt0=Tangent.*[vtU0(:) vtU0(:)];
V0=Vn0+Vt0;


%% 
if ~isempty(Xgrid)
    [X,Y]=meshgrid(Xgrid,Ygrid);
    U=zeros(size(X));
    V=zeros(size(X));
    ncp=prod(size(X)); % we will use linear index for matrices..
    % 
    for icp=1:ncp;  % loop on collocation point
        %         if sqrt(X(icp)^2+Y(icp)^2)<0.05
        %             kbd
        %         end
        nCurrentPointInside=0; % counter of how many time the point has negative z in panel coordinates
        bUseCPValue=0;
        for j=1:m;
            % CONVERT THE COLOCATION POINT TO LOCAL PANEL COORDS
            % collocation point
            xt=X(icp)-pt1(j,1);
            zt=Y(icp)-pt1(j,2);
            % second panel point
            x2t=pt2(j,1)-pt1(j,1);
            z2t=pt2(j,2)-pt1(j,2);
            % collocation point in panel coord
            x=xt.*cos(th(j))+zt.*sin(th(j));
            z=-xt.*sin(th(j))+zt.*cos(th(j));
            % point 2 in panel coordinates (point 1 is 0,0)
            x2=x2t.*cos(th(j))+z2t.*sin(th(j));
            z2=0;

            % Test if inside
            if z<0
                nCurrentPointInside=nCurrentPointInside +1  ;
            end

            %% Principal value if on top of surface
            if abs(z)<10^-3 && x>=0 && x<=x2
                %                 vPrincip=Tangent(j,:)*vtD(j); % TODO I could do something a bit smoother BL style
                %                 U(icp)=U(icp)+vPrincp(1);
                %                 V(icp)=V(icp)+vPrincp(2);

                bUseCPValue=1;
                jCP=j;
%                 disp('breaking')
                break
            end


            % COMPUTE R AND THETA VALUES FOR THE COLOC. PoINT
            r1=sqrt(x.^2+z.^2);
            r2=sqrt((x-x2).^2+z.^2);
            th1=atan2(z,x);
            th2=atan2(z,x-x2);

            %% Doublet velocity
            %             if(ic==j)
            %                 ul=0;
            %                 wl=-1/(2*pi)*(1/x-1/(x-x2)); % Katz 10.33
            %             else;
            ul=1/(2*pi)*(z./(r1.^2)-z./(r2.^2));
            wl=-1/(2*pi)*(x./(r1.^2)-(x-x2)./(r2.^2));
            %             end;
            % back to ref frame
            ud=ul.*cos(-th(j))+wl.*sin(-th(j));
            wd=-ul.*sin(-th(j))+wl.*cos(-th(j));

            %% Source Velocity due to panel j at collocation point ic
            %         if(ic==j)
            %             uls=0;
            %             wls=0.5;   % Katz 10.24
            %         else
            uls=1/(2*pi)*log(r1/r2);  % Katz 10.17
            wls=1/(2*pi)*(th2-th1);   % Katz 10.18
            %         end 
            % return velocity to global ref. frame
            us=uls*cos(-th(j))+wls*sin(-th(j));
            ws=-uls*sin(-th(j))+wls*cos(-th(j));

            U(icp)=U(icp)+Sigmas(j)*us+mu(j)*ud; 
            V(icp)=V(icp)+Sigmas(j)*ws+mu(j)*wd; 
        end;
        if bUseCPValue
            U(icp)=Vall(jCP,1);
            V(icp)=Vall(jCP,2);
        else
            %% Doublet velocity (WAKE PANEL)
            r=sqrt((X(icp)-pt2(m,1)).^2 +(Y(icp)-pt2(m,2)).^2);
            u=1/(2*pi).* (Y(icp)./(r.^2));
            w=-1/(2*pi).*(X(icp)-pt2(m,1))./(r.^2);
            U(icp)=U(icp)+mu(m+1)*u;
            V(icp)=V(icp)+mu(m+1)*w;

            %% RHS velocity  (could be outside the loop, but it's because of Vall)
            U(icp)=U(icp)+U0*cos(al);  
            V(icp)=V(icp)+U0*sin(al);  
        end 
        if nCurrentPointInside==m
%             disp('it happen')
            U(icp)=NaN;
            V(icp)=NaN;
        end
    end % loop on icp
    U(isnan(U))=0; 
    V(isnan(V))=0; 
end




if nargin==0

    figure, quiver(co(:,1),co(:,2),Vall(:,1),Vall(:,2)); hold all, plot(pt1(:,1),pt1(:,2),'k')

    figure,hold all,  plot(pt1(:,1),pt1(:,2),'k')
    quiver(co(2,1),co(2,2),Vall(2,1),Vall(2,2)); 
    quiver(co(1,1),co(1,2),Vall(1,1),Vall(1,2),'k'); 
    quiver(co(m,1),co(m,2),Vall(m,1),Vall(m,2),'k');
    quiver(co(m-1,1),co(m-1,2),Vall(m-1,1),Vall(m-1,2));


    figure,hold all
    plot(x_ref, Cp_ref,'k'),axis ij
    plot(ptvel(1:end-1,1)/chord+off,cpall(1:end-1),'k+'),axis ij
    % plot(pt2(1:end-1,1),cppot,'d'),axis ij
    %     plot(pt2(1:end-1,1),cppot2,'d'),axis ij
    xlim([0 1])
    ylim([-1.8 1])


    %%
    %     figure, hold all, quiver(co(:,1),co(:,2),Normal(:,1),Normal(:,2),'k'); quiver(co(:,1),co(:,2),Tangent(:,1),Tangent(:,2));
    %       figure, hold all, quiver(co(:,1),co(:,2),Normal(:,1).*vnU0,Normal(:,2).*vnU0,'k'); quiver(co(:,1),co(:,2),Tangent(:,1).*vtU0,Tangent(:,2).*vtU0);
    %       figure, hold all, quiver(co(:,1),co(:,2),V0(:,1),V0(:,2),'k');axis equal

    keyboard
end




% ----- Outputs
chord=max(AirfoilPoints(:,1))-min(AirfoilPoints(:,1));
Gamma_tot=g(m+1);
Cl=2*Gamma_tot/(U0*chord); % check that, chord??
Mus=g; % used to be 1:end-2..
AI=aDvel;
% due to the finite difference, the solution is well approximated at panel points
cpall=cpall(1:end-1);
pt=pt2(1:end-1,1);


%% END Script 2DDoubletSource
