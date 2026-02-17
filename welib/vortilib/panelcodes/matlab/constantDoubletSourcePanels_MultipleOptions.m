function [pt,co,cpall,cpdonly,Cl,Gamma_tot,Mu,sig,AI ] = constantDoubletSourcePanels_MultipleOptions( AirfoilPoints,U0,alpha)

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
    %x_ref =x_ref-0.5;

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
bDoubleSource=0; 
bSourceJustOthers=0; 
bRHS_U0=1; 
bRHS_NegVel=0; 
bRHS_NegSig=0; 


%% init 
m= size(ept,1)-1; % 'NUMBER OF PANELS');
n=m+1;
ep=zeros(n,2); 
pt1=zeros(m,2);
pt2=zeros(m,2);
co=zeros(m,2); 
th=zeros(1,m);
aDvel=zeros(n,n);
bDvel=zeros(n,n);
aDpot=zeros(n,n);
rhsvel=zeros(m,1);
rhsU0=zeros(m,1);
rhspot=zeros(m,1);
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

% ESTABLISH COLOCATION POINTS
co=(pt2+pt1)./2;

% ESTABLISH SOURCE STRENGTHS (SIGMA=V DOT N)
for i=1:m;
    sig(i)=U0*(cos(al).*sin(th(i))-sin(al).*cos(th(i)));
    if bDoubleSource
        sigD(i)=2*U0*(cos(al).*sin(th(i))-sin(al).*cos(th(i)));
    else
        sigD(i)=U0*(cos(al).*sin(th(i))-sin(al).*cos(th(i)));
    end
    if bRHS_NegSig
        sigD(i)=-sigD(i);
    end
end; 
Q0t=U0*(cos(al).*cos(th)+sin(al).*sin(th));
Q0n=U0*(cos(al).*sin(th)-sin(al).*cos(th));
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
        % SAVE PANEL LENGTHS
        if(ic==1)
            dl(j)=x2;
        end;

        % COMPUTE R AND THETA VALUES FOR THE COLOC. PoINT
        r1=sqrt(x.^2+z.^2);
        r2=sqrt((x-x2).^2+z.^2);
        th1=atan2(z,x);
        th2=atan2(z,x-x2);
        %% Doublet Potential
        % DOUBLET POTENTIAL % Katz 11.64 (like 11.40 for source velocity)
        if(ic==j)
            aDpot(ic,j)=0.5;
        else;
            aDpot(ic,j)=-1/(2*pi).*(th2-th1); 
        end;
        %% Source Potential ( Katz 11.63)
        % ADD THEM UP TO GIVE THE RHS
        if(ic==j)
            temp_Spot=temp_Spot+sig(j)./pi.*(x.*log(r1));
        else;
            temp_Spot=temp_Spot+sig(j)./(2*pi).*(x.*log(r1)-(x-x2).*log(r2)+z.*(th2-th1));
        end;


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
        if ic~=j
            temp_Svel_n= temp_Svel_n+ sigD(j)*( -us*sin(th(ic))+ws*cos(th(ic)));
            temp_Svel_others_n= temp_Svel_others_n+ sigD(j)*( -us*sin(th(ic))+ws*cos(th(ic)));
        else
            if ~bSourceJustOthers
                temp_Svel_n= temp_Svel_n+ sigD(j)*( -us*sin(th(ic))+ws*cos(th(ic)));
            end
            Svel_n_self(ic)= sigD(j)*( -us*sin(th(ic))+ws*cos(th(ic)));
        end
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
%     if ic<m
%         rhsU0(ic,1)=(Q0n(ic+1)+Q0n(ic))/2;
%     else
%         rhsU0(ic,1)=Q0n(ic);
%     end
    % source influence (minus sign)
    rhsvel(ic,1)=-temp_Svel_n;
    Svel_n_others(ic,1)=temp_Svel_others_n;

    %% Doublet Potential (WAKE PANEL)
    xw=co(ic,1)-pt2(m,1);
    zw=co(ic,2)-pt2(m,2);
    dthw=-atan(zw./xw);
    aDpot(ic,n)=-1/(2*pi).*(dthw);
    %% RHS potential 
    rhspot(ic,1)=temp_Spot; % MANU modif
end 
if bRHS_NegVel
    signrhsvel=-1;
else
    signrhsvel=1;
end
rhs=signrhsvel*rhsvel+bRHS_U0*rhsU0;


% ADD AN EXPLICIT KUTTA CONDITION
if bSimplifyMatrix
    aDvel(:,1)=aDvel(:,1)-aDvel(:,n);
    aDvel(:,n-1)=aDvel(:,n-1)+aDvel(:,n);
    aDvel=aDvel(1:(n-1),1:(n-1));
    aDpot(:,1)=aDpot(:,1)-aDpot(:,n);
    aDpot(:,n-1)=aDpot(:,n-1)+aDpot(:,n);
    aDpot=aDpot(1:(n-1),1:(n-1));

    if bAddExtraLineCondition
         aDvel(n,1:m)=1;
         rhs(n)=0;
    end
else
    rhs(n)=0;
    rhspot(n)=0;
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

% SOLVE FOR THE SOLUTION VECTOR OF DOUBLET STRENGTHS
g=aDvel\rhs;
gpot=aDpot\rhspot;




%% Doublet formulation - Potential part
Q0t=U0*(cos(al).*cos(th)+sin(al).*sin(th));
for i=1:m;
    phi(i)=U0*(co(i,1).*cos(al)+co(i,2).*sin(al))+gpot(i); %<<<<<gpot!!!   also remeber phi_infty=x cos alpha + y sin alpha
    phib(i)=gpot(i);
end;
vel=0;
for i=1:m-1;
    r=(dl(i+1)+dl(i))./2;
    vel2=(phi(i)-phi(i+1))./r;
    vel=(phib(i)-phib(i+1))./r-(Q0t(i+1)+Q0t(i))/2;  % TODO MANU REMBER THAT THIS IS VERY IMPORTANT !!!!! IF You don't use the mean, results are a bit shifted
    cppot(i)=1-vel.^2/U0^2;
    cppot2(i)=1-vel2.^2/U0^2;
end; 



%% Doublet formulation
% CONVERT DOUBLET STRENGTHS INTO TANGENTIAL
% VELOCITIES ALONG THE AIRFOIL SURFACE AND CP'S
% ON EACH PANEL.
vDothers=bDvel(:,(1:(m+1)))*g(1:(m+1)); % !!!!!!!!!! here the m+1 should be important, you need the influence of the wake as well in the total velocity
vDothers=vDothers(1:m)';
vSothers=(bSvel(:,1:m)*sigD(1:m)');
vSothers=vSothers(1:m)';
% for i=1:m;
%     phi(i)=co(i,1).*cos(al)+co(i,2).*sin(al)+g(i);
%     phib(i)=g(i);
% end;
% vel=0;
% for i=1:m-1;
%     r=(dl(i+1)+dl(i))./2;
%     vel=(phi(i)-phi(i+1))./r;
% %     vel=(phib(i)-phib(i+1))./r+Q0t(i);
%      vel=(phib(i)-phib(i+1))./r-(Q0t(i+1)+Q0t(i))/2;  % TODO MANU REMBER THAT THIS IS VERY IMPORTANT !!!!! IF You don't use the mean, results are a bit shifted
%      velDoubletPrincipalb(i)=(phib(i)-phib(i+1))./r; 
%     velsrc=(vSothers(i)+vSothers(i+1))/2;
%     vel=vel;
% %     cppot(i)=1-vel.^2;
% %     cppot2(i)=1-vel2.^2;
%     cpdonly(i)=1-vel.^2;
% end; 


%% Copied again!!!!!!!!

for i=1:m;
    phi(i)=U0*(co(i,1).*cos(al)+co(i,2).*sin(al))+g(i);
    phiU0(i)=U0*(co(i,1).*cos(al)+co(i,2).*sin(al));
end; i=m+1;

for i=1:m-1;
    r(i)=(dl(i+1)+dl(i))./2;
    velD(i)=-(g(i)-g(i+1))./r(i);
    velU0(i)=-(phiU0(i)-phiU0(i+1))./r(i);
    vel(i)=(phi(i)-phi(i+1))./r(i);
    cpdonly(i)=1-vel(i)^2/U0^2;
end;
velD(m)=velD(m-1);
velU0(m)=velU0(m-1);
% load('s2_dump');
% kbd

%% Velocity formulations
ptvel=pt2*0;
for i=1:m;
    if(i~=1&&i~=m)
        % normal case
        r=sqrt((co(i+1,1)-co(i-1,1)).^2  +(co(i+1,2)-co(i-1,2)).^2);
        vDprincipal(i)=(g(i+1)-g(i-1))./r;
%         r=sqrt((co(i+1,1)-co(i,1)).^2  +(co(i+1,2)-co(i,2)).^2);
%         vDprincipal(i)=(g(i+1)-g(i))./r;
%         r=sqrt((co(i+1,1)-co(i,1)).^2  +(co(i+1,2)-co(i,2)).^2);
%         vDprincipal(i)=(g(i+1)-g(i))./r;
        v0t(i)=U0*(cos(al).*cos(th(i-1))+sin(al).*sin(th(i-1))+cos(al).*cos(th(i+1))+sin(al).*sin(th(i+1)))/2 ; 
        ptvel(i,:)=co(i,:);
    elseif(i==1) ;
        r=sqrt((co(2,1)-co(1,1)).^2 +(co(2,2)-co(1,2)).^2);
        vDprincipal(i)=(g(2)-g(1))./r;
        v0t(i)=U0*(cos(al).*cos(th(i))+sin(al).*sin(th(i))+cos(al).*cos(th(i+1))+sin(al).*sin(th(i+1)))/2 ; 
        ptvel(i,:)=(co(2,:)+co(1,:))/2;
    elseif(i==m) ;
        r=sqrt((co(m,1)-co(m-1,1)).^2 +(co(m,2)-co(m-1,2)).^2);
        vDprincipal(i)=(g(m)-g(m-1))./r;
        v0t(i)=U0*(cos(al).*cos(th(m))+sin(al).*sin(th(m))+cos(al).*cos(th(m-1))+sin(al).*sin(th(m-1)))/2 ; 
        ptvel(i,:)=(co(m,:)+co(m-1,:))/2;
    end;
     v0t(i)=U0*(cos(al).*cos(th(i))+sin(al).*sin(th(i)));

end; 
ptvel=pt2;
% vDprincipalinterp=interp1(pt2(:,1),vDprincipal,co(:,1));
vDprincipalinterp=vDprincipal;
% vDprincipalinterp(1:end-1)=-velDoubletPrincipalb; % hack
vDothers2=[ (vDothers(2:end)+vDothers(1:end-1))/2 vDothers(end)  ];
vSothers2=[ (vSothers(2:end)+vSothers(1:end-1))/2 vSothers(end)  ];


% vel=v0t + 1*vDothers + vDprincipalinterp./2 +1*vSothers2;
vel=velU0 + 1*vDothers + velD./2 +1*vSothers;
cpall=1-(vel).^2/U0^2;




if nargin==0
    figure,plot(1:length(GammaDS),GammaDS,'k',1:length(g),g,'+',1:length(gpot),gpot,'.',1:length(GammaD),GammaD,'k--')
    legend('Gamma DS pot','g current','gpot current','Gamma D')
    figure,hold all,plot(rhs,'k-'),plot(rhsU0,'r.'),plot(rhsvel,'o'),plot(Svel_n_others,'g'),plot(Svel_n_self,'b'),plot(aDvel*g,'+'),title('rhs'),legend('rhs','rhsU0','rhsvelS','Svelothers','Svelself','a*g')
    figure,hold all,plot(aDvel*g-rhs),title('rhs'),legend('rhs','a*g')
    figure,hold all
    % plot(co(1:m,1),v0t)
    % plot(co(1:m,1),vSothers)
    % plot(co(1:m,1),vDothers)
    % plot(co(1:m,1),vDprincipal/2)
    plot(v0t)
    plot(vSothers)
    plot(vDothers)
    plot(vDprincipal/2)
    legend('U0','Sources','Doublets Neighbor','Doublet Principal');
    
    figure,hold all
    plot(x_ref,Cp_ref,'k'),axis ij
    plot(pt2(1:end-1,1)/chord+0.5,cpdonly,'+'),axis ij  % this on has only doublet contrib, has to be wrong...
    % plot(pt2(1:end-1,1),cp2(1:end-1),'.'),axis ij
    % plot(co(1:end-1,1),cp2(1:end-1),'.'),axis ij
    plot(ptvel(1:end-1,1)/chord+0.5,cpall(1:end-1),'k+'),axis ij
    % plot(pt2(1:end-1,1),cppot,'d'),axis ij
%     plot(pt2(1:end-1,1),cppot2,'d'),axis ij
%     xlim([0 1])
%     ylim([-10 1])


    keyboard
end




% ----- Outputs
chord=max(AirfoilPoints(:,1))-min(AirfoilPoints(:,1));
Gamma_tot=g(end);
Cl=2*Gamma_tot/(U0*chord); % check that, chord??
Mu=g; % used to be 1:end-2..
AI=aDvel;
% due to the finite difference, the solution is well approximated at panel points
cpdonly=cpdonly;
cpall=cpall(1:end-1);
pt=pt2(1:end-1,1);


%% END Script 2DDoubletSource
