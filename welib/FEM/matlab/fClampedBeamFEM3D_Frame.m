function [Mr,Kr,f,x,Ux,Uy,Uz,Vx,Vy,Vz,MM,KK,sModes,iModes]=fClampedBeamFEM3D_Frame(x0,EIx0,EIy0,EIz0,EA0,m0,A0,E,G,Kt,nel)

% Use FEM to get the modes of a rotating or non-rotating cantilever beam
%
% NOTE: input values can be vectors or scalars.
%       If they are scalars, then a beam with constant properties and of length L=x0 is used;
%       If they are vectors, then linear interpolation is used. The dimension of the inputs does not need to match nel
% 
% INPUTS
%   x0   : (1xn) Span vector of the beam [m]
%   EI0  : (1xn) Elastic Modulus times Second Moment of Area of cross section [Nm2]
%   m0   : (1xn) Mass per length along the beam [kg/m]
%   nel  : Number of elements
%
% OUTPUTS
%   f : (1 x nMl)   Modes Frequencies
%   x : (1 x nel)   Span vector
%   U : (nM x nel)  Modes Shapes Displacements
%   V : (nM x nel)  Modes Shapes Slopes
%   K : (nM x nel)  Modes Shapes Curvature
%
% AUTHOR: E.Branlard 

%% Test Function
if nargin==0
    nel=100; % Number of Elements
    % --- Tube2 - Structural data
    L   = 100                      ;% Beam Length [m]
    G   = 79.3e9                   ;% Shear modulus. Steel: 79.3  [Pa] [N/m^2]
    E   = 210e9                    ;% Young modulus [Pa] [N/m^2]
    D   = 8                        ;
    t   = 0.045                    ;
    A0  = pi*( (D/2)^2 - (D/2-t)^2);
    rho = 7850                     ;
    m0  = rho*A0                   ;
    EIx0 = E*pi/32*(D^4-(D-2*t)^4); % Polar second moment of area [m^4]
    EIy0 = E*pi/64*(D^4-(D-2*t)^4); % Planar second moment of area [m^4]
    EIz0 = EIy0                   ;
    EA0  = E*A0                   ;
    Kt=pi/64*(D^4 - (D-2*t)^4); % Torsion constant [m^4]

    [Mr,Kr,f,x,Ux,Uy,Uz,Vx,Vy,Vz,MM,KK,sModes,iModes]=fClampedBeamFEM3D_Frame(L,EIx0,EIy0,EIz0,EA0,m0,A0,E,G,Kt,nel);

    addpath('../BeamTheory/');
    x_th=linspace(0,L,nel+1);
    [f_th,x_th,U_th] = fUniformBeamTheory('transverse-unloaded-clamped-free',EIy0,rho,A0,L,'x',x_th,'norm','tip_norm');
    %% Comparison of Modes and Frequencies with theory
    vColrs=lines;
    for i=1:6
        figure(i),clf,hold all,box on;
        if isequal(sModes{i}(1:2),'uz')
            plot3(x_th/L,0*x_th,U_th(iModes(i),:),'k.');
            title(sprintf('Bending Z %d',iModes(i)));
            f_ref=f_th(iModes(i));
        elseif isequal(sModes{i}(1:2),'uy')
            plot3(x_th/L,U_th(iModes(i),:),0*x_th,'k.');
            title(sprintf('Bending Y %d',iModes(i)));
            f_ref=f_th(iModes(i));
        end
        if isequal(sModes{i}(1:2),'vx')
            plot( (x'+Ux(:,i))/L ,Vx(:,i)*180/pi,'-','Color',vColrs(i,:));
            title(sprintf('Torsion %d',iModes(i)));
            f_ref=NaN;
        elseif isequal(sModes{i}(1:2),'ux')
            plot(x/L,Ux(:,i),'-','Color',vColrs(1,:));
            plot(x/L,Uy(:,i),'-','Color',vColrs(2,:));
            plot(x/L,Uz(:,i),'-','Color',vColrs(3,:));
            plot(x/L,Vx(:,i)*180/pi,'-','Color',vColrs(4,:));
            f_ref=NaN;
            legend('Ux','Uy','Uz','Vx');
            title(sprintf('Longi %d',iModes(i)));
        else
            plot3( (x'+Ux(:,i))/L, Uy(:,i), Uz(:,i),'-','Color',vColrs(i,:));
            xlim([0 1.1])
            ylim([-1 1])
            zlim([-1 1])
            axis equal
            view(3)
        end
        fprintf('%s - f=%8.3f   -  f=%8.3f   df=%9.5f\n',sModes{i},f(i),f_ref,f_ref-f(i))
    end

    return
end



% --------------------------------------------------------------------------------}
%% --- INPUTS 
% --------------------------------------------------------------------------------{
if length(x0)==1
    % Constant beam properties
    x0   = [0 x0]    ;
    EIx0 = [1 1]*EIx0;
    EIy0 = [1 1]*EIy0;
    EIz0 = [1 1]*EIz0;
    EA0  = [1 1]*EA0 ;
    A0   = [1 1]*A0  ;
    m0   = [1 1]*m0  ;
end




%% Parameters
ndof=6; % Degrees of Freedom per Node
% --------------------------------------------------------------------------------}
%% --- Boundary Conditions (BC)
% --------------------------------------------------------------------------------{
% Clamped-Free BC
% bcdof(1)=1; bcval(1)=0;  % first element deflection u=0
% bcdof(2)=2; bcval(2)=0;  % first element slope v=0

[MM,KK,x]=fBeamMatrices3D_Frame6DOF(x0,EIx0,EIy0,EIz0,EA0,m0,A0,E,G,Kt,nel);


% --- Apply BC
% Clamped-Free:  - Removing 6 first column and row of M and K
Kr=KK;
Mr=MM;
for i=1:ndof
    Kr(1,:)=[]; Kr(:,1) = [];
    Mr(1,:)=[]; Mr(:,1) = [];
end

%% --- Solve EVA
[Q,Lambda]=eig(Kr,Mr);
Omega2=diag(Lambda);
[Omega2,Isort]=sort(Omega2);
f=sqrt(Omega2)/(2*pi);
Q=Q(:,Isort);


%% Mode shapes 
Qz=zeros(6,size(Q,2));
Q=[Qz;Q];
Ux=Q(1:ndof:end,:);
Uy=Q(2:ndof:end,:);
Uz=Q(3:ndof:end,:);
Vx=Q(4:ndof:end,:);
Vy=Q(5:ndof:end,:);
Vz=Q(6:ndof:end,:);


%% Mode detection - Not easy
nModes=size(Q,2);
sModes=cell(1,nModes);
iModes=zeros(1,nModes);
iz=0;
iy=0;
it=0;
ix=0;
for i=1:nModes
    sModes{i}='NA';
    % -- Longi - Late
    UxCum=1/(nel+1)*sum(abs(Q(1:ndof:end,i)));
    UyCum=1/(nel+1)*sum(abs(Q(2:ndof:end,i)));
    UzCum=1/(nel+1)*sum(abs(Q(3:ndof:end,i)));
    VxCum=1/(nel+1)*sum(abs(Q(4:ndof:end,i)));
    [~,iMax]=max([UxCum,UyCum,UzCum,VxCum]);
    switch iMax
        case 1
            % --- Extension/Contraction
            ix=ix+1;
            sModes{i}=['ux' num2str(ix)];
            [~,iMax]=max(abs(Ux(:,i)));
            fact=abs(Ux(iMax,i))*sign(Ux(end,i));
            iModes(i)=ix;
            Q(:,i)=Q(:,i)/fact;
        case {2,3}
            % --- Bending
            if max(abs(Uz(:,i)))>max(abs(Uy(:,i)))
                iz=iz+1;
                sModes{i}=['uz' num2str(iz)];
                [~,iMax]=max(abs(Uz(:,i)));
                fact=Uz(iMax,i);
                iModes(i)=iz;
            else
                iy=iy+1;
                sModes{i}=['uy' num2str(iy)];
                [~,iMax]=max(abs(Uy(:,i)));
                fact=Uy(iMax,i);
                iModes(i)=iy;
            end
            Q(:,i)=Q(:,i)/fact;
        case 4
            % --- Torsion
            it=it+1;
            sModes{i}=['vx' num2str(it)];
            [~,iMax]=max(abs(Vx(:,i)));
            fact=abs(Vx(iMax,i))*sign(Vx(end,i));
            Q(:,i)=Q(:,i)/fact;
            iModes(i)=it;
    end
%     fprintf('f=%8.4f  -  %s   Ux=%12e  Uy=%12e  Uz=%12e  Vx=%12e \n',f(i),sModes{i},UxCum,UyCum,UzCum,VxCum);
end;


%%
Ux=Q(1:ndof:end,:);
Uy=Q(2:ndof:end,:);
Uz=Q(3:ndof:end,:);
Vx=Q(4:ndof:end,:);
Vy=Q(5:ndof:end,:);
Vz=Q(6:ndof:end,:);
