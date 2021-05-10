%% Documentation
% FEA to calculate Natural frequencies and modes of Cantilevered uniform beam
%% Init
clear all;  close all;
addpath('../BeamTheory');
% addpath('C:/Svn/WindLoadTeam/Matlab/Tools/lib_Plot/');

%% Main Params 
nel        = 150; % Number of Elements
bPlot      = 1  ;
nModesPlot = 12 ;
nModesSelect = 15 ;
bNoTension = 0  ; % NOTE: dangerous, results may depend on number of elements;
bNoTorsion = 0  ; % NOTE: dangerous, results may depend on number of elements;
bNoBending = 0  ; % NOTE: dangerous, results may depend on number of elements;

%% Structural Data
% --- Tube
L   = 100                      ;
G   = 79.3e9                   ;% Shear modulus. Steel: 79.3  [Pa] [N/m^2]
E   = 210e9                    ;% Young modulus [Pa] [N/m^2]
D   = 8                        ;
t   = 0.045                    ;
A   = pi*( (D/2)^2 - (D/2-t)^2); % [m^2]
Iy  = pi/64*(D^4-(D-2*t)^4)    ; % [m^4]
rho = 7850                     ;
m   = rho*A                    ;


% --- Torsion 
% Kt: Saint-Venant's torsion constant, or Torsion constant
% Torsion equation:    Ix/3A Mass txddot  + G Kt /L tx = Torsional Moment
Ix = pi/32*(D^4 - (D-2*t)^4); % Polar second moment of area [m^4]
if bNoTorsion
    Kt=1e10;
else
    Kt=pi/64*(D^4 - (D-2*t)^4); % Torsion constant
end

% --- Tension 
if bNoTension
    EA=E*A*10000; % NOTE: HACK to add elongation stiffness
else
    EA=E*A;
end
% --- Tension 
if bNoBending
    Iy=Iy*100000;
end
% --- Bending symmetry
Iz=Iy;

%% Bending and Torsion Theory 
x_th=linspace(0,L,50);
[fb_th,x_th,U_th,V_th,K_th]      = fUniformBeamTheory('transverse-unloaded-clamped-free',E*Iy,rho,A,L,'x',x_th,'norm','tip_norm');
[ft_th,~,Vx_th] = fUniformBeamTheoryTorsion('torsion-unloaded-clamped-free',G,Kt,Ix,rho,A,L,'x',x_th,'norm','tip_norm');
[fl_th,~,Ux_th] = fUniformBeamTheoryLongi  ('longi-unloaded-clamped-free',E,rho,A,L,'x',x_th,'norm','tip_norm');

% function [freq,x,ModesU,ModesV,ModesK,p] = fUniformBeamTheory(Type,EI,rho,A,L,varargin)


% --- FEM
[Mr,Kr,f,x,Ux,Uy,Uz,Vx,Vy,Vz,MM,KK,sModes,iModes]=fClampedBeamFEM3D_Frame(L,E*Ix,E*Iy,E*Iz,EA,m,A,E,G,Kt,nel);

%% Mode shapes - Normalized to unit tip
P    = zeros(size(Ux,1),3);
P(:,1) = x;
%th=x_th*L;
%% Plotting Mode Shapes
close all;
if bPlot
    vColrs=lines;
    for i=1:nModesPlot
        figure(i),clf,hold all,box on;
        if isequal(sModes{i}(1:2),'uz')
            plot3(x_th/L,0*x_th,U_th(iModes(i),:),'k.');
            title(sprintf('Bending Z %d',iModes(i)));
        elseif isequal(sModes{i}(1:2),'uy')
            plot3(x_th/L,U_th(iModes(i),:),0*x_th,'k.');
            title(sprintf('Bending Y %d',iModes(i)));
        end
        if isequal(sModes{i}(1:2),'vx')
            plot(  x_th/L            ,Vx_th(iModes(i),:),'k.');
            plot( (P(:,1)+Ux(:,i))/L ,Vx(:,i),'-','Color',vColrs(i,:));
            title(sprintf('Torsion %d',iModes(i)));
        elseif isequal(sModes{i}(1:2),'ux')
            plot(  x_th/L,Ux_th(iModes(i),:),'k.');
            plot(P(:,1)/L,Ux(:,i),'-','Color',vColrs(1,:));
            plot(P(:,1)/L,Uy(:,i),'-','Color',vColrs(2,:));
            plot(P(:,1)/L,Uz(:,i),'-','Color',vColrs(3,:));
            plot(P(:,1)/L,Vx(:,i)*180/pi,'-','Color',vColrs(4,:));
            legend('Ux','Uy','Uz','Vx');
            title(sprintf('Longi %d',iModes(i)));
        else
            plot3( (P(:,1)+Ux(:,i))/L , P(:,2)+Uy(:,i), P(:,3)+Uz(:,i),'-','Color',vColrs(i,:));
            xlim([0 1.1])
            ylim([-1 1])
            zlim([-1 1])
            axis equal
            view(3)
        end
    end
    fDispatchFigs(1);
%     legend('Theory','FEM')
end



nModesSelect=min(nModesSelect,size(Mr,2));
for i=1:nModesSelect 
    if isequal(sModes{i}(1:2),'uz') || isequal(sModes{i}(1:2),'uy')
        f_ref=fb_th(iModes(i));
    elseif sModes{i}(1:2)=='vx' 
        f_ref=ft_th(iModes(i));
    elseif sModes{i}(1:2)=='ux' 
        f_ref=fl_th(iModes(i));
    else
        f_ref=NaN;
    end;
    vColrs=lines;
    fprintf('%s - f=%8.3f   -  f=%8.3f   df=%9.5f\n',sModes{i},f(i),f_ref,f_ref-f(i))
end

%%
