%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3D plot of induced velocities
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization
Init
InitBEM
load 'WEquilibria10.mat'

%% Parameters for the algorithm
V0=[10];
nturns=10;
deltadeg=5;
dt=pi/180/omega*deltadeg;             %time step [s]
Tstop=2*pi/omega*nturns;           %simulation time  [s]
Vtime=0:dt:(Tstop);
w_guess=-2.5;   %initial guess for induced velocity [m/s]
Model='Constant'; figoff=600;
BigStorage=1;
% Steady parameters
tilt=0;
cone=0;
% Unsteady parameters
Vpitch_of_t=0;
VV0_of_t(1,:)=[0 0 V0];
%% 
yaw=20; Wguess=Wyaw_saved; W0guess=W0yaw_saved;figoff=1212
%yaw=0; Wguess=W_saved; W0guess=W0_saved;figoff=1200
YawModel=1;
Vyaw_of_t=yaw;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
UnsteadyBEM;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% figure 
% hold on
% plot(Vtime,P,'b')



%% Index preparation
m(1)=1;
m(2)=floor((360-120)/deltadeg+1);
m(3)=floor((360-240)/deltadeg+1);
Vpsi=Mpsi(1:floor(360/deltadeg),1);
Ib1=m(1):floor(360/deltadeg)*(nturns-1);
Ib2=Ib1+m(2)-1;
Ib3=Ib1+m(3)-1;
npsi=length(Vpsi);

ne2=ne;

%%%%%%%%
npsi2=floor(npsi*2/3);


Wnn=zeros(npsi2,nB,ne2);
Wnnn=zeros(npsi2,ne2);
Y=zeros(npsi2,ne2);
X=zeros(npsi2,ne2);
Ipsi=(0:(nturns-2))*(360/deltadeg);


j=1;
% averaging
for idPsi=1:npsi2
    for e=1:ne2
        for idB=1:nB
            Wnn(idPsi,idB,e)=mean( MW(Ipsi+m(idB)+(idPsi-1),3,idB,e   ) );
        end
        X(idPsi,e+1)=r(e)*cosd(Vpsi(idPsi));
        Y(idPsi,e+1)=r(e)*sind(Vpsi(idPsi));
        Wnnn(idPsi,e+1)=mean(Wnn(idPsi,1:3,e));
        Supercat(j,1:3)=[X(idPsi,e) Y(idPsi,e) Wnnn(idPsi,e)];
        j=j+1;
    end
end
X(:,1)=9*r(1)/10*cosd(Vpsi(1:npsi2));
Y(:,1)=9*r(1)/10*sind(Vpsi(1:npsi2));
Wnnn(:,1)=0*idPsi;

[Points Faces Colrs]=get3DWindTurbine(yaw/4,3);

%% 
figure(figoff)
hold on
colormap(hot)
cmap = colormap;
colormap(cmap([60 60 60 60 60 60 60 60 60 60 60 60 60 60 60 60 60 60 60 60 60 60 60 60 60 60 60 60 60 60 60 60:-1:1],:))
meshz(X,Y,-Wnnn)
surf(X,Y,-Wnnn,'FaceColor','interp','FaceLighting','phong')
colorbar
xlabel('x_3 [m]')
ylabel('y_3 [m]')
zlabel('Induced velocity normal to the rotor [m/s]')
xlim([-R R])
ylim([-R R])
%zlim([0 max(max(-Wnnn))])

grid on

%%

figure(figoff+1)
surface(Y(:,2:end),X(:,2:end),Wnnn(:,2:end))
colorbar
colormap(hot)
cmap = colormap;
colormap(cmap(1:60,:))
xlim([-R R])
ylim([-R R])
xlabel('y_3 [m]')
ylabel('x_3 [m]')


shading interp;
%colormap(pink);
grid on;
box on;
%axis equal;
%view([35, 30]);
patch('Vertices',Points(:,[2 1 3]),'Faces',Faces,'FaceColor',[0.7 0.7 0.7],'LineWidth',0.1)


%% Plotting
% 
% %% Mongfg
% xrange = linspace(-R,R,20);
% yrange = linspace(-R,R,20);
% [XI YI]=meshgrid(xrange,yrange);
% [XI,YI,ZI] = griddata(Supercat(:,1),Supercat(:,2),Supercat(:,3),XI,YI);
% figure
% surface(XI,YI,ZI);
% 
% %% dsf
% figure
% hold on
% surface(X,Y,Wnnn)
% quiver3(XI,YI,0*ZI,0*XI,0*YI,ZI/3,0.5);
% 
% 
% %% fdkkjgk
% 
% [U,V,W] = surfnorm(XI,YI,ZI);
% figure
% quiver3(XI,YI,ZI,U,V,W,0.5);
% 
% %% jhjh
% 
% clc
% YY=zeros(npsi2,ne2,1);
% XX=zeros(npsi2,ne2,1);
% WWn=zeros(npsi2,ne2,1);
% XX(:,:,1)=X;
% YY(:,:,1)=Y;
% WWn(:,:,1)=Wnnn;
% XX(:,:,2)=X;
% YY(:,:,2)=Y;
% WWn(:,:,2)=Wnnn;
% UU=0*XX;
% 
% 
% 
% %%
% figure
% xrange = linspace(-R,R,8);
% yrange = linspace(-R,R,8);
% zrange = 0;
% [cx cy cz] = meshgrid(xrange,yrange,zrange);
% hcones = coneplot(XX,YY,UU,UU,UU,WWn,cx,cy,cz,5);
% set(hcones,'FaceColor','red','EdgeColor','none')

