function [ eta Ftot Mtot Fdrag Finertia vz Stored_u Stored_dudt dFtot dMtot dFdrag dFinertia]=  fHydroCalcFinal(vt,q,vf,vk,vA,vphases,g,rhow,h,zBot,D,Cm,CD,nz,bWheeler,bFloating,bVzTo0)
% This function does the summation over different sinusoidal waves components and computes: elevation, total load and moment

%% Inputs
% % ------ Time
% vt=0 ; % either scalar or a vector.. but watch out q for now is scalar
% % ------- Position
% q(1)=0;% x0
% q(2)=0;% theta
% q(3)=0;
% q(4)=0;
% % ------- Component related inputs
% vf= ; %: vector of frequencies of the different wave component
% vk= ; % vector of wave number as returned by fgetDisperion
% vA=6 ; % wave amplitudes ( the wave goes from +-A) !!! NOTE
% vphases= ;% wave phases
% % ------- Environment
% g=9.81; % gravity
% rhow = 1025; % kg/m3 %rho water!!! NOTE
% % ------- Cylinder properties
% h=10; : water depth % positive !!!! NOTE
% D=1  ;% cylinder diameter
% zBot=-5; % depth of floater negative!!!! NOTE
% Cm= 1.0; %  !!!!!!!!!!!!!!!!!!!!!!!!!!!! This is Cm not CM, CM=2 Cm=1 NOTE
% CD= 0.60; % Drag coeff
% % ------- Algo
% nz=50;  % number of computation points over the depth
% bWheeler=0; % if true (1) use wheeler stretching
% bFloating=1; % if true (1) the moment is computed at z=0 and only from zBot and not from -h
% bVzTo0=1; % if true (1) the vector vz goes to 0 and not to eta

%% Init
Ncomp=length(vf);
vOmega=2*pi*vf;
%% Allocations
nt=length(vt);
eta=zeros(1,nt);
Ftot=zeros(1,nt);
Mtot=zeros(1,nt);
Stored_u=zeros(1,nt);
Stored_dudt=zeros(1,nt);
Fdrag=zeros(1,nt);
Finertia=zeros(1,nt);


%%  Degree of freedom
x0=q(1); theta=q(2); xDot=q(3); thetaDot=q(4);
x=x0;

if bWheeler && bFloating
    error('Verify formula')
end

%% Computing water elevation -> useless, only for wheeler
for it=1:nt % time loop
    t=vt(it);
    eta(it) =  sum(vA.*cos(vOmega*t - vk*x +vphases));
    
    % Vector over depth
    if bFloating    
        vz=linspace(zBot,eta(it),nz); % for floating
        if bVzTo0
            vz=linspace(zBot,0,nz); % for floating
        end
    else
        vz=linspace(-h,eta(it),nz); % for monopile
        if bVzTo0
            vz=linspace(-h,0,nz); % for monopile
        end
    end



    % Wheeler Correction for one time step, workes for monopile, not floating
    %
    if bWheeler
        vzc=(vz-eta(it))*h/(h+eta(it)); % zcalc slide 17 lecture 2
    else
        vzc=vz;
    end


    % Now, loop on depth to get loads and moments
    u   =zeros(1,length(vz));
    dudt=zeros(1,length(vz));
    if sum(abs(vA))~=0
        for iz=1:nz % loop on depth
            %     x=x0+vz(iz)*theta; % !!!!!!!!!!!!!1 we assume small displacement, so we'll keep
            for ip=1:Ncomp % component loop
                comp_u   =  vOmega(ip)  *(cosh(vk(ip)*(vzc(iz)+h))/sinh(vk(ip)*h)) * vA(ip) * cos( vOmega(ip) *t - vk(ip) * x + vphases(ip) );
                comp_dudt= -vOmega(ip)^2*(cosh(vk(ip)*(vzc(iz)+h))/sinh(vk(ip)*h)) * vA(ip) * sin( vOmega(ip) *t - vk(ip) * x + vphases(ip) );
                u(iz)   =u(iz)   +comp_u;    % summing all components for a given time
                dudt(iz)=dudt(iz)+comp_dudt; % summing all components for a given time
            end
        end
    end
    % Loads in N/m at different z locations
    dFdrag=0.5*rhow*CD*D* abs(u-xDot).*(u-xDot);            % vectorial computation over depth
    dFinertia= rhow * (Cm+1) * pi * (D/2)^2 * dudt; % idem %!!!!!!!!!!!! Watch out between Cm and CM
    dFtot=dFdrag+dFinertia;

    if bFloating
        dMtot=dFtot.*(vzc);  % the moment is computed aroung z=0
        vz_int=vzc; 
    else 
        vz_int=vzc+h;   
        dMtot=dFtot.*(vzc+h);  % the moment is computed from the sea bed, where z=-h
    end

    % % integration over depth 
    Ftot(it)=trapz(vzc,dFtot);
    Fdrag(it)=trapz(vzc,dFdrag);
    Finertia(it)=trapz(vzc,dFinertia);
    Mtot(it)=trapz(vz_int,dMtot);

    % we store the wave velocity at z=0
    Stored_u(it)=u(end);
    Stored_dudt(it)=dudt(end);

end
