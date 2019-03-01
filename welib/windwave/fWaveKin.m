function [ eta Stored_u Stored_dudt]=  fWaveKin(vt,vf,vk,vA,vphases,x,h,zBot,nz,bWheeler,bFloating,bVzTo0)
% This function does the summation over different sinusoidal waves components and computes: elevation, total load and moment

%% Inputs
% % ------ Time
% vt=0 ; % either scalar or a vector.. but watch out q for now is scalar
% % ------- Component related inputs
% vf= ; %: vector of frequencies of the different wave component
% vk= ; % vector of wave number as returned by fgetDisperion
% vA=6 ; % wave amplitudes ( the wave goes from +-A) !!! NOTE
% vphases= ;% wave phases
% % ------- Cylinder properties
% h=10; : water depth % positive !!!! NOTE
% zBot=-5; % depth of floater negative!!!! NOTE
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
    % we store the wave velocity at z=0
    Stored_u(it)=u(end);
    Stored_dudt(it)=dudt(end);

end
