function [fc,Mr,Kr,Fr,Dr,Psic,Omegac,Is,Psi_sm] = fCraigBamptonReduce(Im,nModesCB,MM,KK,F,DD)
% Performs the CraigBampton reduction of a system given some input master dofs index
% and a number of modes
%
% INPUTS
%   Im      : index of master nodes
%   nModesCB: number of CB modes to keep
%   MM, KK  : Mass and stiffness matrix
%
% INPUTS (Optional)
%   nModesCB: number of CB modes to keep
%
% OUTPUTS
%   fc: critical frequency
%   Mr,Kr,Fr,Dr: reduced mass, stiffness, force and damping  matrices
%
% AUTHOR: E. Branlard
%
%
%% Test case
if nargin==0
    L = 100                  ;
    EI= 1.868211939147334e+12;
    Mass = L*8.828201296825122e+03;
    KK=EI/(L^3)*[12      6*L   -12       6*L;
    6*L  4*L^2 -6*L   2*L^2;
    -12     -6*L    12      -6*L;
    6*L  2*L^2 -6*L   4*L^2];
    MM=Mass/420*[156      22*L   54       -13*L;
    22*L  4*L^2  13*L  -3*L^2;
    54       13*L   156      -22*L;
    -13*L -3*L^2 -22*L   4*L^2];
   [fc,Mr,Kr] = fCraigBamptonReduce(3,2,MM,KK);
    %% --- Solve EVA
    [~,Lambda]=eig(Kr,Mr);
    [Omega2,~]=sort(diag(Lambda));
    f= sqrt(Omega2)/(2*pi);
    for i=1:min(8,size(Mr,1))
        fprintf('f%d=%8.3f  Rayleigh Ratio=%.5f\n',i,f(i),  (f(i)/fc)^2 );
    end
    return
end
% --- Optional arguments
if ~exist('F','var') ; F  = []; end;
if ~exist('DD','var'); DD = []; end;
Fr=[];
Dr=[];
fc=NaN;


% Slaves are complementary to masters
Is = setdiff(1:size(MM,1),Im);

% Rearranging - NOTE: master will be first nodes in reduced matrix Mr and Kr
Mmm=MM(Im,Im); Mms=MM(Im,Is); Mss=MM(Is,Is);
Kmm=KK(Im,Im); Kms=KK(Im,Is); Kss=KK(Is,Is);
Kss1Ksm=Kss\(Kms');
Psi_sm = -Kss1Ksm;

% --- Solve EVP for constrained (c) system
[Psic,Lambdac]=eig(Kss,Mss);
% Sorting by ascending frequencies
[Omega2c,Isort]=sort(diag(Lambdac));
Psic    = Psic   (:,Isort);
Lambdac = Lambdac(:,Isort);
if ~isempty(Omega2c)>0
    fc= sqrt(Omega2c(1))/(2*pi);
end
% --- Taking only thefirst few modes
Psic    = Psic   (:,1:nModesCB);
Lambdac = Lambdac(:,1:nModesCB);
Omegac  = sqrt(Omega2c(1:nModesCB));
nm=length(Im);

% --- Using the T matrix:
% T=[eye(nm)  zeros(nm,nModesCB); -Kss1Ksm   Psic];
% MM=[Mmm Mms; Mms' Mss];
% KK=[Kmm Kms; Kms' Kss];
% Mr=T' * MM * T;
% Kr=T' * KK * T;

% --- Using handmade formulae
Mr11=Mmm-(Kss1Ksm')*Mms' - Mms*Kss1Ksm + (Kss1Ksm')*Mss*Kss1Ksm;

Kr11=Kmm-Kms*Kss1Ksm;
Mr12=(Mms-(Kss1Ksm')*Mss)*Psic;

% Building reduced matrix - NOTE: masters are on top
Mr=[Mr11 Mr12; Mr12' eye(nModesCB)];
Kr=[Kr11 zeros(nm,nModesCB); zeros(nModesCB,nm)  Lambdac(1:nModesCB,:)];


if ~isempty(DD)
    error('Not done')
end
if ~isempty(F)
    error('Not done. CB assumes no force in modes')
end


