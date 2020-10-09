function [fc,Mr,Kr,Fr,Dr] = fGuyanReduce(Im,MM,KK,F,DD)
% Performs the Guyan reduction of a system given some input master dofs index

% INPUTS
%
%
% AUTHOR: E. Branlard
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
   [fc,Mr,Kr] = fGuyanReduce([1 3],MM,KK);
    %% --- Solve EVA
    [Q,Lambda]=eig(Kr,Mr);
    [Omega2,Isort]=sort(diag(Lambda));
    Q=Q(:,Isort);
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

% Rearranging
Mmm=MM(Im,Im); Mms=MM(Im,Is); Mss=MM(Is,Is);
Kmm=KK(Im,Im); Kms=KK(Im,Is); Kss=KK(Is,Is);
Kss1Ksm=Kss\(Kms');
% Kss1=Kss^(-1);
% Kss1Ksm=Kss^(-1)*Kms';
% T=[eye(size(Mmm)) ; -Kss1Ksm ];

% MM=[Mmm Mms; Mms' Mss];
% KK=[Kmm Kms; Kms' Kss];
% MM=[MM(Im,Im) MM(Im,Is); MM(Is,Im) MM(Is,Is)]
% KK=[KK(Im,Im) KK(Im,Is); KK(Is,Im) KK(Is,Is)]
% Mr=T' * MM * T;
% Kr=T' * KK * T;

Kr=Kmm-Kms*Kss1Ksm;
Mr=Mmm-(Kss1Ksm')*Mms' - Mms*Kss1Ksm + (Kss1Ksm')*Mss*Kss1Ksm;

if ~isempty(DD)
    Dmm=DD(Im,Im); Dms=DD(Im,Is); Dss=DD(Is,Is);
    Dr=Dmm-(Kss1Ksm')*Dms' - Dms*Kss1Ksm + (Kss1Ksm')*Dss*Kss1Ksm;
%     DD=[Dmm Dms; Dms' Dss];
%     Dr=T' * KK * T;
end
if ~isempty(F)
    Fm=F(Im);
    Fs=F(Is);
    Fr= Fm - Kss1Ksm*Fs;
    %Fr=T' * F;
%     Fr=T' * F;
end


% Solve EVP for constrained system
[~,Lambdac]=eig(Kss,Mss);
[Omega2c]=sort(diag(Lambdac));
if length(Omega2c)>0
    fc= sqrt(Omega2c(1))/(2*pi);
end
