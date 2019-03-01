function [ TargetError ] = fBEMCTFind(WT,Sim,Wind,Algo, Amplitude,CTTarget )
Algo.PrescribedGamma=Algo.PrescribedGamma*Amplitude;
Algo.bPrescribedGamma=1;
% fprintf('Amplitude %4.2f - ',Amplitude);
[ BEM ] = fRunBEM(WT,Sim,Wind,Algo);
TargetError=BEM.CT-CTTarget;
end

