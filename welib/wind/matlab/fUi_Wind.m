function [Uind, Grad]=fUi_Wind(CPs, Time, Wind, Algo)

% CPs np x 3 or np x 2, checks are done by the callee!!!

%% Init
[nCps, ndim]=size(CPs);
Uind=zeros(nCps,ndim);
Grad=zeros(nCps,ndim*ndim);
if Algo.bComputeGrad
else
    Grad=[];
end

%%
if length(findstr('Constant',Wind.Model))>0
    % Constant Wind
    for id=1:ndim
        Uind(:,id)=Wind.V0(id);
    end
else
    warning('To redo')
    keyboard
end
