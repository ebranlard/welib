function [kk]=fBuildGlobalMatrix(kk,k,index)
%  Assembly of element matrices into the system matrix
%
% INPUTS
%     kk - system matrix
%     k  - element matri
%     index - d.o.f. vector associated with an element
% 
edof = length(index);
for i=1:edof
    ii=index(i);
    for j=1:edof
        jj=index(j);
        kk(ii,jj)=kk(ii,jj)+k(i,j);
    end
end

