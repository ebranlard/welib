function [DOFindex]=fElementDOFIndex(iel,nnel,ndof)
% Compute system dofs associated with each element in one- dimensional problem
%
% INPUTS:
%   DOFindex - system dof vector associated with element "iel"
%   iel - element number whose system dofs are to be determined
%   nnel - number of nodes per element
%   ndof - number of dofs per node 
edof  = nnel*ndof            ;
iStart = (iel-1)*(nnel-1)*ndof;
DOFindex=iStart+(1:edof);


