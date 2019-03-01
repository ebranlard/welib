function [sig,ires]=fHawc2_signal(labels,Results,lab,varargin)
I=1:length(labels);

if(isempty(varargin))
    % then it's a signal that do not depend on radius
    string_searched=sprintf('%s',lab);
else
    r=varargin{1};
    string_searched=sprintf('%s, R=%5.1f',lab,r);
%     string_searched=sprintf('%s, R=%5.2f',lab,r); 
end
ires=I(cellfun(@(x)~isempty(x),strfind(labels,string_searched)));
if length(ires)>1
    warning('More than one signal found for the same string searched in hawc2 labels')
    kbd
elseif length(ires)==0
    % label not found, we return NaN
    sig=nan(size(Results,1),1); 
    ires=[];
elseif length(ires)==1
    % OK, thats good
    sig=Results(:,ires);
else
    error('')
end



% Fortran format for radius:  20 FORMAT (A, F5.1)
% Fortran format for radius:  20 FORMAT (A, F5.2)
