function [ s_flat s ] = struct2str( param ,sep)
%STRUCT2STR 
%   Detailed explanation goes here
if nargin==1
    sep=', ';
end

if isempty(param) || ~isstruct(param)
    s_flat='';
    s='';
    return
end
% values=structfun(@(x)num2str(x),param,'UniformOutput',0);
names=fieldnames(param);
% val=num2cell(values);

% couldn't find a good vectorial way to do it... it takes too much time to think of that
s_flat='';
for i=1:length(names)
    s_new=sprintf('%s=%s',names{i},num2str(getfield(param,names{i})));
    s_flat=[s_flat s_new];
    s{i,1}=s_new;
    if i~=length(names)
        s_flat=[s_flat sep];
    end
end

