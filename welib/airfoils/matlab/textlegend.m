function [ h ] = textlegend( txt,Position )
%TEXTLEGEND Summary of this function goes here
%   Detailed explanation goes here
if nargin==1
    Position='North-East';
end

a = axis;
wdth = a(2)-a(1);
ht = a(4)-a(3);
switch Position
    case 'North-East'
        pos = [a(2)-0.05*wdth a(4)-0.05*ht];
        h = text(pos(1),pos(2),txt,'HorizontalAlignment','right','VerticalAlignment','top') % correspond to Nort-East
    case 'North-West'
        pos = [a(1)+0.05*wdth a(4)-0.05*ht];
        h = text(pos(1),pos(2),txt,'HorizontalAlignment','left','VerticalAlignment','top') % correspond to Nort-West
    otherwise
        error('not done')
end

