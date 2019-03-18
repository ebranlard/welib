function [slope]=fAlphaSlope(alpha0,alpha,Cl)
% returns the slope dCl/dalpha using two points around alpha0
[~,i]=min(abs((alpha-alpha0)));
if i==1
    kbd
end
i1=i-1;
i2=i+1;
slope = (Cl(i2)-Cl(i1))/(alpha(i2)-alpha(i1));

