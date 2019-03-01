function [sse Fitted_Curve]=fFitChord(params,x,y)
a=params(1);
b=params(2);
c=params(3);
d=params(4);
e=params(5);
Fitted_Curve=acos( exp(x*a) ).*acos( exp((1-x)*b) ).*sin(c*x+d*x.^2+e);
Error_Vector=Fitted_Curve - y;
% minimize is the sum of squares error
sse=sum(Error_Vector.^2);
