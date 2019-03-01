function [sse Fitted_Curve]=fFitChord(params,x,y)
a=params(1);
b=params(2);
c=params(3);

Fitted_Curve=a./(x+b)+c;
Error_Vector=Fitted_Curve - y;
% minimize is the sum of squares error
sse=sum(Error_Vector.^2);
