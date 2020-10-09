function [Fz]=fTowerTopExcitation(vt)
F0     = 5e6                         ;


Fz   = zeros(size(vt));
b05  = vt<5;
b510 = vt<10 & vt>=5;
b10  = vt>=10;

Fz(b05)  = vt(b05)/5*F0;
Fz(b510) = F0;
Fz(b10)  = 0 ;
