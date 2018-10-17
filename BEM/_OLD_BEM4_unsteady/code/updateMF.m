function [F,M]=updateMF(m_r,L,M1,dm1,dm2,IM,x,v)
    M=[m_r*L+M1 -dm1*sin(x)
        -dm1*sin(x) IM];
    GF1=(v)^2*cos(x)*dm1;
    GF2=cos(x)*dm2;
    F=[GF1 GF2];
 end
