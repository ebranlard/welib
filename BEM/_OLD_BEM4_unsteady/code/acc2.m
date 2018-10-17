function g=acc2(t,x,v,Fext,d,M,K,Fargs)
    F=Fext(t,Fargs);
    %g = ( F -K.*x  -d.*v) * inv(M);
    g = inv(M)*(F-K*x-d*v) ;
end
