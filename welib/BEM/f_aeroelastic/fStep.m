function [y] = fStep( t,t0,t1,tau,A0,A1,k)
y=A0;
if(t>t0 && t<t1) 
    y=A1;
    if(t<t0+k*tau)
        y=A0+(A1-A0)*(1-exp(-(t-t0)/tau));
    end
    if(t>t1-k*tau)
        y=A0+(A1-A0)*(exp(-(t-t1+k*tau)/tau));
    end
end
    

end

