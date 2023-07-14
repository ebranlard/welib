function [ KB ] = fBetzCirculation(l,R,Vr)
    % All parameters have dimensions [m]
    % l is not l_bar=l/R,  sometimes l_bar=1/lambda so l=R/lambda
    KB=(Vr/R).^2./((l/R)^2+(Vr/R).^2) ;
end

