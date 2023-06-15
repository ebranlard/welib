function [ KB ] = fBetzCirculation(l_bar,Vx)
    % This function does not return the circulation properly speaking but the dimensionless circulation -> consider calling it factor
    % NOW ALL PARAMETERS ARE Dimensionless
    % l_bar=l/R,  sometimes l_bar=1/lambda so l=R/lambda
    if max(Vx)>1
        error('Vector should be normalized to one - dimensionless')
    end
    KB=Vx(:)'.^2./(l_bar^2+Vx(:)'.^2) ;
end

