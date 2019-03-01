function [ KP ] = fCirculationPrandtl( lambda,B,r_bar )
    % this is not the circulation properly speaking since it does not have dimension
    lambda_r=lambda*r_bar(:)';
    l_bar=1/lambda;
    F = fTipLossPrandtl(lambda,B,r_bar);
    KB=fCirculationBetz(l_bar,r_bar);
    KP = F(:)'.*KB(:)';
end
