function [ F ] = fTipLossPrandtl( lambda,B,r_bar )
    lambda_r=lambda*r_bar(:)';
    F = 2/pi*acos(exp(-B/2*sqrt(1+lambda^2)*(1-lambda_r/lambda)));
end
