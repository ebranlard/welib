function [ G ] = fTipLossGoldsteinOkulov( l_bar,B,Vx )
    %l_bar=h/(2piR)
%     G = fCirculationGoldstein( l_bar,w,R,B,Vx );
    K = fGoldsteinFactor( l_bar,B,Vx );
    KB= fCirculationBetz(l_bar,Vx);
%     G = G(:)'./KB(:)';
    G = K(:)'./KB(:)';
end
