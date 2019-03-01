function [ K ] = fCirculationGoldstein( l_bar,B,Vx )
    % this is not the circulation properly speaking
    %l_bar=h/(2piR)
%     h=l_bar*2*pi*R;
    [ K ] = fGoldsteinFactor( l_bar,B,Vx );
%     Gamma=h*w/B*K(:)';
end
