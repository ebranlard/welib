function [vk]=fgetDispersion(vf,h,g)
% Solves for the dispersion relation once and for all
% h: water depth >0
% g: gravity 9.81

for ip=1:length(vf) % loop on frequencies
    if h<100
         [ vk(ip) ] = fkSolve( vf(ip),h,g );
    else 
        % Not solving the dispersion relation since we are in deep water
        vOmega=2*pi*vf;
        vk=vOmega.^2/g;
    end
end
end


function [ k ] = fkSolve( f,h,g )
warning off

omega=2*pi*f;
eq=@(k) -omega^2+g*k.*tanh(k*h);
% k=fsolve(eq,(omega^2)/g);
k=fzero(eq,[0 100]); % seems a bit faster


warning on
end
