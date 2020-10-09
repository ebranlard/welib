% --- Sub Function for ode solver
function [M,D,K,F] = fTowerTopMDKR( t, gz, gzp, M, K, D, iTT, fTT, Ke_soil, Isoil)
    % iTT: index of tower top translational DOF
    % fTT: function handle f(t) to get excitation at tower top
    global nit
    nDOF_tot=size(K,1);

    % --- Prescribed force
    [Fz] = fTT(t);

    F=zeros(nDOF_tot,1);
    F(iTT)=Fz; % Tower top DOF get Force

    if ~isempty(Isoil)
        F(Isoil) = F(Isoil) - Ke_soil(3:4,3:4)*gz(Isoil);
    end


    %% --- Output to screen
    nit=nit+1; 
    if mod(nit,100)==1
        fprintf('it=%7d   t=%4.8fs  tip_disp=%9.2e m\n',nit,t,gz(iTT));
        if isnan(gz(iTT))
            error()
        end
    end
end
