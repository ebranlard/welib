function Ydot = fTowerTopYDot( t, Y, M_L, M_U, K, D, iTT, fTT, Ke_soil, Isoil)
    % Y    = [ x    ; xdot  ] ;
    % Ydot = [ xdot ; xddot ] ;
    % xddot = fAcc;
    nDOF=size(K,1);
    gz = Y(     1:nDOF  );
    gzp= Y(nDOF+1:2*nDOF);

    [~,~,~,F] = fTowerTopMDKR( t, gz, gzp, [], K, D, iTT, fTT, Ke_soil, Isoil);

    RHS  = -K*gz -D*gzp +F;
    gzpp = M_U\(M_L\RHS)  ;
    %gzpp = MM\RHS  ;
    Ydot = [gzp ; gzpp ];

end
