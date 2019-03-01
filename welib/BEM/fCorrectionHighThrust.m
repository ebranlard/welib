function [ a CT_loc] = fCorrectionHighThrust(CTCorrection,a, Cn ,phi, a_last,sigma, F,Ftip, CT_loc )
    switch CTCorrection
        case 'Shen'
            ac=1/3;
            Iac=a>ac;
            a(Iac)=(CT_loc(Iac)/4*Ftip(Iac)-Ftip(Iac)*ac^2)./(1-2*ac*Ftip(Iac)) ;
            %%% Glauert correction
        case 'Glauert'
            ac=0.3;
            Iac=a>ac;
            A=sigma.*Cn./sind(phi)^2;
            for i=1:length(a)
                a(i)=fzero(@(aa) -A(i)+aa*(4*F(i)+2*A(i))+aa.^2*(-5*F(i)-A(i))+3*F(i)*aa.^3    ,[0 1]);
            end
            %%% Glauert Exact correction
        case 'GlauertExact'
            ac=0.3;
            error();
            if a>ac
                A=sigma(e)*Cn/sind(phi)^2;
                asolutions=fGlauertSolutions(F,A);
                a=asolutions(whichmin(abs(asolutions-a_last)));
            end
            %%% Glauert correction REQUIRES RELAXATION
        case 'GlauertCT'
            ac=0.3;
            Iac=a>ac;
            fg=0.25*(5-3*a(Iac));
            a(Iac)=CT_loc(Iac)./(4*F(Iac).*(1-fg.*a(Iac)));
            %SperaExact correction
        case 'Spera'
            ac=0.34;
            Iac=a>ac;
            K=4*F(Iac).*(sind(phi(Iac))).^2./(sigma(Iac).*Cn(Iac));
            a(Iac)=0.5*(2+K*(1-2*ac)-sqrt((K*(1-2*ac)+2 ).^2 + 4*(K*ac^2-1)    )  );
            %Spera correction REQUIRES RELAXATION
        case 'SperaCT'
            ac=0.34;
            Iac=a>ac;
            fgs=ac./a(Iac)*(2-ac./a(Iac));
            a=CT_loc(Iac)./(4*F(Iac).*(1-fgs.*a(Iac)));
            %WE handbook correction
        case 'HandbookCT'
            Ict=CT_loc>0.96;
            a(Ic)=1./F(Ic).*(0.143 + sqrt( 0.0203-0.6427 *(0.889-CT_loc(Ic) ) ));
            %Aerodyn correction, REQUIRES RELAXATION
        case 'AeroDyn'
            %         error();
            CT_loc=min(max(-2.0,CT_loc),2);
            if CT_loc>0.96*F    %
                a=0.1432+sqrt(-0.55106+0.6427*CT_loc./F);
            else
                a=0.5*(1-sqrt(1-CT_loc./F));
            end
        case 'Hawc'
            k = [-0.0017077,0.251163,0.0544955,0.0892074]; % glauert correction
            CT_loc = CT_loc./Ftip;
            a = k(4)*CT_loc.^3+k(3)*CT_loc.^2+k(2)*CT_loc+k(1);
        case 'none'
            %% 
            a=a;
        otherwise
            error('High thrust correction unknown')
    end
    if(imag(a)~=0)
        warning('crash');
        %     keyboard
        a=real(a);
    end
end %end  function


function [ as ] = fGlauertSolutions(F,A )
    a1 = (-sqrt(3)*sqrt(-1)/2-1/2)*(sqrt(368*F^3+12*A*F^2+87*A^2*F-8*A^3) ...
        /(2*3^(7/2)*F^(3/2)) ...
        +(-290*F^3+15*A*F^2-24*A^2*F+2*A^3)/(1458*F^3)) ...
        ^(1/3) ...
        +(sqrt(3)*sqrt(-1)/2-1/2)*(-11*F^2-8*A*F+A^2) ...
        /(81*F^2 ...
        *(sqrt(368*F^3+12*A*F^2+87*A^2*F-8*A^3)/(2*3^(7/2)*F^(3/2)) ...
        +(-290*F^3+15*A*F^2-24*A^2*F+2*A^3)/(1458*F^3)) ...
        ^(1/3))+(5*F+A)/(9*F);

    a2 = (sqrt(3)*sqrt(-1)/2-1/2)*(sqrt(368*F^3+12*A*F^2+87*A^2*F-8*A^3) ...
        /(2*3^(7/2)*F^(3/2)) ...
        +(-290*F^3+15*A*F^2-24*A^2*F+2*A^3)/(1458*F^3)) ...
        ^(1/3) ...
        +(-sqrt(3)*sqrt(-1)/2-1/2)*(-11*F^2-8*A*F+A^2) ...
        /(81*F^2 ...
        *(sqrt(368*F^3+12*A*F^2+87*A^2*F-8*A^3)/(2*3^(7/2)*F^(3/2)) ...
        +(-290*F^3+15*A*F^2-24*A^2*F+2*A^3)/(1458*F^3)) ...
        ^(1/3))+(5*F+A)/(9*F);
    a3 = (sqrt((368*F^3+12*A*F^2+87*A^2*F-8*A^3)/F)/(2*3^(7/2)*F) ...
        +(-290*F^3+15*A*F^2-24*A^2*F+2*A^3)/(1458*F^3))^(1/3) ...
        +(-11*F^2-8*A*F+A^2)/ ...
        (81*F^2 *(sqrt((368*F^3+12*A*F^2+87*A^2*F-8*A^3)/F) ...
        /(2*3^(7/2)*F) ...
        +(-290*F^3+15*A*F^2-24*A^2*F+2*A^3) ...
        /(1458*F^3)) ...
        ^(1/3))+(5*F+A)/(9*F);
    as=[real(a1) real(a2) real(a3)];

end

