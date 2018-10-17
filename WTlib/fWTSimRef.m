function [R]= fWTSimRef(WT,Opts,Algo)

% --------------------------------------------------------------------------------
% ---  Computing parametric CP(lambda,pitch) curve and finding optimal lambda_opt
% --------------------------------------------------------------------------------
[Parametric, CPmax, lambda_opt, pitch_opt]= fWTParametricCP(WT,Opts.vLambda,Opts.vPitch,Opts.bPlot,Algo);

% --------------------------------------------------------------------------------
% --- Computing RPM(WS) curve given: lambda_opt, max and min RPM
% --------------------------------------------------------------------------------
[vRPMRef,Regions]= fWTGetRefRPM(WT,Opts.vWS_out,Opts.WS_Startup,lambda_opt,Opts.TipSpeedMax,Opts.OmegaMin,Opts.bPlot);

% --------------------------------------------------------------------------------
% --- Find Pitch given the RPM curve
% --------------------------------------------------------------------------------
%%
vPITCHRef= zeros(1,length(Opts.vWS_out)); % Unknown 
WT.Spec.vSIMRef=[Opts.vWS_out(:) vRPMRef(:) vPITCHRef(:)];

R=fWTFindPitch(WT ,Opts.vWS_out,Opts.Pref,Opts.bPlot,Opts.bOptimBelowPref,Algo);

%% Regions


%% Return results
% returning results by adding field to the previous R, which contain allready the final PowerCurve data
R.vSIMRef       = [Opts.vWS_out(:) vRPMRef(:) R.PITCH(:)];
R.CPlambdaPitch = Parametric ;
R.OmegaMax      = Opts.TipSpeedMax/WT.Rotor.R ;
R.lambda_opt    = lambda_opt;
R.pitch_opt     = pitch_opt ;
R.CPmax         = CPmax     ;
R.Regions       = Regions   ;

%% DEBUG
% kbd

%%
% [ Sim ]  = fInitSim(WT);
% [ Sim ]  = fSetSim( Sim, WT, R.vSIMRef ); 
% [ Wind ] = fInitWind(  ); 
% keyboard
% [ R2 ]=fWTSimulation('BEM',WT,Sim,Wind,Algo);


%% Interpolatoin of Power curve
Opts.vWS_out = R.vSIMRef(:,1) ;
vRPM         = R.vSIMRef(:,2) ;
RCPlp        = R.CPlambdaPitch;

WS     = linspace(Opts.vWS_out(1),Opts.vWS_out(end),100);
Omega  = interp1(Opts.vWS_out,vRPM*2*pi/60,WS)          ;
Pitch  = interp1(Opts.vWS_out,R.PITCH,WS)               ;
P      = interp1(Opts.vWS_out,R.Power,WS)               ;
lambda = Omega*WT.Rotor.R./WS                           ;
CPout  = cell2mat({R.PowerCurveData.CP})                ;
CP     = interp1(Opts.vWS_out,CPout,WS)                 ;

[~,i3]= min(abs(P-Opts.Pref));
Regions.WS3 =WS(i3);


%% Plotting
if(Opts.bPlot)
   

    [~,i3]= min(abs(P-Opts.Pref));
    Regions.WS3 =WS(i3);
   
    wsmin=min(WS);
    wsmax=max(WS);
   
    figure
    hold on
    plot(WS,Omega/R.OmegaMax,'-','LineWidth',2 ,'Color',fColrs(3))
    plot(WS,P/(Opts.Pref/1000),'-','LineWidth',2, 'Color',fColrs(1))
    plot(WS,Pitch/max(Pitch),'-','LineWidth',2 ,'Color',fColrs(2))
    plot(WS,CP/CPmax,'k-','LineWidth',2 )
    plot(WS,lambda/lambda_opt,'-','LineWidth',2,'Color',fColrs(4) )
   
    % Regions
    PLOT_RANGE=[0 max(1.135,max(lambda/lambda_opt))];
    plot([Regions.WS0 Regions.WS0],PLOT_RANGE,'k','LineWidth',2 )
    plot([Regions.WS1 Regions.WS1],PLOT_RANGE,'k','LineWidth',2 )
    plot([Regions.WS2 Regions.WS2],PLOT_RANGE,'k','LineWidth',2 )
    plot([Regions.WS3 Regions.WS3],PLOT_RANGE,'k','LineWidth',2 )

    tpos=1.1;
    if Regions.WS0~=Regions.WS1 
        text((Regions.WS0+Regions.WS1)/2-1,tpos,'R0')
    end
    text((Regions.WS1+Regions.WS2)/2-1,tpos,'R1')
    text((Regions.WS2+Regions.WS3)/2-1,tpos,'R2')
    text((Regions.WS3+Regions.WS4)/2-1,tpos,'R3')

    legend('\Omega','Power','Pitch','Cp','Lambda','Location','North','Orientation','Horizontal')
    xlabel('Wind Speed [m/s]')
    ylabel('\Omega/\Omega_{max} , P/P_{ref} , \theta/\theta_{max}, C_P/C_{P,opt} , \lambda/\lambda_{opt}  [.]')
    %grid on
    box on
    ylim(PLOT_RANGE)
    xlim([wsmin wsmax])
    title('ControllerRegionsDimLess')
   
   
    figure
    hold on
    plot(lambda,CP,'k-','LineWidth',2  )
   
    xlabel('\lambda [.]')
    ylabel('Cp [.]')
    %grid on
    box on
    title('ControllerRegionsCpLambda')
   
   
    figure, hold on,grid on, box on
    CP2=Parametric.Data.CP;
    CP2(CP2<0)=NaN;
    surf(RCPlp.Data.lambda,RCPlp.Data.PITCH,CP2,'FaceAlpha',0.7,'EdgeColor',[0. 0. 0.],'EdgeAlpha',0.5)
    plot3(lambda_opt,pitch_opt,CPmax,'ko','MarkerSize',6,'MarkerFaceColor','k');
    plot3(lambda(i3),Pitch(i3),CP(i3),'ko','MarkerSize',6,'MarkerFaceColor','k');
    plot3(lambda,Pitch,CP+0.001,'k','LineWidth',2);
    xlabel('\lambda [.]')
    ylabel('Pitch [deg]')
    zlabel('C_P [.]')
    axis ij
    view(3)
    dispatchFigs(1)
end



% %%
%
% I = find(beta > 0);
% dpdbeta_FrozenWake=dpdbeta_FrozenWake;
% % temp = polyfit(beta(I(1)-1:end)*180/pi,dpdbeta_FrozenWake(I(1)-1:end),1);
% temp = polyfit(beta(I(1):end)*180/pi,dpdbeta_FrozenWake(I(1):end),1)
% KK = temp(2)/temp(1)
%
% dPdpitch_0 = 0*temp(1)+temp(2) ;%  ([kW/deg])
% dPdpitch_10 = 10*temp(1)+temp(2);
% dPdpitch_20 = 20*temp(1)+temp(2);
% dQdpitch_0 =  dPdpitch_0/par.Omega;  %  ([kW/deg/rad/s])
%
% dQdpitch_0=dQdpitch_0*1000*180/(pi);
%
%
%
%
% %% Inertia
% Rotor.I=355455.2;
% I_r=Rotor.I;
% I_shaft=110529;
% %I_g=45791;
% n_g2=1; % CHECK
% I_t=(I_r + I_shaft );
%
% I_t=I_t;
%
% R=par.sfull(end);
% %% PARAM AERO
% Cp_opt=0.49;
% lambda_opt=8;
% eta_track=0.9;
%
% %% Region 2 control of Omega
% omega_Omega_g=0.1*2*pi; % Response frequency [Hz]  0.1->0.2
% zeta_Omega_g=0.65; % Damping [.]
%
% KPg=2*par.eta*zeta_Omega_g*omega_Omega_g* I_t;
% KIg=par.eta * omega_Omega_g^2 *(1+zeta_Omega_g^2 ) * I_t ;
%
% %% Region 3
% omega_Omega=0.1*2*pi; % Response frequency [Hz]
% zeta_Omega=0.65; % Damping [.]
%
%
% KP=(2*par.eta*zeta_Omega*omega_Omega* I_t -1/par.eta* (-par.Pr/par.Omega^2)  ) / -dQdpitch_0    ;
% %KP_constanttorque=(2*par.eta*zeta_Omega*omega_Omega*(I_r+n_g2*I_g)  ) / -dQdpitch_0
% KI=(omega_Omega^2 *(1+zeta_Omega^2 ) * I_t ) / -dQdpitch_0  ;
%
%
%
%
% %%
% % Cp tracking K factor
% K=eta_track*0.5*par.rho*(pi*R^2)*Cp_opt *R^3/lambda_opt^3  /1000; %[kNm/(rad/s)^2]
%
% disp(sprintf('Cp tracking K: %.2f   [kNm/(rad/s)^2]',K))
% disp(sprintf('Gain scheduling KK: %.2f',KK))
%
% disp(sprintf( [
% '\\begin{itemize}\n'  ...
% '    \\item $k_{P_g} = %.3f$\n' ...
% '    \\item $k_{I_g} = %.3f$\n' ...
% '\\end{itemize}'],KPg,KIg))
%
% disp(sprintf( [
% '\\begin{itemize}\n'  ...
% '    \\item $\\omega_{\\Omega} = %.3f$ [rad/s]\n' ...
% '    \\item $\\zeta_{\\Omega} = %.3f$  [.]\n' ...
% '    \\item $\\parcial{Q}{\\theta} = %.3f$ (obtained from \\autoref{fig:dpdpitch})\n' ...
% '\\end{itemize}'],omega_Omega,zeta_Omega,dQdpitch_0))
%
% disp(sprintf( [
% '\\begin{itemize}\n'  ...
% '    \\item $k_{P} = %.3f$\n' ...
% '    \\item $k_{I} = %.3f$\n' ...
% '\\end{itemize}'],KP,KI))
%
%
%
% disp(sprintf( [
% '            general constant  %f ; KP\n'  ...
% '            general constant  %f ; KI\n' ...
% '            general constant 0.00 ; Channel 12:\n'  ...
% '            general constant  %f ; KPg\n' ...
% '            general constant  %f ; KIg\n' ...
%             ],KP,KI,KPg,KIg))
%
%
%
%
%
%
%
%
%
%
%
%         InitClear
% KK=3.16;
% theta=0:30;
% G=1./(1+theta/KK);
% figure(1)
% plot(theta,G,'k','LineWidth',2)
% text('Interpreter','latex','String','$$G=\frac{1}{1+\frac{\theta_{m}}{K_K}}$$', 'Position',[15 0.8], 'FontSize',14)
% text('Interpreter','latex','String',sprintf('$$K_K\\ =\\ %.1f\\ [^{\\circ}]$$',KK), 'Position',[15 0.65], 'FontSize',14)
% xlabel('Measured pitch angle \theta_m [deg]')
% ylabel('Gain Scheduling G [.]')
% %grid on
% box on
% title('GainScheduling')
% xlim([-1 30])
%
% %%
% omrel=10;
% theta=0:30;
% omratio=0:0.1:2;
% Gnl=1+(omratio-1).^2/(omrel-1)^2;
% figure(2)
% plot(omratio,Gnl,'k','LineWidth',2)
% text('Interpreter','latex','String','$$G_{nl}=1+\frac{(\Omega-\Omega_r)^2}{\Omega_r^2 (\omega_{rel}-1)^2}$$', 'Position',[1 1.07], 'FontSize',14)
% text('Interpreter','latex','String',sprintf('$$\\omega_{rel}\\ =\\ %d\\ [.]$$',omrel), 'Position',[1 1.04], 'FontSize',14)
% xlabel('Rotor speed ratio \Omega/\Omega_r [.]')
% ylabel('Gain Scheduling - Non linear gain G_{nl} [.]')
% %grid on
% box on
% title('GainSchedulingNL')
% xlim([0 2])
% ylim([0.9 1.1])
