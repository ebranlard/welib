function [vRPM,Regions]= fWTGetRefRPM(WT,vWS,WS_Startup,lambda_opt,TipSpeedMax,OmegaMin,bPlot)
% Computes the reference RPM curve at vWS 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Determination of omega in Region 1 and above
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
OmegaMax = TipSpeedMax/WT.Rotor.R          ;

%% Some Wind Speed Regions
% WS below which we dont operate (and above which the RPM may be Omega Min)
Regions.WS0 = WS_Startup;
% WS above which region 1 start (respecting OmegaMin)
Regions.WS1 = max(WS_Startup, OmegaMin*WT.Rotor.R/lambda_opt);
% WS above which region 2 starts (above which the rotational speed is constant)
Regions.WS2 =OmegaMax*WT.Rotor.R/lambda_opt;
% Max Wind speed
Regions.WS4 =max(vWS);


if(Regions.WS2<WS_Startup)
    fprintf('!!! Warning: Optimal Cp tracking impossible lambda_opt too high\n');
end
if(Regions.WS0~=Regions.WS1)
    fprintf('!!! Warning: Optimal Cp tracking not done bewteen %.2f and %.2f\n',Regions.WS0,Regions.WS1);
end




%% Loop on wind speed to find appropriate RPM
vRPM    = zeros(1,length(vWS));
for iWS = 1:length(vWS)
    WS=vWS(iWS);
    if WS<WS_Startup 
        % Not enough wind
        vRPM(iWS)    = 0; 
    elseif WS >= Regions.WS2
        % Constant RPM (Max RPM)
        vRPM(iWS)    = OmegaMax*60/(2*pi); % 
    else
        % Region of varying RPM
        Omega1  = lambda_opt*WS/WT.Rotor.R;
        if Omega1 < OmegaMin
            % Optimal rotational speed too low, we need to restrict it
            Omega1  = OmegaMin;
        else
            % Using Optimal CP tracking
        end
        vRPM(iWS)    = Omega1*60/(2*pi);
    end
end
%%

%%
if bPlot
    figure, hold all
    plot(vWS,vRPM)
    xlabel('WS')
    ylabel('RPM')

    % 
    PLOT_RANGE=[floor(min(vRPM)-1), ceil(max(vRPM)+1)];
    plot([Regions.WS0 Regions.WS0],PLOT_RANGE,'k','LineWidth',2 )
    plot([Regions.WS1 Regions.WS1],PLOT_RANGE,'k','LineWidth',2 )
    plot([Regions.WS2 Regions.WS2],PLOT_RANGE,'k','LineWidth',2 )

    tpos=mean(PLOT_RANGE);
    if Regions.WS0~=Regions.WS1 
        text((Regions.WS0+Regions.WS1)/2-1,tpos,'R0')
    end
    text((Regions.WS1+Regions.WS2)/2-1,tpos,'R1')
    text((Regions.WS2+Regions.WS4)/2-1,tpos,'R2-R3')
    ylim(PLOT_RANGE)
end
