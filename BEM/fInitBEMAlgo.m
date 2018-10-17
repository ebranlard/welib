function BEM =fInitBEMAlgo()

BEM.Kind='';
BEM.CTCorrection='AeroDyn';
BEM.SwirlMethod='AeroDyn';
BEM.NumSect=1;
BEM.bTipLoss=logical(1);
BEM.bHubLoss=logical(0); %!!!
BEM.bTipLossCl=logical(0);
BEM.TipLossMethod='Glauert';
BEM.bYawModel=logical(0);
BEM.bDynaWake=logical(0);  % for unsteady BEM!!!!
BEM.w_guess=-2;
