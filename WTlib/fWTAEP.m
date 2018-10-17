function AEP=fWTAEP() 

hpy=24*365.25;% hours per year
hoursRepartition=wblpdf(Vspeed,6,2)*hpy;
AEP=sum(hoursRepartition.*Pv)/1000/1000; % in MW.h/year



