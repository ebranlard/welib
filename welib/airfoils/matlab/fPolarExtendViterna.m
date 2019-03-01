function [newpolar]=fAirfoilExtendPolarViterna(orgpolar,uniformPointsFlag)
%If uniformPointsFlag is one, then it outputs in a uniform point
%distribution. If 0, then it matching with sent data and maintains original
%points




%%
AspectRatio=25;
maxcd = 1.11+0.018*AspectRatio;
%MaxCD = 1.56;
MaxCD=maxcd;

FlatPlateCD='no';
%%
%find points for matching.
% indCLmax=find(orgpolar(:,2)==max(orgpolar(:,2)));
% indCLmin=find(orgpolar(:,2)==min(orgpolar(:,2)));

[CLmax,indCLmax]=max(orgpolar(:,2));
[CLmin,indCLmin]=min(orgpolar(:,2));
pastPosStall=4;
pastNegStall=2;
%the following extends index for hi and low by pastPosStall and pastNegStall 
if (indCLmax+pastPosStall)<=length(orgpolar(:,1))
    indmatchhi=indCLmax+pastPosStall; 
else
    indmatchhi=length(orgpolar(:,1));
end 
if (indCLmin-pastNegStall)>=1
    indmatchlo=indCLmin-pastNegStall; 
else
    indmatchlo=1;
end 

%the following extends index for hi and low to include angles past neg and
%pos stall by pastPosStall and pastNegStall
if orgpolar(end,1)<orgpolar(indCLmax,1)+pastPosStall
    indmatchhi=length(orgpolar(:,1));
else
    [temp,indmatchhi]=max(orgpolar(orgpolar(:,1)<orgpolar(indCLmax,1)+pastPosStall));
end
if orgpolar(1,1)>orgpolar(indCLmin,1)-pastNegStall
    indmatchlo=1;
else
    [temp,indmatchlo]=min(orgpolar(orgpolar(:,1)>orgpolar(indCLmin,1)-pastNegStall));
end


AlfaHi=orgpolar(indmatchhi,1); 
CLHi=orgpolar(indmatchhi,2); 
CDHi=orgpolar(indmatchhi,3); 
CMHi=orgpolar(indmatchhi,4);

AlfaLo=orgpolar(indmatchlo,1); 
CLLo=orgpolar(indmatchlo,2); 
CDLo=orgpolar(indmatchlo,3); 
CMLo=orgpolar(indmatchlo,4);

%fprintf(1,'\n%s\n',['lower alfa match = ' num2str(AlfaLo) ', upper alfa match = ' num2str(AlfaHi)]);
%% Assemble Vitnera coeffs
SinAlfaHi=sin(AlfaHi*pi/180);  CosAlfaHi=cos(AlfaHi*pi/180);

A2 = (CLHi-MaxCD/2*sin(2*AlfaHi*pi/180))*SinAlfaHi/(CosAlfaHi^2);
B2 = (CDHi-MaxCD*SinAlfaHi^2)/CosAlfaHi; %disp(['B2 = ' num2str(B2)]);

%% CMCoeff
ind=find(orgpolar(:,2)<=0); 
if isempty(ind)
    ind=1;
end
ind=ind(end);
if ind==length(orgpolar(:,2))
    newpolar=orgpolar;
    return;
end
CM0=interp1([orgpolar(ind,2) orgpolar(ind+1,2)],[orgpolar(ind,4) orgpolar(ind+1,4)],[0],'linear','extrap');%CM at zero lift
XM=(-CMHi+CM0)/(CLHi*CosAlfaHi+CDHi*SinAlfaHi);
CMCoeff= (XM-0.25)/(tan((AlfaHi-90)*pi/180));
%% NewTable
if uniformPointsFlag
    newAoAarray=[-180:5:-30 -29:1:29 30:5:180]';
    
else
    newAoAarray=[-180:5:-30 -29:1:29 30:5:180];
    newAoAarray=[newAoAarray(newAoAarray<AlfaLo),orgpolar(orgpolar(:,1)>=AlfaLo&orgpolar(:,1)<=AlfaHi,1)',newAoAarray(newAoAarray>AlfaHi)]';
end


newpolar=zeros(length(newAoAarray),4); 
newpolar(:,1)=newAoAarray;

indstart=find(newpolar(:,1)>=orgpolar(indmatchlo,1)); 
indstart=indstart(1);

indend=find(newpolar(:,1)<=orgpolar(indmatchhi,1));
indend=indend(end);

newpolar(indstart:indend,2)=interp1(orgpolar(:,1),orgpolar(:,2),newpolar(indstart:indend,1),'cubic');
newpolar(indstart:indend,3)=interp1(orgpolar(:,1),orgpolar(:,3),newpolar(indstart:indend,1),'cubic');
newpolar(indstart:indend,4)=interp1(orgpolar(:,1),orgpolar(:,4),newpolar(indstart:indend,1),'cubic');
%% ViternaFill
CLAdj=0.7;%ratio on positive polar
FPAlfa=max([abs(orgpolar(indmatchhi,1)) abs(orgpolar(indmatchlo,1))]);
for i=1:length(newpolar(:,1))
    alfa=newpolar(i,1);
    SinAlfa=sin(alfa*pi/180); CosAlfa=cos(alfa*pi/180);
    if alfa<-180+AlfaHi
        Ang=180+alfa;
        SinAng=sin(Ang*pi/180); CosAng=cos(Ang*pi/180);
        newpolar(i,2)=Ang/AlfaHi*CLHi*CLAdj;
        newpolar(i,3)=MaxCD*SinAng^2+B2*CosAng; %disp(['Alfa = ' num2str(alfa) ', cl = ' num2str(newpolar(i,2)) ', cd = ' num2str(newpolar(i,3))])
    elseif alfa<-90 && alfa>=-180+AlfaHi
        Ang=180+alfa;
        SinAng=sin(Ang*pi/180); CosAng=cos(Ang*pi/180);
        newpolar(i,2)=CLAdj*(MaxCD/2*sin(2*Ang*pi/180)+A2*(CosAng^2)/SinAng);
        newpolar(i,3)=MaxCD*SinAng^2+B2*CosAng;
    elseif alfa<=-max([abs(AlfaHi) abs(AlfaLo)]) && alfa>=-90
        Ang=-alfa;
        SinAng=sin(Ang*pi/180); CosAng=cos(Ang*pi/180);
        newpolar(i,2)=CLAdj*(-MaxCD/2*sin(2*Ang*pi/180)-A2*(CosAng^2)/SinAng);
        newpolar(i,3)=MaxCD*SinAng^2+B2*CosAng;
    elseif alfa>-AlfaHi && alfa<AlfaLo
        newpolar(i,2)=CLAdj*(-CLHi+(alfa+AlfaHi)/(AlfaHi+AlfaLo)*(CLHi+CLLo)); %disp(alfa);
        newpolar(i,3)=CDLo+(-alfa+AlfaLo)*(CDHi-CDLo)/(AlfaHi+AlfaLo);
        %Ang=-alfa;
        %SinAng=sin(Ang*pi/180); CosAng=cos(Ang*pi/180);
        %newpolar(i,2)=CLAdj*(-MaxCD/2*sin(2*Ang*pi/180)-A2*(CosAng^2)/SinAng);
        %newpolar(i,3)=MaxCD*SinAng^2+B2*CosAng;
    elseif alfa<=90 && alfa>AlfaHi
        newpolar(i,2)=MaxCD/2*sin(2*alfa*pi/180)+A2*(CosAlfa^2)/SinAlfa;
        newpolar(i,3)=MaxCD*SinAlfa^2+B2*CosAlfa;
    elseif alfa>90 && alfa<=180-AlfaHi
        Ang=180-alfa;
        SinAng=sin(Ang*pi/180); CosAng=cos(Ang*pi/180);
        newpolar(i,2)=CLAdj*(-MaxCD/2*sin(2*Ang*pi/180)-(A2*CosAng^2)/SinAng);
        newpolar(i,3)=MaxCD*SinAng^2+B2*CosAng;
    elseif alfa>180-AlfaHi
        Ang=alfa-180;
        SinAng=sin(Ang*pi/180); CosAng=cos(Ang*pi/180);
        newpolar(i,2)=Ang/AlfaHi*CLHi*CLAdj;
        newpolar(i,3)=MaxCD*SinAng^2+B2*CosAng;
    end
    %remove any negative CD values
    if newpolar(i,3)<0.03 && abs(alfa)>20, newpolar(i,3)=0.03; end;
end
   
switch lower(FlatPlateCD)
    case 'yes'
        for i=1:length(newpolar(:,1))
            if abs(newpolar(i,1))>FPAlfa, newpolar(i,3)=newpolar(i,2)*tan(newpolar(i,1)*pi/180); end;
        end
end
%% GetCM
for i=1:length(newpolar(:,1))
    alfa=newpolar(i,1);
    if alfa>=AlfaLo && alfa<=AlfaHi, continue; end
    if abs(alfa)<165
        x=CMCoeff*tan((alfa-90)*pi/180)+0.25;
        newpolar(i,4)=CM0-x*(newpolar(i,2)*cos(alfa*pi/180)+newpolar(i,3)*sin(alfa*pi/180));
    end
end
% newpolar((find(newpolar(:,1)==165)),4)=-0.4; newpolar((find(newpolar(:,1)==170)),4)=-0.5;
% newpolar((find(newpolar(:,1)==175)),4)=-0.25; newpolar((find(newpolar(:,1)==180)),4)=0.0;
% newpolar((find(newpolar(:,1)==-165)),4)=0.35; newpolar((find(newpolar(:,1)==-170)),4)=0.4;
% newpolar((find(newpolar(:,1)==-175)),4)=0.2; newpolar((find(newpolar(:,1)==-180)),4)=0.0;

% %% Clean up inside range (redundant after "fill in AoAs" section @ top)
% newpolar(indstart+1:indend-1,2)=interp1(polar(:,1),polar(:,2),newpolar(indstart+1:indend-1,1),'cubic');
% newpolar(indstart+1:indend-1,3)=interp1(polar(:,1),polar(:,3),newpolar(indstart+1:indend-1,1),'cubic');
% newpolar(indstart+1:indend-1,4)=interp1(polar(:,1),polar(:,4),newpolar(indstart+1:indend-1,1),'cubic');        
        
%% Attempt to clean up CM curve
indCMstart=find(newpolar(:,1)<=-50); indCMstart=indCMstart(end);
tmparray1(:,1)=[newpolar(indCMstart-5:indCMstart,1)]; tmparray1(:,2)=[newpolar(indCMstart-5:indCMstart,4)];
tmparray2(:,1)=[newpolar(indstart:indend,1)]; tmparray2(:,2)=[newpolar(indstart:indend,4)];
tmparray=[tmparray1;tmparray2];
for i=1:length(newpolar(:,1))
    alfa=newpolar(i,1);
    if alfa>-50 && alfa < AlfaLo
        newpolar(i,4)=interp1(tmparray(:,1),tmparray(:,2),newpolar(i,1),'cubic');
    end
end
clear tmparray1 tmparray2 tmparray 
indCMstart=find(newpolar(:,1)<=-160); indCMstart=indCMstart(end);
tmparray1(:,1)=-180; tmparray1(:,2)=0;
tmparray2(:,1)=[newpolar(indCMstart:indCMstart+10,1)]; tmparray2(:,2)=[newpolar(indCMstart:indCMstart+10,4)];
tmparray=[tmparray1;tmparray2];
for i=1:length(newpolar(:,1))
    alfa=newpolar(i,1);
    if alfa<-160
        newpolar(i,4)=interp1(tmparray(:,1),tmparray(:,2),newpolar(i,1),'cubic');
    end
end
clear tmparray1 tmparray2 tmparray 
indCMstart=find(newpolar(:,1)<=160); indCMstart=indCMstart(end);
tmparray1(:,1)=180; tmparray1(:,2)=0;
tmparray2(:,1)=[newpolar(indCMstart-10:indCMstart,1)]; tmparray2(:,2)=[newpolar(indCMstart-10:indCMstart,4)];
tmparray=[tmparray2;tmparray1];
for i=1:length(newpolar(:,1))
    alfa=newpolar(i,1);
    if alfa>160
        newpolar(i,4)=interp1(tmparray(:,1),tmparray(:,2),newpolar(i,1),'cubic');
    end
end 
%%
% scrsz = get(0,'ScreenSize'); 
% hh=figure;
% set(hh,'Position',[10 10 scrsz(3)*.85 scrsz(4)*.85])
% subplot('Position',[.06 .53 .425 .40]);
% plot(newpolar(:,1),newpolar(:,2));
% subplot('Position',[.54 .53 .425 .40]);
% plot(newpolar(:,1),newpolar(:,3)); 
% subplot('Position',[.06 .06 .425 .40]); 
% plot(newpolar(:,1),newpolar(:,4));
% subplot('Position',[.54 .06 .425 .40]); 
% plot(newpolar(:,1),newpolar(:,2)./newpolar(:,3));
% subplot('Position',[.06 .53 .425 .40]);
% hold on; plot(orgpolar(:,1),orgpolar(:,2),'r+'); hold off; 
% subplot('Position',[.54 .53 .425 .40]); 
% hold on; plot(orgpolar(:,1),orgpolar(:,3),'r+'); hold off; 
% subplot('Position',[.06 .06 .425 .40]); 
% hold on; plot(orgpolar(:,1),orgpolar(:,4),'r+'); hold off; 
% subplot('Position',[.54 .06 .425 .40]); 
% hold on; plot(orgpolar(:,1),orgpolar(:,2)./orgpolar(:,3),'r+'); hold off;    
end