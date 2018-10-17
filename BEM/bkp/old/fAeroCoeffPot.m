function ClCdCm = fAeroCoeffPot(alpha,Profiles,ProfileSet,rel_thickness,Re,bReInterp,bThicknessInterp,bRough)
% alpha in deg
% find lift dataInit
if(bRough)
    Data=Profiles.DataRough;
    vRe=Profiles.ReRough;
else
    Data=Profiles.Data;
    vRe=Profiles.Re;
end


ClCdCm=zeros(length(ProfileSet(1,:)),3);
temp=zeros(3,3);
for i=1:length(ProfileSet(1,:))
    i1=ProfileSet(2,i);
    i2=ProfileSet(3,i);
    iRe=whichvalue(vRe{i1}(:),Re(i));
    Ialpha=find(Data{i1}{iRe}(:,1)<10 & Data{i1}{iRe}(:,1)>-10); %alpha0 in [-5 5]
    Ialpha0Neg=find(Data{i1}{iRe}(Ialpha,2)<0,1,'last');   
    Ialpha0Pos=Ialpha0Neg+1;
    alpha01=interp1([Data{i1}{iRe}(Ialpha(Ialpha0Neg),2) Data{i1}{iRe}(Ialpha(Ialpha0Pos),2)],[Data{i1}{iRe}(Ialpha(Ialpha0Neg),1) Data{i1}{iRe}(Ialpha(Ialpha0Pos),1)],0);
    
    if ProfileSet(1,i) == 1
        ClCdCm(i,:) = [2*pi*sind(alpha(i)-alpha01) 0 0];
%         [Profiles.thickness_rel(i1) alpha01]
    else
        if ~bThicknessInterp
            ClCdCm(i,:) = [2*pi*sind(alpha(i)-alpha01) 0 0];
        else
            Ialpha=find(Data{i2}{iRe}(:,1)<10 & Data{i2}{iRe}(:,1)>-10); %alpha0 in [-10 10]
            Ialpha0Neg=find(Data{i2}{iRe}(Ialpha,2)<0,1,'last');   
            if(isempty(Ialpha0Neg))
                alpha02=0;%cylinder
            else
                Ialpha0Pos=Ialpha0Neg+1;
                alpha02=interp1([Data{i2}{iRe}(Ialpha(Ialpha0Neg),2) Data{i2}{iRe}(Ialpha(Ialpha0Pos),2)],[Data{i2}{iRe}(Ialpha(Ialpha0Neg),1) Data{i2}{iRe}(Ialpha(Ialpha0Pos),1)],0);
            end
%             [Profiles.thickness_rel(i1) alpha01 alpha02 Profiles.thickness_rel(i2)]
            temp(1,1:3)=[2*pi*sind(alpha(i)-alpha01) 0 0];
            temp(2,1:3)=[2*pi*sind(alpha(i)-alpha02) 0 0];
            % Interpolating for thickness
            ClCdCm(i,:) = interp1(Profiles.thickness_rel(i1:i2),temp(1:2,1:3),rel_thickness(i),'linear');
        end
    end
end
