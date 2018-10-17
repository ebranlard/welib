function ClCdCm = fAeroCoeff(alpha,Profiles,ProfileSet,rel_thickness,Re,bReInterp,bThicknessInterp,bRough)
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
    iRe1inf=whichvalue(vRe{i1}(:),Re(i));
    iRe2inf=whichvalue(vRe{i2}(:),Re(i));
    if(bReInterp)
        %First Profile
        if((length(vRe{i1})>=iRe1inf+1) && vRe{i1}(iRe1inf)<Re(i))
            temp(1,1:3) = interp1(Data{i1}{iRe1inf}(:,1),Data{i1}{iRe1inf}(:,2:4),alpha(i),'linear');    %Reynolds inf
            temp(2,1:3) = interp1(Data{i1}{iRe1inf+1}(:,1),Data{i1}{iRe1inf+1}(:,2:4),alpha(i),'linear');%Reynolds sup
            % Interpolating for Reynolds
            temp(3,1:3)=interp1(vRe{i1}(iRe1inf:iRe1inf+1),temp(1:2,1:3),Re(i),'linear' );
        else
            temp(3,1:3) = interp1(Data{i1}{iRe1inf}(:,1),Data{i1}{iRe1inf}(:,2:4),alpha(i),'linear');    %Reynolds inf
        end
        if ProfileSet(1,i) == 1  
            ClCdCm(i,:) = temp(3,1:3);
        elseif ~bThicknessInterp
            ClCdCm(i,:) = temp(3,1:3);
        else
            %Second Profile
            if((length(vRe{i2})>=iRe2inf+1) && vRe{i2}(iRe2inf)<Re(i))
                temp(1,1:3) = interp1(Data{i2}{iRe2inf}(:,1),Data{i2}{iRe2inf}(:,2:4),alpha(i),'linear');    %Reynolds inf
                temp(2,1:3) = interp1(Data{i2}{iRe2inf+1}(:,1),Data{i2}{iRe2inf+1}(:,2:4),alpha(i),'linear');%Reynolds sup
                % Interpolating for Reynolds
                temp(2,1:3)=interp1(vRe{i2}(iRe2inf:iRe2inf+1),temp(1:2,1:3),Re(i) ,'linear');
            else
                temp(2,1:3) = interp1(Data{i2}{iRe2inf}(:,1),Data{i2}{iRe2inf}(:,2:4),alpha(i),'linear');    %Reynolds inf
            end      
            % Interpolating for thickness
            ClCdCm(i,:) = interp1(Profiles.thickness_rel(i1:i2),temp(3:-1:2,1:3),rel_thickness(i),'linear');
                if(sum(isnan(Profiles.thickness_rel(i1:i2)))+sum(isnan(temp(3:-1:2,1:3)))+sum(isnan(rel_thickness(i)))>0)
                    warning('here');
                    keyboard
                end
        end
    else
        iRe=whichvalue(vRe{i1}(:),Re(i));
        if ProfileSet(1,i) == 1
            ClCdCm(i,:) = interp1(Data{i1}{iRe}(:,1),Data{i1}{iRe}(:,2:4),alpha(i),'linear');
        else
            %temp(1,1:3) = interp1(Data(n1:n2,1),Data(n1:n2,2:4),alpha(i));
            %first profile
            temp(3,1:3) = interp1(Data{i1}{iRe}(:,1),Data{i1}{iRe}(:,2:4),alpha(i),'linear');
            
            if ~bThicknessInterp
                 ClCdCm(i,:)=temp(3,1:3);
            else
                %second profile
                if(length(Data{i2})<iRe)
                    iRe=length(Data{i2});
                end
                temp(2,1:3) = interp1(Data{i2}{iRe}(:,1),Data{i2}{iRe}(:,2:4),alpha(i),'linear');
                % Interpolating for thickness
                ClCdCm(i,:) = interp1(Profiles.thickness_rel(i1:i2),temp(3:-1:2,1:3),rel_thickness(i),'linear');
                if(sum(isnan(Profiles.thickness_rel(i1:i2)))+sum(isnan(temp(3:-1:2,1:3)))+sum(isnan(rel_thickness(i)))>0)
                    warning('here');
                    keyboard
                end
            end
        end
    end
end
% ClCdCm=zeros(length(ProfileSet(1,:)),3);
% for i=1:length(ProfileSet(1,:))
%     if ProfileSet(1,i) == 1
%         n1 = 1+sum(Profiles.ndata(1:ProfileSet(2,i)-1));
%         n2 = 1+sum(Profiles.ndata(1:ProfileSet(2,i)-1))+Profiles.ndata(ProfileSet(2,i))-1;
%         ClCdCm(i,:) = interp1(Data(n1:n2,1),Data(n1:n2,2:4),alpha(i));
%     else
%         n1 = 1+sum(Profiles.ndata(1:ProfileSet(2,i)-1));
%         n2 = 1+sum(Profiles.ndata(1:ProfileSet(2,i)-1))+Profiles.ndata(ProfileSet(2,i))-1;
%         temp(1,1:3) = interp1(Data(n1:n2,1),Data(n1:n2,2:4),alpha(i));
%         n1 = 1+sum(Profiles.ndata(1:ProfileSet(3,i)-1));
%         n2 = 1+sum(Profiles.ndata(1:ProfileSet(3,i)-1))+Profiles.ndata(ProfileSet(3,i))-1;
%         temp(2,1:3) = interp1(Data(n1:n2,1),Data(n1:n2,2:4),alpha(i));
%         ClCdCm(i,:) = interp1(Profiles.thickness_rel(ProfileSet(2:3,i)),temp(1:2,1:3),rel_thickness(i));
%     end
% end
