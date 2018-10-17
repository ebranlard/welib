function ClCdCm = fAeroCoeff(alpha,Profiles,ProfileSet,rel_thickness)
% alpha in deg
% find lift data
ClCdCm=zeros(length(ProfileSet(:,1)),3);
for i=1:length(ProfileSet(:,1))
    if ProfileSet(i,1) == 1
        n1 = 1+sum(Profiles.ThicknessVec(1:ProfileSet(i,2)-1,2));
        n2 = 1+sum(Profiles.ThicknessVec(1:ProfileSet(i,2)-1,2))+Profiles.ThicknessVec(ProfileSet(i,2),2)-1;
        ClCdCm(i,:) = interp1(Profiles.Data(n1:n2,1),Profiles.Data(n1:n2,2:4),alpha(i));
    else
        n1 = 1+sum(Profiles.ThicknessVec(1:ProfileSet(i,2)-1,2));
        n2 = 1+sum(Profiles.ThicknessVec(1:ProfileSet(i,2)-1,2))+Profiles.ThicknessVec(ProfileSet(i,2),2)-1;
        temp(1,1:3) = interp1(Profiles.Data(n1:n2,1),Profiles.Data(n1:n2,2:4),alpha(i));
        n1 = 1+sum(Profiles.ThicknessVec(1:ProfileSet(i,3)-1,2));
        n2 = 1+sum(Profiles.ThicknessVec(1:ProfileSet(i,3)-1,2))+Profiles.ThicknessVec(ProfileSet(i,3),2)-1;
        temp(2,1:3) = interp1(Profiles.Data(n1:n2,1),Profiles.Data(n1:n2,2:4),alpha(i));
        ClCdCm(i,:) = interp1(Profiles.ThicknessVec(ProfileSet(i,2:3),1),temp(1:2,1:3),rel_thickness(i));
    end
end