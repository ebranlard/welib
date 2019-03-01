
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LiftDataFun functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ClCdCm = LiftDataFun(alpha,par,Data);
% finde lift data
for i=1:length(par.ProfileSet(:,1))
    if par.ProfileSet(i,1) == 1
        n1 = 1+sum(Data.ThicknessVec(1:par.ProfileSet(i,2)-1,2));
        n2 = 1+sum(Data.ThicknessVec(1:par.ProfileSet(i,2)-1,2))+Data.ThicknessVec(par.ProfileSet(i,2),2)-1;
        ClCdCm(i,:) = interp1(Data.ProfileData(n1:n2,1),Data.ProfileData(n1:n2,2:4),alpha(i));
    else
        n1 = 1+sum(Data.ThicknessVec(1:par.ProfileSet(i,2)-1,2));
        n2 = 1+sum(Data.ThicknessVec(1:par.ProfileSet(i,2)-1,2))+Data.ThicknessVec(par.ProfileSet(i,2),2)-1;
        temp(1,1:3) = interp1(Data.ProfileData(n1:n2,1),Data.ProfileData(n1:n2,2:4),alpha(i));
        n1 = 1+sum(Data.ThicknessVec(1:par.ProfileSet(i,3)-1,2));
        n2 = 1+sum(Data.ThicknessVec(1:par.ProfileSet(i,3)-1,2))+Data.ThicknessVec(par.ProfileSet(i,3),2)-1;
        temp(2,1:3) = interp1(Data.ProfileData(n1:n2,1),Data.ProfileData(n1:n2,2:4),alpha(i));
        ClCdCm(i,:) = interp1(Data.ThicknessVec(par.ProfileSet(i,2:3),1),temp(1:2,1:3),par.thickness(i));
    end
end