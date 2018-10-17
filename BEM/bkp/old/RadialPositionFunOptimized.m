function [Rotor]=RadialPositionFunOptimized(Data,Rotor);
%% stations without hub length
% find thickness and profile data
if Rotor.rhub > 0
    Rotor.sfull = linspace(Rotor.init.s(1),Rotor.init.s(end),Rotor.Ngrid+1);
    Rotor.s = Rotor.sfull(1:end-1)';
else % skip center point if no hub element
    Rotor.sfull = linspace(Rotor.init.s(1),Rotor.init.s(end),Rotor.Ngrid+2);
    Rotor.s = Rotor.sfull(2:end-1)';
end

for i=1:Rotor.Ngrid
    %rotor.c(i,1) = interp1(Data.AeData(:,1),Data.AeData(:,2),rotor.s(i));
    %rotor.thickness(i,1) = interp1(Data.AeData(:,1),Data.AeData(:,3),rotor.s(i));
    %rotor.theta(i,1) = interp1(Data.PitchAxis(:,4),Data.PitchAxis(:,5),rotor.s(i))*pi/180;
    
    Rotor.opt_rel_thickness(i,1) = interp1(Rotor.s,Rotor.opt_rel_thickness,Rotor.s(i));
    Rotor.opt_abs_thickness(i,1) = interp1(Rotor.init.s,Rotor.init.abs_thickness,Rotor.s(i));
%     temp = find(Data.ThicknessVec(:,1) >= Rotor.abs_thickness(i));
%     if length(temp) == 1
%         Rotor.ProfileSet(i,:) = [1,temp,temp];
%     elseif length(temp) == length(Data.ThicknessVec)
%         Rotor.ProfileSet(i,:) = [1,1,1];
%     else
%         Rotor.ProfileSet(i,:) = [2,temp(1)-1,temp(1)];
%     end
   profilesup=find(Data.ThicknessVec(:,1) >= Rotor.opt_rel_thickness(i));
   profileinf=find(Data.ThicknessVec(:,1) <= Rotor.opt_rel_thickness(i));
   profileinf=max(profileinf);
   if isempty(profilesup)
%        disp(' oué oué')
       profilesup=profileinf;
   end
   profilesup=profilesup(1);
   Rotor.ProfileSet(i,:) = [(profileinf~=profilesup)+1,profileinf,profilesup];
end
%% stations now have hub length
Rotor.s = Rotor.s+Rotor.rhub;
Rotor.sfull = Rotor.sfull+Rotor.rhub;

end

