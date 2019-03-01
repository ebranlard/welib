%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MainFun functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [par res Cp P]=fPower(user,par,Data)

Ngrid = 30;
dbeta = 0.01; % used for computing slope

%% stations without hub length
% find thickness and profile data
if par.r > 0
    par.sfull = linspace(Data.AeData(1,1),Data.AeData(end,1),Ngrid+1)
    par.s = par.sfull(1:end-1)';
else % skip center point if no hub element
    par.sfull = linspace(Data.AeData(1,1),Data.AeData(end,1),Ngrid+2);
    par.s = par.sfull(2:end-1)';
end

for i=1:Ngrid
    par.c(i,1) = interp1(Data.AeData(:,1),Data.AeData(:,2),par.s(i));
    par.thickness(i,1) = interp1(Data.AeData(:,1),Data.AeData(:,3),par.s(i));
    par.theta(i,1) = interp1(Data.PitchAxis(:,4),Data.PitchAxis(:,5),par.s(i))*pi/180;
   profilesup=find(Data.ThicknessVec(:,1) >= par.thickness(i));
   profileinf=find(Data.ThicknessVec(:,1) <=par.thickness(i));
   profileinf=max(profileinf);
   if isempty(profilesup)
%        disp(' oué oué')
       profilesup=profileinf;
   end
   profilesup=profilesup(1);
   par.ProfileSet(i,:) = [(profileinf~=profilesup)+1,profileinf,profilesup];
end

   profilesup=find(Data.ThicknessVec(:,1) >= par.thickness(i));
   profileinf=find(Data.ThicknessVec(:,1) <=par.thickness(i));
   profileinf=max(profileinf);
   if isempty(profilesup)
%        disp(' oué oué')
       profilesup=profileinf;
   end
   profilesup=profilesup(1);
   par.ProfileSet(i,:) = [(profileinf~=profilesup)+1,profileinf,profilesup];

%% stations now have hub length
par.s = par.s+par.r;
par.sfull = par.sfull+par.r;
par.sigma = Data.Nb*par.c./(2*pi*par.s);


% compute power and pitch for each wind speed
beta_temp = 0;
P=zeros(1,length(par.Uvec));
beta=zeros(1,length(par.Uvec));
res.an=zeros(Ngrid,length(par.Uvec));
res.at=zeros(Ngrid,length(par.Uvec));

for i=1:length(par.Uvec)
    par.an = zeros(Ngrid,1);
    par.at = par.an;
    par.Uinf = par.Uvec(i);
    [P(i),par] = AeroFun(beta_temp,par,Data);
    if P(i) > par.Pref
        beta_temp = fzero(@AeroFun,beta_temp+0.1,[],par,Data,1);
        [P(i),par] = AeroFun(beta_temp,par,Data);
    end
    res.an(:,i) = par.an;
    res.at(:,i) = par.at;
    beta(i) = beta_temp;
    disp(['Uinf = ',num2str(par.Uvec(i)),'  P = ',num2str(P(i)),'  pitch = ',num2str(beta(i)*180/pi)])
end
res.beta=beta;
% compute dP/dpitch, frozen wake
par.Induction = 0;
for i=1:length(par.Uvec)
    par.Uinf = par.Uvec(i);
    par.an = res.an(:,i);
    par.at = res.at(:,i);
end

Cp = P./(0.5*par.rho*par.Uvec.^3*pi*par.sfull(end)^2);

end
