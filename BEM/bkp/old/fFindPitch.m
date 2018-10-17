%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MainFun functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [par res]=fFindPitch(user,par,Data)

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
        %        disp(' oué oué')  %
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

% compute dP/dpitch, frozen wake
par.Induction = 0;
for i=1:length(par.Uvec)
    par.Uinf = par.Uvec(i);
    par.an = res.an(:,i);
    par.at = res.at(:,i);
    P1 = AeroFun(beta(i)+dbeta,par,Data);
    P2 = AeroFun(beta(i)-dbeta,par,Data);
    dpdbeta_FrozenWake(i) = (P1-P2)/(2*dbeta*180/pi)*1e-3;
end

% compute kk factor based on frozen wake
I = find(beta > 0);
% temp = polyfit(beta(I(1)-1:end)*180/pi,dpdbeta_FrozenWake(I(1)-1:end),1);
temp = polyfit(beta(I(1):end)*180/pi,dpdbeta_FrozenWake(I(1):end),1);
KK = temp(2)/temp(1);

dPdpitch_0 = 0*temp(1)+temp(2);
dPdpitch_10 = 10*temp(1)+temp(2);
dPdpitch_20 = 20*temp(1)+temp(2);

dQdpitch_0 =  dPdpitch_0/par.Omega;
dQdpitch_10 = dPdpitch_10/par.Omega;
dQdpitch_20 = dPdpitch_20/par.Omega;



Cp = P./(0.5*par.rho*par.Uvec.^3*pi*par.sfull(end)^2);



% write results
fid = fopen('output.txt','w');
fprintf(fid,'%s \n\n',['Based on information from HAWC2 inputfile: ',user.MainPath,user.HtcFileName]);
fprintf(fid,'%s \n',['Aerodynamic file: ',user.MainPath,Data.AeFileName]);
temp = [Data.AeData(:,1)+par.r,Data.AeData(:,2:end),interp1(Data.PitchAxis(:,4),Data.PitchAxis(:,5),Data.AeData(:,1))];
fprintf(fid,'%s \n',['radius chord    thickness profile twist']);
fprintf(fid,'%s \n',['  m      m          %              deg  ']);
fprintf(fid,'%6.2f  %6.2f  %6.1f  %5d  %8.2f \n',temp');
fprintf(fid,'\n%s \n\n',['Profile data file: ',Data.PcFileName]);
fprintf(fid,'%s \n','  U        P       pitch  dP/dTheta        Cp');
fprintf(fid,'%s \n','                         frozen wake');
fprintf(fid,'%s \n',' m/s      kW        deg     kw/deg          - ');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
res.all = [par.Uvec' P'*1e-3 beta'*180/pi dpdbeta_FrozenWake' Cp']';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(fid,'%5.1f  %9.2f  %6.2f  %10.1f %10.3f\n',res.all);
fprintf(fid,'\n  KK = %3.2f deg,  based on frozen wake',KK);
fprintf(fid,'\n  dP/dTheta|0deg  = %3.2f kW/deg',dPdpitch_0);
fprintf(fid,'\n  dP/dTheta|10deg = %3.2f kW/deg',dPdpitch_10);
fprintf(fid,'\n  dP/dTheta|20deg = %3.2f kW/deg',dPdpitch_20);
fprintf(fid,'\n\n  dQ/dTheta|0deg  = %3.2f kW/deg',dQdpitch_0);
fprintf(fid,'\n  dQ/dTheta|10deg = %3.2f kW/deg',dQdpitch_10);
fprintf(fid,'\n  dQ/dTheta|20deg = %3.2f kW/deg',dQdpitch_20);


fclose(fid);

save('data/Workspace.mat')
end
