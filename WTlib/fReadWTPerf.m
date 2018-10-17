function []=fReadWTPerf(FileName)
global Algo Profiles Rotor Aero Simulation Controller Nacelle Environment 
fid=fopen(FileName);
% Header
fgetl(fid);
fgetl(fid);
fgetl(fid);
% Input config
fgetl(fid);
A=cell2mat(textscan(fid,'%c %*[^\n]\n',3));
if(~istrue(A(2))) error('Non Dimensional Data Not Supported'); end
if(~istrue(A(3))) error('non Metric Not Supported'); end
% Model config
fgetl(fid);
A=cell2mat(textscan(fid,'%f %*[^\n]\n',4));
Algo.NumSect=A(1);
Algo.nbIt=A(2);
Algo.aTol=A(3);
Algo.swTol=A(4);
% Algo config
fgetl(fid);
A=cell2mat(textscan(fid,'%c %*[^\n]\n',8));
Algo.TipLoss=istrue(A(1));
Algo.HubLoss=istrue(A(2));
Algo.Swirl=istrue(A(3));
Algo.SkewWake=istrue(A(4));
Algo.AdvanceBrakeState=istrue(A(5));
if(istrue(A(6)))
    Algo.correction='Unknown';
else
    Algo.correction='AeroDyn';
end
Algo.AIDrag=istrue(A(7));
Algo.TIDrag=istrue(A(8));
% Turbine Data
fgetl(fid);
A=cell2mat(textscan(fid,'%f %*[^\n]\n',8));
Rotor.nB=A(1);
Rotor.R=A(2);
Rotor.rhub=A(3);
Rotor.BladeLength=Rotor.R-Rotor.rhub;
Rotor.cone=A(4);
Nacelle.tilt=A(5);
Controller.yaw=A(6);
Rotor.HubHeight=A(7);
Rotor.ne=A(8);
% Blade Data
fgetl(fid);
A=cell2mat(textscan(fid,'%f %f %f %f %*[^\n]\n',Rotor.ne));
Rotor.r=A(:,1);
Rotor.twist=A(:,2);
Rotor.chord=A(:,3);
Rotor.ProfileSet=ones(3,Rotor.ne);
Rotor.ProfileSet(2,:)=A(:,4);
Rotor.ProfileSet(3,:)=A(:,4);
% Aerodynamic data
fgetl(fid);
A=cell2mat(textscan(fid,'%f %*[^\n]\n',3));
Environment.rho=A(1);
Environment.KinVisc=A(2);
Aero.Wind.nu=A(3);
Algo.UseCm=istrue(cell2mat(textscan(fid,'%c %*[^\n]\n',1)));
Profiles.n=cell2mat(textscan(fid,'%d %*[^\n]\n',1));
A=textscan(fid,'%s %*[^\n]\n',Profiles.n);
A=A{1};
Profiles.Files=A;
% IO Data
fgetl(fid);
A=cell2mat(textscan(fid,'%c %*[^\n]\n',4));
Simulation.KiloOut=istrue(A(2));
if(istrue(A(4))) error('Input TSR not Supported'); end
A=textscan(fid,'%s %*[^\n]\n',1);
A=A{1};
if(~isequal(A{1},'"mps"')) error('Wind Speed will be in mps!'); end
% Combined case
fgetl(fid);
A=cell2mat(textscan(fid,'%d %*[^\n]\n',1)); 
Simulation.CombinedCase.n=A(1);
fgetl(fid);
if(Simulation.CombinedCase.n>0)
    A=cell2mat(textscan(fid,'%f,%f,%f,%f %*[^\n]\n',Simulation.CombinedCase.n)); 
    Simulation.CombinedCase.Cases=A;
end
% Parametric Analysis
fgetl(fid);
A=cell2mat(textscan(fid,'%d %*[^\n]\n',3)); % row col tab
A=cell2mat(textscan(fid,'%c %*[^\n]\n',5)); % output
A=cell2mat(textscan(fid,'%f,%f,%f %*[^\n]\n',3));
Simulation.Parametric.Pitch=A(1,:);
Simulation.Parametric.Omega=A(2,:)/60*2*pi;
Simulation.Parametric.WS=A(3,:);
fclose(fid);

%% Read Airfoils
MainPath=dirname(FileName);
Profiles.Data=cell(1,length(Profiles.Files));
Profiles.Re=cell(1,length(Profiles.Files));
for i=1:length(Profiles.Files)
    s = regexprep(Profiles.Files(i), '"', '');
    fid=fopen([MainPath,s{1}]);
    A=textscan(fid,'%s %*[^\n]\n',1);
    if(isempty(cell2mat(regexpi(A{1},'aerodyn'))))
        error(['Not AeroDyn File: ',s{1} ])
    else
        fgetl(fid);
        fgetl(fid);
        nTab=cell2mat(textscan(fid,'%d %*[^\n]\n',1)) ;
        Profiles.Data{i}=cell(1,nTab);
        Profiles.Re{i}=zeros(1,nTab);
        for j=1:nTab
            Re=cell2mat(textscan(fid,'%f %*[^\n]\n',1));
            cell2mat(textscan(fid,'%f %*[^\n]\n',7)) ; 
            A=cell2mat(textscan(fid,'%f %f %f %f %*[^\n]\n'));
            A(isnan(A))=-999;
            Profiles.Data{i}{j}=A;
            Profiles.Re{i}(j)=Re;
            fgetl(fid);
        end
    end
    
    fclose(fid);
end
