function fExportWTPerf(WT,Simulation,Wind,Algo,DirPath)
if(isequal(DirPath,''))
    DirPath='data/export_wt_perf/';
end
mkdir(DirPath);
mkdir([DirPath,'airfoils/']);
fid=fopen([DirPath,Simulation.Name,'.wtp'],'w');
fprintf(fid,'-----  WT_Perf Input File ---------------------\n');
fprintf(fid,'%s\n',Simulation.Name);
fprintf(fid,'Exported by Manu\n');
fprintf(fid,'-----  Input Configuration  ---------------------\n');
fprintf(fid,'True\t\t Echo:\n');
fprintf(fid,'True\t\t DimenInp:\n');
fprintf(fid,'True\t\t Metric:\n');
fprintf(fid,'-----  Model Configuration  ---------------------\n');
fprintf(fid,'%d\t\t\t\t\t NumSect:\n',Algo.NumSect);
fprintf(fid,'%d\t\t\t\t MaxIter:\n',Algo.nbIt);
fprintf(fid,'%e\t\t ATol:\n',Algo.aTol);
fprintf(fid,'%e\t\t SWTol:\n',Algo.swTol);
fprintf(fid,'-----  Algo Configuration  ---------------------\n');
fprintf(fid,'%s\t\t TipLoss:\n',totrue(Algo.bTipLoss));
fprintf(fid,'%s\t\t HubLoss:\n',totrue(Algo.bHubLoss));
fprintf(fid,'%s\t\t Swirl:\n',totrue(Algo.bSwirl));
fprintf(fid,'%s\t\t SkewWake:\n',totrue(Algo.bSkewWake));
fprintf(fid,'%s\t\t AdvBrake:\n',totrue(Algo.bAdvanceBrakeState));
fprintf(fid,'%s\t\t IndProp:\n',totrue(0));       % !!! IndProp always to false
fprintf(fid,'%s\t\t AIDrag:\n',totrue(Algo.bAIDrag));
fprintf(fid,'%s\t\t TIDrag:\n',totrue(Algo.bTIDrag));
fprintf(fid,'-----  Turbine Data  ---------------------\n');
fprintf(fid,'%d\t\t NumBlade:\n',WT.Rotor.nB);
fprintf(fid,'%f\t\t RotorRad:\n',WT.Rotor.R);
fprintf(fid,'%f\t\t HubRad:\n',WT.Rotor.rhub);
fprintf(fid,'%f\t\t PreCone:\n',WT.Rotor.cone);
fprintf(fid,'%f\t\t Tilt:\n',Nacelle.tilt);
fprintf(fid,'%f\t\t Yaw:\n',Controller.yaw);
fprintf(fid,'%f\t\t HubHt:\n',WT.Rotor.HubHeight);
fprintf(fid,'%d\t\t NumSeg:\n',WT.Rotor.ne);
fprintf(fid,'   RElm   Twist    Chord  AFfile  PrntElem\n');
M=zeros(WT.Rotor.ne,4);
M(:,1)=WT.Rotor.r;
M(:,2)=WT.Rotor.twist;
M(:,3)=WT.Rotor.chord;
M(:,4)=WT.Rotor.ProfileSet(2,:); % to be improved?
fprintf(fid,'%f\t%f\t%f\t%d\tTrue\n',M');
fprintf(fid,'-----  Aero Data  ---------------------\n');
fprintf(fid,'%f\t\t\t Rho:\n',WT.DefaultEnvironment.rho);
fprintf(fid,'%e\t\t KinVisc:\n',WT.DefaultEnvironment.KinVisc);
fprintf(fid,'%f\t\t\t ShearExp:\n',Wind.nu);
fprintf(fid,'%s\t\t\t\t UseCm:\n',totrue(Algo.UseCm)); 
fprintf(fid,'%d\t\t\t\t\t NumAF:\n',Profiles.n);
for i=1:Profiles.n
    FileNameRel=['airfoils/',Simulation.Name,'_aero_',num2str(i),'.dat'];
    FileName=[DirPath,'airfoils/',Simulation.Name,'_aero_',num2str(i),'.dat'];
    fprintf('%s\n',FileName);
    fprintf(fid,'"%s"\n',FileNameRel);
    % More on creating airfoils!!!
    fidb=fopen(FileName,'w');
    fprintf(fidb,'AeroDyn airfoil file.\n');
    fprintf(fidb,'Exported by Manu\n');
    fprintf(fidb,'Simulation %s\n',Simulation.Name);
    nTab=length(Profiles.Re{i});
    fprintf(fidb,'%d\tNumber of airfoil tables\n',nTab);
    for iRe=1:nTab
        fprintf(fidb,'%f\tReynolds\n',Profiles.Re{i}(iRe));
        warning('Add a line here for version alpha of WTperf with Control Settings')
        fprintf(fidb,'-1 Stall Angle\n');
        fprintf(fidb,'-1 Zero Cn angle\n');
        fprintf(fidb,'-1 Cn slope for 0 lift\n');
        fprintf(fidb,'-1 Cn extrapolated\n');
        fprintf(fidb,'-1 Cn at stall value\n');
        fprintf(fidb,'-1 Angle of attack minimum Cd\n');
        fprintf(fidb,'-1 Minimum Cd\n');
        fprintf(fidb,'%f\t%f\t%f\t%f\n',Profiles.Data{i}{iRe}(:,1:4)');
        fprintf(fidb,'EOT\n');
    end
    fclose(fidb);
end
fprintf(fid,'-----  Output Config Data  ---------------------\n');
fprintf(fid,'True\t\t TabDel:\n'); %tabdel
fprintf(fid,'%s\t\t KFact:\n',totrue(Simulation.KiloOut));
fprintf(fid,'True\t\t WriteBED:\n');%blade element data
fprintf(fid,'False\t\t InputTSR:\n');%tsr
fprintf(fid,'"mps"\t\t SpdUnits:\n');%tsr
fprintf(fid,'-----  Combined Case  ---------------------\n');
nCases=size(Simulation.CombinedCases,1);
fprintf(fid,'%d\t\t NumCases:\n',nCases);
fprintf(fid,'WS or TSR   RotSpd   Pitch \n');
if(nCases>0)
    Cases=Simulation.CombinedCases(:,1:3).*[1 60/(2*pi) 1];
    fprintf(fid,'%f %f %f\n',Cases'); 
end
fprintf(fid,'-----  Parametric Analysis  ---------------------\n');
fprintf(fid,'3\t\t\t ParRow: 1rpm 2 pitch 3speed\n');
fprintf(fid,'2\t\t\t ParCol:\n');
fprintf(fid,'1\t\t\t ParTab:\n');
fprintf(fid,'True\t\t OutPwr:\n');
fprintf(fid,'True\t\t OutCp:\n');
fprintf(fid,'True\t\t OutTrq:\n');
fprintf(fid,'True\t\t Outflp:\n');
fprintf(fid,'True\t\t OutThr:\n');
p=Simulation.ParametricPitch;
P=[min(p) max(p) (max(p)-min(p))/length(p);
p=Simulation.ParametricRPM;
RPM=[min(p) max(p) (max(p)-min(p))/length(p);
p=Simulation.ParametricWS;
WS=[min(p) max(p) (max(p)-min(p))/length(p);
fprintf(fid,'%f,%f,%f\t\t PitchSt, PitEnd, PitDel:\n',P); 
fprintf(fid,'%f,%f,%f\t\t OmgSt, OmgEnd, OmgDel:\n',RPM);  % oups I forgor if it's Omega or RPM!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!check previous versions!!
fprintf(fid,'%f,%f,%f\t\t SpdSt, SpdEnd, SpedDel:\n',WS); 
fclose(fid);
disp(['Data exported to :',DirPath])
