function [Res,Flags,Sty,Codes,LW,Colrs]=fLoadCodesResults(Sims,Params)
%% Loading
Res={}; Flags={}; Sty={}; Codes={}; LW={};Colrs={};

for i=1:length(Sims)
    if (length(Sims{1})==5) 
        error('Sims is 5')
    end
    if (length(Sims{1})==6) 
        error('Sims is 6')
    end
%     try
        iBl     = 2 ; 
        iSty    = 3 ; 
        iColr   = 4 ; 
        iLW     = 5 ; 
        iFlag   = 6 ; 
        iFolder = 7 ; 
        fprintf('Loading:%s\n',Sims{i}{iFolder})
        bAppend=true;
        Blade=fLoadBlade(Sims{i}{iBl});
        switch Sims{i}{1}
            case 'AL-CFD' 
                Res{end+1}= fLoadCFD_AL(Sims{i}{iFolder},'Mexico',1,Params.nB,Blade);
            case '3D-CFD' 
                Res{end+1}= fLoadMeas(Sims{i}{iFolder},Blade);
            case 'VC'
                Res{end+1}= fLoadHawc2(Sims{i}{iFolder},sprintf('%.2f_y%.1f_i4',Params.U0,Params.yaw_hawc),Blade); 
            case 'VCout'
                Res{end+1}=fLoadVC(Sims{i}{iFolder},Blade,Params.U0,Params.yaw_hawc);
            case 'BEM'
                Res{end+1}= fLoadHawc2(Sims{i}{iFolder},sprintf('%.2f_y%.1f_i1',Params.U0,Params.yaw_hawc),Blade);
            case 'MyBEM'
                A= load(Sims{i}{iFolder});
                Res{end+1}=A.MyBEM;
            case 'SHEN'
                Res{end+1}=fLoadShen(Sims{i}{iFolder},Blade);
        end
        if (bAppend)
            Codes{end+1} = Sims{i}{1} ; 
            Sty{end+1}   = Sims{i}{iSty} ; 
            LW{end+1} =    Sims{i}{iLW} ; 
            if(length(Sims{i}{iFlag})==0)
                Flags{end+1} = Sims{i}{iFolder}(end-09:end) ; 
            else
                Flags{end+1} = Sims{i}{iFlag} ; 
            end
            Colrs{end+1} = Sims{i}{iColr} ; 
            Res{end}.folder=Sims{i}{iFolder};
        end
%     catch
%         %
%         warning('Error loading')
%     end
end



%% Postpro or results
for i=1:length(Res)
    Res{i}.Ir=whichvalue(Res{i}.vr,Params.r_ref);
    Res{i}.Uit=zeros(length(Res{i}.psiT),length(Res{i}.vr));
    Res{i}.Uin =zeros(length(Res{i}.psiT),length(Res{i}.vr));
    if(isfield(Res{i},'Ut')) 
        for ir=1:length(Res{i}.vr)
            Res{i}.Uit(:,ir)=Res{i}.Ut(Res{i}.IT,ir)-cosd(Res{i}.psiT)*Params.U0*sind(Params.yaw);
            Res{i}.Uin(:,ir)=Res{i}.Un(Res{i}.IT,ir)-Params.U0*cosd(Params.yaw);
        end
    end
    if(isfield(Res{i},'Power')) 
        Res{i}.CP=Res{i}.Power/Params.P0;
    end
    if(isfield(Res{i},'CP')) 
        Res{i}.CP_last=mean(Res{i}.CP(Res{i}.IT));
    end
    if(isfield(Res{i},'Thrust')) 
        Res{i}.CT=Res{i}.Thrust/Params.T0;
    end
    if(isfield(Res{i},'CT')) 
        Res{i}.CT_last=mean(Res{i}.CT(Res{i}.IT));
    end
end
% Meas= fLoadMeas('data/MexicoData/loads2.1.dat',Blade); 
% Meas.Ir=whichvalue(Meas.vr,r_ref);
