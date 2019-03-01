function [AeData WT] = fReadXbladeParamFile(FileName,WT)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read file for AE data
%-------------------------------------------------------------------------
fid = fopen(FileName);
if fid == -1
    disp('  ')
    disp('==============================================================')
    disp(['file "',FileName,'" could not be found'])
    disp('--------------------------------------------------------------')
    return
end
tline = fgets(fid);
%%% Parameters  Has to be made more general!!!!
A=cell2mat(textscan(fid,'%*s %f %*[^\n]\n',13));
WT.Spec.Pmax=A(1);
WT.DefaultEnvironment.rho=A(2);
WT.Rotor.nB=A(3);
WT.Rotor.R=A(4);
WT.Rotor.HubHeight=A(5);
%A(6)  mic Tower Dist??
WT.Nacelle.tilt=A(7);
WT.Rotor.cone=A(8);
WT.Spec.OmegaMax=A(9);
WT.Spec.OmegaMin=A(10);

%Algo.TipLoss=A(11);
% corr3D
% skip induction


while tline ~= -1
    tline = fgets(fid);    
    if ~isempty(strfind(tline,'aeMatrix'))
        % now we can read Ae data
        AeData=cell2mat(textscan(fid,'%f %f %f %f %f %f %f %*[^\n]\n'));
    end
end
%%% ... read rpm and pitch ... one day

%%% Post pro
WT.Rotor.BladeLength=AeData(end,1)-AeData(1,1);
WT.Rotor.rhub=AeData(1,1); %%% I nedd to decide, check WTperf

if(length(unique(AeData(:,7)))>1)
    disp('!Warning: Script does not support several Pc Sets, first one will be used')
end
fclose(fid);


