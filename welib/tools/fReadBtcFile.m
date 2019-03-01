function Data = fReadBtcFile(BtcFileName);
MainPath=dirname(BtcFileName);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read BTC file for pitch axis info
%-------------------------------------------------------------------------
fid = fopen(BtcFileName);
if fid == -1
    disp('  ')
    disp('==============================================================')
    disp(['file "',BtcFileName,'" could not be found'])
    disp('--------------------------------------------------------------')
    return
end
tline = fgets(fid);
iP=1;
while tline ~= -1
    tline = fgets(fid);
    if ~isempty(strfind(tline,'PATH DATA_DEFAULT'))
        Data.Path=strtrim(sscanf(tline,'%*s DATA_DEFAULT %s;'));
    end
    if ~isempty(strfind(tline,'FILENAME AERODYNAMIC_LAYOUT_FILE'))
        Data.AeFileName=strtrim(sscanf(tline,'%*s AERODYNAMIC_LAYOUT_FILE %s;'));
    end
    if ~isempty(strfind(tline,'FILENAME PROFILE_COEFFICIENT_FILE'))
        Data.PcFileName=strtrim(sscanf(tline,'%*s PROFILE_COEFFICIENT_FILE %s;'));
    end
    if ~isempty(strfind(tline,'NODE BLADE_ALL'))
        Data.PitchAxis(iP,:)=sscanf(tline,'%*s BLADE_ALL %f %f %f %f ;')';
        iP=iP+1;
    end
    if ~isempty(strfind(tline,'AERODYN BLADES_LAYOUT'))
        Data.AeSet=sscanf(tline,'%*s BLADES_LAYOUT %d');
    end   
end
fclose(fid);
 