function [ Data Profiles ] = fReadBtcFileAll(BtcFileName)

Data = fReadBtcFile(BtcFileName);
PathSplitChar='\';
if isempty(regexp(Data.Path,'[\\]'))
    PathSplitChar='/';
end
if isempty(regexp(Data.AeFileName,'[\\]|[//]'))
    Data.AeFileName=sprintf('%s%cBH_Aerodyn_Ae.%s',Data.Path,PathSplitChar,Data.AeFileName);
end
if isempty(regexp(Data.PcFileName,'[\\]|[//]'))
    Data.PcFileName=sprintf('%s%cBH_Aerodyn_Pc.%s',Data.Path,PathSplitChar,Data.PcFileName);
end
Data.AeData = fReadAeFile(Data.AeFileName,Data.AeSet,7);
PcSet=Data.AeData(1,7);
Profiles = fReadPcFile(Data.PcFileName,PcSet);

