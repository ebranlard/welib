function [ Data Profiles ] = fReadHtcFileAll( HtcFileName,BladeBodyName)

Data = fReadHtcFile(HtcFileName,BladeBodyName);
Data.AeData = fReadAeFile(Data.AeFileName,Data.AeSet(1),4);
PcSet=Data.AeData(1,4);
Profiles = fReadPcFile(Data.PcFileName,PcSet);

