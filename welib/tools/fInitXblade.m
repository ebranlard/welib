function [WT]=fInitXblade(Files,WT,Opts)
%File order is important!!!!
if(~isfield(Opts,'Extended'))
    Opts.Extended=0;
end

WT.Sources.ParamFile=Files{1};
WT.Sources.PcFile=Files{2};
WT.Sources.Extended=Opts.Extended;
%% Scanning files
[AeData WT]=fReadXbladeParamFile(Files{1},WT);  %!!! Modifies a lot of other global variables
PcSet=AeData(1,7);
if(Opts.Extended)
    [Profiles]=fReadPcFileExtended(Files{2},PcSet);
else
    [Profiles]=fReadPcFile(Files{2},PcSet);
end
%% Radial parameters,XBlade, they already include rhub
WT.Profiles=Profiles;

WT.Sources.Rotor.r             = AeData(:,1);  
WT.Sources.Rotor.chord         = AeData(:,2);
WT.Sources.Rotor.thickness_rel = AeData(:,3);
WT.Sources.Rotor.twist         = AeData(:,4);
