if ~exist('Wind','var') || ~isfield(Wind,'bOptionsSet')
    %------------------------------------
    wind_SetDefault
    %------------------------------------
end
if ~Wind.bOptionsSet
    %% Loading Current Options
    [Wind]=fReadMasterFileSection(CommonData.MasterFile,'Wind',Wind);
end
