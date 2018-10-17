function WT=fInitHawc(Files,WT,Opts)


%% Scanning files
Ihtc=whichfile(Files, '(\.htc)$');
Ispec=whichfile(Files, '(Spec\.dat)$');
Ibtc=whichfile(Files, '(\.btc)$');
Ipc =whichfile(Files, 'pc');
Iae =whichfile(Files, '(\.ae)$');

format='';
if(isempty(Ibtc))
    format='hawc';
else
    format='bhawc';
end

%% Read spec if present
if(~isempty(Ispec))
%     disp('Reading Spec file')
    WT=fReadSpec(WT,Files{Ispec});
end




%% Now for the hawc files
if(isempty(Ipc)||isempty(Iae))
%     disp('Reading full body file')
    if(isequal(format,'hawc'))
        [ Data WT.Profiles ] = fReadHtcFileAll(Files{Ihtc},'blade1');
    elseif(isequal(format,'bhawc'))
        [ Data WT.Profiles ] = fReadBtcFileAll(Files{Ibtc});
    else
        error('Provide htc or btc file')
    end
else
    nCol=4;
    if(isequal(format,'hawc'))
        Data = fReadHtcFile(Files{Ihtc},'blade1');
    elseif(isequal(format,'bhawc'))
        Data = fReadBtcFile(Files{Ibtc});
        nCol=7;
    else
        error('Provide htc or btc file')
    end
    if(~isempty(Ipc) && ~isempty(Iae))
        if(~isfield(Data,'AeSet'))
            Data.AeSet=1;
        end
        Data.AeData = fReadAeFile(Files{Iae},Data.AeSet(1),nCol);  %chord and thickness
        PcSet=Data.AeData(1,nCol); % replaced nCol by 7
        WT.Profiles = fReadPcFile(Files{Ipc},PcSet);
    else
        error('Aero files missing')
    end
end

%% HAWC and bHAWC have AeData that does not contain the hub radius
if(WT.Rotor.rhub==0)
    warning('Hub radius is 0.. Usually hawc files do not contain hub radius')
    WT.Rotor.rhub=0;
end

WT.Sources.Rotor.r             = Data.AeData(:,1)+WT.Rotor.rhub;  
WT.Sources.Rotor.chord         = Data.AeData(:,2);
WT.Sources.Rotor.thickness_rel = Data.AeData(:,3);

if(isequal(format,'hawc'))
    % Hawc convention we use the twist from the Htc file
    Stations=Data.PitchAxis;
    twist=Stations(:,5);
    rtwist=Stations(:,4)+WT.Rotor.rhub;
else
    % BHAWC convention, we use the twist from the Ae file
    twist=Data.AeData(:,4);
    rtwist=WT.Rotor.r;
end
% Dealing with sign
if(mean(twist)<0) sign=-1; else sign =1; end

% Dealing with problem of interpolation
if max(rtwist)<max(WT.Sources.Rotor.r)
    disp('! For twist interpolation, last blade section in htc file should be at bigger the last of the ae file. Im replacing it....')
    rtwist(end)=max(WT.Sources.Rotor.r);
end
if min(rtwist)>min(WT.Sources.Rotor.r)
    disp('! For twist interpolation, first blade section in htc file should be at smaller (usually 0) than the one in the ae file. Im replacing it....')
    rtwist(1)=WT.Sources.Rotor.r(1);
end

%  Interpolating twist
WT.Sources.Rotor.twist     = interp1(rtwist,twist,WT.Sources.Rotor.r)*sign;    
