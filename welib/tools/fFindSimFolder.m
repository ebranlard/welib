function [ Out ] = fFindSimFolder( DataPath,sim,WS,RPM,PITCH,distinct )
Out.DataPath=DataPath;

%% finding simulation path
listing = dir(Out.DataPath);
listing=listing([listing(:).isdir]);
if(length(listing)==2)
    error('Check Out data path');
end

s=regexpi({listing.name},sim);
I=1;found=0;
while I<=length(s)
    if(length(s{I})==1)
        found=1; break;
    end
    I=I+1;
end
if(~found)
    warning(sprintf('Simulation %s not found in folder %s',sim,DataPath))
else
    Out.WTPath=sprintf('%s/',listing(I).name);

    %% Looking at cases present
    listing = dir([Out.DataPath Out.WTPath]);
    listing=listing([listing(:).isdir]);
    s=regexpi({listing.name},sprintf('Vw%04.1f_%04.1frpm',WS,RPM));
    I=1;found=[];
    while I<=length(s)
        if(length(s{I})==1)
            found=[found I];
        end
        I=I+1;
    end
    if(isempty(found))
        %error(sprintf('Case not found in %s',[Out.DataPath Out.SimPath]))
        return
    end
    % TO BE IMPROVE TO INCLUDE PITCH CHECK
    if(length(found)>1)
        % using distinction
        listing=listing(found);
        s=regexpi({listing.name},distinct);
        I=1;found=0;
        while I<=length(s)
            if(length(s{I})==1)
                found=1; break;
            end
            I=I+1;
        end
        if(~found)
            error(sprintf('Case not found in %s',[Out.DataPath Out.WTPath]))
        end
        Out.SimPath=sprintf('%s%s/',Out.WTPath,listing(I).name);
    else
        Out.SimPath=sprintf('%s%s/',Out.WTPath,listing(found(1)).name);
    end
end
end

