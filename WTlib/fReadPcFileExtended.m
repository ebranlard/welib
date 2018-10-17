function Profiles = fReadPcFileExtended(PcFileName,PcSet)
if nargin==1
    PcSet = 1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read profile Profiles file for aero info
%-------------------------------------------------------------------------
fid = fopen(PcFileName);
if fid == -1
    disp('  ')
    disp('==============================================================')
    disp(['file "',PcFileName,'" could not be found'])
    disp('--------------------------------------------------------------')
    return
end
NrSet = fscanf(fid,'%d',1);
fgets(fid);
if NrSet < PcSet
    disp('  ')
    disp('==============================================================')
    disp(['Not enough ae sets in file "',PcFileName,'"'])
    disp('--------------------------------------------------------------')
    return
end
for i=1:PcSet-1
    NrSubSet = fscanf(fid,'%d',1);
    fgets(fid);
    for j=1:NrSubSet;
        temp = fscanf(fid,'%d',2);fgets(fid);
        for k=1:temp(2)
            fgets(fid);
        end
    end
end
NrSubSet = fscanf(fid,'%d',1);
Label=fgets(fid);
% fprintf('PC set #%d - Label: %s\n', PcSet, strtrim(Label))
for i=1:NrSubSet
    temp = fscanf(fid,'%f',3);   
    Headers=fgets(fid);
    ncolumns=textscan(Headers,'%s','delimiter',' ','MultipleDelimsAsOne',1);
    ncolumns=length(ncolumns{:});
    TempVec(i,:) = temp(3:-1:2);
    ndatasets=ncolumns/6;
    n=temp(2)*ncolumns;
    bigpolar=textscan(fid,'%f',n);
    bigpolar=reshape(bigpolar{:},[ncolumns n/ncolumns])';
    k=1;
    kk=1;
    for j=1:ndatasets
        if(bigpolar(1,2+(j-1)*6)==9 || ndatasets==1)            
             Profiles.Data{i}{k}=bigpolar(:,3+(j-1)*6:6+(j-1)*6 );
             Profiles.Re{i}(k)=bigpolar(1,1+(j-1)*6)/(10^6); % in million!!!
            k=k+1;
        end
        if(bigpolar(1,2+(j-1)*6)==1 || ndatasets==1)            
             Profiles.DataRough{i}{kk}=bigpolar(:,3+(j-1)*6:6+(j-1)*6 );
             Profiles.ReRough{i}(kk)=bigpolar(1,1+(j-1)*6)/(10^6); % in million!!!
            kk=kk+1;
        end        
    end
end
fclose(fid);
Profiles.Files = PcFileName;
Profiles.n=length(TempVec(:,1));
Profiles.thickness_rel=TempVec(:,1)';
Profiles.ndata=TempVec(:,2)';
Profiles.N=[1 9]; % to be improved
