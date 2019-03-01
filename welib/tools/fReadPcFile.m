function Profiles= fReadPcFile(PcFileName,PcSet)
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
fprintf('PC set #%d - Label: %s\n', PcSet, strtrim(Label))
k=0;
for i=1:NrSubSet
    temp = fscanf(fid,'%f',3);    fgets(fid);
    TempVec(i,:) = temp(3:-1:2);
%     for j=1:temp(2)
%         k=k+1;
%         Profiles.Data(k,:) = fscanf(fid,'%f',4);fgets(fid);
%     end
    Profiles.Data{i}{1} = cell2mat(textscan(fid,'%f %f %f %f %*[^\n]\n',temp(2)));    
    Profiles.Re{i}(1)=1;
end
fclose(fid);
% Profiles.AeFileName = AeFileName;
Profiles.Files = PcFileName;
Profiles.n=length(TempVec(:,1));
Profiles.thickness_rel=TempVec(:,1)';
Profiles.ndata=floor(TempVec(:,2)');
%Profiles.I99=Profiles.rel_thickness<100;
