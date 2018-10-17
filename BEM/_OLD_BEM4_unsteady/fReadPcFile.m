function Profiles = fReadPcFile(PcFileName,PcSet)
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
NrSet = fscanf(fid,'%d',1);fgets(fid);
if NrSet < PcSet
    disp('  ')
    disp('==============================================================')
    disp(['Not enough ae sets in file "',PcFileName,'"'])
    disp('--------------------------------------------------------------')
    return
end
for i=1:PcSet-1
    NrSubSet = fscanf(fid,'%d',1);fgets(fid);
    for j=1:NrSubSet;
        temp = fscanf(fid,'%d',2);fgets(fid);
        for k=1:temp(2)
            fgets(fid);
        end
    end
end
NrSubSet = fscanf(fid,'%d',1);fgets(fid);
k=0;
for i=1:NrSubSet
    temp = fscanf(fid,'%d',3);fgets(fid);
    Profiles.ThicknessVec(i,:) = temp(3:-1:2);
    for j=1:temp(2)
        k=k+1;
        Profiles.Data(k,:) = fscanf(fid,'%f',4);fgets(fid);
    end
end
fclose(fid);
% Profiles.AeFileName = AeFileName;
Profiles.PcFileName = PcFileName;
Profiles.n=length(Profiles.ThicknessVec(:,1));
Profiles.rel_thickness=Profiles.ThicknessVec(:,1)';
Profiles.ndata=Profiles.ThicknessVec(:,2)';
Profiles.I99=Profiles.rel_thickness<100;
