function AeData = fReadAeFile(AeFileName,AeSet)
if nargin==1
    AeSet = 1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read Ae file for aero info
%-------------------------------------------------------------------------
fid = fopen(AeFileName);
if fid == -1
    disp('  ')
    disp('==============================================================')
    disp(['file "',AeFileName,'" could not be found'])
    disp('--------------------------------------------------------------')
    return
end
NrSet = fscanf(fid,'%d',1);fgets(fid);
if NrSet < AeSet(1)
    disp('  ')
    disp('==============================================================')
    disp(['Not enough ae sets in file "',AeFileName,'"'])
    disp('--------------------------------------------------------------')
    return
end
for i=1:AeSet(1)-1
    temp = fscanf(fid,'%d',2);fgets(fid);
    for j=1:temp(2);
        fgets(fid);
    end
end
temp = fscanf(fid,'%d',2);fgets(fid);
for i=1:temp(2)
    AeData(i,:) = fscanf(fid,'%f',4);
end
PcSet = AeData(1,4);
fclose(fid);


