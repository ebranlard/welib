function [ Data ] = ReadHtcFileAll( user)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read HTC file for pitch axis info
%-------------------------------------------------------------------------
fid = fopen([user.MainPath,user.HtcFileName]);
if fid == -1
    Flag = 0;
    disp('  ')
    disp('==============================================================')
    disp(['file "',user.HtcFileName,'" could not be found'])
    disp('--------------------------------------------------------------')
    return
end
tline = fgets(fid);
while tline ~= -1
    tline = fgets(fid);
    % read pitch axis data
    if isempty(strfind(tline,'name'))+isempty(strfind(tline,user.BladeMainBodyName)) == 0
        while isempty(strfind(tline,'begin')) | isempty(strfind(tline,'c2_def'))
            tline = fgets(fid);
        end
        tline = fgets(fid);
        I1 = strfind(tline,'nsec')+4;
        I2 = strfind(tline,';')-1;
        nsec = str2num(tline(I1:I2));
        for i=1:nsec
            tline = fgets(fid);
            I1 = strfind(tline,'sec')+3;
            I2 = strfind(tline,';')-1;
            Data.PitchAxis(i,:) = str2num(tline(I1:I2));
        end
    end
    %read aerodynamic data
    if isempty(strfind(tline,'begin'))+isempty(strfind(tline,'aero')) == 0
        if isempty(strfind(tline,'aerod')) == 1
            while isempty(strfind(tline,'end'))+isempty(strfind(tline,'aero')) ~= 0
                tline = fgets(fid);
                if strfind(tline,'nblades')
                    I1 = strfind(tline,'nblades')+7;
                    I2 = strfind(tline,';')-1;
                    Data.Nb = str2num(tline(I1:I2));
                end
                if strfind(tline,'ae_filename')
                    I1 = strfind(tline,'ae_filename')+11;
                    I2 = strfind(tline,';')-1;
                    I3 = strfind(tline(I1:I2),' ');
                    i = 1;
                    while i == I3(i);
                        i=i+1;
                    end
                    temp = tline(I1+i:I2);
                    AeFileName = strtrim([user.MainPath,temp]);
                end
                if strfind(tline,'pc_filename')
                    I1 = strfind(tline,'pc_filename')+11;
                    I2 = strfind(tline,';')-1;
                    I3 = strfind(tline(I1:I2),' ');
                    i = 1;
                    while i == I3(i);
                        i=i+1;
                    end
                    temp = tline(I1+i:I2);
                    PcFileName = strtrim([user.MainPath,temp]);
                end
                if strfind(tline,'ae_sets')
                    I1 = strfind(tline,'ae_sets')+7;
                    I2 = strfind(tline,';')-1;
                    AeSet = str2num(tline(I1:I2));
                end
            end
        end
    end
end
fclose(fid);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read Ae file for aero info
%-------------------------------------------------------------------------
fid = fopen(AeFileName);
if fid == -1
    Flag = 0;
    disp('  ')
    disp('==============================================================')
    disp(['file "',AeFileName,'" could not be found'])
    disp('--------------------------------------------------------------')
    return
end
NrSet = fscanf(fid,'%d',1);fgets(fid);
if NrSet < AeSet(1)
    Flag = 0;
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
    Data.AeData(i,:) = fscanf(fid,'%f',4);
end
PcSet = Data.AeData(1,4);
fclose(fid);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read profile data file for aero info
%-------------------------------------------------------------------------
fid = fopen(PcFileName);
if fid == -1
    Flag = 0;
    disp('  ')
    disp('==============================================================')
    disp(['file "',PcFileName,'" could not be found'])
    disp('--------------------------------------------------------------')
    return
end
NrSet = fscanf(fid,'%d',1);fgets(fid);
if NrSet < PcSet
    Flag = 0;
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
    Data.ThicknessVec(i,:) = temp(3:-1:2);
    for j=1:temp(2)
        k=k+1;
        Data.ProfileData(k,:) = fscanf(fid,'%f',4);fgets(fid);
    end
end
fclose(fid);
Data.AeFileName = AeFileName;
Data.PcFileName = PcFileName;
