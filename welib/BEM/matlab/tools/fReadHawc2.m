function [sig,Freq,Time,Flag,Binary,label_sensor,ylabels,description_sensor] = fReadHawc2(FileName)
% Reads HAWC2 binary and ascii results file
% ------------------------------------------------------------------------
% [sig,Freq,Time,Flag,Binary] = ReadHAWC2(FileName);
% filename should be without extension
% e.g. FileName = '..\res\test'
% output:
% sig:      all data in file
% Freq:     sampling frequency
% Time:     total simulation time in file
% Flag:     0 if file do not exist, 1 if reading scucceed,
%           2 if some data but not for full time span
% Binary:   1 if fils is in binary format, 0 if ascii
% ------------------------------------------------------------------------
% Made by Bjarne Skovmose Kalles�e
% 24/9-2008
% modified 18/12-2008: file error handelig added
% modified 26/2-2009: error msg for empty date file
% ------------------------------------------------------------------------
%% reading scale factors from *.sel file
fid = fopen([FileName,'.sel'], 'r');
if fid == -1
    Flag = 0;
    disp('  ')
    disp('==============================================================')
    disp(['file "',FileName,'" could not be found'])
    disp('--------------------------------------------------------------')
    sig = [];Freq = 0;Time=0;Binary=[];
    return
end
fgets(fid); fgets(fid); fgets(fid); % skip 3 lines
fgets(fid); fgets(fid); fgets(fid); % skip 3 lines
fgets(fid); fgets(fid); % skip 2 lines
tline = fscanf(fid,'%f'); % read line number 9
N = tline(1);
Nch = tline(2);
Time = tline(3);
FileType = fscanf(fid,'%s');
fclose(fid);
Freq = N/Time;

%% Binary data
if FileType(1:6) == 'BINARY'
    Binary = 1;
    fid = fopen([FileName,'.dat'], 'r');
    sig = fread(fid,[N,Nch],'int16');
    if isempty(sig)
        Flag = 0;
        disp('  ')
        disp('==============================================================')
        disp(['no data in file "',FileName,'.dat"'])
        disp('--------------------------------------------------------------')
        sig = [];Freq = 0;Time=0;Binary=[];
        return
    end
    ScaleFactor = dlmread([FileName,'.sel'],'',[9+Nch+5 0 9+2*Nch+4 0]);
    sig = sig*diag(ScaleFactor);
    fclose(fid);
    Flag = 1;
    
    %% ascii data
elseif FileType(1:5) == 'ASCII'
    Binary = 0;
    sig = dlmread([FileName,'.dat']);
    if isempty(sig)
        Flag = 0;
        disp('  ')
        disp('==============================================================')
        disp(['no data in file "',FileName,'.dat"'])
        disp('--------------------------------------------------------------')
        sig = [];Freq = 0;Time=0;Binary=[];
        return
    end
    if length(sig(:,1)) < N
        Time = length(sig(:,1))/Freq;
        Flag = 2;
        disp('  ')
        disp('==============================================================')
        disp(['not full data in file "',FileName,'.dat"'])
        disp('--------------------------------------------------------------')
%         return
    end
    Flag = 1;
    
    %% error
else
    disp('unknown file type')
    Flag = 0;
end



%% ADD Ewan Machefaux to get labels in an array for giving title to plot
fid = fopen([FileName,'.sel'], 'r');
fgets(fid); fgets(fid); fgets(fid); % skip 3 lines
fgets(fid); fgets(fid); fgets(fid); % skip 3 lines
fgets(fid); fgets(fid); fgets(fid); % skip 3 lines
fgets(fid); fgets(fid); fgets(fid); % skip 3 lines

for i=1:length(sig(1,:))
    %    label_sensor{i,:} = fscanf(fid,'%s'); % read
    A = textscan(fid,'%f\t %27c\t %s %[^\n]$\n'); % read line number i + offset
    fgets(fid);
    label_sensor{i,:}=A{2};
    ylabels{i,:}=A{3};
    description_sensor{i,:}=A{4};
end
fclose(fid);




