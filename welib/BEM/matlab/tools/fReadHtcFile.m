function Data = fReadHtcFile(HtcFileName,BladeBodyName);
MainPath=dirname(HtcFileName);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read HTC file for pitch axis info
%-------------------------------------------------------------------------
fid = fopen(HtcFileName);
if fid == -1
    disp('  ')
    disp('==============================================================')
    disp(['file "',HtcFileName,'" could not be found'])
    disp('--------------------------------------------------------------')
    return
end
tline = fgets(fid);
while tline ~= -1
    tline = fgets(fid);
    % read pitch axis data
    if isempty(strfind(tline,'name'))+isempty(strfind(tline,BladeBodyName)) == 0
        while isempty(strfind(tline,'begin')) | isempty(strfind(tline,'c2_def'))
            tline = fgets(fid);
        end
        tline = fgets(fid);
        I1 = strfind(tline,'nsec')+4;
        I2 = strfind(tline,';')-1;
        nsec = str2num(tline(I1:I2));
        if(isempty(nsec))
            error('Comments are not well supported by script fReadHtcFile.')
        end

        
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
            % if BUG check for lines with begin and aero in the file
            % for instante begin output_at_time aero makes this loop crash
            while isempty(strfind(tline,'end'))+isempty(strfind(tline,'aero')) ~= 0
                tline = fgets(fid);
                if strfind(tline,'nblades')
                    I1 = strfind(tline,'nblades')+7;
                    I2 = strfind(tline,';')-1;
                    Data.Nb = str2num(tline(I1:I2));
                end
                if strfind(tline,'ae_filename')
                    % I really need to program this with regexp, this is nasty
                    % If there is no space before the ; the stuff below will crash
                    % So, adding a space before the ;...
                    tline = strrep(tline,';',' ;');

                    I1 = strfind(tline,'ae_filename')+11;
                    I2 = strfind(tline,';')-1;
                    I3 = strfind(tline(I1:I2),' ');
                    i = 1;
                    while i == I3(i);
                        i=i+1;
                    end
                    temp = tline(I1+i:I2);
                    Data.AeFileName = strtrim([MainPath,temp]);
                end
                if strfind(tline,'pc_filename')
                    % So, adding a space before the ;...
                    tline = strrep(tline,';',' ;');
                    %
                    I1 = strfind(tline,'pc_filename')+11;
                    I2 = strfind(tline,';')-1;
                    I3 = strfind(tline(I1:I2),' ');
                    i = 1;
                    while i == I3(i);
                        i=i+1;
                    end
                    temp = tline(I1+i:I2);
                    Data.PcFileName = strtrim([MainPath,temp]);
                end
                if strfind(tline,'ae_sets')
                    I1 = strfind(tline,'ae_sets')+7;
                    I2 = strfind(tline,';')-1;
                    Data.AeSet = str2num(tline(I1:I2));
                end
            end
        end
    end
end
fclose(fid);
 
