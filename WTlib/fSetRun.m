function [ Run ] = fSetRun(varargin)
if(nargin==1)
    WS=varargin{1}(1);
    RPM=varargin{1}(2);
    PITCH=varargin{1}(3);
    if(length(varargin{1})==3)
        YAW=PITCH*0;
    else
        YAW=varargin{1}(4);
    end
else
    WS=varargin{1};
    RPM=varargin{2};
    PITCH=varargin{3};
    if(nargin==4)
        YAW=0;
    end
end

Run.Name=sprintf('%04.1fmps_%04.1frpm_%05.2fpitch_%04.1fyaw',WS,RPM,PITCH,YAW);    
Run.UniqueName='';
Run.WS=WS;
Run.RPM=RPM; 
Run.PITCH=PITCH;
Run.YAW=YAW;
Run.Date=datevec(now);
Run.sDate=sprintf('%04d-%02d-%02d_%02d%02d',Run.Date(1:end-1));
Run.UniqueName=Run.sDate;
end

