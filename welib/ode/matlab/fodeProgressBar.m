function status=fodeProgressBar(t,y,flag)
    persistent tf tstart;

    % test function
    if nargin ==0
        clc;
        fodeProgressBar(10,[],'init');
        fodeProgressBar(1,[],[]);     pause(0.1);
        fodeProgressBar(5,[],[]);     pause(0.1);
        fodeProgressBar(10,[],[]);    pause(0.1);
        fodeProgressBar(10,[],'done'); pause(0.1);
        return
    end


    if isempty(flag)
        % Integration steps
        ts=mean(t);
        progress=100*ts/tf;
        textprogressbar(progress);
        status=0;
    else
        switch flag
            case 'init'     % Initializing progress bar
                tstart=tic;
                tf=max(t);
                textprogressbar('ODE integration: ');
            case 'done'     % Finishing status function
                tf=[];
                fprintf([' time: ' num2str(toc(tstart)) 's']);
                tstart=[];
                textprogressbar(''); % Reinit and prints CR
                status=1;
            otherwise
                error('Unknown flag %s',flag);
        end
    end
end


function textprogressbar(c)
    % This function creates a text progress bar. It should be called with a 
    % STRING argument to initialize and terminate. Otherwise the number correspoding 
    % to progress in % should be supplied.
    % INPUTS:   C   Either: Text string to initialize or terminate 
    %                       Percentage number to show progress 
    % OUTPUTS:  N/A
    % Example:  Please refer to demo_textprogressbar.m

    % Author: Paul Proteus (e-mail: proteus.paul (at) yahoo (dot) com)
    % Version: 1.0
    % Changes tracker:  29.06.2010  - First version

    % Inspired by: http://blogs.mathworks.com/loren/2007/08/01/monitoring-progress-of-a-calculation/

    %% Initialization
    persistent strCR;           %   Carriage return pesistent variable

    % Vizualization parameters
    strPercentageLength = 7;   %   Length of percentage string (must be >5)
    strDotsMaximum      = 12;   %   The total number of dots in a progress bar

    %% Main 

    if isempty(strCR) && ~ischar(c)
        % Progress bar must be initialized with a string
        error('The text progress must be initialized with a string');
    elseif isempty(strCR) && ischar(c)
        % Progress bar - initialization
        fprintf('%s',c);
        strCR = -1;
    elseif ~isempty(strCR) && ischar(c)
        % Progress bar  - termination
        strCR = [];  
        fprintf([c '\n']);
    elseif isnumeric(c)
        % Progress bar - normal progress
        c = floor(c);
        percentageOut = [num2str(c) '%%'];
        percentageOut = [percentageOut repmat(' ',1,strPercentageLength-length(percentageOut)-1)];
        nDots = floor(c/100*strDotsMaximum);
        dotOut = ['[' repmat('.',1,nDots) repmat(' ',1,strDotsMaximum-nDots) ']'];
        strOut = [percentageOut dotOut];

        % Print it on the screen
        if strCR == -1
            % Don't do carriage return during first run
            fprintf(strOut);
        else
            % Do it during all the other runs
            fprintf([strCR strOut]);
        end

        % Update carriage return
        strCR = repmat('\b',1,length(strOut)-1);

    else
        % Any other unexpected input
        error('Unsupported argument type');
    end
end
