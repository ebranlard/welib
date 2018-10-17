classdef TimeManager < handle
    %TIMEMANAGER Summary of this class goes here
    %   Detailed explanation goes here
    
     properties (SetAccess=private)
        time;
        tmax;
        dt;
        t0;
        dt_step;
        dt_substep;
        itime;
        ntmax;
        current_step_time;
        previous_step_time;
        previous_substep_time;
    end
    properties
    end
    methods 
        % constructor 
        function obj=TimeManager(tmax,dt);
           obj.itime = -1   ; 
           obj.dt=dt;
           obj.setTmax(tmax);
           obj.setStepTime(0) ; 
           obj.t0=0 ;  % One day
        end

        function setTmax(obj,t_max)
               obj.tmax=t_max;
               obj.ntmax=floor(t_max/obj.dt)+1;
        end
        function setDt(obj,dt)
               obj.dt=dt;
               obj.ntmax=floor(obj.tmax/obj.dt)+1;
        end

        function setSubStepTime(obj,t_in)
               obj.previous_substep_time=obj.time; 
               obj.dt_substep=t_in-obj.time; 
               obj.time=t_in; 
        end

        function setStepTime(obj,t_in)
               obj.previous_substep_time = obj.time ; 
               obj.previous_step_time = obj.current_step_time ; 
               obj.itime              = obj.itime+1           ; 
               obj.current_step_time  = t_in                  ; 
               obj.time               = t_in                  ; 
               obj.dt_step=obj.current_step_time-obj.previous_step_time; 
        end

        % Uses dt to increment time
        function bContinue=IncrementTime(obj);
               new_time=obj.current_step_time+obj.dt;             
               if obj.itime== obj.ntmax-1
                   bContinue=false;
               else
                   bContinue=true;
                   obj.setStepTime(new_time);
               end
%                if new_time>obj.tmax
%                    bContinue=false;
%                else
%                    bContinue=true;
%                    obj.setStepTime(new_time);
%                end
        end

        % Resetting clock to t0
        function resetClock(obj)
            obj.itime=-1;
            obj.time=obj.t0;
            obj.current_step_time=obj.t0;
            obj.setStepTime(obj.t0);
        end

        function txtProgressBar(obj)
            fprintf('%s%02d %%\n',char([8 8 8 8 8]),floor(100*obj.itime/(obj.ntmax-1)));
        end

        function printStep(obj)
            fprintf('--------- Time step %d/%d, time %.2f        \n',obj.itime,obj.ntmax-1,obj.time);
        end

        function t=vt(o)
            t=o.t0+o.dt:o.dt:o.tmax;
        end

    end
    
end

