function Time = fInitTime(dt,tmax,t0)
if nargin==0
    dt=0; tmax=0; t0=0;
elseif nargin==1
    tmax=0;t0=0;
elseif nargin==2
    t0=0;
end


Time=TimeManager(tmax,dt);
% Time.dt=dt;
% Time.tmax=tmax;
% Time.t0=t0;
% Time.t=t0;
% Time.itime=0;

    
