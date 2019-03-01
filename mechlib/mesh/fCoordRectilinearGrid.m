function [ic dc]=fCoordRectilinearGrid(x0,vx)
    dc=0;
    ic=fbinary_search(vx,x0); % ic can be -1
    if(ic~=-1 & ic<length(vx))
        dc=(x0-vx(ic))/(vx(ic+1)-vx(ic));
    end
