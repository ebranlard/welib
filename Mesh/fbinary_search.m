function [i_inf]=fbinary_search(x,x0)
    % x a sorted vector (typically a grid vector)
    % Performs binary search and return the largest index such that x(i) <= x0
    if(nargin==0); unit_test(); return; end

    i_inf=1;
    i_sup=numel(x);

    % Safety test
    if(x0<x(1))
        i_inf=-1;
        return
    end
    if(x0>=x(end))
        i_inf=i_sup;
        return
    end


    % We loop until we narrow down to one index
    while (i_inf+1<i_sup)
        mid=(floor((i_inf+i_sup)/2));
        if (x(mid)<=x0)
            i_inf=mid;
        else
            i_sup=mid;
        end
    end
end

function unit_test()
    fbinary_search(1:10,1)   == 1
    fbinary_search(1:10,10)  == 10
    fbinary_search(1:10,11)  == 10
    fbinary_search(1:10,0)   == -1
    fbinary_search(1:10,5.1) == 5
    fbinary_search(1:10,5)   == 5
end
