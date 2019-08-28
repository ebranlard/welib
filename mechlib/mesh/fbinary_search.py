import numpy as np
    
def fbinary_search(x = None,x0 = None): 
    # x a sorted vector (typically a grid vector)
# Performs binary search and return the largest index such that x(i) <= x0
    if (len(varargin) == 0):
        unit_test()
        return i_inf
    
    i_inf = 1
    i_sup = np.asarray(x).size
    # Safety test
    if (x0 < x(1)):
        i_inf = - 1
        return i_inf
    
    if (x0 >= x(end())):
        i_inf = i_sup
        return i_inf
    
    # We loop until we narrow down to one index
    while (i_inf + 1 < i_sup):

        mid = (int(np.floor((i_inf + i_sup) / 2)))
        if (x(mid) <= x0):
            i_inf = mid
        else:
            i_sup = mid

    
    return i_inf
    
    
def unit_test(): 
    fbinary_search(np.arange(1,10+1),1) == 1
    fbinary_search(np.arange(1,10+1),10) == 10
    fbinary_search(np.arange(1,10+1),11) == 10
    fbinary_search(np.arange(1,10+1),0) == - 1
    fbinary_search(np.arange(1,10+1),5.1) == 5
    fbinary_search(np.arange(1,10+1),5) == 5
    return
    
    return i_inf