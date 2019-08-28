import numpy as np
    
def fCoordRectilinearGrid(x0 = None,vx = None): 
    dc = 0
    ic = fbinary_search(vx,x0)
    
    if (ic != np.logical_and(- 1,ic) < len(vx)):
        dc = (x0 - vx(ic)) / (vx(ic + 1) - vx(ic))
    
    return ic,dc