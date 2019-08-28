import numpy as np
    
def fCoordRegularGrid(x0 = None,xBox = None,dBox = None): 
    C = (x0 - xBox) / dBox
    
    ic = int(np.floor(C)) + 1
    
    dc = C + 1 - ic
    
    return ic,dc