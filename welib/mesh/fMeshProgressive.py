import numpy as np
    
def fMeshProgressive(x_range = None,x_high = None,h_high = None,p = None): 
    x_high = np.arange(x_high(1),x_high(2)+h_high,h_high)
    x_sup = fMeshProgressiveExtent(x_high(end()),x_range(2),h_high,p)
    x_inf = fMeshProgressiveExtent(x_high(1),x_range(1),h_high,p)
    x = np.array([x_inf,x_high,x_sup])
    return x