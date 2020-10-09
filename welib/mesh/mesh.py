import numpy as np
##

class Mesh():
    # constructor
    def __init__(o,v1 = None,v2 = None,v3 = None): 
        o.values=None
        if v2 is None and v3 is None:
            # 1D
            o.nDim = 1
            o.v1 = v1
            o.v2 = []
            o.v3 = []
            o.dCell     = np.array([v1[1]-v1[0]])
            o.xMesh_min = np.array([v1(1)])
            o.xMesh_max = np.array([v1(end())])
            o.n = np.array([len(v1)])
            o.bRegular = np.mean(np.diff(v1)) == v1(2) - v1(1)
            if not o.bRegular :
                o.dCell = []
        elif v3 is None:
            # 2D
            o.nDim = 2
            o.v1 = v1
            o.v2 = v2
            o.v3 = []
            o.dCell     = np.array([v1[1]-v1[0], v2[1]-v2[0]])
            o.xMesh_min = np.array([v1[0],v2[0]])
            o.xMesh_max = np.array([v1[-1],v2[-1]])
            o.n         = np.array([len(v1),len(v2)])
            o.bRegular = sum(np.abs(np.diff(v1,2)))< 1e-12 and sum(np.abs(np.diff(v2,2)))< 1e-12
            if not o.bRegular :
                o.dCell = []
        else:
            # 3D
            o.nDim = 3
            o.v1 = v1
            o.v2 = v2
            o.v3 = v3
            o.dCell     = np.array([ v1[1]-v1[0], v2[1]-v2[0], v3[1]-v3[0] ])
            o.xMesh_min = np.array([v1[0],v2[0],v3[0]])
            o.xMesh_max = np.array([v1[-1],v2[-1],v3[-1]])
            o.n         = np.array([len(v1),len(v2),len(v3)])
            o.bRegular =  sum(np.abs(np.diff(v1,2))) < 1e-12 and sum(np.abs(np.diff(v2,2)))< 1e-12 and sum(np.abs(np.diff(v3,2)))< 1e-12 
            if not o.bRegular :
                o.dCell = []
        
        o.dv1 = np.diff(o.v1)
        o.dv2 = np.diff(o.v2)
        o.dv3 = np.diff(o.v3)
        o.nGridPoints = np.prod(o.n)
        o.nCells = np.prod(o.n - 1)
        o.bGridPointsComputed      = False
        o.bCellCentersComputed     = False
        o.bCellVolumesComputed     = False
        o.bFakeCellVolumesComputed = False

    def __repr__(o):
        s=''
        s += 'nDim                    : {}\n'.format(o.nDim)
        s += 'dCell                   : {}\n'.format(o.dCell                   )
        s += 'xMesh_min               : {}\n'.format(o.xMesh_min               )
        s += 'xMesh_max               : {}\n'.format(o.xMesh_max               )
        s += 'n                       : {}\n'.format(o.n                       )
        s += 'bRegular                : {}\n'.format(o.bRegular                )
        s += 'nGridPoints             : {}\n'.format(o.nGridPoints             )
        s += 'nCells                  : {}\n'.format(o.nCells                  )
        s += 'bGridPointsComputed     : {}\n'.format(o.bGridPointsComputed     )
        s += 'bCellCentersComputed    : {}\n'.format(o.bCellCentersComputed    )
        s += 'bCellVolumesComputed    : {}\n'.format(o.bCellVolumesComputed    )
        s += 'bFakeCellVolumesComputed: {}\n'.format(o.bFakeCellVolumesComputed)
        s += 'v1                      : {}\n'.format(o.v1                      )
        s += 'v2                      : {}\n'.format(o.v2                      )
        s += 'v3                      : {}\n'.format(o.v3                      )
        s += 'dv1                     : {}\n'.format(o.dv1                     )
        s += 'dv2                     : {}\n'.format(o.dv2                     )
        s += 'dv3                     : {}\n'.format(o.dv3                     )
        return s


        
    def getCellVolume(o): 
        if o.bRegular:
            Volume = np.prod(o.dCell)
        else:
            raise Exception('Cell volume not defined')
        
        return Volume
        
        
    def getCellVolumes(o): 
        if (not o.bCellVolumesComputed ):
            o.CellVolumes = np.zeros((o.nCells,1))
            if o.bRegular:
                o.CellVolumes[:] = np.prod(o.dCell)
            else:
                if o.nDim == 1:
                    o.CellVolumes[:] = o.dv1
                elif o.nDim == 2:
                    p = 0
                    for i in np.arange(o.n[0]):
                        for j in np.arange(o.n[1]):
                            o.CellVolumes[p] = o.dv1[i] * o.dv2[j]
                            p = p + 1
                elif o.nDim == 3:
                    p = 0
                    for i in np.arange(o.n[0]):
                        for j in np.arange(o.n[1]):
                            for k in np.arange(o.n[2]):
                                o.CellVolumes[p] = o.dv1[i] * o.dv2[j] * o.dv3[k]
                                p = p + 1
            o.bCellVolumesComputed = True
        
        Volumes = o.CellVolumes
        return Volumes
        
        
    def getFakeCellVolumes(o): 
        if (not o.bFakeCellVolumesComputed ):
            o.FakeCellVolumes = np.zeros((o.nGridPoints,1))
            if o.bRegular:
                o.FakeCellVolumes = np.prod(o.dCell)
            else:
                if o.nDim == 1:
                    o.FakeCellVolumes = (np.array([0,o.dv1]) + np.array([o.dv1,0])) / 2
                elif o.nDim == 2:
                    X1 = (np.array([0,o.dv1]) + np.array([o.dv1,0])) / 2
                    X2 = (np.array([0,o.dv2]) + np.array([o.dv2,0])) / 2
                    p = 0
                    for i in np.arange(o.n[0]):
                        for j in np.arange(o.n[1]):
                            o.FakeCellVolumes[p] = X1[i] * X2[j]
                            p = p + 1
                    Zc = []
                elif o.nDim == 2:
                    X1 = (np.array([0,o.dv1]) + np.array([o.dv1,0])) / 2
                    X2 = (np.array([0,o.dv2]) + np.array([o.dv2,0])) / 2
                    X3 = (np.array([0,o.dv3]) + np.array([o.dv3,0])) / 3
                    p = 0
                    for i in np.arange(o.n[0]):
                        for j in np.arange(o.n[1]):
                            for k in np.arange(o.n[2]):
                                o.FakeCellVolumes[p] = X1[i] * X2[j] * X3[k]
                                p = p + 1
            o.bFakeCellVolumesComputed = True
        
        Volumes = o.FakeCellVolumes
        return Volumes
        
        
    def getFlatCellCenters(o): 
        if not o.bCellCentersComputed :
            CC = np.zeros((o.nCells,o.nDim))
            if o.nDim == 1:
                X1 = o.v1[:-1] + o.dv1 / 2
                CC[:,0] = X1
            elif o.nDim == 2:
                X1 = o.v1[:-1] + o.dv1 / 2
                X2 = o.v2[:-1] + o.dv2 / 2
                p = 0
                for i in np.arange(o.n[0]-1):
                    for j in np.arange(o.n[1]-1):
                        CC[p,0] = X1[i]
                        CC[p,1] = X2[j]
                        p = p + 1
                Zc = []
            elif o.nDim == 3:
                X1 = o.v1[:-1] + o.dv1 / 2
                X2 = o.v2[:-1] + o.dv2 / 2
                X3 = o.v3[:-1] + o.dv3 / 2
                p = 0
                for i in np.arange(o.n[0]-1):
                    for j in np.arange(o.n[1]-1):
                        for k in np.arange(o.n[2]-1):
                            CC[p,0] = X1[i]
                            CC[p,1] = X2[j]
                            CC[p,2] = X3[k]
                            p = p + 1
            o.bCellCentersComputed = True
            o.CC=CC
        return o.CC
        
        
    def getFlatGridPoints(o = None): 
        if not o.bGridPointsComputed :
            GP = np.zeros((o.nGridPoints,o.nDim))
            if o.nDim == 1:
                GP[:,0]= o.v1
            elif o.nDim == 2:
                p = 0
                for i in np.arange(o.n[0]):
                    for j in np.arange(o.n[1]):
                        GP[p,0] = o.v1[i]
                        GP[p,1] = o.v2[j]
                        p = p + 1
                Zc = []
            elif o.nDim == 3:
                p = 0
                for i in np.arange(o.n[0]):
                    for j in np.arange(o.n[1]):
                        for k in np.arange(o.n[2]):
                            p = p + 1
                            GP[p,0] = o.v1[i]
                            GP[p,1] = o.v2[j]
                            GP[p,2] = o.v3[k]
            o.bGridPointsComputed = True
            o.GP=GP
        return o.GP
        
        
    def flattenValues(o,m = None): 
        # Takes Mesh values m(1:nd,1:n1,1:n2,1:n3)
        nGridPoints = np.prod(o.n)
        s = m.shape
        if (len(s) == o.nDim):
            m = np.reshape(m, tuple(np.array([1,np.transpose(s)])), order="F")
            nd = 1
        else:
            nd = m.shape[0]
        
        nInput = np.prod(m.shape) / nd
        if (nInput != nGridPoints):
            raise Exception('Wrong mesh input')
        
        v_flat = np.zeros((nd,nGridPoints))
        if o.nDim == 1:
            v_flat = m[:,:]
        elif o.nDim == 2:
            p = 0
            for i in np.arange(o.n[0]):
                for j in np.arange(o.n[1]):
                    p = p + 1
                    v_flat[:,p] = m[:,i,j]
        elif o.nDim ==3:
            p = 0
            for i in np.arange(o.n[0]):
                for j in np.arange(o.n[1]):
                    for k in np.arange(o.n[2]):
                        p = p + 1
                        v_flat[:,p] = m[:,i,j]
        
        return v_flat
        
