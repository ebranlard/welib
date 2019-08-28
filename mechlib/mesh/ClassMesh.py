import numpy as np
##

class ClassMesh():
    

        # constructor
        
    def __init__(o = None,v1 = None,v2 = None,v3 = None): 
        o.bGridPointsComputed = False
        o.bCellCentersComputed = False
        o.bCellVolumesComputed = False
        o.bFakeCellVolumesComputed = False
        if len(varargin) == 1:
            o.nDim = 1
            o.v1 = v1
            o.v2 = []
            o.v3 = []
            o.dCell[1-1] = v1(2) - v1(1)
            o.xMesh_min = np.array([v1(1)])
            o.xMesh_max = np.array([v1(end())])
            o.n = np.array([len(v1)])
            o.bRegular = mean(diff(v1)) == v1(2) - v1(1)
            if not o.bRegular :
                o.dCell = []
        else:
            if len(varargin) == 2:
                o.nDim = 2
                o.v1 = v1
                o.v2 = v2
                o.v3 = []
                o.dCell[1-1] = v1(2) - v1(1)
                o.dCell[2-1] = v2(2) - v2(1)
                o.xMesh_min = np.array([v1(1),v2(1)])
                o.xMesh_max = np.array([v1(end()),v2(end())])
                o.n = np.array([len(v1),len(v2)])
                o.bRegular = sum(np.abs(diff(v1,2))) < np.logical_and(1e-12,sum(np.abs(diff(v2,2)))) < 1e-12
                if not o.bRegular :
                    o.dCell = []
            else:
                if len(varargin) == 3:
                    o.nDim = 3
                    o.v1 = v1
                    o.v2 = v2
                    o.v3 = v3
                    o.dCell[1-1] = v1(2) - v1(1)
                    o.dCell[2-1] = v2(2) - v2(1)
                    o.dCell[3-1] = v3(2) - v3(1)
                    o.xMesh_min = np.array([v1(1),v2(1),v3(1)])
                    o.xMesh_max = np.array([v1(end()),v2(end()),v3(end())])
                    o.n = np.array([len(v1),len(v2),len(v3)])
                    o.bRegular = sum(np.abs(diff(v1,2))) < np.logical_and(1e-12,sum(np.abs(diff(v2,2)))) < np.logical_and(1e-12,sum(np.abs(diff(v1,2)))) < 1e-12
                    if not o.bRegular :
                        o.dCell = []
        
        o.dv1 = diff(o.v1)
        o.dv2 = diff(o.v2)
        o.dv3 = diff(o.v3)
        o.nGridPoints = prod(o.n)
        o.nCells = prod(o.n - 1)
        o.bGridPointsComputed = False
        o.bCellCentersComputed = False
        o.bCellVolumesComputed = False
        o.bFakeCellVolumesComputed = False
        return
        
        
    def getCellVolume(o = None): 
        if o.bRegular:
            Volume = prod(o.dCell)
        else:
            raise Exception('Cell volume not defined')
        
        return Volume
        
        
    def getCellVolumes(o = None): 
        if (not o.bCellVolumesComputed ):
            o.CellVolumes = np.zeros((o.nCells,1))
            if o.bRegular:
                o.CellVolumes = prod(o.dCell)
            else:
                if (1) == (o.nDim):
                    o.CellVolumes = o.dv1
                else:
                    if (2) == (o.nDim):
                        p = 0
                        for i in np.arange(1,o.n(1) - 1+1).reshape(-1):
                            for j in np.arange(1,o.n(2) - 1+1).reshape(-1):
                                p = p + 1
                                o.CellVolumes[p-1] = o.dv1(i) * o.dv2(j)
                    else:
                        if (3) == (o.nDim):
                            p = 0
                            for i in np.arange(1,o.n(1) - 1+1).reshape(-1):
                                for j in np.arange(1,o.n(2) - 1+1).reshape(-1):
                                    for k in np.arange(1,o.n(3) - 1+1).reshape(-1):
                                        p = p + 1
                                        o.CellVolumes[p-1] = o.dv1(i) * o.dv2(j) * o.dv3(k)
            o.bCellVolumesComputed = True
        
        Volumes = o.CellVolumes
        return Volumes
        
        
    def getFakeCellVolumes(o = None): 
        if (not o.bFakeCellVolumesComputed ):
            o.FakeCellVolumes = np.zeros((o.nGridPoints,1))
            if o.bRegular:
                o.FakeCellVolumes = prod(o.dCell)
            else:
                if (1) == (o.nDim):
                    o.FakeCellVolumes = (np.array([0,o.dv1]) + np.array([o.dv1,0])) / 2
                else:
                    if (2) == (o.nDim):
                        X1 = (np.array([0,o.dv1]) + np.array([o.dv1,0])) / 2
                        X2 = (np.array([0,o.dv2]) + np.array([o.dv2,0])) / 2
                        p = 0
                        for i in np.arange(1,o.n(1)+1).reshape(-1):
                            for j in np.arange(1,o.n(2)+1).reshape(-1):
                                p = p + 1
                                o.FakeCellVolumes[p-1] = X1(i) * X2(j)
                        Zc = []
                    else:
                        if (3) == (o.nDim):
                            X1 = (np.array([0,o.dv1]) + np.array([o.dv1,0])) / 2
                            X2 = (np.array([0,o.dv2]) + np.array([o.dv2,0])) / 2
                            X3 = (np.array([0,o.dv3]) + np.array([o.dv3,0])) / 3
                            p = 0
                            for i in np.arange(1,o.n(1)+1).reshape(-1):
                                for j in np.arange(1,o.n(2)+1).reshape(-1):
                                    for k in np.arange(1,o.n(3)+1).reshape(-1):
                                        p = p + 1
                                        o.FakeCellVolumes[p-1] = X1(i) * X2(j) * X3(k)
            o.bFakeCellVolumesComputed = True
        
        Volumes = o.FakeCellVolumes
        return Volumes
        
        
    def getFlatCellCenters(o = None): 
        if not o.bCellCentersComputed :
            Xc = np.zeros((o.nCells,1))
            Yc = np.zeros((o.nCells,1))
            Zc = np.zeros((o.nCells,1))
            if (1) == (o.nDim):
                X1 = o.v1(np.arange(1,end() - 1+1)) + o.dv1 / 2
                Xc = X1
                Yc = []
                Zc = []
            else:
                if (2) == (o.nDim):
                    X1 = o.v1(np.arange(1,end() - 1+1)) + o.dv1 / 2
                    X2 = o.v2(np.arange(1,end() - 1+1)) + o.dv2 / 2
                    p = 0
                    for i in np.arange(1,o.n(1) - 1+1).reshape(-1):
                        for j in np.arange(1,o.n(2) - 1+1).reshape(-1):
                            p = p + 1
                            Xc[p-1] = X1(i)
                            Yc[p-1] = X2(j)
                    Zc = []
                else:
                    if (3) == (o.nDim):
                        X1 = o.v1(np.arange(1,end() - 1+1)) + o.dv1 / 2
                        X2 = o.v2(np.arange(1,end() - 1+1)) + o.dv2 / 2
                        X3 = o.v3(np.arange(1,end() - 1+1)) + o.dv3 / 2
                        p = 0
                        for i in np.arange(1,o.n(1) - 1+1).reshape(-1):
                            for j in np.arange(1,o.n(2) - 1+1).reshape(-1):
                                for k in np.arange(1,o.n(3) - 1+1).reshape(-1):
                                    p = p + 1
                                    Xc[p-1] = X1(i)
                                    Yc[p-1] = X2(j)
                                    Zc[p-1] = X3(k)
            o.CC = np.array([Xc,Yc,Zc])
            o.bCellCentersComputed = True
        
        CC = o.CC
        return CC
        
        
    def getFlatGridPoints(o = None): 
        if not o.bGridPointsComputed :
            Xc = np.zeros((o.nGridPoints,1))
            Yc = np.zeros((o.nGridPoints,1))
            Zc = np.zeros((o.nGridPoints,1))
            if (1) == (o.nDim):
                Xc = o.v1
                Yc = []
                Zc = []
            else:
                if (2) == (o.nDim):
                    p = 0
                    for i in np.arange(1,o.n(1)+1).reshape(-1):
                        for j in np.arange(1,o.n(2)+1).reshape(-1):
                            p = p + 1
                            Xc[p-1] = o.v1(i)
                            Yc[p-1] = o.v2(j)
                    Zc = []
                else:
                    if (3) == (o.nDim):
                        p = 0
                        for i in np.arange(1,o.n(1)+1).reshape(-1):
                            for j in np.arange(1,o.n(2)+1).reshape(-1):
                                for k in np.arange(1,o.n(3)+1).reshape(-1):
                                    p = p + 1
                                    Xc[p-1] = o.v1(i)
                                    Yc[p-1] = o.v2(j)
                                    Zc[p-1] = o.v3(k)
            o.GP = np.array([Xc,Yc,Zc])
            o.bGridPointsComputed = True
        
        GP = o.GP
        return GP
        
        
    def flattenValues(o = None,m = None): 
        # Takes Mesh values m(1:nd,1:n1,1:n2,1:n3)
        nGridPoints = prod(o.n)
        s = m.shape
        if (len(s) == o.nDim):
            m = np.reshape(m, tuple(np.array([1,np.transpose(s)])), order="F")
            nd = 1
        else:
            nd = m.shape[1-1]
        
        nInput = prod(m.shape) / nd
        if (nInput != nGridPoints):
            warn('Wrong mesh input')
            kbd
        
        v_flat = np.zeros((nd,nGridPoints))
        if (1) == (o.nDim):
            v_flat = m(:,:)
        else:
            if (2) == (o.nDim):
                p = 0
                for i in np.arange(1,o.n(1)+1).reshape(-1):
                    for j in np.arange(1,o.n(2)+1).reshape(-1):
                        p = p + 1
                        v_flat[:,p-1] = m(:,i,j)
            else:
                if (3) == (o.nDim):
                    p = 0
                    for i in np.arange(1,o.n(1)+1).reshape(-1):
                        for j in np.arange(1,o.n(2)+1).reshape(-1):
                            for k in np.arange(1,o.n(3)+1).reshape(-1):
                                p = p + 1
                                v_flat[:,p-1] = m(:,i,j)
        
        return v_flat
        