import numpy as np

class VTK_Misc:
    def __init__(self):
        self.file      = None
        self.nPoints   = 0
        self.nData     = 0
        self.bFileOpen = False

def vtk_new_ascii_file(filename, label, mvtk):
    mvtk.file      = open(filename, 'w')
    mvtk.bFileOpen = True
    mvtk.nData     = -1
    mvtk.nPoints   = 0
    mvtk.file.write('# vtk DataFile Version 2.0\n')
    mvtk.file.write(f'{label}\n')
    mvtk.file.write('ASCII\n\n')
    return True

def vtk_close_file(mvtk):
    if mvtk.bFileOpen and mvtk.file is not None:
        mvtk.file.close()
        mvtk.bFileOpen = False

def vtk_dataset_polydata(points, mvtk, bladeFrame=False):
    pts = points
    if pts.shape[0] == 3 and pts.shape[1] != 3:
        pts = pts.T
    nPoints         = pts.shape[0]
    mvtk.nPoints    = nPoints
    mvtk.file.write('DATASET POLYDATA\n')
    mvtk.file.write(f'POINTS {nPoints} double\n')
    for i in range(nPoints):
        mvtk.file.write(f'{pts[i,0]:.8e} {pts[i,1]:.8e} {pts[i,2]:.8e}\n')
    mvtk.file.write('\n')

def vtk_lines(lines, mvtk):
    L = lines
    if L.shape[0] == 2 and L.shape[1] != 2:
        L = L.T
    nLines       = L.shape[0]
    mvtk.nData   = nLines
    mvtk.file.write(f'LINES {nLines} {3*nLines}\n')
    for i in range(nLines):
        mvtk.file.write(f'2 {L[i,0]} {L[i,1]}\n')
    mvtk.file.write('\n')

def vtk_quad(Q, mvtk):
    Qarr = Q
    if Qarr.shape[0] == 4 and Qarr.shape[1] != 4:
        Qarr = Qarr.T
    nQuads      = Qarr.shape[0]
    mvtk.nData  = nQuads
    mvtk.file.write(f'POLYGONS {nQuads} {5*nQuads}\n')
    for i in range(nQuads):
        mvtk.file.write(f'4 {Qarr[i,0]} {Qarr[i,1]} {Qarr[i,2]} {Qarr[i,3]}\n')
    mvtk.file.write('\n')

def vtk_cell_data_init(mvtk):
    mvtk.file.write(f'CELL_DATA {mvtk.nData}\n')

def vtk_cell_data_scalar(D, sname, mvtk):
    mvtk.file.write(f'SCALARS {sname} double\n')
    mvtk.file.write('LOOKUP_TABLE default\n')
    Dflat = D.ravel()
    for val in Dflat:
        mvtk.file.write(f'{val:.8e}\n')

def WrVTK_Segments(filename, mvtk, SegPoints, SegConnct, SegGamma, SegEpsilon=None, bladeFrame=False):
    """
    Write VTK file for vortex segments.
    SegPoints: (n,3) or (3,n) array of segment points
    SegConnct: (n,2) or (2,n) array of segment connectivity (VTK 0-based)
    SegGamma:  (n,) array of segment strengths
    SegEpsilon: (n,) array (optional, not used here)
    bladeFrame: bool, not used here
    """
    if vtk_new_ascii_file(filename, 'Sgmt', mvtk):
        vtk_dataset_polydata(SegPoints, mvtk, bladeFrame)
        vtk_lines(SegConnct, mvtk)
        vtk_cell_data_init(mvtk)
        vtk_cell_data_scalar(SegGamma, 'Gamma', mvtk)
        vtk_close_file(mvtk)

def LatticeToPanlConnectivity(LatticePoints):
    """
    Converts a 3D lattice of points (nSpan, nDepth, 3)
    into a flat array of points and a connectivity array for quads.
    Returns:
        Points: (nPoints, 3)
        Connectivity: (nQuads, 4)  # VTK quad connectivity (0-based)
    """
    # Accept both (3, nSpan, nDepth) and (nSpan, nDepth, 3)
    #if LatticePoints.shape[0] == 3:
    #    pts = np.moveaxis(LatticePoints, 0, -1)  # (nSpan, nDepth, 3)
    #else:
    #    pts = LatticePoints
    pts = LatticePoints
    nSpan, nDepth, _ = pts.shape
    print('nSpan', nSpan, 'nDepth', nDepth)
    Points = pts.reshape(-1, 3)
    Connectivity = []
    for i in range(nSpan-1):
        for j in range(nDepth-1):
            n0 = i*nDepth + j
            n1 = (i+1)*nDepth + j
            n2 = (i+1)*nDepth + (j+1)
            n3 = i*nDepth + (j+1)
            Connectivity.append([n0, n1, n2, n3])
    Connectivity = np.array(Connectivity, dtype=int)
    return Connectivity, Points

def WrVTK_Lattice(filename, mvtk, LatticePoints, LatticeGamma, LatticeData3d=None, bladeFrame=False):
    """
    Write VTK file for a lattice of panels (quads).
    LatticePoints: (3, nSpan, nDepth) or (nSpan, nDepth, 3)
    LatticeGamma:  (nSpan-1, nDepth-1) or (nSpan, nDepth)
    LatticeData3d: (n, nSpan-1, nDepth-1) or (n, nSpan, nDepth) (optional, not used here)
    bladeFrame: bool, not used here
    """
    Connectivity, Points = LatticeToPanlConnectivity(LatticePoints)
    if vtk_new_ascii_file(filename, '', mvtk):
        vtk_dataset_polydata(Points, mvtk, bladeFrame)
        vtk_quad(Connectivity, mvtk)
        vtk_cell_data_init(mvtk)
        vtk_cell_data_scalar(LatticeGamma, 'Gamma', mvtk)