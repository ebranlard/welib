import numpy as np
import unittest
    
def fCoordRegularGrid(x0,xBox,dBox):
    C = (x0 - xBox) / dBox
    ic = int(np.floor(C)) + 1
    dc = C + 1 - ic
    return ic-1,dc
    
def interp_coeff_mp4(i2,dx2,nx): 
    """
    NOTE:
      input and output index from 0 to nx-1 !!!
    """
    i2=i2+1 # TODO, waiting for script to be updated
    # Index of other cells
    i1 = i2 - 1
    i3 = i2 + 1
    i4 = i2 + 2
    # Normalized distance to other cells
    dx1 = dx2 + 1.0
    dx3 = 1.0 - dx2
    dx4 = 2.0 - dx2
    # M-prime-4 kernel
    ax1 = 0.5 * (2.0 - dx1) ** 2 * (1.0 - dx1)
    ax2 = 1.0 - 2.5 * dx2 ** 2 + 1.5 * dx2 ** 3
    ax3 = 1.0 - 2.5 * dx3 ** 2 + 1.5 * dx3 ** 3
    ax4 = 0.5 * (2.0 - dx4) ** 2 * (1.0 - dx4)
    if i2==nx-1:
        # My hack for one sided right
        #         ax3 = 2*(ax3)+ax4;
        #         ax2 = (ax2+ax4);
        #         i4=1; ax4=0;# arbitrary, i1 will not be used
        i1 ,i2 ,i3 ,i4  = 1 ,1 ,1 ,1
        ax1,ax2,ax3,ax4 = 0.,0.,0.,0.
    elif i2 == 1:
        # # My hack for one sided left
        #         ax2 = 2*(ax2)+ax1;
        #         ax3 = (ax3+ax1);
        #         i1=1; ax1=0; # arbitrary, i1 will not be used
        i1 ,i2 ,i3 ,i4  = 1 ,1 ,1 ,1
        ax1,ax2,ax3,ax4 = 0.,0.,0.,0.
    elif i2 < 1:
        # Should not happen
        i1 ,i2 ,i3 ,i4  = 1 ,1 ,1 ,1
        ax1,ax2,ax3,ax4 = 0.,0.,0.,0.
    elif (i2 > nx - 1):
        # Should not happen
        i1 ,i2 ,i3 ,i4  = 1 ,1 ,1 ,1
        ax1,ax2,ax3,ax4 = 0.,0.,0.,0.
    elif i1 <= 0 or i2 <= 0:
        # Might happen if on grid
        i1 ,i2 ,i3 ,i4  = 1 ,1 ,1 ,1
        ax1,ax2,ax3,ax4 = 0.,0.,0.,0.
    elif i4 > nx or i3 > nx:
        # Might happen if on grid
        i1 ,i2 ,i3 ,i4  = 1 ,1 ,1 ,1
        ax1,ax2,ax3,ax4 = 0.,0.,0.,0.
    return ax1,ax2,ax3,ax4,i1-1,i2-1,i3-1,i4-1 # NOTE: returns -1

def interp_coeff_lambda3(i2,dx2,nx):
    """
    NOTE:
      input and output index from 0 to nx-1 !!!
    """
    i2=i2+1 # TODO, waiting for script to be updated
    # Find index of other cells
    i1 = i2 - 1
    i3 = i2 + 1
    i4 = i2 + 2
    # Find normalised distance to other cells
    dx1 = dx2 + 1.0
    dx3 = 1.0 - dx2
    dx4 = 2.0 - dx2
    # lambda 3 kernel
    ax1 = 1.0 / 6.0 * (1.0 - dx1) * (2.0 - dx1) * (3.0 - dx1)
    ax2 = 1.0 / 2.0 * (1 - dx2 ** 2) * (2 - dx2)
    ax3 = 1.0 / 2.0 * (1 - dx3 ** 2) * (2 - dx3)
    ax4 = 1.0 / 6.0 * (1.0 - dx4) * (2.0 - dx4) * (3.0 - dx4)
    if i2==nx-1:
        i1 ,i2 ,i3 ,i4  = 1 ,1 ,1 ,1
        ax1,ax2,ax3,ax4 = 0.,0.,0.,0.
    elif i2 == 1:
        i1 ,i2 ,i3 ,i4  = 1 ,1 ,1 ,1
        ax1,ax2,ax3,ax4 = 0.,0.,0.,0.
    elif i2 < 1:
        # Should not happen
        i1 ,i2 ,i3 ,i4  = 1 ,1 ,1 ,1
        ax1,ax2,ax3,ax4 = 0.,0.,0.,0.
    elif (i2 > nx - 1):
        # Should not happen
        i1 ,i2 ,i3 ,i4  = 1 ,1 ,1 ,1
        ax1,ax2,ax3,ax4 = 0.,0.,0.,0.
    elif i1 <= 0 or i2 <= 0:
        # Might happen if on grid
        i1 ,i2 ,i3 ,i4  = 1 ,1 ,1 ,1
        ax1,ax2,ax3,ax4 = 0.,0.,0.,0.
    elif i4 > nx or i3 > nx:
        # Might happen if on grid
        i1 ,i2 ,i3 ,i4  = 1 ,1 ,1 ,1
        ax1,ax2,ax3,ax4 = 0.,0.,0.,0.
    return ax1,ax2,ax3,ax4,i1-1,i2-1,i3-1,i4-1


def interp_p2m_mp4_1d(Part,nPart,part_values,nval,v1,n1,bRegular=None): 
    """
    xBox: Origin of the grid
    dBox: Vector containing the length of a grid cell in the three direction
    part       (npart, 1:ndim)
    part_values(npart, 1:nval)

    M4' kernel needs four points
             xp
     i1   i2   i3   i4
      |    | .  |    |
    """
    Part = np.transpose(Part) # TODO
    part_values = np.transpose(part_values)
    
    if (part_values.shape[0] != nval):
        pass
    
    if (Part.shape[1] != nPart):
        raise Exception('Part has wrong size')
    
    if (n1 < 4):
        print('No guarantee with tiny grid ')
    
    mesh = np.zeros((nval,n1,n2))
    ic   = np.array([0,0])
    dc   = np.zeros(2)
    xBox = np.zeros(2)
    dBox = np.zeros(2)
    xBox[0] = v1[0]
    xBox[1] = v2[0]
    dBox[0] = v1[1] - v1[0]
    dBox[1] = v2[1] - v2[0]

    if bRegular == 1:
        fIndexX= lambda x: fCoordRegularGrid(x,xBox[0],dBox[0])
    else:
        fIndexX = lambda x: fCoordRectilinearGrid(x,v1)

    for i in np.arange(nPart):
        # Coordinates and distance in grid index space (to nearest left grid point)
        ic[0],dc[0] = fIndexX(Part[0,i])
        # Getting the M'4 kernel coefficients
        ax1,ax2,ax3,ax4,i1,i2,i3,i4 = interp_coeff_mp4(ic[0],dc[0],n1)
        # Project to field
        mesh[:nval,i1] +=  ax1 * part_values[:nval,i]
        mesh[:nval,i2] +=  ax2 * part_values[:nval,i]
        mesh[:nval,i3] +=  ax3 * part_values[:nval,i]
        mesh[:nval,i4] +=  ax4 * part_values[:nval,i]
    
    return mesh

def interp_p2m_lambda3_2d(Part,nPart,part_values,nval,v1,v2,n1,n2,bRegular=None): 
    """
    xBox: Origin of the grid
    dBox: Vector containing the length of a grid cell in the three direction
    part       (npart, 1:ndim)
    part_values(npart, 1:nval)

    M4' kernel needs four points
             xp
     i1   i2   i3   i4
      |    | .  |    |
    """
    Part = np.transpose(Part) # TODO
    part_values = np.transpose(part_values)
    
    if (part_values.shape[0] != nval):
        pass
    
    if (Part.shape[1] != nPart):
        raise Exception('Part has wrong size')
    
    if (n1 < 4):
        print('No guarantee with tiny grid ')
    
    if (n2 < 4):
        print('No guarantee with tiny grid ')
    
    mesh = np.zeros((nval,n1,n2))
    ic   = np.array([0,0])
    dc   = np.zeros(2)
    xBox = np.zeros(2)
    dBox = np.zeros(2)
    xBox[0] = v1[0]
    xBox[1] = v2[0]
    dBox[0] = v1[1] - v1[0]
    dBox[1] = v2[1] - v2[0]

    if bRegular == 1:
        fIndexX= lambda x: fCoordRegularGrid(x,xBox[0],dBox[0])
        fIndexY= lambda y: fCoordRegularGrid(y,xBox[1],dBox[1])
    else:
        fIndexX = lambda x: fCoordRectilinearGrid(x,v1)
        fIndexY = lambda x: fCoordRectilinearGrid(y,v2)

    for i in np.arange(nPart):
        # Coordinates and distance in grid index space (to nearest left grid point)
        ic[0],dc[0] = fIndexX(Part[0,i])
        ic[1],dc[1] = fIndexY(Part[1,i])
        # Getting the M'4 kernel coefficients
        ax1,ax2,ax3,ax4,i1,i2,i3,i4 = interp_coeff_lambda3(ic[0],dc[0],n1)
        ay1,ay2,ay3,ay4,j1,j2,j3,j4 = interp_coeff_lambda3(ic[1],dc[1],n2)
        # Combine 1D kernels into the 2D kernel and apply to the field
        mesh[:nval,i1,j1] +=  ax1 * ay1 * part_values[:nval,i]
        mesh[:nval,i2,j1] +=  ax2 * ay1 * part_values[:nval,i]
        mesh[:nval,i3,j1] +=  ax3 * ay1 * part_values[:nval,i]
        mesh[:nval,i4,j1] +=  ax4 * ay1 * part_values[:nval,i]
        mesh[:nval,i1,j2] +=  ax1 * ay2 * part_values[:nval,i]
        mesh[:nval,i2,j2] +=  ax2 * ay2 * part_values[:nval,i]
        mesh[:nval,i3,j2] +=  ax3 * ay2 * part_values[:nval,i]
        mesh[:nval,i4,j2] +=  ax4 * ay2 * part_values[:nval,i]
        mesh[:nval,i1,j3] +=  ax1 * ay3 * part_values[:nval,i]
        mesh[:nval,i2,j3] +=  ax2 * ay3 * part_values[:nval,i]
        mesh[:nval,i3,j3] +=  ax3 * ay3 * part_values[:nval,i]
        mesh[:nval,i4,j3] +=  ax4 * ay3 * part_values[:nval,i]
        mesh[:nval,i1,j4] +=  ax1 * ay4 * part_values[:nval,i]
        mesh[:nval,i2,j4] +=  ax2 * ay4 * part_values[:nval,i]
        mesh[:nval,i3,j4] +=  ax3 * ay4 * part_values[:nval,i]
        mesh[:nval,i4,j4] +=  ax4 * ay4 * part_values[:nval,i]
    
    return mesh


def interp_p2m_mp4_2d(Part,nPart,part_values,nval,v1,v2,n1,n2,bRegular=None): 
    """
    xBox: Origin of the grid
    dBox: Vector containing the length of a grid cell in the three direction
    part       (npart, 1:ndim)
    part_values(npart, 1:nval)

    M4' kernel needs four points
             xp
     i1   i2   i3   i4
      |    | .  |    |
    """
    Part = np.transpose(Part) # TODO
    part_values = np.transpose(part_values)
    
    if (part_values.shape[0] != nval):
        pass
    
    if (Part.shape[1] != nPart):
        raise Exception('Part has wrong size')
    
    if (n1 < 4):
        print('No guarantee with tiny grid ')
    
    if (n2 < 4):
        print('No guarantee with tiny grid ')
    
    mesh = np.zeros((nval,n1,n2))
    ic   = np.array([0,0])
    dc   = np.zeros(2)
    xBox = np.zeros(2)
    dBox = np.zeros(2)
    xBox[0] = v1[0]
    xBox[1] = v2[0]
    dBox[0] = v1[1] - v1[0]
    dBox[1] = v2[1] - v2[0]

    if bRegular == 1:
        fIndexX= lambda x: fCoordRegularGrid(x,xBox[0],dBox[0])
        fIndexY= lambda y: fCoordRegularGrid(y,xBox[1],dBox[1])
    else:
        fIndexX = lambda x: fCoordRectilinearGrid(x,v1)
        fIndexY = lambda x: fCoordRectilinearGrid(y,v2)

    for i in np.arange(nPart):
        # Coordinates and distance in grid index space (to nearest left grid point)
        ic[0],dc[0] = fIndexX(Part[0,i])
        ic[1],dc[1] = fIndexY(Part[1,i])
        # Getting the M'4 kernel coefficients
        ax1,ax2,ax3,ax4,i1,i2,i3,i4 = interp_coeff_mp4(ic[0],dc[0],n1)
        ay1,ay2,ay3,ay4,j1,j2,j3,j4 = interp_coeff_mp4(ic[1],dc[1],n2)
        # Combine 1D kernels into the 2D kernel and apply to the field
        mesh[:nval,i1,j1] +=  ax1 * ay1 * part_values[:nval,i]
        mesh[:nval,i2,j1] +=  ax2 * ay1 * part_values[:nval,i]
        mesh[:nval,i3,j1] +=  ax3 * ay1 * part_values[:nval,i]
        mesh[:nval,i4,j1] +=  ax4 * ay1 * part_values[:nval,i]
        mesh[:nval,i1,j2] +=  ax1 * ay2 * part_values[:nval,i]
        mesh[:nval,i2,j2] +=  ax2 * ay2 * part_values[:nval,i]
        mesh[:nval,i3,j2] +=  ax3 * ay2 * part_values[:nval,i]
        mesh[:nval,i4,j2] +=  ax4 * ay2 * part_values[:nval,i]
        mesh[:nval,i1,j3] +=  ax1 * ay3 * part_values[:nval,i]
        mesh[:nval,i2,j3] +=  ax2 * ay3 * part_values[:nval,i]
        mesh[:nval,i3,j3] +=  ax3 * ay3 * part_values[:nval,i]
        mesh[:nval,i4,j3] +=  ax4 * ay3 * part_values[:nval,i]
        mesh[:nval,i1,j4] +=  ax1 * ay4 * part_values[:nval,i]
        mesh[:nval,i2,j4] +=  ax2 * ay4 * part_values[:nval,i]
        mesh[:nval,i3,j4] +=  ax3 * ay4 * part_values[:nval,i]
        mesh[:nval,i4,j4] +=  ax4 * ay4 * part_values[:nval,i]
    
    return mesh


def interp_m2p_mp4_2d(Part,nPart,mesh,nval,v1,v2,n1,n2,bRegular = None): 
    """
    xBox: Origin of the grid
    dBox: Vector containing the length of a grid cell in the three direction
    part       (npart, 1:ndim)
    part_values(npart, 1:nval)
    mesh(1:nval,n1,n2)
    M4' kernel needs four points
            xp
    i1   i2   i3   i4
     |    | .  |    |
    """
    Part = np.transpose(Part)

    if (mesh.shape[0] != nval):
        raise Exception('Mesh has wrong size')
    
    if (Part.shape[1] != nPart):
        print('Part has wrong size')
    
    part_values = np.zeros((nval,nPart))

    ic   = np.array([0,0])
    dc   = np.zeros(2)
    xBox = np.zeros(2)
    dBox = np.zeros(2)
    xBox[0] = v1[0]
    xBox[1] = v2[0]
    dBox[0] = v1[1] - v1[0]
    dBox[1] = v2[1] - v2[0]

    if bRegular == 1:
        fIndexX= lambda x: fCoordRegularGrid(x,xBox[0],dBox[0])
        fIndexY= lambda y: fCoordRegularGrid(y,xBox[1],dBox[1])
    else:
        fIndexX = lambda x: fCoordRectilinearGrid(x,v1)
        fIndexY = lambda x: fCoordRectilinearGrid(y,v2)

    for i in np.arange(nPart):
        # Coordinates and distance in grid index space (to nearest left grid point)
        ic[0],dc[0] = fIndexX(Part[0,i])
        ic[1],dc[1] = fIndexY(Part[1,i])
        # Getting the M'4 kernel coefficients
        ax1,ax2,ax3,ax4,i1,i2,i3,i4 = interp_coeff_mp4(ic[0],dc[0],n1)
        ay1,ay2,ay3,ay4,j1,j2,j3,j4 = interp_coeff_mp4(ic[1],dc[1],n2)
        # Combine 1D kernels into the 2D kernel and apply to the part
        part_values[:nval,i] = ax1 * ay1 * mesh[:nval,i1,j1] +\
                               ax2 * ay1 * mesh[:nval,i2,j1] +\
                               ax3 * ay1 * mesh[:nval,i3,j1] +\
                               ax4 * ay1 * mesh[:nval,i4,j1] +\
                               ax1 * ay2 * mesh[:nval,i1,j2] +\
                               ax2 * ay2 * mesh[:nval,i2,j2] +\
                               ax3 * ay2 * mesh[:nval,i3,j2] +\
                               ax4 * ay2 * mesh[:nval,i4,j2] +\
                               ax1 * ay3 * mesh[:nval,i1,j3] +\
                               ax2 * ay3 * mesh[:nval,i2,j3] +\
                               ax3 * ay3 * mesh[:nval,i3,j3] +\
                               ax4 * ay3 * mesh[:nval,i4,j3] +\
                               ax1 * ay4 * mesh[:nval,i1,j4] +\
                               ax2 * ay4 * mesh[:nval,i2,j4] +\
                               ax3 * ay4 * mesh[:nval,i3,j4] +\
                               ax4 * ay4 * mesh[:nval,i4,j4]

    part_values = np.transpose(part_values)
    return part_values

def interp_m2p_lambda3_2d(Part,nPart,mesh,nval,v1,v2,n1,n2,bRegular = None): 
    """
    xBox: Origin of the grid
    dBox: Vector containing the length of a grid cell in the three direction
    part       (npart, 1:ndim)
    part_values(npart, 1:nval)
    mesh(1:nval,n1,n2)
    """
    Part = np.transpose(Part)

    if (mesh.shape[0] != nval):
        raise Exception('Mesh has wrong size')
    
    if (Part.shape[1] != nPart):
        print('Part has wrong size')
    
    part_values = np.zeros((nval,nPart))

    ic   = np.array([0,0])
    dc   = np.zeros(2)
    xBox = np.zeros(2)
    dBox = np.zeros(2)
    xBox[0] = v1[0]
    xBox[1] = v2[0]
    dBox[0] = v1[1] - v1[0]
    dBox[1] = v2[1] - v2[0]

    if bRegular == 1:
        fIndexX= lambda x: fCoordRegularGrid(x,xBox[0],dBox[0])
        fIndexY= lambda y: fCoordRegularGrid(y,xBox[1],dBox[1])
    else:
        fIndexX = lambda x: fCoordRectilinearGrid(x,v1)
        fIndexY = lambda x: fCoordRectilinearGrid(y,v2)

    for i in np.arange(nPart):
        # Coordinates and distance in grid index space (to nearest left grid point)
        ic[0],dc[0] = fIndexX(Part[0,i])
        ic[1],dc[1] = fIndexY(Part[1,i])
        # Getting the M'4 kernel coefficients
        ax1,ax2,ax3,ax4,i1,i2,i3,i4 = interp_coeff_lambda3(ic[0],dc[0],n1)
        ay1,ay2,ay3,ay4,j1,j2,j3,j4 = interp_coeff_lambda3(ic[1],dc[1],n2)
        # Combine 1D kernels into the 2D kernel and apply to the part
        part_values[:nval,i] = ax1 * ay1 * mesh[:nval,i1,j1] +\
                               ax2 * ay1 * mesh[:nval,i2,j1] +\
                               ax3 * ay1 * mesh[:nval,i3,j1] +\
                               ax4 * ay1 * mesh[:nval,i4,j1] +\
                               ax1 * ay2 * mesh[:nval,i1,j2] +\
                               ax2 * ay2 * mesh[:nval,i2,j2] +\
                               ax3 * ay2 * mesh[:nval,i3,j2] +\
                               ax4 * ay2 * mesh[:nval,i4,j2] +\
                               ax1 * ay3 * mesh[:nval,i1,j3] +\
                               ax2 * ay3 * mesh[:nval,i2,j3] +\
                               ax3 * ay3 * mesh[:nval,i3,j3] +\
                               ax4 * ay3 * mesh[:nval,i4,j3] +\
                               ax1 * ay4 * mesh[:nval,i1,j4] +\
                               ax2 * ay4 * mesh[:nval,i2,j4] +\
                               ax3 * ay4 * mesh[:nval,i3,j4] +\
                               ax4 * ay4 * mesh[:nval,i4,j4]

    part_values = np.transpose(part_values)
    return part_values


# --------------------------------------------------------------------------------}
# --- High level functions 
# --------------------------------------------------------------------------------{
def interp_p2m(Part,nPart,n,v_p,nDim,v1,v2=None, v3=None, kernel='mp4', bRegular=None): 
    nVal = v_p.shape[1]
    if kernel=='mp4':
        if nDim==2:
            MeshValues = interp_p2m_mp4_2d(Part,nPart,v_p,nVal,v1,v2,n[0],n[1],bRegular)
        else:
            raise NotImplementedError('3d')
    elif kernel=='lambda3':
        if nDim==2:
            raise NotImplementedError('p2m lambda3')
            MeshValues = interp_p2m_lambda3_2d(Part,nPart,v_p,nVal,v1,v2,n[0],n[1],bRegular)
        else:
            raise NotImplementedError('3d')
    else:
        raise Exception('Unknown Interpolation kernel')
    
    return MeshValues

def interp_m2p(PartP,nPart,n,MeshValues,nDim,v1,v2=None,v3=None,kernel='mp4',bRegular=True):
    nVal = MeshValues.shape[0]
    if kernel=='mp4':
        if nDim==2:
            v_p = interp_m2p_mp4_2d(PartP,nPart,MeshValues,nVal,v1,v2,n[0],n[1],bRegular)
        else:
            raise NotImplementedError('3d')
    elif kernel=='lambda3':
        if nDim==2:
            v_p = interp_m2p_lambda3_2d(PartP,nPart,MeshValues,nVal,v1,v2,n[0],n[1],bRegular)
        else:
            raise NotImplementedError('3d')
    else:
        raise Exception('Unknown Interpolation kernel')
    
    return v_p



# --------------------------------------------------------------------------------}
# --- High level function to be moved to particle class!!! TODO 
# --------------------------------------------------------------------------------{
def part2mesh(Part,mesh,kernel='mp4'): 
    """ Project particles into a grid 
    NOTE: comes from fParticleProjection
    """
    if mesh.nDim == 2:
        v_p = np.zeros((Part.nPart,2))
        v_p[:,0] = Part.Intensity
        v_p[:,1] = Part.Volume
    elif mesh.nDim == 2:
        v_p = np.zeros((Part.nPart,4))
        v_p[:,:3] = Part.Intensity
        v_p[:, 3] = Part.Volume
    
    mesh.values = interp_p2m(Part.P, Part.nPart, mesh.n, v_p, mesh.nDim, v1=mesh.v1, v2=mesh.v2, v3=mesh.v3, kernel=kernel, bRegular=mesh.bRegular)
    return mesh.values, v_p





# --------------------------------------------------------------------------------}
# --- TESTS
# --------------------------------------------------------------------------------{
class TestProjection(unittest.TestCase):

    def test_coord_regular(self):
        np.testing.assert_almost_equal(fCoordRegularGrid(1.1,0,2), (1-1,0.55))

    def test_mp4(self):
        np.testing.assert_almost_equal(interp_coeff_mp4(3-1,0.1,5),(-0.040500, 0.97650, 0.068500, -0.0045000, 2-1, 3-1, 4-1, 5-1))

    def test_lambda3(self):
        np.testing.assert_almost_equal(interp_coeff_lambda3(3-1,0.1,5), (-0.028500, 0.94050, 0.10450, -0.016500,2-1 ,3-1 ,4-1 ,5-1))

    def test_p2m_mp4_lambda3(self):
        # --- Init
        nPart            = 3
        nDim             = 2
        nVal             = 2
        Part             = np.zeros((nPart,nDim))
        Part[0,:]        = np.array([0,0])
        Part[1,:]        = np.array([0.5,0.5])
        Part[2,:]        = np.array([1.5,0.5])
        part_values      = np.zeros((nPart,nVal))
        part_values[0,:] = np.array([1,2])
        part_values[1,:] = np.array([1,2])
        part_values[2,:] = np.array([2,3])
        v1               = np.array([0 ,1,2,3 ,4])
        v2               = np.array([-2,0,2,4])
        n1               = len(v1)
        n2               = len(v2)
        bRegular         = True
        # --- Test mp4
        mesh0 = np.array([
        [   0.0087891,   -0.1083984,   -0.0283203,    0.0029297],
        [  -0.0791016,    0.9755859,    0.2548828,   -0.0263672],
        [  -0.0791016,    0.9755859,    0.2548828,   -0.0263672],
        [   0.0087891,   -0.1083984,   -0.0283203,    0.0029297],
        [   0.0000000,    0.0000000,    0.0000000,    0.0000000]])
        mesh1 = np.array([
        [   0.0131836,   -0.1625977,   -0.0424805,    0.0043945],
        [  -0.1186523,    1.4633789,    0.3823242,   -0.0395508],
        [  -0.1186523,    1.4633789,    0.3823242,   -0.0395508],
        [   0.0131836,   -0.1625977,   -0.0424805,    0.0043945],
        [   0.0000000,    0.0000000,    0.0000000,    0.0000000]])
        mesh= interp_p2m_mp4_2d(Part,nPart,part_values,nVal,v1,v2,n1,n2,bRegular)
        np.testing.assert_almost_equal(mesh[0,:,:], mesh0)
        np.testing.assert_almost_equal(mesh[1,:,:], mesh1)
        # --- Test lambda3
        mesh0 = np.array([
        [   0.0068359,   -0.1025391,   -0.0341797,    0.0048828],
        [  -0.0615234,    0.9228516,    0.3076172,   -0.0439453],
        [  -0.0615234,    0.9228516,    0.3076172,   -0.0439453],
        [   0.0068359,   -0.1025391,   -0.0341797,    0.0048828],
        [   0.0000000,    0.0000000,    0.0000000,    0.0000000]])
        mesh1 = np.array([
        [   0.0102539,   -0.1538086,   -0.0512695,    0.0073242],
        [  -0.0922852,    1.3842773,    0.4614258,   -0.0659180],
        [  -0.0922852,    1.3842773,    0.4614258,   -0.0659180],
        [   0.0102539,   -0.1538086,   -0.0512695,    0.0073242],
        [   0.0000000,    0.0000000,    0.0000000,    0.0000000]])
        mesh= interp_p2m_lambda3_2d(Part,nPart,part_values,nVal,v1,v2,n1,n2,bRegular)
        np.testing.assert_almost_equal(mesh[0,:,:], mesh0)
        np.testing.assert_almost_equal(mesh[1,:,:], mesh1)
                       
    def test_m2p_mp4_lambda3(self):
        # --- Init
        nPart = 4;
        nDim  = 2;
        nVal  = 3;
        Part      = np.zeros((nPart,nDim));
        Part[0,:] = [0,0]            ;
        Part[1,:] = [1.5,0.5]        ;
        Part[2,:] = [2.0,1.0]        ;
        Part[3,:] = [4,4]            ;
        v1 = np.array([0 ,1,2,3 ,4]);
        v2 = np.array([-2,0,2,4]);
        n1=len(v1);
        n2=len(v2)
        [X,Y]=np.meshgrid(v1,v2);
        mesh=np.zeros((nVal,n1,n2));
        mesh[0,:,:]=1+X.T*0;
        mesh[1,:,:]=2+Y.T*0;
        mesh[2,:,:]=X.T * Y.T;
        bRegular=True
        v_p_ref = np.array([
        [   0.0000000,    0.0000000,    0.0000000],
        [   1.0000000,    2.0000000,    0.7500000],
        [   1.0000000,    2.0000000,    2.0000000],
        [   0.0000000,    0.0000000,   -0.0000000]])
        # --- Test mp4
        v_p= interp_m2p_mp4_2d(Part,nPart,mesh,nVal,v1,v2,n1,n2,bRegular)
        np.testing.assert_almost_equal(v_p, v_p_ref)
        # --- Test lambda3
        v_p= interp_m2p_lambda3_2d(Part,nPart,mesh,nVal,v1,v2,n1,n2,bRegular)
        np.testing.assert_almost_equal(v_p, v_p_ref)

    def test_m2p_p2m(self):
        # test mesh2p and then p2m should give the same
        # NOTE: still sound boundary effects
        # TODO: Make this test low level (no Part and mesh dependencies!)
        from welib.vortilib.particles.initialization import init_from_mesh
        from welib.vortilib.elements.InviscidVortexPatch import ivp_omega
        from welib.mesh import mesh
        #  Setting up mesh
        nx = 7
        ny = 5
        vx = np.linspace(-1.0,1.0,nx)
        vy = np.linspace(-1.0,1.0,ny)
        X,Y = np.meshgrid(vx,vy)
        Omega=np.zeros(X.shape)
        Z=ivp_omega(X,Y,3,False)
        # 
        mesh = mesh.Mesh(vy,vx)
        mesh.values=np.zeros((1,ny,nx))
        mesh.values[0,:,:]=Z
        MeshValues_in = mesh.values
        # --- Mesh2Particles
        #Part = init_from_mesh(mesh, location='CellCenter')
        Part = init_from_mesh(mesh, location='GridPoint')
        # --- Particles2Mesh
        MeshValues,v_p = part2mesh(Part,mesh)
        # --- test
        #import matplotlib.pyplot as plt
        #from mpl_toolkits import mplot3d
        #plt.figure()
        #ax = plt.axes(projection='3d')
        #ax.plot_surface(X[:-2,:-2],Y[:-2,:-2],(MeshValues_in[0,:-2,:-2]-MeshValues[0,:-2,:-2]), edgecolor='none',alpha=1,antialiased=False)
        #plt.show()
        np.testing.assert_almost_equal(MeshValues[0,:-2,:-2],MeshValues_in[0,:-2,:-2]) 

        pass

if __name__=='__main__':
    unittest.main()
