""" 
Quick summary of functions:

vdel_* : (A.nabla)B
curl_* : curl(B)

intgrad: (grad f) -> f

"""

import numpy as np
#import numpy.matlib
import scipy.sparse as sparse
import scipy.sparse.linalg

def vdel_cylindrical_axi(A, B, vr, vz, shape_in=('r','z')):
    """ Compute (A.nabla)B  assuming axi-symmetric cylindrical coordiates 
    INPUTS:
        A = (Ar, Az) with Ar (nr x nz)
        B = (Br, Bz) with Ar (nr x nz)
    OUTPUTS:
       (A.nabla)B
    """
    Ar, Az = A
    Br, Bz = B
    if shape_in==('r','z'):
        dBr_dr, dBr_dz = np.gradient(Br, vr, vz)
        dBz_dr, dBz_dz = np.gradient(Bz, vr, vz)
    elif shape_in==('z','r'):
        dBr_dz, dBr_dr = np.gradient(Br, vz, vr)
        dBz_dz, dBz_dr = np.gradient(Bz, vz, vr)
    return (
            Ar*dBr_dr + Az*dBr_dz
            , Ar*dBz_dr + Az*dBz_dz
            )
def curl_cylindrical_axi(A, vr, vz, shape_in=('r','z')):
    """ Compute curl(A) assuming axi-symmetric cylindrical coordiates 
    """
    Ar, Az = A
    if shape_in==('r','z'):
        _, dAr_dz = np.gradient(Ar, vr, vz)
        dAz_dr, _ = np.gradient(Az, vr, vz)
    elif shape_in==('z','r'):
        pass
        #
    return dAr_dz - dAz_dr # along theta

    
def intgrad(gf, nodes, constant=0, tol=1e-16):
    """ 
    Integrates a gradient field (inverse of np.gradient) in 1D, 2D or 3D.
    Finds f such that grad(f) = gf.
    A second order finite difference method is used, setup as A.f = gf.

    INPUTS:
     - gf: (tuple of) gradient fields in each directions:
         in 1D: 
           gf = (fx) of fx  : fx is the vector of values df/dx of length nx
         in 2D: 
           gf = (fx,fy)     : fx=df/dx, fy=df/dy, of shape (ny x nx)
         in 3D: 
           gf = (fx,fy,fz)  : fx=df/dx, fy=df/dy, fz=df/fz, of shape (ny x nx, x nz)

     - nodes: (tuple of) vectors defining the position of the nodes in each directions
         in 1D: 
           nodes=(x) or x  : x of length nx 
         in 2D: 
           nodes=(x, y)    : x of length nx, y of length ny
         in 3D: 
           nodes=(x, y, z) : x of length nx, y of length ny, z of length nz

     - constant: Constant use to offset the integrated field
                This value will be prescribed at the lowest corner (xmin, ymin, zmin)

     - tol: tolerance used in the least sqaure inversion of the finite difference system

    OUTPUTS:
     - f : integrated gradient field such that grad(f) = gf
           The integration constant if such that f(xmin) = constant

    AUTHOR: E. Branlard

    Inspired by Matlab script from John D'Errico:
        https://www.mathworks.com/matlabcentral/fileexchange/9734-inverse-integrated-gradient

    """
    # --- Default arguments
    if not isinstance(gf, tuple):
        gf=(gf,)
    if not isinstance(nodes, tuple):
        nodes=(nodes,)

    nDim = len(gf)
    dims = gf[0].shape
    # --- Safety checks
    for i,fi in enumerate(gf):
        if len(fi.shape)!=nDim:
            raise Exception('Inconsistent dimensions. Gradient field number {} has shape {} but should have maximum {} dimensions'.format(i, fi.shape, nDim))
        if fi.shape!=dims:
            raise Exception('Gradient fields should have the same shapes. Gradient field number {} has shape {} but the first field has shape {}'.format(i, fi.shape, dims))
    if any(np.array(dims)<2):
        raise Exception('Fields must have at least 2 values in each direction')
    if len(nodes)!=nDim:
        raise Exception('Nodes argument should have the same length as field argument')

    delta = [0]*nDim
    gradf = [0]*nDim
    for i,x in enumerate(nodes):
        dx       = np.diff(x)
        delta[i] = np.concatenate( (dx.ravel(), [dx[-1]]) ) # NOTE we duplicate the last dx for convenience
        gradf[i] = gf[i].ravel(order='F') 

    # --- For now on we stop the general dimension until we decide on an order
    if nDim==1:
        nx,ny,nz = dims[0], 1, 1
    elif nDim==2:
        ny,nx,nz = dims[0], dims[1], 1       # NOTE order, due to mesh grid!
    elif nDim==3:
        ny,nx,nz = dims[0], dims[1], dims[2] # NOTE order, due to mesh grid!
    else:
        raise Exception('Only implemented in 1D, 2D and 3D for now')

    # --- Setting up finite difference system
    # The system will be stored in a sparse matrix A at the end
    neq     = nDim*nx*ny*nz   # number of difference equations
    nCoeffs = 2               # number of coefficients used in finite diff
    row  = np.zeros((neq, nCoeffs)) # row indices of sparse matrix 
    col  = np.zeros((neq, nCoeffs)) # col indices of sparse matrix 
    dat  = np.zeros((neq, nCoeffs)) # values of sparse matrix 
    rhs  = np.zeros((neq))
    L = 0 # Cumulative index

    def diff(iField, kind, di, ix, iy, iz, L):
        """ helper function to store finite diff coefficients and positions in A matrix """
        fxi = gradf[iField]  # field of interest, e.g. fx, fy or fz
        dxi = delta[iField]  # diffs of mesh vector in direction of interest, e.g. diff(x)
        if kind=='central' and any(np.array([len(ix),len(iy),len(iz)]))==0:
            return L  # No central finite diff if not enough indices
        # Getting flat indices
        iX,iY,iZ = np.meshgrid(ix,iy,iz)
        iC   = (iY + iX*ny + iZ*ny*nx).ravel() 
        ixi = (iX,iY,iZ)[iField].ravel() # flat indices for the direction of interest only
        m = len(iC)
        if   kind=='forward':
            coeff = 1.0 / dxi[ixi]
            col[ L:L+m, : ] = np.column_stack( (iC   , iC+di)  )
        elif kind=='backwrd':
            coeff = 1.0 / dxi[ixi]
            col[ L:L+m, : ] = np.column_stack( (iC-di, iC   )  )
        elif kind=='central':
            coeff = 1.0 / (dxi[ixi-1]+dxi[ixi])
            col[ L:L+m, : ] = np.column_stack( (iC-di, iC+di)  )
        row[ L:L+m, : ] = np.column_stack( [np.arange(L,L+m)]*2 )
        dat[ L:L+m, : ] = np.column_stack( (coeff*-1, coeff*1)  )
        rhs[ L:L+m    ] = fxi[iC]
        return L+m

    if nDim>0: # --- dfx/dx finite differences
        L = diff(0, 'forward', ny     , [0   ]           , np.arange(0,ny)  , np.arange(0,nz)  , L)
        L = diff(0, 'backwrd', ny     , [nx-1]           , np.arange(0,ny)  , np.arange(0,nz)  , L)
        L = diff(0, 'central', ny     , np.arange(1,nx-1), np.arange(0,ny)  , np.arange(0,nz)  , L)
    if nDim>1: # --- dfy/dy finite differences
        L = diff(1, 'forward', 1      ,  np.arange(0,nx) , [0]              , np.arange(0,nz)  , L)
        L = diff(1, 'backwrd', 1      ,  np.arange(0,nx) , [ny-1]           , np.arange(0,nz)  , L)
        L = diff(1, 'central', 1      ,  np.arange(0,nx) , np.arange(1,ny-1), np.arange(0,nz)  , L)
    if nDim>2: # --- dfz/dz finite differences
        L = diff(2, 'forward', (nx*ny), np.arange(0,nx)  , np.arange(0,ny)  , [0]              , L)
        L = diff(2, 'backwrd', (nx*ny), np.arange(0,nx)  , np.arange(0,ny)  , [nz-1]           , L)
        L = diff(2, 'central', (nx*ny), np.arange(0,nx)  , np.arange(0,ny)  , np.arange(1,nz-1), L)

    # --- Ravelling for scipy sparse
    dat = dat.ravel()
    row = row.ravel()
    col = col.ravel()

    # --- First column of A matrix (needed to eliminate constant of integration)
    b   =  col==0
    bNot = np.logical_not(b)
    dat1 = dat[b]
    row1 = row[b]
    col1 = col[b]
    A0 = sparse.coo_matrix( (dat1, (row1, col1)), shape=(neq, 1 )).toarray().ravel()

    # --- Remaining columns (sparse matrix)
    # The system is reduced by one column (used for the integration constant) to make it of full rank 
    dat2 = dat[bNot]
    row2 = row[bNot]
    col2 = col[bNot]-1
    A    = sparse.coo_matrix( (dat2, (row2, col2)), shape=(neq, nx*ny*nz-1))

    # --- Eliminate the first column using the constant of integration 
    rhs = rhs - A0 * constant

    # --- Solve the final system of equations
    f = sparse.linalg.lsqr(A, rhs, atol= tol, btol=tol)
    f = f[0]
    if nDim==1:
        f = np.concatenate(([constant],f))
    elif nDim==2:
        f = np.reshape(np.concatenate(([constant],f)),(ny,nx), order='F')
    elif nDim==3:
        f = np.reshape(np.concatenate(([constant],f)),(ny,nx,nz), order='F')

    return f



if __name__=='__main__':
    import matplotlib.pyplot as plt

    np.set_printoptions(edgeitems=3,infstr='inf', linewidth=175, nanstr='NA', precision=5,
     suppress=False, threshold=1000, formatter=None)
    # --- 1D
    dx=0.01
    xp = np.arange(0,1+dx/2,dx)
    f  = np.exp(xp)
    fx = np.exp(xp)
#     fhat = intgrad(fx,xp, constant = 1)
    #fx = np.gradient(f,xp)
#     fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
#     fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
#     ax.plot(xp, fx   ,'k-' , label='')
#     ax.plot(xp, fhat ,'o'  , label='')
#     ax.set_xlabel('')
#     ax.set_ylabel('')
#     ax.legend()
#     ax.tick_params(direction='in')
#     plt.show()

    # --- 2D
    dx=0.01
    dy=0.02
    xp     = np.arange(0,1+dx/2,dx)
    yp     = np.arange(0,1+dx/2,dy)
    # NOTE: mesh grid has shape: len(y) x len(x), i.e. first dim is y
    x,y    = np.meshgrid(xp,yp)       
    f      = np.exp(x+y) + np.sin((x-2*y)*3)
    fy,fx  = np.gradient(f, yp, xp)      
    fhat = intgrad((fx,fy),(xp,yp),constant = 1)

#     # --- 3D
#     xp = np.linspace(0,1,4)
#     yp = np.linspace(0,1,5)
#     zp = np.linspace(0,1,3)
#     [x,y,z] = np.meshgrid(xp,yp,zp);
#     f = np.exp(x+y+z) + np.sin((x-2*y+3*z)*3);
#     [fy,fx,fz]=np.gradient(f,yp,xp,zp)
#     fhat = intgrad((fx,fy,fz),(xp,yp,zp), constant = 1)
#     tic,fhat = intgrad3(fx,fy,fz,xp,yp,zp,1);toc
#     %   Time required was 51.4336 seconds
# 
#     std(f(:)-fhat(:))
#     from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import
#     import matplotlib.pyplot as plt
#     fig =plt.figure()
#     ax = fig.gca(projection='3d')
# #     surf = ax.plot_surface(x, y, f   , linewidth=0, antialiased=False, color='r')
# #     surf = ax.plot_surface(x, y, fhat+0.1, linewidth=0, antialiased=False, color='b')
#     surf = ax.plot_surface(x, y, f-fhat, linewidth=0, antialiased=False, color='b')
#     plt.show()
    print(  np.std(f.ravel() - fhat.ravel())  )
    print(  np.mean(np.abs(f.ravel() - fhat.ravel()))  )





    pass

