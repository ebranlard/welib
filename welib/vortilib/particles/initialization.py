import numpy as np
from .particles import Particles
from .projection import interp_m2p


def init_from_mesh_and_fun(fOmega, mesh, location='CellCenter'):
    """ Initialize a set of particles from mesh volumes, mesh nodes, and a vortiity function """
    nDim = mesh.nDim
    if location=='CellCenter':
        # Putting particles at the middle of the cell
        cell_volumes = mesh.getCellVolumes()
        cell_Centers = mesh.getFlatCellCenters()
        nPart = cell_Centers.shape[0]
        Part = Particles(nPart=nPart, nDim=mesh.nDim)
        Part.setP(cell_Centers)
        # Computing Omega from analytical function and setting particle intensities
        if nDim==2:
            Omega = fOmega(Part.P[:,0],Part.P[:,1])
            Omega=Omega.reshape(Part.Intensity.shape)
            cell_volumes=cell_volumes.reshape(Omega.shape)
        elif nDim==2:
            Omega = fOmega(Part.P[:,0],Part.P[:,1],Part.P[:,2])
        else:
            raise
        Part.setIntensity(Omega * cell_volumes)
        Part.setVolume(cell_volumes)
    else:
        raise NotImplementedError()

    return Part

#     if ('mesh_GridPoints_OmegaAnalytical') == (InitParams.Method):
#         mesh = InitParams.mesh
#         # Putting particles at the middle of the cell
#         cell_volumes = mesh.getCellVolumes()
#         PP = mesh.getFlatGridPoints()
#         nPart = PP.shape[1-1]
#         Part.reset(nPart)
#         Part.setP(PP)
#         # Computing Omega from analytical function and setting particle intensities
#         if Algo.nDim == 2:
#             Omega = InitParams.fOmega(Part.P(:,1),Part.P(:,2))
#         else:
#             Omega = InitParams.fOmega(Part.P(:,1),Part.P(:,2),Part.P(:,3))
#         Part.setIntensity(np.multiply(Omega,cell_volumes))
#         Part.setVolume(cell_volumes)


# if ('mesh_CellCenterQuasiRandom_OmegaAnalytical') == (InitParams.Method):
#    mesh = InitParams.mesh
#    if not mesh.bRegular :
#        raise Exception('Random init with non regular grid not suported')
#    # Putting particles at the middle of the cell
#    cell_volume = mesh.getCellVolume()
#    cell_Centers = mesh.getFlatCellCenters()
#    nPart = cell_Centers.shape[1-1]
#    # rand is between [0 1] we transform it between -1/2h 1/2h
#    RandP = np.zeros((nPart,Algo.nDim))
#    for id in np.arange(1,Algo.nDim+1).reshape(-1):
#        RandP[:,id-1] = (np.random.rand(nPart,1) - 1 / 2) * mesh.dCell(id)
#    # Setting Particles position
#    Part.reset(nPart)
#    Part.setP(cell_Centers + RandP)
#    # Computing Omega from analytical function and setting particle intensities
#    if Algo.nDim == 2:
#        Omega = InitParams.fOmega(Part.P(:,1),Part.P(:,2))
#    else:
#        Omega = InitParams.fOmega(Part.P(:,1),Part.P(:,2),Part.P(:,3))
#    Part.setIntensity(Omega * cell_volume)
#    Part.setVolume(cell_volume)

def init_from_mesh(mesh, location='CellCenter', kernel='mp4'):
    if location=='CellCenter':
        # Putting particles at the middle of the cell
        cell_volumes = mesh.getCellVolumes()
        cell_Centers = mesh.getFlatCellCenters()
        nPart = cell_Centers.shape[0]
        Part = Particles(nPart=nPart, nDim=mesh.nDim)
        Part.setP(cell_Centers)
        v_p = interp_m2p(Part.P,Part.nPart,mesh.n,mesh.values,mesh.nDim,kernel=kernel,v1=mesh.v1,v2=mesh.v2,v3=mesh.v3,bRegular=mesh.bRegular)
        if mesh.nDim == 2:
            Part.setIntensity(v_p[:,0])
            cell_volumes=cell_volumes.ravel()
        else:
            Part.setIntensity(v_p[:,0:3])
        Part.setVolume(cell_volumes)
    elif location=='GridPoint':
        # Putting particles at the grid points
        cell_volumes = mesh.getFakeCellVolumes()
        grid_points  = mesh.getFlatGridPoints()
        nPart = grid_points.shape[0]
        Part = Particles(nPart=nPart, nDim=mesh.nDim)
        Part.setP(grid_points)
        v_p = interp_m2p(Part.P,Part.nPart,mesh.n,mesh.values,mesh.nDim,kernel=kernel,v1=mesh.v1,v2=mesh.v2,v3=mesh.v3,bRegular=mesh.bRegular)
         #    meshVal_flat = mesh.flattenValues(InitParams.meshValues)
        if mesh.nDim == 2:
             Part.setIntensity(v_p[:,0])
             cell_volumes=cell_volumes.ravel()
        else:
             Part.setIntensity(v_p[:,0:3])
        Part.setVolume(cell_volumes)
#        meshIntensity = np.transpose(meshVal_flat(1,:))
#        meshVol = np.transpose(meshVal_flat(2,:))
#    else:
#        meshIntensity = np.transpose(meshVal_flat(np.arange(1,3+1),:))
#        meshVol = np.transpose(meshVal_flat(2,:))
#        #meshVol(meshVol<=0)=cell_volumes(meshVol<=0);
#        meshOmega = meshIntensity / meshVol
#        #             kbd
#        Part.setIntensity(meshOmega * cell_volumes))
#        #             Part.setIntensity(meshIntensity); # alpha_p=omega*V
#        #             v_p= interp_m2p_mex(Part.P,Part.nPart,mesh.xmesh_min,mesh.dCell,mesh.n, InitParams.meshValues, Algo.nDim, Algo.InterpolationKernel);
#        #Part.setIntensity(v_p(:,1)); #  Gamma=omega*A
#             else:
#                 log_error('Not implemented')
# # 
# #     return Part,Field

    else:
        raise NotImplementedError()
    return Part


# if ('mesh_Random_OmegaAnalytical') == (InitParams.Method):
#     mesh = InitParams.mesh
#     # finding spatial extent and total number
#     nPart = 1
#     TotalVolume = 1
#     extent = np.zeros((1,Algo.nDim))
#     xmin = np.zeros((1,Algo.nDim))
#     for id in np.arange(1,Algo.nDim+1).reshape(-1):
#         nPart = nPart * mesh.n(id)
#         extent[id-1] = mesh.xmesh_max(id) - mesh.xmesh_min(id)
#         xmin[id-1] = mesh.xmesh_min(id)
#         TotalVolume = TotalVolume * extent(id)
#     part_volume = TotalVolume / nPart
#     # rand is between [0 1] we transform it between  [0 extent] + xmin
#     RandP = np.zeros((nPart,Algo.nDim))
#     for id in np.arange(1,Algo.nDim+1).reshape(-1):
#         RandP[:,id-1] = np.random.rand(nPart,1) * extent(id) + xmin(id)
#     #  Setting Particle Positions
#     Part.reset(nPart)
#     Part.setP(RandP)
#     # Computing Omega from analytical function and setting particle intensities
#     if Algo.nDim == 2:
#         Omega = InitParams.fOmega(Part.P(:,1),Part.P(:,2))
#     else:
#         Omega = InitParams.fOmega(Part.P(:,1),Part.P(:,2),Part.P(:,3))
#     Part.setIntensity(Omega * part_volume)
#     Part.setVolume(part_volume)
