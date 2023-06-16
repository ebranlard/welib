import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# Local 
class PointMesh:
    def __init__(self, nPoints, name='mesh', Connectivity=None, RefPoint=None):
        """
        """
        self.Position        = np.zeros((nPoints, 3))
        self.TranslationDisp = np.zeros((nPoints, 3))
        self.TranslationVel  = np.zeros((nPoints, 3))
        self.TranslationAcc  = np.zeros((nPoints, 3))
        self.RefOrientation  = np.zeros((nPoints, 3,3))
        self.Orientation     = np.zeros((nPoints, 3,3))
        self.RotationVel     = np.zeros((nPoints, 3))
        self.RotationAcc     = np.zeros((nPoints, 3))
        self.Force           = np.zeros((nPoints, 3))
        self.Moment          = np.zeros((nPoints, 3))
        self.Connectivity = Connectivity
        self.name         = name

        self.RefPoint = RefPoint

        self.RefOrientation[:,:,:] = np.eye(3)
        self.Orientation   [:,:,:] = np.eye(3)

    @property
    def nNodes(self): 
        return self.Position.shape[0]

    def rigidBodyMotion(self, RefPoint=None, u=(0,0,0), u_dot=(0,0,0), u_ddot=(0,0,0), theta=(0,0,0), omega=(0,0,0), omega_dot=(0,0,0), R_b2g=None, q=None, qd=None, qdd=None, rot_type='smallRot_OF'):
        """
        Apply a rigid body motion to the mesh:
        u         : translation     of reference point
        u_dot     : translation vel of reference point
        u_ddot    : translation acc of reference point
        R_b2g     : rotation matrix from body 2 global
        omega     : rotational vel of body
        omega_dot : rotational acc of body
        RefPoint: Define the reference point (on the undisplaced mesh) where the Rigid Body Motion is to be applied, very important parameter
        """
        from welib.yams.rotations import rotMat
        # --- Default arguments
        if RefPoint is None:
            RefPoint = self.RefPoint
            if RefPoint is None: 
                raise Exception('Provide a reference point for rigid body motion')

        if qd is None:
            qd=(0,0,0,0,0,0)
        if qdd is None:
            qdd=(0,0,0,0,0,0)
        if q is not None:
            u         = np.array([q[0]   , q[1]   , q[2]])
            u_dot     = np.array([qd[0]  , qd[1]  , qd[2]])
            u_ddot    = np.array([qdd[0] , qdd[1] , qdd[2]])
            theta     = np.array([q[3]   , q[4]   , q[5]])
            omega     = np.array([qd[3]  , qd[4]  , qd[5]]) # NOTE  , small angle approx here
            omega_dot = np.array([qdd[3] , qdd[4] , qdd[5]]) # NOTE , small angle approx here


        if R_b2g is None:
#             R_b2g = BodyXYZ_A(theta[0], theta[1], theta[2])# matrix body 2 global, order XYZ
            R_b2g = rotMat(theta, rot_type)
            #if rot_type=='smallRot_OF':
            #    R_b2g = smallRot_OF(theta[0], theta[1], theta[2]).T # TO MATCH OPENFAST !!!
            #elif rot_type=='smallRot':
            #    R_b2g = smallRot_A(theta[0], theta[1], theta[2])
            #else:
            #    raise Exception('Rotation type not supported: {}'.format(rot_type))
        # --- Sanitation    
        u         = np.asarray(u)
        u_dot     = np.asarray(u_dot)
        u_ddot    = np.asarray(u_ddot)
        theta     = np.asarray(theta)
        omega     = np.asarray(omega)
        omega_dot = np.asarray(omega_dot)
        # --- Motion
        r_AB0  = (self.Position[:,:] - RefPoint)
        r_AB   = (R_b2g.dot(r_AB0.T)).T
        om_x_r = (np.cross(omega, r_AB))
        self.TranslationDisp[:,:] = u + (r_AB - r_AB0)
        self.TranslationVel [:,:] = u_dot  + om_x_r
        self.TranslationAcc [:,:] = u_ddot + np.cross(omega_dot, r_AB) + np.cross(omega, om_x_r)
        self.Orientation [:,:,:]  = R_b2g.T # TODO TODO TODO Shouldn't we also have RefOrientation here?
        self.RotationVel [:,:]    = omega
        self.RotationAcc [:,:]    = omega_dot

    def perturbNode(self, iNode, q):
        from welib.yams.rotations import BodyXYZ_A, smallRot_OF, smallRot_A
        u         = np.array([q[0]   , q[1]   , q[2]])
        theta     = np.array([q[3]   , q[4]   , q[5]])
        R_b2g = smallRot_OF(theta[0], theta[1], theta[2]).T # TO MATCH OPENFAST !!!
        # --- Motion
        self.TranslationDisp[iNode,:] += u 
        self.Orientation  [iNode,:,:]  = R_b2g.T # TODO?

    def mapLoadsToPoint(self, P):
        """ Map Force and Moment fields to a given point"""
        P = np.asarray(P)
        F = np.zeros(3)
        M = np.zeros(3)
        r_PP0 = self.Position[:,:] + self.TranslationDisp[:,:] - P
        dM = np.cross(r_PP0, self.Force[:,:])
        F = np.sum(self.Force[:,:], axis=0) 
        M = np.sum(self.Moment+dM , axis=0) 
        return F, M

    def transferMotion2IdenticalMesh(self, target):
        target.TranslationDisp = self.TranslationDisp.copy()
        target.TranslationVel  = self.TranslationVel .copy()
        target.TranslationAcc  = self.TranslationAcc .copy()
        target.Orientation     = self.Orientation    .copy()
        target.RotationVel     = self.RotationVel    .copy()
        target.RotationAcc     = self.RotationAcc    .copy()

    def backupValues(self):
        """ Backup mesh motion values, useful at an operating point """
        self.TranslationDisp_0 = self.TranslationDisp.copy()
        self.TranslationVel_0  = self.TranslationVel .copy()
        self.TranslationAcc_0  = self.TranslationAcc .copy()
        self.Orientation_0     = self.Orientation    .copy()
        self.RotationVel_0     = self.RotationVel    .copy()
        self.RotationAcc_0     = self.RotationAcc    .copy()

    def restoreValues(self):
        """ Resstore stored values by backupValues() """
        self.TranslationDisp = self.TranslationDisp_0.copy()
        self.TranslationVel  = self.TranslationVel_0.copy()
        self.TranslationAcc  = self.TranslationAcc_0.copy()
        self.Orientation     = self.Orientation_0.copy()
        self.RotationVel     = self.RotationVel_0.copy()
        self.RotationAcc     = self.RotationAcc_0.copy()

    def packLoads(self, FthenM=True):
        if FthenM:
            # Force for all nodes + Moment for all nodes
            M = np.vstack((self.Force,self.Moment))
        else:
            # Force+Moment for each node
            M=np.zeros((self.nNodes*2, 3))
            for i in range(self.nNodes):
                M[i*3+0: i*3+3] = self.Force[i,:]
                M[i*3+3: i*3+6] = self.Moment[i,:]
        return M.flatten()

    # --------------------------------------------------------------------------------}
    # --- IO  
    # --------------------------------------------------------------------------------{
    def toConnectedObject(self, diameters=5):
        from welib.FEM.graph import ConnectedObject
        if self.Connectivity is None:
            #  TODO set a connectivity with one point per element..
            raise NotImplementedError('Cannot convert to connected object when Connectivity is None')

        Props={'shape':'cylinder','type':0, 'Diam':diameters} # TODO
        Nodes        = self.Position
        Connectivity = self.Connectivity
        self.obj = ConnectedObject(self.name, Nodes, Connectivity, Props)
        return self.obj

    def printDebug(self):
        for j in range(self.nNodes):
            print('Node {:6d}'.format(j+1))
            print('pos  {:12.6f}{:12.6f}{:12.6f}'.format(*self.Position       [j,:]))
            print('disp {:12.6f}{:12.6f}{:12.6f}'.format(*self.TranslationDisp[j,:]))
            print('Tvel {:12.6f}{:12.6f}{:12.6f}'.format(*self.TranslationVel [j,:]))
            print('Rvel {:12.6f}{:12.6f}{:12.6f}'.format(*self.RotationVel    [j,:]))
            print('Tacc {:12.6f}{:12.6f}{:12.6f}'.format(*self.TranslationAcc [j,:]))
            print('Racc {:12.6f}{:12.6f}{:12.6f}'.format(*self.RotationAcc    [j,:]))
            print('ROri {:12.6f}{:12.6f}{:12.6f}'.format(*self.RefOrientation [j,0,:]))
            print('ROri {:12.6f}{:12.6f}{:12.6f}'.format(*self.RefOrientation [j,1,:]))
            print('ROri {:12.6f}{:12.6f}{:12.6f}'.format(*self.RefOrientation [j,2,:]))
            print('Orie {:12.6f}{:12.6f}{:12.6f}'.format(*self.Orientation    [j,0,:]))
            print('Orie {:12.6f}{:12.6f}{:12.6f}'.format(*self.Orientation    [j,1,:]))
            print('Orie {:12.6f}{:12.6f}{:12.6f}'.format(*self.Orientation    [j,2,:]))



# --------------------------------------------------------------------------------}
# --- Storage 
# --------------------------------------------------------------------------------{
class MeshStorage():
    def __init__(self, mesh, time, name='TS1'):
        nt = len(time)
        nPoints = mesh.nNodes
        self.vTranslationDisp = np.zeros((nt, nPoints, 3))
        self.vTranslationVel  = np.zeros((nt, nPoints, 3))
        self.vTranslationAcc  = np.zeros((nt, nPoints, 3))
        self.vOrientation     = np.zeros((nt, nPoints, 3,3))
        self.vRotationVel     = np.zeros((nt, nPoints, 3))
        self.vRotationAcc     = np.zeros((nt, nPoints, 3))
        self.vForce           = np.zeros((nt, nPoints, 3))
        self.vMoment          = np.zeros((nt, nPoints, 3))
        self.name='TS1'
        self.time=time

    def store(self, mesh, it):
        # TODO reshape all this
        self.vTranslationDisp[it, :, :3] = mesh.TranslationDisp
        self.vTranslationVel [it, :, :3] = mesh.TranslationVel
        self.vTranslationAcc [it, :, :3] = mesh.TranslationAcc
        #self.vOrientation     = np.zeros((nt, nPoints, 3,3)) # TODO reshape..
        self.vRotationVel [it, :, :3]    = mesh.RotationVel
        self.vRotationAcc [it, :, :3]    = mesh.RotationAcc
        self.vForce       [it, :, :3]    = mesh.Force
        self.vMoment      [it, :, :3]    = mesh.Moment

    def toConnectedObject(self, mesh):
        obj = mesh.toConnectedObject()
        obj.addTimeSeries(self.time, displ=self.vTranslationDisp, name=self.name)
        return obj

    def toJSON3D(self, mesh, jsonfile):
        from welib.plot.json3d import JSON3DFile
        obj = self.toConnectedObject(mesh)
        js = JSON3DFile()
        js.addObject(obj)
        js.write(jsonfile)






