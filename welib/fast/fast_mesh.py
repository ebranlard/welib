import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# Local 
class PointMesh:
    def __init__(self, nPoints, name='mesh', RefPoint=None):
        # NOTE: for now, using Fortran convention for memory order
        self.Position         = np.zeros((3, nPoints))
        self.TranslationDisp  = np.zeros((3, nPoints))
        self.TranslationVel   = np.zeros((3, nPoints))
        self.TranslationAcc   = np.zeros((3, nPoints))
        self.RefOrientation = np.zeros((3,3, nPoints))
        self.Orientation    = np.zeros((3,3, nPoints))
        self.RotationVel = np.zeros((3, nPoints))
        self.RotationAcc = np.zeros((3, nPoints))
        self.Force  = np.zeros((3, nPoints))
        self.Moment = np.zeros((3, nPoints))
        self.Connectivity=None
        self.RefPoint = RefPoint
        self.name = name

        for j in range(self.nNodes):
            self.RefOrientation[:,:,j] = np.eye(3)
            self.Orientation   [:,:,j] = np.eye(3)

    @property
    def nNodes(self): 
        return self.Position.shape[1]

    def rigidBodyMotion(self, u=(0,0,0), u_dot=(0,0,0), u_ddot=(0,0,0), theta=(0,0,0), omega=(0,0,0), omega_dot=(0,0,0), R_b2g=None, q=None, qd=None, qdd=None):
        """
        Apply a rigid body motion to the mesh:
        u         : translation     of reference point
        u_dot     : translation vel of reference point
        u_ddot    : translation acc of reference point
        R_b2g     : rotation matrix from body 2 global
        omega     : rotational vel of body
        omega_dot : rotational acc of body
        """
        from welib.yams.rotations import BodyXYZ_A, SmallRot_DCM
        # --- Default arguments
        if self.RefPoint is None:
            self.RefPoint = self.Position[:,0]

        if q is not None:
            u         = np.array([q[0], q[1], q[2]])
            u_dot     = np.array([q[0], q[1], q[2]])
            u_ddot    = np.array([q[0], q[1], q[2]])
            theta     = np.array([q[3], q[4], q[5]])
            omega     = np.array([q[3], q[4], q[5]]) # NOTE, small angle approx here
            omega_dot = np.array([q[3], q[4], q[5]]) # NOTE, small angle approx here


        if R_b2g is None:
#             R_b2g = BodyXYZ_A(theta[0], theta[1], theta[2])# matrix body 2 global, order XYZ
            R_b2g = SmallRot_DCM(theta[0], theta[1], theta[2]).T # TO MATCH OPENFAST !!!
        # --- Sanitation    
        u         = np.asarray(u)
        u_dot     = np.asarray(u_dot)
        u_ddot    = np.asarray(u_ddot)
        theta     = np.asarray(theta)
        omega     = np.asarray(omega)
        omega_dot = np.asarray(omega_dot)
        # --- Motion

        for j in range(self.nNodes):
            r_AB0 = self.Position[:,j] - self.RefPoint
            r_AB  = R_b2g.dot(r_AB0)
            om_x_r = np.cross(omega, r_AB)
            self.TranslationDisp[:,j] = u + (r_AB - r_AB0)
            self.TranslationVel [:,j] = u_dot  + om_x_r
            self.TranslationAcc [:,j] = u_ddot + np.cross(omega_dot, r_AB) + np.cross(omega, om_x_r)
            self.Orientation [:,:,j]  = R_b2g.T
            self.RotationVel [:,j]    = omega
            self.RotationAcc [:,j]    = omega_dot

    def mapLoadsToPoint(self,P, R=None):
        """ Map Force and Moment fields to a given point"""
        P = np.asarray(P)
        F = np.zeros(3)
        M = np.zeros(3)
        if R is None:
            R=np.eye(3)
        for j in range(self.nNodes):
            F0 = self.Force[:,j]
            M0 = self.Moment[:,j]
            P0 = self.Position[:,j] + self.TranslationDisp[:,j]
            r = P0-P
            dM = np.cross(r, F0)
            F+=F0
            Mloc = R.dot(M0+dM)
#             if j==0:
#                 print('P    {:16.3f}{:16.3f}{:16.3f}'.format(*P))
#                 print('Pos  {:16.3f}{:16.3f}{:16.3f}'.format(*self.Position[:,j]))
#                 print('TD   {:16.3f}{:16.3f}{:16.3f}'.format(*self.TranslationDisp[:,j]))
#                 print('P0   {:16.3f}{:16.3f}{:16.3f}'.format(*P0))
#                 print('r    {:16.3f}{:16.3f}{:16.3f}'.format(*r))
#                 print('F0   {:16.3f}{:16.3f}{:16.3f}'.format(*F0))
#                 print('M0   {:16.3f}{:16.3f}{:16.3f}'.format(*M0))
#                 print('dM   {:16.3f}{:16.3f}{:16.3f}'.format(*dM))
#                 print('Mloc {:16.3f}{:16.3f}{:16.3f}'.format(*Mloc))
            M+=R.dot(M0 + dM)
            #M+=M0 
        return F, M

    def transferMotion2IdenticalMesh(self, target):
        target.TranslationDisp = self.TranslationDisp
        target.TranslationVel  = self.TranslationVel 
        target.TranslationAcc  = self.TranslationAcc 
        target.Orientation     = self.Orientation    
        target.RotationVel     = self.RotationVel    
        target.RotationAcc     = self.RotationAcc    

    # --------------------------------------------------------------------------------}
    # --- IO  
    # --------------------------------------------------------------------------------{
    def toConnectedObject(self, diameters=5):
        from welib.FEM.graph import ConnectedObject
        if self.Connectivity is None:
            #  TODO set a connectivity with one point per element..
            raise NotImplementedError('Cannot convert to connected object when Connectivity is None')

        Props={'shape':'cylinder','type':0, 'Diam':diameters} # TODO
        Nodes        = self.Position.T
        Connectivity = self.Connectivity
        self.obj = ConnectedObject(self.name, Nodes, Connectivity, Props)
        return self.obj

    def printDebug(self):
        for j in range(self.nNodes):
            print('Node {:6d}'.format(j+1))
            print('pos  {:12.6f}{:12.6f}{:12.6f}'.format(*self.Position       [:,j]))
            print('disp {:12.6f}{:12.6f}{:12.6f}'.format(*self.TranslationDisp[:,j]))
            print('Tvel {:12.6f}{:12.6f}{:12.6f}'.format(*self.TranslationVel [:,j]))
            print('Rvel {:12.6f}{:12.6f}{:12.6f}'.format(*self.RotationVel    [:,j]))
            print('Tacc {:12.6f}{:12.6f}{:12.6f}'.format(*self.TranslationAcc [:,j]))
            print('Racc {:12.6f}{:12.6f}{:12.6f}'.format(*self.RotationAcc    [:,j]))
            print('ROri {:12.6f}{:12.6f}{:12.6f}'.format(*self.RefOrientation [0,:,j]))
            print('ROri {:12.6f}{:12.6f}{:12.6f}'.format(*self.RefOrientation [1,:,j]))
            print('ROri {:12.6f}{:12.6f}{:12.6f}'.format(*self.RefOrientation [2,:,j]))
            print('Orie {:12.6f}{:12.6f}{:12.6f}'.format(*self.Orientation    [0,:,j]))
            print('Orie {:12.6f}{:12.6f}{:12.6f}'.format(*self.Orientation    [1,:,j]))
            print('Orie {:12.6f}{:12.6f}{:12.6f}'.format(*self.Orientation    [2,:,j]))



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
        self.vTranslationDisp[it, :, :3] = mesh.TranslationDisp.T
        self.vTranslationVel [it, :, :3] = mesh.TranslationVel.T
        self.vTranslationAcc [it, :, :3] = mesh.TranslationAcc.T
        #self.vOrientation     = np.zeros((nt, nPoints, 3,3)) # TODO reshape..
        self.vRotationVel [it, :, :3]    = mesh.RotationVel.T
        self.vRotationAcc [it, :, :3]    = mesh.RotationAcc.T
        self.vForce       [it, :, :3]    = mesh.Force.T
        self.vMoment      [it, :, :3]    = mesh.Moment.T

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






