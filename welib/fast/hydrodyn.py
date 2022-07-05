import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
# Local 
from welib.tools.clean_exceptions import *
from welib.weio import FASTInputFile

# Submodules
from welib.fast.hydrodyn_morison import Morison
from welib.fast.hydrodyn_waves import Waves


class HydroDyn:

    def __init__(self, filename=None, hdData=None):
        """ 
        INPUTS:
          - filename: optional, name of HydroDyn input file, or OpenFAST or HydroDyn Driver input file
                  NOTE: if an fst or dvr file is provided, init is called
        """
        # Main object data
        self.File        = None
        self.outFilename = None
        self.p           = {}
        self.m           = {}
        self.u           = {}
        self.y           = {}
        self._graph      = None

        # Read HydroDyn file (and optionally Driver file)
        if filename is not None:
            File = FASTInputFile(filename)
            if 'WaveMod' in File.keys(): # User provided a HydroDyn input file
                self.File = File
                self.outFilename = os.path.splitext(filename)[0]+'.pyHD.outb'
            else:
                if 'HydroFile' in File.keys(): # User provided and OpenFAST inputs file
                    hdFilename = os.path.join(os.path.dirname(filename), File['HydroFile'].replace('"','') )
                    outFilenames= [filename.replace('.fst',ext) for ext in ['.outb','.out'] if os.path.exists(filename.replace('.fst',ext))]
                    if len(outFilenames)==0:
                        self.outFilename = filename.replace('.fst','.outb')
                    else:
                        self.outFilename = outFilenames[0]

                elif 'HDInputFile' in File.keys(): # User provided a hydrodyn driver input file
                    hdFilename = os.path.join(os.path.dirname(filename), File['HDInputFile'].replace('"','') )
                    self.outFilename = filename.replace('.dvr','.HD.out')
                else:
                    raise Exception()

                self.File = FASTInputFile(hdFilename)

                # Calling init since we know important environmental conditions
                self.init(Gravity = File['Gravity'], WtrDens=File['WtrDens'], WtrDpth=File['WtrDpth'])


        elif hdData is not None:
            self.File = hdData


    def __repr__(self):
        s='<{} object>:\n'.format(type(self).__name__)
        s+='|properties:\n'
        s+='|- File: (input file data)\n'
        s+='|methods:\n'
        s+='|- init\n'
        return s

    # --------------------------------------------------------------------------------}
    # --- Functions for general FEM model (jacket, flexible floaters)
    # --------------------------------------------------------------------------------{
    def init(self, Gravity = 9.81, WtrDens='1025', WtrDpth=0, MSL2SWL=0):
        """
        Initialize HydroDyn model 

        gravity: position of transition point
        """
        u=dict() # guess
        y=dict() # guess


        f = self.File
        # Handle "default" value in input file
        try:
            WtrDpth=float(f['WtrDpth']) # If default
        except:
            if WtrDpth is None:
                raise Exception('Provide WtrDepth if default in file')
        try:
            MSL2SWL=float(f['MSL2SWL']) # If default
        except:
            if MSL2SWL is None:
                raise Exception('Provide MSL2SWL if default in file')
        try:
            WtrDens=float(f['WtrDens']) # If default
        except:
            if WtrDens is None:
                raise Exception('Provide WtrDens if default in file')

        self.p['WtrDpth'] = WtrDpth + MSL2SWL
        self.p['WtrDens'] = WtrDens 
        self.p['Gravity'] = Gravity

        graph = self.graph

        # --- Morison
        # NOTE: graph will be copied
        # Division occurs in constructor
        self.morison = Morison(graph=self.graph, File=self.File, WtrDpth=self.p['WtrDpth'], MSL2SWL=MSL2SWL)
        # TODO there is  mess with MSL2SWL
        Nodes = self.morison.NodesBeforeSwap

        # --- Transfer of nodes 
        Waves_WaveKin_Nodes = Nodes
        Current_NodesZ      = Nodes[:,2]

        # --- Waves Inits
        self.Waves = Waves(File=self.File, WtrDpth=self.p['WtrDpth'], MSL2SWL=MSL2SWL)
        self.Waves.init(Gravity=self.p['Gravity'])

        # --- WvStretch_Init in HydroDyn.f90
        NStepWave   = self.Waves.p['NStepWave']
        nodeInWater = np.zeros( (NStepWave, len(Waves_WaveKin_Nodes)    ))
        WaveDynP    = np.zeros( (NStepWave, len(Waves_WaveKin_Nodes)    ))
        WaveVel     = np.zeros( (NStepWave, len(Waves_WaveKin_Nodes), 3 ))
        WaveAcc     = np.zeros( (NStepWave, len(Waves_WaveKin_Nodes), 3 ))
        WtrDpth   = self.p['WtrDpth']
        if f['WaveStMod']==0:
            for j, p in enumerate(Waves_WaveKin_Nodes):
                if p[2] < -WtrDpth or p[2] >0:
                    pass # all is zero
                else:
                    nodeInWater[:,j] = 1 # for all time steps
        else:
            raise NotImplementedError()
        #             CASE ( 1 )                 ! Vertical stretching.
        #                ! Vertical stretching says that the wave kinematics above the mean sea level
        #                !   equal the wave kinematics at the mean sea level.  The wave kinematics
        #                !   below the mean sea level are left unchanged:
        #                IF (   ( WaveKinzi(J) < -WtrDpth ) .OR. ( WaveKinzi(J) > WaveElev(I,J) ) ) THEN   ! .TRUE. if the elevation of the point defined by WaveKinzi(J) lies below the seabed or above the instantaneous wave elevation (exclusive)
        #                   WaveDynP   (I,J  )  = 0.0
        #                   WaveVel    (I,J,:)  = 0.0
        #                   WaveAcc    (I,J,:)  = 0.0
        #                   nodeInWater(I,J  )  = 0
        #                ELSE 
        #                   nodeInWater(I,J  )  = 1
        #                   IF   ( WaveKinzi(J) >= 0.0_ReKi ) THEN
        #                      ! Set the wave kinematics to the kinematics at mean sea level for locations above MSL, but below the wave elevation.
        #                      WaveDynP   (I,J  )  = WaveDynP0  (I,J  )
        #                      WaveVel    (I,J,:)  = WaveVel0   (I,J,:)
        #                      WaveAcc    (I,J,:)  = WaveAcc0   (I,J,:)
        #                   END IF
        #                   ! Otherwise, do nothing because the kinematics have already be set correctly via the various Waves modules
        #                END IF
        #             CASE ( 2 )                 ! Extrapolation stretching.
        #             ! Extrapolation stretching uses a linear Taylor expansion of the wave
        #             !   kinematics (and their partial derivatives with respect to z) at the mean
        #             !   sea level to find the wave kinematics above the mean sea level.  The
        #             !   wave kinematics below the mean sea level are left unchanged:
        #                IF (   ( WaveKinzi(J) < -WtrDpth ) .OR. ( WaveKinzi(J) > WaveElev(I,J) ) ) THEN   ! .TRUE. if the elevation of the point defined by WaveKinzi(J) lies below the seabed or above the instantaneous wave elevation (exclusive)
        #                   WaveDynP   (I,J  )  = 0.0
        #                   WaveVel    (I,J,:)  = 0.0
        #                   WaveAcc    (I,J,:)  = 0.0
        #                   nodeInWater(I,J  )  = 0
        #                ELSE 
        #                   nodeInWater(I,J  )  = 1
        #                   wavekinzloc = WaveKinzi(J)
        #                   WavePVel0loc = WavePVel0   (I,J,1)
        #                   IF   ( WaveKinzi(J) >= 0.0_ReKi ) THEN
        #                      ! Set the wave kinematics to the kinematics at mean sea level for locations above MSL, but below the wave elevation.
        #                      WaveDynP   (I,J  )  = WaveDynP0  (I,J  ) + WaveKinzi(J)*WavePDynP0  (I,J  )
        #                      WaveVel    (I,J,:)  = WaveVel0   (I,J,:) + WaveKinzi(J)*WavePVel0   (I,J,:)
        #                      WaveAcc    (I,J,:)  = WaveAcc0   (I,J,:) + WaveKinzi(J)*WavePAcc0   (I,J,:)
        #                   END IF
        #                   ! Otherwise, do nothing because the kinematics have already be set correctly via the various Waves modules
        #    ! Set the ending timestep to the same as the first timestep
        WaveDynP[NStepWave-1,:  ]  = WaveDynP [0,:  ]
        WaveVel [NStepWave-1,:,:]  = WaveVel  [0,:,:]
        WaveAcc [NStepWave-1,:,:]  = WaveAcc  [0,:,:]

        # --- Morison Init
        initData = {}
        initData['nodeInWater'] = nodeInWater
        initData['WaveVel']     = WaveVel
        initData['WaveAcc']     = WaveAcc
        initData['WaveDynP']    = WaveDynP
        initData['WaveTime']    = self.Waves.p['WaveTime']
        initData['Gravity']     = self.p['Gravity']
        initData['WtrDens']     = self.p['WtrDens']
        #print('WaveTime',self.Waves.p['WaveTime'])
        uMor, yMor = self.morison.init(initData)
        u['Morison'] = uMor
        y['Morison'] = yMor
        self.u = u 
        self.y = y 
        return u, y

    
    def calcOutput(self, t, x=None, xd=None, xo=None, u=None, y=None, optsM=None):
        yMor = self.morison.calcOutput(t, x=x, xd=xd, xo=xo, u=u['Morison'], y=y['Morison'], opts=optsM)
        y['Morison'] = yMor

        return y



    def linearize_RigidMotion2Loads(self, q0=None, qd0=None, qdd0=None, dq=None, dqd=None, dqdd=None, 
            RefPointMotion=None, RefPointMapping=None,
            moveWithOP=True,
            optsM=None, saveFile=None,
            around=None
            ):
        """ 
        Linearize the outputs force at the reference point with respect to
           motions of the reference point (assuming rigid body motion of the structure)
           Return: M=df/dqdd C=df/dqd K=df/dq 

        -q0, qd0, qdd0: 6-array of rigid DOFs positions, velocities and accelerations at platform ref 
        -dq, dqd, dqdd: 6-array of perturbations around q0, qd0 and qdd0

        - RefPoint: point used to defined rigid body motion, important parameter

        -optsM: options for Morison Calculation

        NOTE: only Morison for now
        """
        from welib.system.linearization import numerical_jacobian
        print('Computing HydroDyn linearized model, about: ({:.3f} {:.3f} {:.3f} {:.3f} {:.3f} {:.3f})'.format(*q0))
        uMesh = self.u['Morison']['Mesh']
        if RefPointMotion is None:
            RefPointMotion = uMesh.RefPoint

        if RefPointMapping is None:
            # TODO add options to map to a point in body coordinates (including rotations)
            RefPointMapping=np.array([0,0,0])
        print('- RefPoint for Motion : ',RefPointMotion)
        print('- RefPoint for Mapping: ',RefPointMapping)

        # --- Define a function that returns outputs
        def fh(q,qd,qdd,p=None):
            # Rigid body motion of the mesh
            self.u['Morison']['Mesh'].rigidBodyMotion(q=q, qd=qd, qdd=qdd, RefPoint=RefPointMotion)
            # Calculate hydrodynamic loads at every nodes 
            self.calcOutput(t=0, u=self.u, y=self.y, optsM=p)
            # Compute integral loads (force&moment) at the reference point (translated but not rotated)
            fh=np.zeros(6)
            if moveWithOP:
                MappingPoint = RefPointMapping + np.array([q[0],q[1],q[2]])
            else:
                MappingPoint = RefPointMapping
            fh[:3], fh[3:] = self.y['Morison']['Mesh'].mapLoadsToPoint(MappingPoint)
            return fh

        # --- Operating point and perturbation sizes
        if q0 is None:
            q0  = np.zeros(6)
        if qd0 is None:
            qd0  = np.zeros(6)
        if qdd0 is None:
            qdd0 = np.zeros(6)
        if dq is None:
            dq   = [0.01]*3 + [0.01]*3
        if dqd is None:
            dqd  = [0.01]*3 + [0.01]*3
        if dqdd is None:
            dqdd = [0.1]*3  + [0.1]*3

        # --- Linearization
        f0 = fh(q0,qd0,qdd0, optsM)
        Kall = -numerical_jacobian(fh, (q0,qd0,qdd0), 0, dq  , optsM)
        Call = -numerical_jacobian(fh, (q0,qd0,qdd0), 1, dqd , optsM)
        Mall = -numerical_jacobian(fh, (q0,qd0,qdd0), 2, dqdd, optsM)

        if around is not None:
            Kall=np.around(Kall,around)
            Call=np.around(Call,around)
            Mall=np.around(Mall,around)
            f0  =np.around(f0,around)


        return Mall, Call, Kall, f0

   
    @property
    def graph(self):
        import copy
        if self._graph is None:
            self._graph = self.File.toGraph()
        return copy.deepcopy(self._graph)


    def elementDivisions(self, e):
        n1, n2 = e.nodes
#         if e.data['Pot'] is False:
#             numDiv = np.ceil(e.length/e.data['DivSize']).astype(int)
#             dl = e.length/numDiv
#             SubNodesPositions = np.zeros((numDiv-1,3))
#             for j in range(numDiv-1):
#                 s = (j+1)/numDiv
#                 SubNodesPositions[j,:] = n1.point * (1-s) + n2.point * s
#                 nodeCount+=1
#             Positions = np.vstack( (n1.point, SubNodesPositions, n2.point ) )
# 
#             N=numDiv
#         else:
#             Positions=np.vstack((n1.point, n2.point))
#             dl=e.length
#             N=1
#         member['MGdensity']=np.zeros(N+1) # TODO
#         member['tMG'] = np.zeros(N+1)   # TODO
#         prop1 = e.nodeProps[0]  # NOTE: t&D are not stored in nodes since they are member dependent
#         prop2 = e.nodeProps[1] 
#         t             = np.linspace(prop1.data['t']  , prop2.data['t']  , N+1)
#         member['R']   = np.linspace(prop1.data['D']/2, prop2.data['D']/2, N+1)
#         member['RMG'] = member['R']+member['tMG']
#         member['Rin'] = member['R']-t


    def memberBuoyancyForce(self, e) :
        if e.data['Pot']:
            return None
        # Member is NOT modeled with Potential Flow Theory
        # --------------------------------------------------------------------------------}
        # ---Buoyancy loads
        # --------------------------------------------------------------------------------{
        # sides: Sections 3.1 and 3.2 ------------------------
        if z1 < 0:  # if segment is at least partially submerged ...
            if z1*z2 <= 0: # special calculation if the slice is partially submerged
                # Check that this is not the 1st element of the member
                if i==0:
                    raise Exception('The lowest element of a Morison member has become partially submerged!  This is not allowed.  Please review your model and create a discretization such that even with displacements, the lowest element of a member does not become partially submerged.')
                h0 = -z1/cosPhi             # distances along element centerline from point 1 to the waterplane
                if abs(dRdl_mg)< 0.0001:      # untapered cylinder case
                    Vs =    np.pi*r1*r1*h0   # volume of total submerged portion
                    if Vs ==0:
                        cx = 0.0  # Avoid singularity, but continue to provide the correct solution
                    else:
                        cr = 0.25*r1*r1*tanPhi/h0
                        cl = 0.5*h0 + 0.125*r1*r1*tanPhi*tanPhi/h0
                        cx = cr*cosPhi + cl*sinPhi
                else: # inclined tapered cylinder case (note I've renamed r0 to rh here##)
                    # NOTE: a0 and b0 always appear as a0b0, never separately.
                    rh   = r1 + h0*dRdl_mg    # radius of element at point where its centerline crosses the waterplane
                    C_1  = 1.0 - dRdl_mg**2 * tanPhi**2
                    # waterplane ellipse shape
                    b0   = rh/np.sqrt(C_1)
                    a0   = rh/((C_1)*cosPhi)             # simplified from what's in ConicalCalcs.ipynb
                    a0b0 = a0*b0
                    C_2  = a0b0*rh*cosPhi - r1**3
                    cl   = -(-0.75*a0b0*rh**2*cosPhi + 0.75*r1**4*C_1 + r1*C_1*C_2) / (dRdl_mg*C_1*C_2)
                    cr   = (0.75*a0b0*dRdl_mg*rh**2*sinPhi)/(C_1*C_2)
                    cx   = cr*cosPhi + cl*sinPhi 
                    Vs   = np.pi*(a0b0*rh*cosPhi - r1**3)/(3.0*dRdl_mg)       
                    # End per plan equations
                    #===================
                pwr = 3
                alpha    = (1.0-mem['alpha'][i])*z1**pwr/(-mem['alpha'][i]*z2**pwr + (1.0-mem['alpha'][i])*z1**pwr)
                Fb  = Vs*p['WtrDens']*g       #buoyant force
                Fr  = -Fb*sinPhi     #radial component of buoyant force
                Fl  = Fb*cosPhi      #axial component of buoyant force
                Moment = -Fb*cx      #This was matt's code        #moment induced about the center of the cylinder's bottom face
                # calculate (imaginary) bottom plate forces/moment to subtract from displacement-based values
                Fl  = Fl  + p['WtrDens']*g*z1* np.pi *r1*r1        
                Moment  = Moment  + p['WtrDens']*g* sinPhi * np.pi/4.0*r1**4       
                # reduce taper-based moment to remove (not double count) radial force distribution to each node 
                Moment  = Moment + Fr*(1.0-alpha)*dl
                F_B1, F_B2 = DistributeElementLoads(Fl, Fr, Moment, sinPhi, cosPhi, sinBeta, cosBeta, alpha)
                #print('Case 1')
                #print('>>> FB_1',F_B1[:3])
                #print('>>> FB_2',F_B2[:3])
                memLoads['F_B'][:, i]   += F_B1  # alpha
                memLoads['F_B'][:, i-1] += F_B2  # 1-alpha
            else: # normal, fully submerged case
                Fl = -2.0*np.pi*dRdl_mg*p['WtrDens']*g*dl*( z1*r1 + 0.5*(z1*dRdl_mg + r1*cosPhi)*dl + 1.0/3.0*(dRdl_mg*cosPhi*dl*dl) )   # from CylinderCalculationsR1.ipynb
                Fr = -np.pi*p['WtrDens']*g*dl*(r1*r1 + dRdl_mg*r1*dl + (dRdl_mg**2*dl**2)/3.0)*sinPhi                          # from CylinderCalculationsR1.ipynb
                Moment = -np.pi*dl*g*p['WtrDens']*(3.0*dl**3*dRdl_mg**4 + 3.0*dl**3*dRdl_mg**2 + 12.0*dl**2*dRdl_mg**3*r1 + 8.0*dl**2*dRdl_mg*r1 + 18.0*dl*dRdl_mg**2*r1*r1 + 6.0*dl*r1*r1 + 12.0*dRdl_mg*r1**3)*sinPhi/12.0   # latest from CylinderCalculationsR1.ipynb

                # precomputed as mem['alpha[i] ... alpha0 = (r1*r1 + 2*r1*r2 + 3*r2**2)/4/(r1*r1 + r1*r2 + r2**2)
                z1d = -min(0.0,z1)
                z2d = -min(0.0,z2)
                pwr = 3
                alpha = mem['alpha'][i]*z2d**pwr/(mem['alpha'][i]*z2d**pwr+(1-mem['alpha'][i])*z1d**pwr)
                # reduce moment to remove (not double count) radial force distribution to each node
                Moment = Moment - Fr*alpha*dl
                F_B1, F_B2 = DistributeElementLoads(Fl, Fr, Moment, sinPhi, cosPhi, sinBeta, cosBeta, alpha)
                #print('Case 2')
                #print('>>> FB_1',F_B1[:3])
                #print('>>> FB_2',F_B2[:3])
                memLoads['F_B'][:,i+1] += F_B1  # alpha
                memLoads['F_B'][:, i]  += F_B2  # 1-alpha



    def memberVolumeSubmerged(self, e, useDiv=False):
        from welib.hydro.tools import tapered_cylinder_geom
        prop1 = e.nodeProps[0]  # NOTE: t&D are not stored in nodes since they are member dependent
        prop2 = e.nodeProps[1] 
        Za = e.nodes[0].point[2]
        Zb = e.nodes[1].point[2]
        if Za>=0 and Zb>=0:
            return 0   # Fully above water
        elif Za<0 and Zb<0: # Fully submerged 
            return self.memberVolumeStructure(e, useDiv=useDiv)
        elif Za < -self.p['WtrDpth'] or Zb < -self.p['WtrDpth']:
            raise NotImplementedError()
        # Partially submerged, interpolated to "0"
        if not useDiv:
            Z0  = np.array([Za, Zb])
            tMG0= np.array([0,0]) # TODO
            R0  = np.array([prop1.data['D']/2, prop2.data['D']/2])
            # Stopping at 0
            Z   = np.array([np.min(Z0), 0])
            tMG = np.interp(Z , Z0, tMG0)
            R   = np.interp(Z , Z0, R0)
            RMG = R + tMG
            l=e.length * np.abs(Z[1]-Z[0])/ np.abs(Z0[1]-Z0[0])
            # get V and CV for marine growth displacement
            Vouter, cVouter = tapered_cylinder_geom(RMG[0], RMG[1], l)
        else:
            raise NotImplementedError()
        return Vouter

    def memberVolumeStructure(self, e, useDiv=False):
        from welib.hydro.tools import tapered_cylinder_geom
        prop1 = e.nodeProps[0]  # NOTE: t&D are not stored in nodes since they are member dependent
        prop2 = e.nodeProps[1] 
        Za = e.nodes[0].point[2]
        Zb = e.nodes[1].point[2]
        if Za < -self.p['WtrDpth'] or Zb < -self.p['WtrDpth']:
            raise NotImplementedError()
        if not useDiv:
            tMG = np.array([0,0]) # TODO
            R   = np.array([prop1.data['D']/2, prop2.data['D']/2])
            RMG = R + tMG
            # get V and CV for marine growth displacement
            Vouter, cVouter = tapered_cylinder_geom(RMG[0], RMG[1], e.length)
        else:
            raise NotImplementedError()
        return Vouter



    def VolumeStructure(self, method='NoDiv'):
        if method=='Morison':
            return self.morison.VolumeStructure
        return np.sum([self.memberVolumeStructure(e, useDiv=method=='Div') for e in self.graph.Elements])

    def VolumeSubmerged(self, method='NoDiv'):
        if method=='Morison':
            return self.morison.VolumeSubmerged
        return np.sum([self.memberVolumeSubmerged(e, useDiv=method=='Div') for e in self.graph.Elements])


    # --------------------------------------------------------------------------------}
    # --- IO/Converters
    # --------------------------------------------------------------------------------{
    def writeSummary(self, filename):
        with open(filename, 'w') as fid:
            self.morison.writeSummary(fid=fid)

    def toYAML(self, filename):
        if self._FEM is None:
            raise Exception('Call `initFEM()` before calling `toYAML`')
        subdyntoYAMLSum(self._FEM, filename, more = self.File['OutAll'])


    def toYAMSData(self, shapes=[0,4], main_axis='z'):
        """ 
        Convert to Data needed by YAMS
        """
        from welib.mesh.gradient import gradient_regular
        p=dict()
        return p


if __name__ == '__main__':
    import sys
    if len(sys.argv)>=1:
        filename=sys.argv[1]
    else:
        #filename='../../data/SparNoRNA/SparNoRNA_HD_RefH.dat'
        filename='_SparNoRNA_HD_RefH.dat'
        filename='_HD_T.dat'
        filename='_HD_T2.dat'
    import welib


    import  welib.weio

    driverfilename=filename.replace('.dat','.dvr')
    dvr=welib.weio.FASTInputFile(driverfilename)
    Gravity = dvr['Gravity']
    WtrDens = dvr['WtrDens']
    WtrDpth = dvr['WtrDpth']


#     hd = welib.weio.FASTInputFile(filename)
# #     hd.write('Out.dat')
#     graph = hd.toGraph()
#     print(graph)

    hd = HydroDyn(filename)
    hd.init(Gravity=Gravity, WtrDens=WtrDens, WtrDpth=WtrDpth)
#     hd.MorisonPositions
    hd.writeSummary(filename.replace('.dat','.HD_python.sum'))
#     print(hd.graph)

    EM = hd.morison.graph.Elements[13]
    E  = hd.graph.Elements[13]
    EM = hd.morison.graph.Elements[0]
    E  = hd.graph.Elements[0]
    print(E)
    print(EM)
    print(EM.MorisonData.keys())
    print()
    print(EM.MorisonData['R'])
    print(EM.MorisonData['Vouter']    , hd.memberVolumeStructure(E, useDiv=False))
    print(EM.MorisonData['Vsubmerged'], hd.memberVolumeSubmerged(E, useDiv=False))

    for e,em in zip(hd.graph.Elements, hd.morison.graph.Elements):
        if abs(em.MorisonData['Vouter']-hd.memberVolumeStructure(e))>1e-8:
           print('ID',e.ID,'V',em.MorisonData['Vouter'], hd.memberVolumeStructure(e)  )
    for e,em in zip(hd.graph.Elements, hd.morison.graph.Elements):
        if abs(em.MorisonData['Vsubmerged']-hd.memberVolumeSubmerged(e) )>1e-8:
            print('ID',e.ID,'V',em.MorisonData['Vsubmerged'], hd.memberVolumeSubmerged(e),   )

#     graph.divideElements(3)
#     print(graph)
#     import numpy as np
#     import matplotlib.pyplot as plt
#     from matplotlib import collections  as mc
#     from mpl_toolkits.mplot3d import Axes3D
#     fig = plt.figure()
#     ax = fig.add_subplot(1,2,1,projection='3d')
#     lines=graph.toLines(output='coord')
#     for l in lines:
#     #     ax.add_line(l)
#         ax.plot(l[:,0],l[:,1],l[:,2])

#     ax.autoscale()
    # ax.set_xlim([-40,40])
    # ax.set_ylim([-40,40])
    # ax.set_zlim([-40,40])
    # ax.margins(0.1)
#     plt.show()

# 
# if __name__ == '__main__':
#     pass
