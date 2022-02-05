import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# Local 
from welib.tools.clean_exceptions import *
from welib.weio import FASTInputFile
from welib.fast.hydrodyn_morison import Morison


class HydroDyn:

    def __init__(self, filename=None, hdData=None):

        self._graph=None
        self.File=None

        # Read SubDyn file
        if filename is not None:
            self.File = FASTInputFile(filename)
        elif hdData is not None:
            self.File = hdData

        # Internal
        self.p={}
        self.m={}
        self._graph=None

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
    def init(self, gravity = 9.81, WtrDens='1025', WtrDpth=0, MSL2SWL=0):
        """
        Initialize HydroDyn model 

        gravity: position of transition point
        """
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

        self.p['WtrDpth'] = WtrDpth + MSL2SWL

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
        self.p['NStepWave'] = np.ceil(f['WaveTMax']/f['WaveDT'] ).astype(int)

        # --- WvStretch_Init in HydroDyn.f90
        nodeInWater = np.zeros( (self.p['NStepWave'], len(Waves_WaveKin_Nodes)    ))
        WaveDynP    = np.zeros( (self.p['NStepWave'], len(Waves_WaveKin_Nodes)    ))
        WaveVel     = np.zeros( (self.p['NStepWave'], len(Waves_WaveKin_Nodes), 3 ))
        WaveAcc     = np.zeros( (self.p['NStepWave'], len(Waves_WaveKin_Nodes), 3 ))
        NStepWave = self.p['NStepWave']
        WtrDpth = self.p['WtrDpth']
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
        #                END IF
        #    ! Set the ending timestep to the same as the first timestep
        WaveDynP[NStepWave-1,:  ]  = WaveDynP [0,:  ]
        WaveVel [NStepWave-1,:,:]  = WaveVel  [0,:,:]
        WaveAcc [NStepWave-1,:,:]  = WaveAcc  [0,:,:]

        # --- Morison Init
        initData={'nodeInWater':nodeInWater}
        self.morison.init(initData)
   
    @property
    def graph(self):
        import copy
        if self._graph is None:
            self._graph = self.File.toGraph()
        return copy.deepcopy(self._graph)


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

#     hd = welib.weio.FASTInputFile(filename)
# #     hd.write('Out.dat')
#     graph = hd.toGraph()
#     print(graph)

    hd = HydroDyn(filename)
    hd.init()
#     hd.MorisonPositions
    hd.writeSummary(filename.replace('.dat','.HD_python.sum'))
#     print(hd.graph)


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
