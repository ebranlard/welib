# --- General
import unittest
import numpy as np
# --- Local
from welib.vortilib.elements.VortexCylinder import vc_tang_u, vcs_tang_u
from welib.wiz.WindTurbine import WindTurbine

class WindFarm(list):
    """ 
    Simple class to ease list of wind turbines manipulations
    """
    def __init__(self,WTs=[],name=''):
        self.name=name
        return super(WindFarm, self).__init__(WTs)

    def append(self,WT):
        if not isinstance(WT,WindTurbine):
            raise Exception('Can only append `WindTurbine` objects')
        return super(WindFarm, self).append(WT)

    def velocity_field(self, X, Y, Z, no_wake=False, U0=None, **kwargs):
        """ 
        Returns the induction field from all wind turbine in the wind farm. 
        If U0 is provided, this free stream is added to the induction field
        If no_wake is True, the induced velocity in the (self-)wakes is set to 0

        INPUTS:
            X, Y, Z: Control points in global coordinates where the flow is to be computed. 
            no_wake: boolean, if true: the induced velocity in the wake is set to 0. 
                     Typically set to true when combining with wake models.
            U0     : vector of size 3, in global coordinates, representing the free stream to be added
                     If None, no velocity field is added to the induction field
        """
        X=np.asarray(X)
        Y=np.asarray(Y)
        Z=np.asarray(Z)
        ux, uy, uz = np.zeros(X.shape), np.zeros(X.shape), np.zeros(X.shape)

        for WT in self:
            ux1,uy1,uz1 = WT.compute_u(X, Y, Z, no_wake=no_wake, only_ind=True, **kwargs)
            ux+=ux1
            uy+=uy1
            uz+=uz1

        if U0 is not None:
            ux+=U0[0]
            uy+=U0[1]
            uz+=U0[2]

        return ux, uy, uz


    @property
    def nWT(self):
        return len(self)

    @property
    def layout(self):
        return np.array([np.array(WT.r_hub) for WT in self]).reshape((self.nWT,3))

    @layout.setter
    def layout(self, layout):
        print('setter')
        layout = np.asarray(layout)
        if layout.shape != (self.nWT,3):
            raise Exception('Wrong layout shape. Layout should be {}'.format((self.nWT,3)))
        for iWT,WT in enumerate(self):
            r_hub=layout[iWT,:]
            print('r_hub',r_hub)
            WT.update_position(r_hub)

    def tostring(self,short=True):
        s ='class WindFarm({}):\n'.format(self.name)
        s+=' - Number of turbines: {:d} \n'.format(len(self))
        s+=' - layout: \n{}\n'.format(self.layout)
        return s

    def __repr__(self):
        return self.tostring(short=False)





def gridLayout(nx,xSpacing,nz,zSpacing=None,hub_height=0,mirror=False):
    """ Returned list of turbine positions on a grid layout
          y is vertical, positive upward
          x is lateral
          z is longi, positive along the wake
      If mirror is true, mirrored turbines are placed at y=-hubheight
    """
    if mirror:
        ny=2
    else:
        ny=1
    nWT = nz * nx * ny
    xWT = np.zeros((nWT))
    yWT = np.zeros((nWT))
    zWT = np.zeros((nWT))
    k   = 0
    for i in np.arange(1,nz+1).reshape(-1):
        for j in np.arange(1,nx+1).reshape(-1):
            k = k + 1
            xWT[k-1] = (j - int(np.floor(nx / 2)) - 1) * xSpacing
            zWT[k-1] = (i - 1) * zSpacing
            yWT[k-1] = hub_height
    if mirror:
        for i in np.arange(1,nz+1).reshape(-1):
            for j in np.arange(1,nx+1).reshape(-1):
                k = k + 1
                xWT[k-1] =  (j - int(np.floor(nx / 2)) - 1) * xSpacing
                zWT[k-1] =  (i - 1) * zSpacing
                yWT[k-1] = - hub_height
    return xWT,yWT,zWT


def windfarm_gridlayout_CTconst(Xcp,Ycp,Zcp,R,CT,U0,nxWT,xWTSpacing,nzWT,zWTSpacing,hub_height=0,mirror=False): 

    # Wind farm layout
    xWT,yWT,zWT= gridLayout(nxWT,xWTSpacing,nzWT,zWTSpacing,hub_height=hub_height,mirror=mirror)

    # Approximate a-CT-gamma relation
    a = 0.5 * (1 - np.sqrt(1 - CT))
    gamma_t = - 2 * U0 * a

    # All turbine are identical with constant CT/gamma 
    nWT      = len(xWT)
    nr       = 1
    vR       = np.zeros((nWT,nr)) + R
    vgamma_t = np.zeros((nWT,nr)) + gamma_t

    return vcs_tang_u(Xcp,Ycp,Zcp,vgamma_t,vR,xWT,yWT,zWT,epsilon=0)








# --------------------------------------------------------------------------------}
# --- TEST 
# --------------------------------------------------------------------------------{
class TestWindFarm(unittest.TestCase):

    def test_WF_layout(self):
        # test funciton that generates a grid layout, with possible mirror
        x,y,z=gridLayout(nx=2,xSpacing=2,nz=2,zSpacing=3,hub_height=2,mirror=True)
        np.testing.assert_almost_equal(x,[-2, 0,-2, 0, -2,  0, -2,  0])
        np.testing.assert_almost_equal(y,[ 2, 2, 2, 2, -2, -2, -2, -2])
        np.testing.assert_almost_equal(z,[ 0, 0, 3, 3,  0,  0,  3,  3])

        #from mpl_toolkits.mplot3d import Axes3D
        #import matplotlib.pyplot as plt
        #fig=plt.figure()
        #ax= fig.add_subplot(111,projection='3d')
        #ax.plot(z,x,y,'+')
        #ax.set_xlabel('z')
        #ax.set_ylabel('x')
        #ax.set_zlabel('y')
        #plt.show()

if __name__ == "__main__":
    unittest.main()
