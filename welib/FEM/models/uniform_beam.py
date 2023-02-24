import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# Local 
from welib.FEM.fem_elements import *
from welib.FEM.fem_model import *
from welib.FEM.utils import xNodesInputToArray



class FEMUniformBeam3d(FEMModel):
    def __init__(self, xNodes, E, G, rho, A, Ixx=None, Iyy=None, Izz=None, Kt=None, 
            elementType='frame3d', nel=None, main_axis='x', gravity=None,
            BC='clamped-free', M_root=None, M_tip=None, K_root=None, K_tip=None # Boundary conditions
            ):
        """ 
        INPUTS:
          - xNodes: define beam length, beam spanwise positions or beam nodes, either:
                -  (scalar) Beam length, for uniform beam [m]
                -  (1xn) Span vector of the beam (for straight beams) [m]
                -  (2xn) Nodes positions x,z along the beam for 2d beam [m]
                -  (3xn) Nodes positions x,y,z along the beam for 3d beam [m]
          - E   : (scalar) Elastic (Young) modulus
          - G   : (scalar) Shear modulus. Approximation: G = E/2/(1+0.3) [N/m^2]
          - rho  : Density [kg/m^3]
          - A    : Beam cross section area along the beam, at nodes [m^2]
          - Ixx  : Second Moment of Area of cross section, at nodes [m^4]
          - Iyy  : Second Moment of Area of cross section, at nodes [m^4]
          - Izz  : Second Moment of Area of cross section, at nodes [m^4]
          - Kt  : (n) Torsion constant, at nodes [m^4]
          - elementType: specify the element type to use along the beam: 
                 'timoshenko3d'
                 'frame3d'
                 'frame3dlin'
          - nel  : Number of elements. If provided Structural propeties and nodes will be interpolated to match nel. 
                   Otherwise, the length of xNodes determines the discretization
          - BC: string defining boundary condition:
                -'clamped-free': clamped at root, free at tip
                -'free-free': free at root, free at tip
          - M_root/tip: (6x6) additional rigid body mass matrix at beam ends
          - K_root/tip: (6x6) additional stiffness matrix at beam ends

          - main_axis: 'x' or 'z', defines the main direction of the beam (if xNodes is scalar)
                       and the main definition of the elements.
                       Mostly relevant for straight beams.
        """ 

        Nodes    = self.setNodes(xNodes, main_axis, nel)
        # --- Create Elements
        Elements=[]
        for i in np.arange(len(Nodes)-1):
            if elementType=='frame3d':
                elem = UniformFrame3dElement(ID=i, E=E, G=G, rho=rho, A=A, Ixx=Ixx, Iyy=Iyy, Izz=Izz, Kv=Kt, nodes=[Nodes[i], Nodes[i+1]], main_axis=main_axis)
            else:
                raise Exception()
            Elements.append(elem)
        # --- Initialize parent class
        FEMModel.__init__(self, Nodes=Nodes, Elements=Elements, mainElementType=elementType, gravity=gravity)

        self.init(BC=BC)

    def init(self, BC):
        # Assemble 
        self.assembly(gravity=None)
        self.applyInternalConstraints()

        # --- Boundary conditions
        if BC=='clamped-free':
            self.Nodes[0].setFixed()
            self.Nodes[-1].setInterfaceBC(6*[idDOF_Leader])
        else:
            raise NotImplementedError()

        self.partition()
        self.applyFixedBC() #, IFixed=None):

        self.eig(normQ='byMax')
        self.setModes(nModesFEM=30)

    def setNodes(self, xNodes, main_axis, nel):
        xNodes, xNodes0, s_span0, interp_needed = xNodesInputToArray(xNodes, main_axis=main_axis, nel=nel)

        nNodes =xNodes.shape[1]
        # --- Create Nodes
        Nodes=[]
        for i in np.arange(nNodes):
            Nodes.append(FEMNode(ID=i, x=xNodes[0,i], y=xNodes[1,i], z=xNodes[2,i])) # DOF=?
        return Nodes


class FEMUniformCylinder3d(FEMUniformBeam3d):
    def __init__(self, xNodes, E, G, rho, D, t, 
            elementType='frame3d', nel=None, main_axis='x', gravity=None,
            BC='clamped-free', M_root=None, M_tip=None, K_root=None, K_tip=None # Boundary conditions
            ):

        Nodes    = self.setNodes(xNodes, main_axis, nel)
        # --- Create Elements
        Elements=[]
        for i in np.arange(len(Nodes)-1):
            if elementType=='frame3d':
                elem = UniformFrame3dElement(ID=i, E=E, G=G, rho=rho, D=D, t=t, nodes=[Nodes[i], Nodes[i+1]], main_axis=main_axis)
                #elem = CylinderFrame3dElement(ID, E, G, rho, D, t=None, nodeIDs=None, nodes=None, main_axis='z', **kwargs):
            elif elementType=='euler-bernoulli3d':
                elem = CylinderTimoshenko3dElement(ID=i, E=E, G=G, rho=rho, D=D, t=t, nodes=[Nodes[i], Nodes[i+1]], main_axis=main_axis, shear=False)
            elif elementType=='timoshenko3d':
                elem = CylinderTimoshenko3dElement(ID=i, E=E, G=G, rho=rho, D=D, t=t, nodes=[Nodes[i], Nodes[i+1]], main_axis=main_axis, shear=True)
            else:
                raise Exception()
            Elements.append(elem)
        # --- Initialize grand parent! class...
        FEMModel.__init__(self, Nodes=Nodes, Elements=Elements, mainElementType=elementType, gravity=gravity)

        self.init(BC=BC)


