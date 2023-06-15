""" 
Define Finite element classes to easily compute stiffness and mass matrices.

An element needs to implement two main methods:
    - elem.Ke(): returns the element stiffness matrix in global coordinates (unless specified)
    - elem.Me(): returns the element mass matrix in global coordinates (unless specified)

Element present in this package:

- FEMElement: template class
|
|-> SubDynBeam3dElement: generic element for SubDyn cylindrical elements

    |->SubDynBeam3dElement: generic element for SubDyn cylindrical elements

    |->SubDynFrame3dElement: SubDyn Euler-Benoulli beam elements (untapered)
    |      defined using: rho, D, t, E, G at both extremeties, but at the element level,
    |      not at node level! That's a SubDyn perk. Two values need to be specified 
    |      for each element

    |->SubDynTimoshenko3dElement: SubDyn Timoshenko beam elements (untapered)
    |      idem 

    |->SubDynCable3dElement: SubDyn cable elements



"""
import numpy as np


from welib.FEM.utils import DCM
# TODO get rid of these:
from welib.FEM.graph import Node as GraphNode
from welib.FEM.graph import Element as GraphElement


# Following the convention of SubDyn
idJointCantilever = 1
idJointUniversal  = 2
idJointPin        = 3
idJointBall       = 4
# NOTE: TypeID is mostly used for rigid links constraints and for graph visualization..
idMemberBeam      = 1
idMemberCable     = 2
idMemberRigid     = 3

idDOF_Fixed    =  0 # Fixed DOF BC
idDOF_Internal = 10 # Free/Internal DOF 
idDOF_Leader   = 20 # Leader DOF


# --------------------------------------------------------------------------------}
# --- NODE: TODO use it. TODO convert to simple dict for ease of usage/portability 
# --------------------------------------------------------------------------------{
class FEMNode(GraphNode):
    def __init__(self, ID, x, y, z=0, Type=idJointCantilever, DOFs=None, **kwargs):
        kwargs['Type'] = Type
        kwargs['DOFs'] = [] if DOFs is None else DOFs
        super(FEMNode,self).__init__(ID, x, y, z, **kwargs)

    def __repr__(self):
        s='<FNode{:4d}> x:{:7.2f} y:{:7.2f} z:{:7.2f}, {:}'.format(self.ID, self.x, self.y, self.z, self.data)
        return s

    def setInterfaceBC(self, I):
        """ 
        I is an array of the size DOFs of values idDOF_Fixed, idDOF_Internal, idDOF_Leader
        """
        self.data['IBC'] = I
        if 'DOFs' in self.data:
            if len(I)!=len(self.data['DOFs']):
                raise Exception('Interface BC should be of the same length as DOFs')

    def setReactionBC(self, I):
        """ 
        I is an array of the size DOFs of values idDOF_Fixed, idDOF_Internal, idDOF_Leader
        """
        self.data['RBC'] = I
        if 'DOFs' in self.data:
            if len(I)!=len(self.data['DOFs']):
                raise Exception('Interface BC should be of the same length as DOFs')

    def setFixed(self):
        if 'DOFs' in self.data:
            self.data['RBC'] = [idDOF_Fixed]*len(self.data['DOFs']) # TODO
        else:
            self.data['RBC'] = idDOF_Fixed # TODO


# --------------------------------------------------------------------------------}
# --- Generic FEM Elements
# --------------------------------------------------------------------------------{
class FEMElement(GraphElement):
    def __init__(self, ID, nodeIDs=None, nodes=None, propset=None, propIDs=None, properties=None, **kwargs):
        super(FEMElement, self).__init__(ID, nodeIDs, nodes, propset, propIDs, properties, **kwargs)
        self.data['Type'] = 'Generic'
        self.data['TypeID'] = None

    def Ke(self):
        raise NotImplementedError()

    def Me(self):
        raise NotImplementedError()

    def __repr__(self):
        s='<FElem{:4d}> nodeIDs: {} {}'.format(self.ID, self.nodeIDs, self.data)
        if self.propIDs is not None:
            s+=' {'+'propIDs:{} propset:{}'.format(self.propIDs, self.propset)+'}'
        if self.nodes is not None:
            s+=' l={:.2f}'.format(self.length)
        return s

class FEM2NodesElement(FEMElement):
    """ 
    Generic FEM element with two nodes. 
    """
    def __init__(self, ID, main_axis, nodeIDs=None, nodes=None, **kwargs):
        super(FEM2NodesElement, self).__init__(ID, nodeIDs=nodeIDs, nodes=nodes, **kwargs)
        self.data['Type'] = 'Generic2NodesElement'
        self.data['axis'] = main_axis

    @property
    def DCM(self):
        n1, n2 = self.nodes[0], self.nodes[1]
        return DCM((n1.x,n1.y,n1.z),(n2.x, n2.y, n2.z), main_axis=self.data['axis'])

class FEMIsotropic2NodesElement(FEM2NodesElement):
    """ 
    Generic FEM element with two nodes and isotropic material properties E, G, rho
    """
    def __init__(self, ID, E, rho, main_axis, G=None, nodeIDs=None, nodes=None, **kwargs):
        """ 
        INPUTS:
         - E : Young's (elastic) modulus [N/m2]
         - G : Shear modulus. For an isotropic material G = E/2(nu+1) with nu the Poission's ratio [N/m^2]
         - rho: Material density [kg/m^3] 
         """
        super(FEMIsotropic2NodesElement, self).__init__(ID, main_axis=main_axis, nodeIDs=nodeIDs, nodes=nodes, **kwargs)
        self.data['Type'] = 'GenericIsotropic2NodesElement'
        self.data['E']   = E
        self.data['G']   = G
        self.data['rho'] = rho
    @property
    def rho(self): return self.data['rho']
    @property
    def E(self): return self.data['E']
    @property
    def G(self): return self.data['G']

# --------------------------------------------------------------------------------}
# --- 2D Beam-like elements 
# --------------------------------------------------------------------------------{
class Frame2dElement(FEMIsotropic2NodesElement):
    def __init__(self, ID, E, rho, I, A, nodeIDs=None, nodes=None, main_axis='z', **kwargs):
        """ 
        Frame 2D element
        Nodal DOF   : (ux uy theta)

        INPUTS:
         - E  : Young's (elastic) modulus [N/m2]
         - rho: Material density [kg/m^3] 
         - I  : Planar second moment of area, 
               If main_axis='x':   Iy=\iint z^2 dy dz [m4] 
               If main_axis='z':   Iy=\iint x^2 dx dy [m4] 
        """
        super(Frame2dElement,self).__init__(ID, E=E, rho=rho, G=None, main_axis=main_axis, nodeIDs=nodeIDs, nodes=nodes, **kwargs)
        self.data['Type'] = 'Frame2dElement'
        self.data['I'] = I
        self.data['A'] = A

    @property
    def mass(self): return self.data['rho'] * self.data['A'] * self.length # [kg]

    def Ke(e, local=False):
        from .frame2d import frame2d_KeMe
        R = None if local else e.DCM
        main_axis = e.data['axis']
        EA = e.data['E'] * e.data['A']
        EI = e.data['E'] * e.data['I']
        # TODO sort out theta, DCM, main_axis
        return frame2d_KeMe(EA, EI, e.length, e.mass, T=0, theta=None, MMFormulation='consistent')[0]
        #return timoshenko_Ke(e.length, e.area, I[0], I[1], I[2], e.kappa,  e.E, e.G,  main_axis=main_axis) # NOTE: shear False for 

    def Me(e, local=False):
        from .frame2d import frame2d_KeMe
        R = None if local else e.DCM
        main_axis = e.data['axis']
        EA = e.data['E'] * e.data['A']
        EI = e.data['E'] * e.data['I']
        # TODO sort out theta, DCM, main_axis
        return frame2d_KeMe(EA, EI, e.length, e.mass, T=0, theta=None, MMFormulation='consistent')[1]


class Beam2dElement(FEMIsotropic2NodesElement):
    def __init__(self, ID, E, rho, I, A, nodeIDs=None, nodes=None, main_axis='z', **kwargs):
        """ 
        Beam 2D element
        Nodal DOF:    (u theta)

        INPUTS:
         - E  : Young's (elastic) modulus [N/m2]
         - rho: Material density [kg/m^3] 
         - I  : Planar second moment of area, 
               If main_axis='x':   Iy=\iint z^2 dy dz [m4] 
               If main_axis='z':   Iy=\iint x^2 dx dy [m4] 
        """
        super(Frame2dElement,self).__init__(ID, E=E, rho=rho, G=None, main_axis=main_axis, nodeIDs=nodeIDs, nodes=nodes, **kwargs)
        self.data['Type'] = 'Beam2dElement'
        self.data['I'] = I
        self.data['A'] = A

    @property
    def mass(self): return self.data['rho'] * self.data['A'] * self.length # [kg]

    def Ke(e, local=False):
        from .beam2d import beam2d_KeMe
        R = None if local else e.DCM
        main_axis = e.data['axis']
        EI = e.data['E'] * e.data['I']
        # TODO sort out theta, DCM, main_axis
        return frame2d_KeMe(EA, EI, e.length, e.mass, T=0, theta=None, MMFormulation='consistent')[0]
        #return timoshenko_Ke(e.length, e.area, I[0], I[1], I[2], e.kappa,  e.E, e.G,  main_axis=main_axis) # NOTE: shear False for 

    def Me(e, local=False):
        from .beam2d import beam2d_KeMe
        R = None if local else e.DCM
        main_axis = e.data['axis']
        EI = e.data['E'] * e.data['I']
        # TODO sort out theta, DCM, main_axis
        return beam2d_KeMe(EI, e.length, e.mass, T=0, theta=None, MMFormulation='consistent' )[1]


# --------------------------------------------------------------------------------}
# --- 3D Beam-like elements
# --------------------------------------------------------------------------------{
class UniformFrame3dElement(FEMIsotropic2NodesElement):
    def __init__(self, ID, E, G, rho, A, Ixx, Iyy, Izz, Kv, nodeIDs=None, nodes=None, main_axis='z', **kwargs):
        super(UniformFrame3dElement,self).__init__(ID, E=E, rho=rho, G=G, main_axis=main_axis, nodeIDs=nodeIDs, nodes=nodes, **kwargs)
        self.data['Type'] = 'UniformFrame3dElement'
        self.data['A']    = A
        self.data['EA']   = E*A
        self.data['Ixx']  = Ixx
        self.data['Iyy']  = Iyy
        self.data['Izz']  = Izz
        self.data['Kv']   = Kv

    @property
    def mass(self): return self.data['rho'] * self.data['A'] * self.length # [kg]
    @property
    def inertias(self): return (self.data['Ixx'], self.data['Iyy'], self.data['Izz'])

    def Fe_g(e, g,local=False):
        """ Element force due to gravity """
        from .timoshenko import timoshenko_Fe_g
        R = np.eye(3) if local else e.DCM
        return timoshenko_Fe_g(e.length, e.data['A'], e.rho, g, R=R, main_axis=e.data['axis'])

    def Fe_o(e, main_axis='z', local=False):
        """ Element force due to other sources (e.g. pretension cable) """
        return np.zeros(12)

    def Ke(e, local=False, T=0):
        from .frame3d import frame3d_KeMe
        R = None if local else e.DCM
        main_axis = e.data['axis']
        E   = e.data['E']
        EIx = E * e.data['Ixx']
        EIy = E * e.data['Iyy']
        EIz = E * e.data['Izz']
        A   = e.data['A']
        G   = e.data['G']
        Kv  = e.data['Kv']
        return frame3d_KeMe(E, G, Kv, E*A, EIx, EIy, EIz, e.length, A, e.mass, T=T, R=R, main_axis=main_axis) [0]

    def Me(e, local=False, T=0):
        from .frame3d import frame3d_KeMe
        R = None if local else e.DCM
        main_axis = e.data['axis']
        E   = e.data['E']
        EIx = E * e.data['Ixx']
        EIy = E * e.data['Iyy']
        EIz = E * e.data['Izz']
        A   = e.data['A']
        G   = e.data['G']
        Kv  = e.data['Kv']
        return frame3d_KeMe(E, G, Kv, E*A, EIx, EIy, EIz, e.length, A, e.mass, T=T, R=R, main_axis=main_axis) [1]

# TODO
# def frame3dlin_KeMe(E,G,Kv1,Kv2,A1,A2,Iy1,Iy2,Iz1,Iz2,L,me1,me2,R=None):
# def timoshenko_KeMe(L, A, Ixx, Iyy, Jzz, kappa, E, G, rho, R=None, main_axis='z') # <<< NOTE: for cylinder only for now


class CylinderBeam3dElement(FEMIsotropic2NodesElement):
    def __init__(self, ID, E, G, rho, D, t=None, nodeIDs=None, nodes=None, main_axis='z', **kwargs):
        """ 
        Untapered cylinder, constant properties.

        INPUTS:
         - D: outer diameter
         - t: wall thickness
        """
        super(CylinderBeam3dElement,self).__init__(ID, E=E, rho=rho, G=G, main_axis=main_axis, nodeIDs=nodeIDs, nodes=nodes, **kwargs)
        self.data['Type'] = 'CylinderBeam3dElement'
        self.data['TypeID'] = idMemberBeam
        self.data['shape'] = 'cylinder'
        self.data['D'] = D
        if t is None:
            t = D/2
        self.data['t'] = t
        self.data['axis'] = main_axis
        self.shear=True
    @property
    def r1(self): return self.data['D']/2 # Outer diameter
    @property
    def r2(self): return self.data['D']/2-self.data['t'] # Inner diameter
    @property
    def area(self): return np.pi*(self.r1**2- self.r2**2)
    @property
    def inertias(self):
        if self.data['axis']=='z':
            Ixx = 0.25*np.pi*(self.r1**4-self.r2**4)
            Iyy = Ixx
            Izz = 2.0*Ixx # Polar
        else:
            Izz = 0.25*np.pi*(self.r1**4-self.r2**4)
            Iyy = Izz
            Ixx = 2.0*Izz # Polar
        return Ixx, Iyy, Izz

    @property
    def D(self): return 2*self.r1

    @property
    def kappa(self): raise NotImplementedError()

    def Fe_g(e, g, main_axis='z', local=False):
        """ Element force due to gravity """
        from .timoshenko import timoshenko_Fe_g
        R = np.eye(3) if local else e.DCM
        return timoshenko_Fe_g(e.length, e.area, e.rho, g, R=R, main_axis=main_axis)

    def Fe_o(e, main_axis='z', local=False):
        """ Element force due to other sources (e.g. pretension cable) """
        return np.zeros(12)

    def Ke(e, local=False):
        from .timoshenko import timoshenko_Ke
        I = e.inertias
        R = None if local else e.DCM
        main_axis = e.data['axis']
        return timoshenko_Ke(e.length, e.area, I[0], I[1], I[2], e.kappa,  e.E, e.G, shear=e.shear, R=R, main_axis=main_axis) # NOTE: shear False for 

    def Me(e, local=False):
        from .timoshenko import timoshenko_Me
        I = e.inertias
        R = None if local else e.DCM
        main_axis = e.data['axis']
        return timoshenko_Me(e.length, e.area, I[0], I[1], I[2], e.rho, R=R, main_axis=main_axis)

class CylinderFrame3dElement(CylinderBeam3dElement):
    def __init__(self, ID, E, G, rho, D, t=None, nodeIDs=None, nodes=None, main_axis='z', **kwargs):
        """ 
        Untapered cylinder, constant properties. Frame3d
        """
        super(CylinderFrame3dElement,self).__init__(ID, E=E, G=G, rho=rho, D=D, t=t, main_axis=main_axis, nodeIDs=nodeIDs, nodes=nodes, **kwargs)
        self.data['Type'] = 'CylinderFrame3dElement'
        self.shear=False

    @property
    def kappa(self): return 0  # shear coefficients are zero for Euler-Bernoulli

class CylinderTimoshenko3dElement(CylinderBeam3dElement):
    def __init__(self, ID, E, G, rho, D, t=None, shear=True, nodeIDs=None, nodes=None, main_axis='z', propset=None, propIDs=None, properties=None, **kwargs):
        """ 
        Untapered cylinder, constant properties. Timoshenko3d
        """
        super(CylinderTimoshenko3dElement,self).__init__(ID, E=E, G=G, rho=rho, D=D, t=t, main_axis=main_axis, nodeIDs=nodeIDs, nodes=nodes, **kwargs)
        self.data['Type'] = 'CylinderTimoshenko3dElement'
        self.shear=shear

    @property
    def kappa(self): 
        r1 = self.data['D']/2
        t  = self.data['t']
        nu      = self.E/(2*self.G) - 1
        D_outer = 2 * r1              # average (outer) diameter
        D_inner = D_outer - 2*t       # remove 2x thickness to get inner diameter
        ratioSq = ( D_inner / D_outer)**2
        kappa =   ( 6 * (1 + nu) **2 * (1 + ratioSq)**2 )/( ( 1 + ratioSq )**2 * ( 7 + 14*nu + 8*nu**2 ) + 4 * ratioSq * ( 5 + 10*nu + 4 *nu**2 ) )
        return kappa

class SubDynBeam3dElement(FEMElement):
    def __init__(self, ID, nodeIDs, nodes=None, propset=None, propIDs=None, properties=None, **kwargs):
        super(SubDynBeam3dElement,self).__init__(ID, nodeIDs, nodes, propset, propIDs, properties, **kwargs)
        self.data['Type'] = 'SubDynBeam3d'
        self.data['TypeID'] = idMemberBeam
        self.data['shape'] = 'cylinder'

    def __repr__(self):
        s='<{:s}Elem{:4d}> nodeIDs: {} {}'.format(self.data['Type'], self.ID, self.nodeIDs, self.data)
        if self.propIDs is not None:
            s+=' {'+'propIDs:{} propset:{}'.format(self.propIDs, self.propset)+'}'
        if self.nodes is not None:
            s+=' l={:.2f}'.format(self.length)
        return s

    @property
    def DCM(self):
        n1, n2 = self.nodes[0], self.nodes[1]
        return DCM((n1.x,n1.y,n1.z),(n2.x, n2.y, n2.z), main_axis='z')

    @property
    def t(self):
        n1, n2 = self.nodeProps[0], self.nodeProps[1] # Using SubDyn convention
        t1, t2 = n1.data['t'], n2.data['t']
        return (t1+t2)/2

    @property
    def r1(self):
        n1, n2 = self.nodeProps[0], self.nodeProps[1] # Using SubDyn convention
        D1, D2 = n1.data['D'], n2.data['D']
        return (D1+D2)/4

    @property
    def r2(self):
        n1, n2 = self.nodeProps[0], self.nodeProps[1] # Using SubDyn convention
        D1, D2 = n1.data['D'], n2.data['D']
        t  = self.t
        if t==0:
            r2 = 0
        else:
            r2 = self.r1 - t
        return r2

    @property
    def area(self): return np.pi*(self.r1**2- self.r2**2)

    @property
    def inertias(self):
        Ixx = 0.25*np.pi*(self.r1**4-self.r2**4)
        Iyy = Ixx
        Jzz = 2.0*Ixx
        return Ixx, Iyy, Jzz

    @property
    def D(self): return 2*self.r1

    @property
    def rho(self): return self.nodeProps[0].data['rho'] # same convention as SubDyn returns density of first node
    @property
    def E(self): return self.nodeProps[0].data['E'] # same convention as SubDyn returns density of first node
    @property
    def G(self): return self.nodeProps[0].data['G'] # same convention as SubDyn returns density of first node

    @property
    def T0(self):
        return -9.990000E+36 # pretension, only for cable


    def Fe_g(e, g, main_axis='z', local=False):
        """ Element force due to gravity """
        from .timoshenko import timoshenko_Fe_g
        R = np.eye(3) if local else e.DCM
        return timoshenko_Fe_g(e.length, e.area, e.rho, g, R=R, main_axis=main_axis)

    def Fe_o(e, main_axis='z', local=False):
        """ Element force due to other sources (e.g. pretension cable) """
        return np.zeros(12)

# --------------------------------------------------------------------------------}
# --- Frame3D
# --------------------------------------------------------------------------------{
class SubDynFrame3dElement(SubDynBeam3dElement):
    def __init__(self, ID, nodeIDs, nodes=None, propset=None, propIDs=None, properties=None, **kwargs):
        super(SubDynFrame3dElement,self).__init__(ID, nodeIDs, nodes, propset, propIDs, properties, **kwargs)
        self.data['Type'] = 'SubDynFrame3d'
        self.data['TypeID'] = idMemberBeam
        self.data['shape'] = 'cylinder'

    @property
    def kappa(self): return 0  # shear coefficients are zero for Euler-Bernoulli

    def Ke(e, main_axis='z', local=False):
        from .frame3d import frame3d_KeMe # TODO TODO
        from .timoshenko import timoshenko_Ke
        I = e.inertias
        R = None if local else e.DCM
        return timoshenko_Ke(e.length, e.area, I[0], I[1], I[2], e.kappa,  e.E, e.G, shear=False, R=R, main_axis=main_axis) # NOTE: shear False for 

    def Me(e, main_axis='z', local=False):
        from .frame3d import frame3d_KeMe # TODO TODO
        from .timoshenko import timoshenko_Me
        I = e.inertias
        R = None if local else e.DCM
        return timoshenko_Me(e.length, e.area, I[0], I[1], I[2], e.rho, R=R, main_axis=main_axis)

# --------------------------------------------------------------------------------}
# --- Timoshenko 3D 
# --------------------------------------------------------------------------------{
class SubDynTimoshenko3dElement(SubDynBeam3dElement):
    def __init__(self, ID, nodeIDs, nodes=None, propset=None, propIDs=None, properties=None, **kwargs):
        super(SubDynTimoshenko3dElement,self).__init__(ID, nodeIDs, nodes, propset, propIDs, properties, **kwargs)
        self.data['Type'] = 'SubDynTimoshenko3d'
        self.data['TypeID'] = idMemberBeam
        self.data['shape'] = 'cylinder'

    @property
    def kappa(self): 
        n1, n2 = self.nodeProps[0], self.nodeProps[1] # Using SubDyn convention
        D1 = n1.data['D']
        D2 = n2.data['D']
        t1 = n1.data['t']
        t2 = n2.data['t']
        r1 = (D1+D2)/4
        t  = (t1+t2)/2
        nu      = self.E/(2*self.G) - 1
        D_outer = 2 * r1              # average (outer) diameter
        D_inner = D_outer - 2*t       # remove 2x thickness to get inner diameter
        ratioSq = ( D_inner / D_outer)**2
        kappa =   ( 6 * (1 + nu) **2 * (1 + ratioSq)**2 )/( ( 1 + ratioSq )**2 * ( 7 + 14*nu + 8*nu**2 ) + 4 * ratioSq * ( 5 + 10*nu + 4 *nu**2 ) )
        return kappa

    def Ke(e, main_axis='z', local=False):
        from .timoshenko import timoshenko_Ke
        I = e.inertias
        R = None if local else e.DCM
        return timoshenko_Ke(e.length, e.area, I[0], I[1], I[2], e.kappa,  e.E, e.G, shear=True, R=R, main_axis=main_axis)

    def Me(e, main_axis='z', local=False):
        from .timoshenko import timoshenko_Me
        I = e.inertias
        R = None if local else e.DCM
        return timoshenko_Me(e.length, e.area, I[0], I[1], I[2], e.rho, R=R, main_axis=main_axis)

# --------------------------------------------------------------------------------}
# --- Cable3D 
# --------------------------------------------------------------------------------{
class SubDynCable3dElement(SubDynBeam3dElement):
    def __init__(self, ID, nodeIDs, nodes=None, propset=None, propIDs=None, properties=None, **kwargs):
        super(SubDynCable3dElement,self).__init__(ID, nodeIDs, nodes, propset, propIDs, properties, **kwargs)
        self.data['Type'] = 'SubDynCable3d'
        self.data['TypeID'] = idMemberCable
        self.data['shape'] = 'cylinder'

    # For compatbility with other beam-like structures
    @property
    def area(self): return 1
    @property
    def G(self): return -9.99e+36   # for yaml only
    @property
    def inertias(self): return (-9.99e+36,-9.99e+36,-9.99e+36) # for yaml only
    @property
    def kappa(self): return -9.99e+36# for yaml only
    @property
    def E(self): return self.EA/self.area

    @property
    def EA(self): return self.nodeProps[0].data['EA'] # using first node, SubDyn convention

    @property
    def L0(e): 
        """ rest length for which pretension would be zero"""
        Eps0 = e.T0/(e.EA)
        return  e.length/(1+Eps0)

    @property
    def T0(self): 
        return self.nodeProps[0].data['T0'] # same convention as SubDyn returns density of first node

    def Ke(e, main_axis='z', local=False):
        from .cable import cable3d_Ke
        I = e.inertias
        R = None if local else e.DCM


        return cable3d_Ke(e.length, e.area, e.E, e.T0, R=R, main_axis=main_axis)

    def Me(e, main_axis='z', local=False):
        from .cable import cable3d_Me
        R = None if local else e.DCM
        return cable3d_Me(e.L0, e.area, e.rho, R=R, main_axis=main_axis) # NOTE: we use L0 for the mass

    def Fe_o(e, main_axis='z', local=False):
        from .cable import cable3d_Fe_T0
        R = None if local else e.DCM
        return cable3d_Fe_T0(e.T0, R=R, main_axis=main_axis)


# --------------------------------------------------------------------------------}
# --- Rigid
# --------------------------------------------------------------------------------{
class SubDynRigid3dElement(SubDynBeam3dElement):
    def __init__(self, ID, nodeIDs, nodes=None, propset=None, propIDs=None, properties=None, **kwargs):
        super(SubDynRigid3dElement,self).__init__(ID, nodeIDs, nodes, propset, propIDs, properties, **kwargs)
        self.data['Type'] = 'SubDynRigid3d'
        self.data['TypeID'] = idMemberRigid
        self.data['shape'] = 'cylinder'

    @property 
    def area(self): return 1  # arbitrary, but used for mass
    @property 
    def D(self): return min(np.sqrt(1/np.pi)*4, self.length*0.05)  # arbitrary, for plotting only
    @property
    def inertias(self): return (-9.99e+36,-9.99e+36,-9.99e+36) # for yaml only
    @property
    def E(self): return -9.99e+36   # for yaml only
    @property
    def G(self): return -9.99e+36   # for yaml only
    @property
    def kappa(self): return -9.99e+36   # for yaml only

    def Ke(e, main_axis='z', local=False):
      return np.zeros((12,12))


    def Me(e, main_axis='z', local=False):
        from .cable import cable3d_Me # NOTE: we use the same as cable
        if e.rho==0:
            return np.zeros((12,12))
        else:
            R = None if local else e.DCM
            return cable3d_Me(e.length, e.area, e.rho, R=R, main_axis=main_axis)
