""" 
Define Finite element classes to easily compute stiffness and mass matrices.

An needs to implement two main methods:
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

# --------------------------------------------------------------------------------}
# --- Generic FEM Element
# --------------------------------------------------------------------------------{
class FEMElement(GraphElement):
    def __init__(self, ID, nodeIDs, nodes=None, propset=None, propIDs=None, properties=None, **kwargs):
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

# --------------------------------------------------------------------------------}
# --- Beam-like elements
# --------------------------------------------------------------------------------{
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
