import os
import numpy as np

from sympy.physics.mechanics import dynamicsymbols
from sympy.physics.mechanics.functions import find_dynamicsymbols
from sympy import diff, Matrix, trigsimp
from sympy import python

from .yams_sympy_tools import smallAngleApprox, cleantex, subs_no_diff , cleanPy, myjacobian
from .yams_kane import YAMSKanesMethod
from .yams_sympy import YAMSFlexibleBody
from welib.tools.tictoc import Timer
from collections import OrderedDict
# --------------------------------------------------------------------------------}
# ---  
# --------------------------------------------------------------------------------{
class YAMSModel(object):
    """ 
    Wrapper class to generate the equation of motions for an assembly of bodies
      - solve the equations 
      - perform small angle approximation
      - export to latex
      - TODO: generate python code
    
    """
    def __init__(self, name=''):
        self.name        = name
        self.ref         = None
        self.coordinates = None
        self.speeds      = None
        self.kdeqsSubs   = None
        self.bodies      = None
        self.kane        = None
        self.g_vect      = None
        self._sa_forcing     = None
        self._sa_mass_matrix = None
        self._sa_EOM         = None
        self.opts        = None
        self.body_loads  = None
        self.var         = [] # Independent variables
        self.smallAngleUsed=None

        self.M=None # Non-linear Mass matrix
        self.F=None # Non-linear Forcing
        self.M0=None # linear Mass matrix
        self.K0=None # linear Stiffness matrix
        self.C0=None # linear Damping matrix
        self.B0=None # linear Forcing

    def __repr__(self):
        s='<{} object "{}" with attributes:>\n'.format(type(self).__name__,self.name)
        s+=' - coordinates:       {}\n'.format(self.coordinates)
        s+=' - speeds:            {}\n'.format(self.speeds)
        s+=' - kdeqsSubs:         {}\n'.format(self.kdeqsSubs)
        s+=' - var:               {}\n'.format(self.var)
        s+=' - smallAngleUsed   : {}\n'.format(self.smallAngleUsed)
        s+=' - number of bodies : {}\n'.format(len(self.bodies))
        s+=' - opts             : {}\n'.format(self.opts)
        s+=' * loads            : {}\n'.format(self.loads)
        return s
    
    @property
    def kdeqs(self):
        return [s-v for (s,v) in self.kdeqsSubs]

    @property
    def loads(self):
        return [f[1] for f in self.body_loads]

    @property
    def EOM(self, Mform='symbolic'):
        if self.kane is None:
            self.kaneEquations(Mform=Mform)

        return (self.fr+self.frstar).subs(self.kdeqsSubs)

    @property
    def coordinates_speed(self):
        return [q.diff(dynamicsymbols._t) for q in self.coordinates]

    @property
    def coordinates_acc(self):
        return [qd.diff(dynamicsymbols._t) for qd in self.coordinates_speed]


    def kaneEquations(self, Mform='symbolic', addGravity=True):
        """ 
        Compute equation of motions using Kane's method
        Mform: form to use for mass matrix
        addGravity: include gravity for elastic bodies
        """
        for sa in ['ref', 'coordinates', 'speeds','kdeqs','bodies','loads']:
            if getattr(self,sa) is None:
                raise Exception('Attribute {} needs to be set before calling `kane` method'.format(sa))

        with Timer('Kane step1',True,silent=True):
            self.kane = YAMSKanesMethod(self.ref.frame, self.coordinates, self.speeds, self.kdeqs)

        # --- Expensive kane step
        with Timer('Kane step 2',True,silent=True):
            #(use  Mform ='symbolic' or 'TaylorExpanded'), Mform='symbolic'
            self.fr, self.frstar  = self.kane.kanes_equations(self.bodies, self.loads, Mform=Mform, addGravity=addGravity, g_vect=self.g_vect)
        self.kane.fr     = self.fr
        self.kane.frstar = self.frstar

    def smallAngleApprox(self, angle_list, extraSubs=None):
        """ 
        Apply small angle approximation to forcing and mass matrix
        """
        extraSubs = [] if extraSubs is None else extraSubs
        # Forcing
        with Timer('Small angle approx. forcing',True,silent=True):
            if self._sa_forcing is None:
                self._sa_forcing=self.kane.forcing
            self._sa_forcing = self._sa_forcing.subs(self.kdeqsSubs).subs(extraSubs)
            self._sa_forcing = smallAngleApprox(self._sa_forcing, angle_list)
        with Timer('Small angle approx. forcing simplify',True,silent=True):
            self._sa_forcing.simplify()
        # Mass matrix
        with Timer('Small angle approx. mass matrix',True,silent=True):
            if self._sa_mass_matrix is None:
                self._sa_mass_matrix=self.kane.mass_matrix
            self._sa_mass_matrix = self._sa_mass_matrix.subs(self.kdeqsSubs).subs(extraSubs)
            self._sa_mass_matrix = smallAngleApprox(self._sa_mass_matrix, angle_list)
        with Timer('Small angle approx. mass matrix simplify',True,silent=True):
            self._sa_mass_matrix.simplify()

        self.smallAngleUsed=angle_list
        

    def smallAngleApproxEOM(self, angle_list, extraSubs=None):
        """ 
        Apply small angle approximation to equation of motion H(x,xd,xdd,..)=0
        """
        extraSubs = [] if extraSubs is None else extraSubs
        EOM=self.EOM
        with Timer('Small angle approx. EOM',True,silent=True):
            EOM = EOM.subs(extraSubs)
            EOM = smallAngleApprox(EOM, angle_list).subs(extraSubs)
        #with Timer('Small angle approx. EOM simplify',True):
        #    EOM.simplify()
        self._sa_EOM = EOM
       
    def smallAngleLinearize(self, op_point=None, noAcc=True, noVel=False, extraSubs=None):
        """ 
        Linearize the equations with small angle approximations
        """
        op_point  = [] if op_point is None else op_point
        extraSubs = [] if extraSubs is None else extraSubs
        if self._sa_EOM is None:
            raise Exception('Run smallAngleApproxEOM first')
        M,C,K,B = self._linearize(op_point=op_point, EOM=self._sa_EOM, noAcc=noAcc, noVel=noVel, extraSubs=extraSubs)

        self._sa_M = M
        self._sa_C = C
        self._sa_K = K
        self._sa_B = B

    def linearize(self, op_point=None, noAcc=True, noVel=False, extraSubs=None):
        """ 
        Linearize the "non" linear equations
        NOTE: EOM has kdeqsSubs in it
        """
        op_point  = [] if op_point is None else op_point
        extraSubs = [] if extraSubs is None else extraSubs
        M,C,K,B = self._linearize(op_point=op_point, EOM=self.EOM, noAcc=noAcc, noVel=noVel, extraSubs=extraSubs)

        self.M0 = M
        self.C0 = C
        self.K0 = K
        self.B0 = B

    def _linearize(self, op_point=None, EOM=None, noAcc=True, noVel=False, extraSubs=None):
        op_point  = [] if op_point is None else op_point
        extraSubs = [] if extraSubs is None else extraSubs
        if EOM is None:
            EOM=self.EOM

        with Timer('Linearization',True,silent=True):
            # NOTE: order important
            op_point=[]
            if noAcc: 
                op_point+=[(qdd,0) for qdd in self.coordinates_acc]
                op_point+=[(diff(qd,dynamicsymbols._t),0) for qd in self.speeds]
            if noVel: 
                op_point+=[(qd,0) for qd in self.coordinates_speed]
                op_point+=[(qd,0) for qd in self.speeds]
        
            # ---  Using kane to get "r"...
            linearizer = self.kane.to_linearizer()
            self.var = self.kane._lin_r
            II = np.argsort([str(v) for v in self.var])
            self.var = list(np.array(self.var)[II])
            # KEEP ME Alternative
            #M, A, B = linearizer.linearize(op_point=op_point ) #o

            M =-EOM.jacobian(self.coordinates_acc  ).subs(extraSubs).subs(op_point)
            C =-EOM.jacobian(self.coordinates_speed).subs(extraSubs).subs(op_point)
            K =-EOM.jacobian(self.coordinates      ).subs(extraSubs).subs(op_point)
            if len(self.var)>0:
                B = EOM.jacobian(self.var              ).subs(extraSubs).subs(op_point)
            else:
                B=Matrix([])


        return M,C,K,B

    def to_EOM(self):
        """ return a class to easily manipulate the equations of motion in place"""
        EOM = self.EOM.subs(self.kdeqsSubs).doit()

        replaceDict=OrderedDict()
        for b in self.bodies:
            if isinstance(b, YAMSFlexibleBody):
                b.replaceDict(replaceDict)

        return EquationsOfMotionQ(EOM, self.coordinates, self.name, replaceDict)

    def saveTex(self, prefix='', suffix='', folder='./', extraSubs=None, header=True, extraHeader=None, variables=['MM','FF','M','C','K','B','MMsa','FFsa','Msa','Csa','Ksa','Bsa','body_details'], doSimplify=False, velSubs=[(0,0)]):
        """ 
        Save forcing and mass matrix to latex file
        """
        extraSubs = [] if extraSubs is None else extraSubs
        name=prefix+self.name
        if len(suffix)>0:
            name=name+'_'+suffix.strip('_')

        filename=os.path.join(folder,name+'.tex')
        with Timer('Latex to {}'.format(filename),True,silent=True):
            with open(filename,'w') as f:
                if header:
                    f.write('Model: {}, \n'.format(self.name.replace('_','\_')))
                    f.write('Degrees of freedom: ${}$, \n'.format(cleantex(self.coordinates)))
                    try:
                        f.write('Small angles:       ${}$\\\\ \n'.format(cleantex(self.smallAngleUsed)))
                    except:
                        pass
                    f.write('Free vars:       ${}$\\\\ \n'.format(cleantex(self.var)))
                if extraHeader:
                    f.write('\\clearpage\n{}\\\\\n'.format(extraHeader))
                if 'F' in variables:
                    FF = self.kane.forcing.subs(self.kdeqsSubs)
                    toTex(f, FF, label='Forcing', fullPage=True, extraSubs=extraSubs, velSubs=velSubs, doSimplify=doSimplify)
                if 'M' in variables:
                    MM = self.kane.mass_matrix.subs(self.kdeqsSubs) # NOTE this is wrong
                    toTex(f, MM, label='Mass matrix', fullPage=True, extraSubs=extraSubs, velSubs=velSubs, doSimplify=doSimplify)
                if 'M0' in variables:
                    toTex(f, self.M0, label='Linearized mass matrix', fullPage=True, extraSubs=extraSubs, velSubs=velSubs, doSimplify=doSimplify)
                if 'C0' in variables:
                    toTex(f, self.C0, label='Linearized damping matrix', fullPage=True, extraSubs=extraSubs, velSubs=velSubs, doSimplify=doSimplify)
                if 'K0' in variables:
                    toTex(f, self.K0, label='Linearized stiffness matrix', fullPage=True, extraSubs=extraSubs, velSubs=velSubs, doSimplify=doSimplify)
                if 'B0' in variables:
                    toTex(f, self.B0, label='Linearized forcing matrix', fullPage=True, extraSubs=extraSubs, velSubs=velSubs, doSimplify=doSimplify)
                if 'FFsa' in variables:
                    FF=subs_no_diff(self._sa_forcing,extraSubs)
                    toTex(f, FF, label='Forcing small angle', fullPage=True, extraSubs=extraSubs, velSubs=velSubs, doSimplify=doSimplify)
                if 'MMsa' in variables:
                    toTex(f, self._sa_mass_matrix, label='Mass matrix small angle', fullPage=True, extraSubs=extraSubs, velSubs=velSubs, doSimplify=doSimplify)
                if 'Msa' in variables:
                    toTex(f, self._sa_M, label='Linearized mass matrix small angle', fullPage=True, extraSubs=extraSubs, velSubs=velSubs, doSimplify=doSimplify)
                if 'Csa' in variables:
                    toTex(f, self._sa_C, label='Linearized damping matrix small angle', fullPage=True, extraSubs=extraSubs, velSubs=velSubs, doSimplify=doSimplify)
                if 'Ksa' in variables:
                    toTex(f, self._sa_K, label='Linearized stiffness matrix small angle', fullPage=True, extraSubs=extraSubs, velSubs=velSubs, doSimplify=doSimplify)
                if 'body_details' in variables:
                    for b in self.bodies:
                        print('b.name',b.name)
                        print('Jv',b._Jv_vect)
                        print('Jo',b._Jo_vect)
                        toTex(f, b._vel, label='Body {} vel'.format(b.name), fullPage=False, extraSubs=extraSubs, velSubs=velSubs, doSimplify=doSimplify)
                        toTex(f, b._acc, label='Body {} acc'.format(b.name), fullPage=False, extraSubs=extraSubs, velSubs=velSubs, doSimplify=doSimplify)
                        toTex(f, b._omega, label='Body {} omega'.format(b.name), fullPage=False, extraSubs=extraSubs, velSubs=velSubs, doSimplify=doSimplify)
                        toTex(f, b._inertial_force, label='Body {} inertialforce'.format(b.name), fullPage=False, extraSubs=extraSubs, velSubs=velSubs, doSimplify=doSimplify)
                        toTex(f, b._inertial_torque, label='Body {} inertialtorque'.format(b.name), fullPage=False, extraSubs=extraSubs, velSubs=velSubs, doSimplify=doSimplify)
                        try:
                            toTex(f, b.inertial_elast, label='Body {} inertialelast'.format(b.name), fullPage=False, extraSubs=extraSubs, velSubs=velSubs, doSimplify=doSimplify)
                        except:
                            pass
#                         try:
                        toTex(f, b._Jv_vect, label='Body {} JvVect'.format(b.name), fullPage=False, extraSubs=extraSubs, velSubs=velSubs, doSimplify=doSimplify)
                        toTex(f, b._Jo_vect, label='Body {} JoVect'.format(b.name), fullPage=False, extraSubs=extraSubs, velSubs=velSubs, doSimplify=doSimplify)
#                         except:
#                             pass
                    for (jac,force) in self.kane._fr_products:
                        if jac!=0:
                            toTex(f, jac, label='frproducts jac', fullPage=False, extraSubs=extraSubs, velSubs=velSubs, doSimplify=doSimplify)
                            toTex(f, force, label='frproducts force', fullPage=False, extraSubs=extraSubs, velSubs=velSubs, doSimplify=doSimplify)


    # --------------------------------------------------------------------------------}
    # --- Export to Python script 
    # --------------------------------------------------------------------------------{
    def savePython(self, prefix='', suffix='', folder='./', extraSubs=None, variables=['MM','FF','MMsa','FFsa','M','C','K','B','Msa','Csa','Ksa','Bsa'], replaceDict=None, doSimplify=False, velSubs=[(0,0)]):
        """ 
        Save forcing, mass matrix and linear model to python package
        """
        extraSubs = [] if extraSubs is None else extraSubs
        name=prefix+self.name
        if len(suffix)>0:
            name=name+'_'+suffix.strip('_')

        filename=os.path.join(folder,name+'.py')


        # --- replace  dict
        if replaceDict is None: 
            replaceDict=OrderedDict()
        for b in self.bodies:
            if isinstance(b, YAMSFlexibleBody):
                b.replaceDict(replaceDict)
        #print(replaceDict)

        with Timer('Python to {}'.format(filename),True,silent=True):
            with open(filename,'w') as f:
                f.write('"""\n')
                f.write('{}\n'.format(self.__repr__()))
                f.write('"""\n')
                f.write('import numpy as np\n')
                f.write('from numpy import cos, sin\n')

                infoToPy(f, self.name, self.coordinates, self.var)

                if 'FF' in variables:
                    forcing = self.kane.forcing.subs(self.kdeqsSubs)
                    forcingToPy(f, self.coordinates, forcing, replaceDict=replaceDict, extraSubs=extraSubs, velSubs=velSubs, doSimplify=doSimplify)
                if 'MM' in variables:
                    MM = self.kane.mass_matrix.subs(self.kdeqsSubs)
                    MMToPy(f, self.coordinates, MM, replaceDict=replaceDict, extraSubs=extraSubs, velSubs=velSubs, doSimplify=doSimplify)
                if 'M' in variables:
                    #M0ToPy(f, self.q, self.M0, replaceDict=replaceDict, extraSubs=extraSubs, velSubs=velSubs, doSimplify=doSimplify)
                    M0ToPy(f, self.coordinates, self.M0, replaceDict=replaceDict, extraSubs=extraSubs, velSubs=velSubs, doSimplify=doSimplify)
                if 'C' in variables:
                    C0ToPy(f, self.coordinates, self.C0, replaceDict=replaceDict, extraSubs=extraSubs, velSubs=velSubs, doSimplify=doSimplify)
                if 'K' in variables:
                    K0ToPy(f, self.coordinates, self.K0, replaceDict=replaceDict, extraSubs=extraSubs, velSubs=velSubs, doSimplify=doSimplify)
                if 'B' in variables:
                    B0ToPy(f, self.coordinates, self.B0, self.var, replaceDict=replaceDict, extraSubs=extraSubs, velSubs=velSubs, doSimplify=doSimplify)

                if 'FFsa' in variables:
                    FF = subs_no_diff(self._sa_forcing,extraSubs)
                    FF = FF.subs(velSubs)
                    if doSimplify:
                        FF=FF.simplify()
                    s, params, inputs, sdofs  = cleanPy(FF, varname='FF', dofs = self.coordinates, indent=4, replDict=replaceDict)
                    f.write('def forcing_sa(t,q=None,qd=None,p=None,u=None,z=None):\n')
                    f.write('    """ Non linear forcing with small angle approximation\n')
                    f.write('    q:  degrees of freedom, array-like: {}\n'.format(sdofs))
                    f.write('    qd: dof velocities, array-like\n')
                    f.write('    p:  parameters, dictionary with keys: {}\n'.format(params))
                    f.write('    u:  inputs, dictionary with keys: {}\n'.format(inputs))
                    f.write('           where each values is a function of time\n')
                    f.write('    """\n')
                    f.write('    if z is not None:\n')
                    f.write('        q  = z[0:int(len(z)/2)] \n')
                    f.write('        qd = z[int(len(z)/2): ] \n')
                    f.write(s)
                    f.write('    return FF\n\n')

                if 'MMsa' in variables:
                    MM=subs_no_diff(self._sa_mass_matrix,extraSubs)
                    MM = MM.subs(velSubs)
                    if doSimplify:
                        MM.simplify()
                    s, params, inputs, sdofs  = cleanPy(MM, varname='MM', dofs = self.coordinates, indent=4, replDict=replaceDict)
                    f.write('def mass_matrix_sa(q=None,p=None,z=None):\n')
                    f.write('    """ Non linear mass matrix with small angle approximation\n')
                    f.write('     q:  degrees of freedom, array-like: {}\n'.format(sdofs))
                    f.write('     p:  parameters, dictionary with keys: {}\n'.format(params))
                    f.write('    """\n')
                    f.write('    if z is not None:\n')
                    f.write('        q  = z[0:int(len(z)/2)] \n')
                    f.write(s)
                    f.write('    return MM\n\n')

                if 'Msa' in variables:
                    MM=subs_no_diff(self._sa_M,extraSubs)
                    MM = MM.subs(velSubs)
                    if doSimplify:
                        MM.simplify()
                    s, params, inputs, sdofs  = cleanPy(MM, varname='MM', dofs = self.coordinates, indent=4, replDict=replaceDict, noTimeDep=True)
                    f.write('def M_lin_sa(q=None,p=None,z=None):\n')
                    f.write('    """ Linear mass matrix with small angle approximation\n')
                    f.write('    q:  degrees of freedom at operating point, array-like: {}\n'.format(sdofs))
                    f.write('    p:  parameters, dictionary with keys: {}\n'.format(params))
                    f.write('    """\n')
                    f.write('    if z is not None:\n')
                    f.write('        q  = z[0:int(len(z)/2)] \n')
                    f.write(s)
                    f.write('    return MM\n\n')

                if 'Csa' in variables:
                    MM=subs_no_diff(self._sa_C,extraSubs)
                    MM = MM.subs(velSubs)
                    if doSimplify:
                        MM.simplify()
                    s, params, inputs, sdofs  = cleanPy(MM, varname='CC', dofs = self.coordinates, indent=4, replDict=replaceDict, noTimeDep=True)
                    f.write('def C_lin_sa(q=None,qd=None,p=None,u=None,z=None):\n')
                    f.write('    """ Linear damping matrix with small angle approximation\n')
                    f.write('    q:  degrees of freedom at operating point, array-like: {}\n'.format(sdofs))
                    f.write('    qd: dof velocities at operating point, array-like\n')
                    f.write('    p:  parameters, dictionary with keys: {}\n'.format(params))
                    f.write('    u:  inputs at operating point, dictionary with keys: {}\n'.format(inputs))
                    f.write('           where each values is a constant!\n')
                    f.write('    """\n')
                    f.write('    if z is not None:\n')
                    f.write('        q  = z[0:int(len(z)/2)] \n')
                    f.write('        qd = z[int(len(z)/2): ] \n')
                    f.write(s)
                    f.write('    return CC\n\n')

                if 'Ksa' in variables:
                    MM=subs_no_diff(self._sa_K,extraSubs)
                    MM = MM.subs(velSubs)
                    if doSimplify:
                        MM.simplify()
                    s, params, inputs, sdofs  = cleanPy(MM, varname='KK', dofs = self.coordinates, indent=4, replDict=replaceDict, noTimeDep=True)
                    f.write('def K_lin_sa(q=None,qd=None,p=None,u=None,z=None):\n')
                    f.write('    """ Linear stiffness matrix with small angle approximation\n')
                    f.write('    q:  degrees of freedom, array-like: {}\n'.format(sdofs))
                    f.write('    qd: dof velocities, array-like\n')
                    f.write('    p:  parameters, dictionary with keys: {}\n'.format(params))
                    f.write('    u:  inputs at operating point, dictionary with keys: {}\n'.format(inputs))
                    f.write('           where each values is a constant!\n')
                    f.write('    """\n')
                    f.write('    if z is not None:\n')
                    f.write('        q  = z[0:int(len(z)/2)] \n')
                    f.write('        qd = z[int(len(z)/2): ] \n')
                    f.write(s)
                    f.write('    return KK\n\n')

                if 'Bsa' in variables:
                    MM=subs_no_diff(self._sa_B,extraSubs)
                    MM = MM.subs(velSubs)
                    if doSimplify:
                        MM.simplify()
                    s, params, inputs, sdofs  = cleanPy(MM, varname='BB', dofs = self.coordinates, indent=4, replDict=replaceDict, noTimeDep=True)
                    f.write('def B_lin_sa(q=None,qd=None,p=None,u=None):\n')
                    f.write('    """ Linear mass matrix with small angle approximation\n')
                    f.write('    q:  degrees of freedom at operating point, array-like: {}\n'.format(sdofs))
                    f.write('    qd: dof velocities at operating point, array-like\n')
                    f.write('    p:  parameters, dictionary with keys: {}\n'.format(params))
                    f.write('    u:  inputs at operating point, dictionary with keys: {}\n'.format(inputs))
                    f.write('           where each values is a constant!\n')
                    f.write('    The columns of B correspond to:   {}\\\\ \n'.format(self.var))
                    f.write('    """\n')
                    f.write(s)
                    f.write('    return BB\n\n')


    def save(self, filename):
        """ Save to pickle """
        import pickle
        with open(filename,'wb') as f:
            pickle.dump(repr(self._sa_forcing),f)

    @staticmethod
    def load(filename):
        """ Save to pickle """
        import pickle
        from sympy.parsing.sympy_parser import parse_expr
        with open(filename,'rb') as f:
            #return parse_expr(pickle.load(f))
            return pickle.load(f)

# --------------------------------------------------------------------------------}
# --- Equation of motions 
# --------------------------------------------------------------------------------{
class EquationsOfMotionQ(object):
    """
    Class to handle equations of motion (EOM), assumed to be of the form:
       e(q,qdot,qddot,p,u) = 0  = F - M qdd
    with e, a vector of dimension n

    The following operations can be applied to the EOM:
      - small angle approximation
      - substitution
      - simplifications
    All operations occur in place, and thus modify the equations
    After thee operations, once can:
      - Split into Mass matrix and right hand side forcing
      - linearize the equations
    Export are possible:
      - to python 
      - to latex
    """

    def __init__(self, EOM, q, name, bodyReplaceDict):
        self.EOM = EOM
        self.q   = q
        # Non linear system
        self.M=None
        self.F=None
        # Linearized matrices
        self.M0=None
        self.K0=None
        self.C0=None
        self.B0=None
        self.input_vars=None
        self.smallAngleUsed=[]

        self.name=name
        self.bodyReplaceDict=bodyReplaceDict

    def subs(self, subs_list):
        """ Apply substitutions to equations of motion """
        self.EOM = self.EOM.subs(subs_list)

    def simplify(self):
        """ Simplify equations of motion """
        with Timer('Simplify',True,silent=True):
            self.EOM.simplify()

    def trigsimp(self):
        """ Trigonometric simplifications of  equations of motion """
        self.EOM = trigsimp(self.EOM)

    def expand(self):
        """ Trigonometric simplifications of  equations of motion """
        self.EOM = self.EOM.expand()

    def smallAngleApprox(self, angle_list, order=1):
        """ 
        Apply small angle approximation to EOM
        """
        with Timer('Small angle approx',True,silent=True):
            self.EOM = smallAngleApprox(self.EOM, angle_list, order=order)
        self.smallAngleUsed+=angle_list

    def mass_forcing_form(self, extraSubs=None):
        """ Extract Mass Matrix and RHS from EOM """
        extraSubs = [] if extraSubs is None else extraSubs
        qd  = [diff(c,dynamicsymbols._t) for c in self.q]
        qdd = [diff(diff(c,dynamicsymbols._t),dynamicsymbols._t) for c in self.q]
        self.M = - myjacobian(self.EOM, qdd)  # mass matrix is jacobian wrt qdd
        self.F = (self.M * Matrix(qdd) + self.EOM).expand() # remainder
        self.F = self.F.subs([(qddi,0) for qddi in qdd  ]) # safety
        self.F = self.F.expand()

    def linearize(self, op_point=None, noAcc=True, noVel=False, extraSubs=None):
        """ """
        op_point  = [] if op_point is None else op_point
        extraSubs = [] if extraSubs is None else extraSubs
        self.M0,self.C0,self.K0,self.B0, self.input_vars = linearizeQ(self.EOM, self.q, op_point=op_point, noAcc=noAcc, noVel=noVel, extraSubs=extraSubs)

    def saveTex(self, prefix='', suffix='', folder='./', extraSubs=None, header=True, extraHeader=None, variables=['M','F','M0','C0','K0','B0'], doSimplify=False, velSubs=[(0,0)]):
        """ 
        Save EOM to a latex file
        """
        extraSubs = [] if extraSubs is None else extraSubs
        name=prefix+self.name
        if len(suffix)>0:
            name=name+'_'+suffix.strip('_')

        filename=os.path.join(folder,name+'.tex')
        with Timer('Latex to {}'.format(filename),True,silent=True):
            with open(filename,'w') as f:
                if header:
                    f.write('Model: {}, \n'.format(self.name.replace('_','\_')))
                    f.write('Degrees of freedom: ${}$, \n'.format(cleantex(self.q)))
                    try:
                        f.write('Small angles:       ${}$\\\\ \n'.format(cleantex(self.smallAngleUsed)))
                    except:
                        pass
                    f.write('Free vars:       ${}$\\\\ \n'.format(cleantex(self.input_vars)))

                if extraHeader:
                    f.write('\\clearpage\n{}\\\\\n'.format(extraHeader))

                if 'F' in variables:
                    toTex(f, self.F, label='Forcing', fullPage=True, extraSubs=extraSubs, velSubs=velSubs, doSimplify=doSimplify)
                if 'M' in variables:
                    toTex(f, self.M, label='Mass matrix', fullPage=True, extraSubs=extraSubs, velSubs=velSubs, doSimplify=doSimplify)
                if 'M0' in variables:
                    toTex(f, self.M0, label='Linearized mass matrix', fullPage=True, extraSubs=extraSubs, velSubs=velSubs, doSimplify=doSimplify)
                if 'C0' in variables:
                    toTex(f, self.C0, label='Linearized damping matrix', fullPage=True, extraSubs=extraSubs, velSubs=velSubs, doSimplify=doSimplify)
                if 'K0' in variables:
                    toTex(f, self.K0, label='Linearized stiffness matrix', fullPage=True, extraSubs=extraSubs, velSubs=velSubs, doSimplify=doSimplify)
                if 'B0' in variables:
                    toTex(f, self.B0, label='Linearized forcing matrix', fullPage=True, extraSubs=extraSubs, velSubs=velSubs, doSimplify=doSimplify)
    
    def savePython(self, prefix='', suffix='', folder='./', extraSubs=None, replaceDict=None, doSimplify=False, velSubs=[(0,0)]):
        """ Export EOM to a python file (will be rewritten) """
        extraSubs = [] if extraSubs is None else extraSubs
        name=prefix+self.name
        if len(suffix)>0:
            name=name+'_'+suffix.strip('_')
        filename=os.path.join(folder,name+'.py')
        # --- replace  dict
        if replaceDict is None: 
            replaceDict=OrderedDict()
        replaceDict.update(self.bodyReplaceDict)
        #print(replaceDict)

        with Timer('Python to {}'.format(filename),True,silent=True):
            with open(filename,'w') as f:
                f.write('"""\n')
                f.write('Equations of motion\n')
                f.write('model name: {}\n'.format(self.name))
                f.write('"""\n')
                f.write('import numpy as np\n')
                f.write('from numpy import cos, sin\n')

                infoToPy(f, self.name, self.q, self.input_vars)

                if self.M is not None:
                    forcingToPy(f, self.q, self.F, replaceDict=replaceDict, extraSubs=extraSubs, velSubs=velSubs, doSimplify=doSimplify)
                    MMToPy(f, self.q, self.M, replaceDict=replaceDict, extraSubs=extraSubs, velSubs=velSubs, doSimplify=doSimplify)

                if self.M0 is not None:
                    M0ToPy(f, self.q, self.M0, replaceDict=replaceDict, extraSubs=extraSubs, velSubs=velSubs, doSimplify=doSimplify)
                    C0ToPy(f, self.q, self.C0, replaceDict=replaceDict, extraSubs=extraSubs, velSubs=velSubs, doSimplify=doSimplify)
                    K0ToPy(f, self.q, self.K0, replaceDict=replaceDict, extraSubs=extraSubs, velSubs=velSubs, doSimplify=doSimplify)
                    B0ToPy(f, self.q, self.B0, self.input_vars, replaceDict=replaceDict, extraSubs=extraSubs, velSubs=velSubs, doSimplify=doSimplify)

        # TODO
        #MM_s = insertQQd(MM, coordinates)
        #RHS_s = insertQQd(RHS, coordinates)


# --------------------------------------------------------------------------------}
# --- Main helper functions  
# --------------------------------------------------------------------------------{
def linearizeQ(EOM, q, op_point=None, noAcc=True, noVel=False, extraSubs=None):
    op_point  = [] if op_point is None else op_point
    extraSubs = [] if extraSubs is None else extraSubs
    qd  = [qi.diff(dynamicsymbols._t) for qi in q]
    qdd = [qdi.diff(dynamicsymbols._t) for qdi in qd]

    with Timer('Linearization',True,silent=True):
        # NOTE: order important
        op_point=[]
        if noAcc: 
            op_point+=[(qddi,0) for qddi in qdd]
        if noVel: 
            op_point+=[(qdi,0) for qdi in qd]
    
        # --- Inputs are dynamic symbols that are not coordinates
        dyn_symbols = find_dynamicsymbols(EOM)
        all_q       = q + qd + qdd
        inputs      = [s for s in dyn_symbols if s not in all_q]
        II = np.argsort([str(v) for v in inputs]) # sorting for consistency
        inputs = list(np.array(inputs)[II])
        # KEEP ME Alternative
        #M, A, B = linearizer.linearize(op_point=op_point ) #o

        M =-EOM.jacobian(qdd).subs(extraSubs).subs(op_point)
        C =-EOM.jacobian(qd ).subs(extraSubs).subs(op_point)
        K =-EOM.jacobian(q  ).subs(extraSubs).subs(op_point)
        if len(inputs)>0:
            B = EOM.jacobian(inputs).subs(extraSubs).subs(op_point)
        else:
            B=Matrix([])
    return M,C,K,B, inputs


def forcingToPy(f, q, forcing, replaceDict=None, extraSubs=None, velSubs=[(0,0)], doSimplify=False):
    extraSubs = [] if extraSubs is None else extraSubs
    forcing = subs_no_diff(forcing, extraSubs)
    forcing = forcing.subs(velSubs)
    if doSimplify:
        forcing=trigsimp(forcing)

    s, params, inputs, sdofs  = cleanPy(forcing, varname='FF', dofs = q, indent=4, replDict=replaceDict)
    f.write('def forcing(t,q=None,qd=None,p=None,u=None,z=None):\n')
    f.write('    """ Non linear mass forcing \n')
    f.write('    q:  degrees of freedom, array-like: {}\n'.format(sdofs))
    f.write('    qd: dof velocities, array-like\n')
    f.write('    p:  parameters, dictionary with keys: {}\n'.format(params))
    f.write('    u:  inputs, dictionary with keys: {}\n'.format(inputs))
    f.write('           where each values is a function of time\n')
    f.write('    """\n')
    f.write('    if z is not None:\n')
    f.write('        q  = z[0:int(len(z)/2)] \n')
    f.write('        qd = z[int(len(z)/2): ] \n')
    f.write(s)
    f.write('    return FF\n\n')


def MMToPy(f, q, MM, replaceDict=None, extraSubs=None, velSubs=[(0,0)], doSimplify=False):
    extraSubs = [] if extraSubs is None else extraSubs
    MM = subs_no_diff(MM,extraSubs)
    MM = MM.subs(velSubs)
    if doSimplify:
        MM=trigsimp(MM)
        #MM.simplify()
    s, params, inputs, sdofs  = cleanPy(MM, varname='MM', dofs = q, indent=4, replDict=replaceDict)
    f.write('def mass_matrix(q=None,p=None,z=None):\n')
    f.write('    """ Non linear mass matrix \n')
    f.write('     q:  degrees of freedom, array-like: {}\n'.format(sdofs))
    f.write('     p:  parameters, dictionary with keys: {}\n'.format(params))
    f.write('    """\n')
    f.write('    if z is not None:\n')
    f.write('        q  = z[0:int(len(z)/2)] \n')
    f.write(s)
    f.write('    return MM\n\n')


def M0ToPy(f, q, MM, replaceDict=None, extraSubs=None, velSubs=[(0,0)], doSimplify=False):
    extraSubs = [] if extraSubs is None else extraSubs
    MM=subs_no_diff(MM,extraSubs)
    MM = MM.subs(velSubs)
    if doSimplify:
        MM.simplify()
    s, params, inputs, sdofs  = cleanPy(MM, varname='MM', dofs = q, indent=4, replDict=replaceDict, noTimeDep=True)
    f.write('def M_lin(q=None,p=None,z=None):\n')
    f.write('    """ Linear mass matrix \n')
    f.write('    q:  degrees of freedom at operating point, array-like: {}\n'.format(sdofs))
    f.write('    p:  parameters, dictionary with keys: {}\n'.format(params))
    f.write('    """\n')
    f.write('    if z is not None:\n')
    f.write('        q  = z[0:int(len(z)/2)] \n')
    f.write(s)
    f.write('    return MM\n\n')

def C0ToPy(f, q, MM, replaceDict=None, extraSubs=None, velSubs=[(0,0)], doSimplify=False):
    extraSubs = [] if extraSubs is None else extraSubs
    MM=subs_no_diff(MM,extraSubs)
    MM = MM.subs(velSubs)
    if doSimplify:
        MM.simplify()
    s, params, inputs, sdofs  = cleanPy(MM, varname='CC', dofs = q, indent=4, replDict=replaceDict, noTimeDep=True)
    f.write('def C_lin(q=None,qd=None,p=None,u=None,z=None):\n')
    f.write('    """ Linear damping matrix \n')
    f.write('    q:  degrees of freedom at operating point, array-like: {}\n'.format(sdofs))
    f.write('    qd: dof velocities at operating point, array-like\n')
    f.write('    p:  parameters, dictionary with keys: {}\n'.format(params))
    f.write('    u:  inputs at operating point, dictionary with keys: {}\n'.format(inputs))
    f.write('           where each values is a constant!\n')
    f.write('    """\n')
    f.write('    if z is not None:\n')
    f.write('        q  = z[0:int(len(z)/2)] \n')
    f.write('        qd = z[int(len(z)/2): ] \n')
    f.write(s)
    f.write('    return CC\n\n')

def K0ToPy(f, q, MM, replaceDict=None, extraSubs=None, velSubs=[(0,0)], doSimplify=False):
    extraSubs = [] if extraSubs is None else extraSubs
    MM=subs_no_diff(MM,extraSubs)
    MM = MM.subs(velSubs)
    if doSimplify:
        MM.simplify()
    s, params, inputs, sdofs  = cleanPy(MM, varname='KK', dofs = q, indent=4, replDict=replaceDict, noTimeDep=True)
    f.write('def K_lin(q=None,qd=None,p=None,u=None,z=None):\n')
    f.write('    """ Linear stiffness matrix \n')
    f.write('    q:  degrees of freedom, array-like: {}\n'.format(sdofs))
    f.write('    qd: dof velocities, array-like\n')
    f.write('    p:  parameters, dictionary with keys: {}\n'.format(params))
    f.write('    u:  inputs at operating point, dictionary with keys: {}\n'.format(inputs))
    f.write('           where each values is a constant!\n')
    f.write('    """\n')
    f.write('    if z is not None:\n')
    f.write('        q  = z[0:int(len(z)/2)] \n')
    f.write('        qd = z[int(len(z)/2): ] \n')
    f.write(s)
    f.write('    return KK\n\n')

def B0ToPy(f, q, MM, input_vars, replaceDict=None, extraSubs=None, velSubs=[(0,0)], doSimplify=False):
    extraSubs = [] if extraSubs is None else extraSubs
    MM=subs_no_diff(MM,extraSubs)
    MM = MM.subs(velSubs)
    if doSimplify:
        MM.simplify()
    s, params, inputs, sdofs  = cleanPy(MM, varname='BB', dofs = q, indent=4, replDict=replaceDict, noTimeDep=True)
    f.write('def B_lin(q=None,qd=None,p=None,u=None):\n')
    f.write('    """ Linear mass matrix \n')
    f.write('    q:  degrees of freedom at operating point, array-like: {}\n'.format(sdofs))
    f.write('    qd: dof velocities at operating point, array-like\n')
    f.write('    p:  parameters, dictionary with keys: {}\n'.format(params))
    f.write('    u:  inputs at operating point, dictionary with keys: {}\n'.format(inputs))
    f.write('           where each values is a constant!\n')
    f.write('    The columns of B correspond to:   {}\\\\ \n'.format(input_vars))
    f.write('    """\n')
    f.write(s)
    f.write('    return BB\n\n')

def infoToPy(f, name, q, u):
    f.write('def info():\n')
    f.write('    """ Return information about current model present in this package """\n')
    f.write('    I=dict()\n')
    f.write('    I[\'name\']=\'{}\'\n'.format(name))
    f.write('    I[\'nq\']={}\n'.format(len(q)))
    f.write('    I[\'nu\']={}\n'.format(len(u)))
    f.write('    I[\'sq\']=[{}]\n'.format(','.join(['\''+repr(qi).replace('(t)','')+'\'' for qi in q])))
    f.write('    I[\'su\']=[{}]\n'.format(','.join(['\''+repr(ui).replace('(t)','')+'\'' for ui in u])))
    f.write('    return I\n\n')


def toTex(f, FF, label='', fullPage=False, extraSubs=None, velSubs=[(0,0)], doSimplify=False):
    extraSubs = [] if extraSubs is None else extraSubs
    if isinstance(FF,list):
        for i,ff in enumerate(FF):
            FF[i] = subs_no_diff(ff,extraSubs)
            FF[i] = FF[i].subs(velSubs)
    else:
        FF = subs_no_diff(FF,extraSubs)
        FF = FF.subs(velSubs)
    if doSimplify:
        FF=trigsimp(FF)
    if len(label)>0:
        f.write('{}:\n'.format(label))
    f.write('\\begin{equation*}\n')
    if fullPage:
        f.write('\\resizebox{\\textwidth}{!}{$\n')
    if isinstance(FF,list):
        for i,ff in enumerate(FF):
            f.write('a[{:}]='.format(i)+cleantex(ff)+'\n')
    else:
        f.write(cleantex(FF))
    if fullPage:
        f.write('$}\n')
    f.write('\\end{equation*}\n')

