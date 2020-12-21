import os
from sympy.physics.mechanics import dynamicsymbols
from sympy import diff

from .yams_sympy_tools import smallAngleApprox, cleantex, subs_no_diff , cleanPy
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
        self._sa_forcing     = None
        self._sa_mass_matrix = None
        self._sa_EOM         = None
        self.opts        = None
        self.body_loads  = None
        self.var         = [] # Independent variables
    
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


    def kaneEquations(self, Mform='symbolic'):
        for sa in ['ref', 'coordinates', 'speeds','kdeqs','bodies','loads']:
            if getattr(self,sa) is None:
                raise Exception('Attribute {} needs to be set before calling `kane` method'.format(sa))

        with Timer('Kane step1',True):
            self.kane = YAMSKanesMethod(self.ref.frame, self.coordinates, self.speeds, self.kdeqs)

        # --- Expensive kane step
        with Timer('Kane step 2',True):
            #(use  Mform ='symbolic' or 'TaylorExpanded'), Mform='symbolic'
            self.fr, self.frstar  = self.kane.kanes_equations(self.bodies, self.loads, Mform=Mform)
        self.kane.fr     = self.fr
        self.kane.frstar = self.frstar

    def smallAngleApprox(self, angle_list, extraSubs=[]):
        """ 
        Apply small angle approximation to forcing and mass matrix
        """
        # Forcing
        with Timer('Small angle approx. forcing',True):
            if self._sa_forcing is None:
                self._sa_forcing=self.kane.forcing
            self._sa_forcing = self._sa_forcing.subs(self.kdeqsSubs).subs(extraSubs)
            self._sa_forcing = smallAngleApprox(self._sa_forcing, angle_list)
        with Timer('Small angle approx. forcing simplify',True):
            self._sa_forcing.simplify()
        # Mass matrix
        with Timer('Small angle approx. mass matrix',True):
            if self._sa_mass_matrix is None:
                self._sa_mass_matrix=self.kane.mass_matrix
            self._sa_mass_matrix = self._sa_mass_matrix.subs(self.kdeqsSubs).subs(extraSubs)
            self._sa_mass_matrix = smallAngleApprox(self._sa_mass_matrix, angle_list)
        with Timer('Small angle approx. mass matrix simplify',True):
            self._sa_mass_matrix.simplify()

        self.smallAngleUsed=angle_list
        

    def smallAngleApproxEOM(self, angle_list, extraSubs=[]):
        """ 
        Apply small angle approximation to equation of motion H(x,xd,xdd,..)=0
        """
        EOM=self.EOM
        with Timer('Small angle approx. EOM',True):
            EOM = EOM.subs(extraSubs)
            EOM = smallAngleApprox(EOM, angle_list).subs(extraSubs)
        #with Timer('Small angle approx. EOM simplify',True):
        #    EOM.simplify()
        self._sa_EOM = EOM
       
    def smallAngleLinearize(self, op_point=[], noAcc=True, noVel=False, extraSubs=[]):
        """ 
        Linearize the equations with small angle approximations
        """
        if self._sa_EOM is None:
            raise Exception('Run smallAngleApproxEOM first')
        M,C,K,B = self.linearize(op_point=op_point, EOM=self._sa_EOM, noAcc=noAcc, noVel=noVel, extraSubs=extraSubs)

        self._sa_M = M
        self._sa_C = C
        self._sa_K = K
        self._sa_B = B

    def linearize(self, op_point=[], EOM=None, noAcc=True, noVel=False, extraSubs=[]):
        if EOM is None:
            EOM=self.EOM

        with Timer('Linearization',True):
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
            # KEEP ME Alternative
            #M, A, B = linearizer.linearize(op_point=op_point ) #o

            M =-EOM.jacobian(self.coordinates_acc  ).subs(extraSubs).subs(op_point)
            C =-EOM.jacobian(self.coordinates_speed).subs(extraSubs).subs(op_point)
            K =-EOM.jacobian(self.coordinates      ).subs(extraSubs).subs(op_point)
            B = EOM.jacobian(self.var              ).subs(extraSubs).subs(op_point)

        return M,C,K,B



    def smallAngleSaveTex(self, prefix='', folder='./', extraSubs=[], header=True, extraHeader=None, variables=['MM','FF','M','C','K','B']):
        """ 
        Save forcing and mass matrix to latex file
        """

        if len(prefix)>0:
            prefix=self.name+'_'+prefix.strip('_')
        else:
            prefix=self.name


        filename=os.path.join(folder,prefix+'.tex')
        with Timer('Latex to {}'.format(filename),True):
            with open(filename,'w') as f:
                if header:
                    f.write('Model: {}, \n'.format(self.name.replace('_','\_')))
                    f.write('Degrees of freedom: ${}$, \n'.format(cleantex(self.coordinates)))
                    f.write('Small angles:       ${}$\\\\ \n'.format(cleantex(self.smallAngleUsed)))
                    f.write('Free vars:       ${}$\\\\ \n'.format(cleantex(self.var)))
                    #if self.opts is not None:
                    #    f.write('Options: \\begin{lstlisting}\n')
                    #    f.write('{}\n'.format(self.opts))
                    #    f.write('\\end{lstlisting}\n')

                if extraHeader:
                    f.write('\\clearpage\n{}\\\\\n'.format(extraHeader))

                if 'FF' in variables:
                    FF=subs_no_diff(self._sa_forcing,extraSubs)
                    #FF.simplify()

                    f.write('Forcing vector:\n')
                    f.write('\\begin{equation*}\n')
                    f.write('\\resizebox{\\textwidth}{!}{$\n')
                    f.write(cleantex(FF))
                    f.write('$}\n')
                    f.write('\\end{equation*}\n')

                if 'MM' in variables:
                    MM=subs_no_diff(self._sa_mass_matrix,extraSubs)
                    #MM.simplify()
                    f.write('Mass matrix:\n')
                    f.write('\\begin{equation*}\n')
                    f.write('\\resizebox{\\textwidth}{!}{$\n')
                    f.write(cleantex(MM))
                    f.write('$}\n')
                    f.write('\\end{equation*}\n')

                if 'M' in variables:
                    MM=subs_no_diff(self._sa_M,extraSubs)
                    #MM.simplify()
                    f.write('Linearized mass matrix:\n')
                    f.write('\\begin{equation*}\n')
                    f.write('\\resizebox{\\textwidth}{!}{$\n')
                    f.write(cleantex(MM))
                    f.write('$}\n')
                    f.write('\\end{equation*}\n')

                if 'C' in variables:
                    MM=subs_no_diff(self._sa_C,extraSubs)
                    #MM.simplify()
                    stretch=True
                    if MM==MM*0:
                        stretch=False
                    f.write('Linearized damping matrix:\n')
                    f.write('\\begin{equation*}\n')
                    if stretch:
                        f.write('\\resizebox{\\textwidth}{!}{$\n')
                    f.write(cleantex(MM))
                    if stretch:
                        f.write('$}\n')
                    f.write('\\end{equation*}\n')

                if 'K' in variables:
                    MM=subs_no_diff(self._sa_K,extraSubs)
                    #MM.simplify()
                    f.write('Linearized stiffness matrix:\n')
                    f.write('\\begin{equation*}\n')
                    f.write('\\resizebox{\\textwidth}{!}{$\n')
                    f.write(cleantex(MM))
                    f.write('$}\n')
                    f.write('\\end{equation*}\n')


                if 'B' in variables:
                    MM=subs_no_diff(self._sa_B,extraSubs)
                    #MM.simplify()
                    f.write('Linearized forcing matrix:\n')
                    f.write('\\begin{equation*}\n')
                    f.write('\\resizebox{\\textwidth}{!}{$\n')
                    f.write(cleantex(MM))
                    f.write('$}\n')
                    f.write('\\end{equation*}\n')

    # --------------------------------------------------------------------------------}
    # --- Export to Python script 
    # --------------------------------------------------------------------------------{
    def smallAngleSavePython(self, prefix='', folder='./', extraSubs=[], variables=['MM','FF','M','C','K','B'], replaceDict=None):
        """ 
        Save forcing, mass matrix and linear model to python package
        """
        from sympy import python
        import sympy.printing.pycode as pycode

        if len(prefix)>0:
            prefix=self.name+'_'+prefix.strip('_')
        else:
            prefix=self.name
        filename=os.path.join(folder,prefix+'.py')


        # --- replace  dict
        if replaceDict is None: 
            replaceDict=OrderedDict()
        for b in self.bodies:
            if isinstance(b, YAMSFlexibleBody):
                b.replaceDict(replaceDict)
        #print(replaceDict)

        with Timer('Python to {}'.format(filename),True):
            with open(filename,'w') as f:
                f.write('"""\n')
                f.write('Model: {}, \n'.format(self.name.replace('_','\_')))
                f.write('Degrees of freedom: ${}$, \n'.format(cleantex(self.coordinates)))
                f.write('Small angles:       ${}$\\\\ \n'.format(cleantex(self.smallAngleUsed)))
                f.write('Free vars:       ${}$\\\\ \n'.format(cleantex(self.var)))
                f.write('"""\n')
                f.write('import numpy as np\n')
                f.write('from numpy import cos, sin\n')


                if 'FF' in variables:
                    FF = subs_no_diff(self._sa_forcing,extraSubs)
                    s, params, inputs, sdofs  = cleanPy(FF, varname='FF', dofs = self.coordinates, indent=4, replDict=replaceDict)
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

                if 'MM' in variables:
                    MM=subs_no_diff(self._sa_mass_matrix,extraSubs)
                    s, params, inputs, sdofs  = cleanPy(MM, varname='MM', dofs = self.coordinates, indent=4, replDict=replaceDict)
                    f.write('def mass_matrix(q=None,p=None,z=None):\n')
                    f.write('    """ Non linear mass matrix \n')
                    f.write('     q:  degrees of freedom, array-like: {}\n'.format(sdofs))
                    f.write('     p:  parameters, dictionary with keys: {}\n'.format(params))
                    f.write('    """\n')
                    f.write('    if z is not None:\n')
                    f.write('        q  = z[0:int(len(z)/2)] \n')
                    f.write(s)
                    f.write('    return MM\n\n')

                if 'M' in variables:
                    MM=subs_no_diff(self._sa_M,extraSubs)
                    s, params, inputs, sdofs  = cleanPy(MM, varname='MM', dofs = self.coordinates, indent=4, replDict=replaceDict)
                    f.write('def M_lin(q=None,p=None,z=None):\n')
                    f.write('    """ Linear mass matrix \n')
                    f.write('    q:  degrees of freedom at operating point, array-like: {}\n'.format(sdofs))
                    f.write('    p:  parameters, dictionary with keys: {}\n'.format(params))
                    f.write('    """\n')
                    f.write('    if z is not None:\n')
                    f.write('        q  = z[0:int(len(z)/2)] \n')
                    f.write(s)
                    f.write('    return MM\n\n')

                if 'C' in variables:
                    MM=subs_no_diff(self._sa_C,extraSubs)
                    s, params, inputs, sdofs  = cleanPy(MM, varname='CC', dofs = self.coordinates, indent=4, replDict=replaceDict)
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

                if 'K' in variables:
                    MM=subs_no_diff(self._sa_K,extraSubs)
                    s, params, inputs, sdofs  = cleanPy(MM, varname='KK', dofs = self.coordinates, indent=4, replDict=replaceDict)
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

                if 'B' in variables:
                    MM=subs_no_diff(self._sa_B,extraSubs)
                    s, params, inputs, sdofs  = cleanPy(MM, varname='BB', dofs = self.coordinates, indent=4, replDict=replaceDict)
                    f.write('def B_lin(q=None,qd=None,p=None,u=None):\n')
                    f.write('    """ Linear mass matrix \n')
                    f.write('    q:  degrees of freedom at operating point, array-like: {}\n'.format(sdofs))
                    f.write('    qd: dof velocities at operating point, array-like\n')
                    f.write('    p:  parameters, dictionary with keys: {}\n'.format(params))
                    f.write('    u:  inputs at operating point, dictionary with keys: {}\n'.format(inputs))
                    f.write('           where each values is a constant!\n')
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

