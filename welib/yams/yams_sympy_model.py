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
    def __init__(self, name='', ref=None, bodies=None, coordinates=None, speeds=None, kdeqsSubs=None,
            body_loads=None, g_vect=None, opts=None):
        """ 
        INPUTS:
         - ref: reference frame body, YAMSInertialBody
         - bodies: list of YAMSFlexibleBody or YAMSRigidBody
         - coordinates: list of degrees of freedom
         - speeds: list of speeds
         - kdeqsSubs: relationships between speeds and coordinates derivatives
         - body_loads: list of loads. A load is a tuple:
                                (body, (point, force) )   or (body, (frame, moment))
         - g_vect: gravity vector g_vect  = -gravity * ref.frame.z


        """
        self.name        = name
        self.ref         = ref
        self.coordinates = coordinates
        self.speeds      = speeds
        self.kdeqsSubs   = kdeqsSubs
        self.bodies      = bodies
        self.body_loads  = body_loads
        self.g_vect      = g_vect
        self.opts        = opts
        # Generated / Internal data
        self.kane        = None
        self._sa_forcing     = None
        self._sa_mass_matrix = None
        self._sa_M           = None
        self._sa_C           = None
        self._sa_K           = None
        self._sa_B           = None
        self._sa_EOM         = None
        self.var         = [] # Independent variables
        self.smallAnglesUsed=[]

        self.M      = None # Non-linear Mass matrix
        self.F      = None # Non-linear Forcing
        self.M0     = None # linear Mass matrix
        self.K0     = None # linear Stiffness matrix
        self.C0     = None # linear Damping matrix
        self.B0     = None # linear Forcing
        self.Points = [] # Points of interest
        self.PointsFrames = [] # Frames for Points of interest 
        self.PointsMotions = [] # Points of interest

    def __repr__(self):
        s='<{} object "{}" with attributes:>\n'.format(type(self).__name__,self.name)
        s+=' - coordinates:       {}\n'.format(self.coordinates)
        s+=' - speeds:            {}\n'.format(self.speeds)
        s+=' - kdeqsSubs:         {}\n'.format(self.kdeqsSubs)
        s+=' - var:               {}\n'.format(self.var)
        s+=' - smallAnglesUsed  : {}\n'.format(self.smallAnglesUsed)
        s+=' - number of bodies : {}\n'.format(len(self.bodies))
        s+=' - opts             : {}\n'.format(self.opts)
        s+=' * loads            : {}\n'.format(self.loads)
        return s

    @property
    def q(self):
        return self.coordinates
    
    @property
    def kdeqs(self):
        return [s-v for (s,v) in self.kdeqsSubs]

    @property
    def loads(self):
        return [f[1] for f in self.body_loads]

    def EOM(self, Mform='symbolic', extraSubs=None):
        if self.kane is None:
            self.kaneEquations(Mform=Mform)
        if extraSubs is not None:
            return (self.fr+self.frstar).subs(self.kdeqsSubs).subs(extraSubs)
        else:
            return (self.fr+self.frstar).subs(self.kdeqsSubs)

    @property
    def coordinates_speed(self):
        return [q.diff(dynamicsymbols._t) for q in self.coordinates]

    @property
    def coordinates_acc(self):
        return [qd.diff(dynamicsymbols._t) for qd in self.coordinates_speed]

    def addForce(self, body, point, force):
        """ """
        self.body_loads.append((body, (point, force)))

    def addMoment(self, body, frame, moment):
        """ """
        self.body_loads.append((body, (frame, moment)))

    def addPoint(self, P, frame=None):
        """ add a point of interest """
        if frame is None:
            frame =self.ref.frame
        self.Points.append(P)
        self.PointsFrames.append(frame)


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

    def smallAngleApprox(self, angle_list, extraSubs=None, order=1):
        """ 
        Apply small angle approximation to forcing and mass matrix
        NOTE: can be called multiple times with different angle list (cumulative effect)
        """
        extraSubs = [] if extraSubs is None else extraSubs
        # Forcing
        with Timer('Small angle approx. forcing',True,silent=True):
            if self._sa_forcing is None:
                self._sa_forcing=self.kane.forcing
            self._sa_forcing = self._sa_forcing.subs(self.kdeqsSubs).subs(extraSubs)
            self._sa_forcing = smallAngleApprox(self._sa_forcing, angle_list, order=order)
        with Timer('Small angle approx. forcing simplify',True,silent=True):
            self._sa_forcing.simplify()
        # Mass matrix
        with Timer('Small angle approx. mass matrix',True,silent=True):
            if self._sa_mass_matrix is None:
                self._sa_mass_matrix=self.kane.mass_matrix
            self._sa_mass_matrix = self._sa_mass_matrix.subs(self.kdeqsSubs).subs(extraSubs)
            self._sa_mass_matrix = smallAngleApprox(self._sa_mass_matrix, angle_list, order=order)
        with Timer('Small angle approx. mass matrix simplify',True,silent=True):
            self._sa_mass_matrix.simplify()

        self.smallAnglesUsed+=angle_list
        

    def smallAngleApproxEOM(self, angle_list, extraSubs=None, order=1):
        """ 
        Apply small angle approximation to equation of motion H(x,xd,xdd,..)=0
        """
        extraSubs = [] if extraSubs is None else extraSubs
        EOM=self.EOM()
        with Timer('Small angle approx. EOM',True,silent=True):
            EOM = EOM.subs(extraSubs)
            EOM = smallAngleApprox(EOM, angle_list, order=order).subs(extraSubs)
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
        (YAMS Model)
        Linearize the "non" linear equations
        NOTE: EOM has kdeqsSubs in it
        """
        op_point  = [] if op_point is None else op_point
        extraSubs = [] if extraSubs is None else extraSubs
        M,C,K,B = self._linearize(op_point=op_point, EOM=self.EOM(), noAcc=noAcc, noVel=noVel, extraSubs=extraSubs)

        self.M0 = M
        self.C0 = C
        self.K0 = K
        self.B0 = B

    def _linearize(self, op_point=None, EOM=None, noAcc=True, noVel=False, extraSubs=None):
        op_point  = [] if op_point is None else op_point
        extraSubs = [] if extraSubs is None else extraSubs
        if EOM is None:
            EOM=self.EOM()

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

    def to_EOM(self, extraSubs=None):
        """ return a class to easily manipulate the equations of motion in place"""
        EOM = self.EOM().subs(self.kdeqsSubs).doit()
        if extraSubs is not None:
            EOM = EOM.subs(extraSubs)

        bodyReplaceDict=OrderedDict()
        for b in self.bodies:
            if isinstance(b, YAMSFlexibleBody):
                b.replaceDict(bodyReplaceDict)

        return EquationsOfMotionQ(EOM, self.coordinates, self.name, bodyReplaceDict)

    def computeBodiesMotion(self):
        ref = self.ref
        for b in self.bodies:
            O = b.origin
            acc =  O.acc(ref.frame)
            vel =  O.vel(ref.frame)
            pos =  O.pos_from(ref.origin).express(ref.frame)
            print(pos)
            print(vel)
            print(acc)


    def computePointsMotion(self, noPosInJac=False):
        ref = self.ref
        qd  = [diff(c,dynamicsymbols._t) for c in self.q]
        qdd = [diff(diff(c,dynamicsymbols._t),dynamicsymbols._t) for c in self.q]
        self.PointsMotions=[]
        for P, pf in zip(self.Points, self.PointsFrames):
            # Express pos, vel, acc in point frame
            pos =  P.pos_from(ref.origin).to_matrix(pf)
            vel =  P.vel(ref.frame).to_matrix(pf)
            acc =  P.acc(ref.frame).to_matrix(pf)

            noAcc = [(qddi,0) for qddi in qdd]
            if noPosInJac:
                noPos = [(qi,0) for qi in self.q]
                accJac = subs_no_diff(acc,noPos)
            else:
                accJac = acc

            Ma= myjacobian(accJac, qdd)
            Ca= myjacobian(accJac.subs(noAcc), qd)
            Ka= myjacobian(accJac.subs(noAcc), self.q)

            self.PointsMotions.append((pos,vel,acc,Ma,Ca,Ka))
        return acc,accJac,Ma,Ca,Ka

    def exportPackage(self, path='', extraSubs=None, smallAngles=None, linearize=True, replaceDict=None, pathtex=None, fullPage=True, silentTimer=True):
        """ 
        Export to a python package
        - replaceDict: replace dictionary for python export
        - path: path to package. If none is provided, export to model.name
        - smallAngle: list of tuples, each item being a tuble with a list of angles, and order of lin:
                [ ([angle1, angle2], 1) ,  ([vc], 2) ]
        """
        # --- Compute Kane's equations of motion
        if self.kane is None:
            with Timer('Kanes equations', silent=silentTimer):
                self.kaneEquations(Mform='TaylorExpanded')

        # --- Use EOM
        with Timer('EOM', silent=silentTimer):
            EOM = self.to_EOM(extraSubs=extraSubs)

        #  --- Small angle approximations
        if smallAngles is not None:
            with Timer('Small angle approx', silent=silentTimer):
                for sa_list_order in smallAngles:
                    sa_list = sa_list_order[0]
                    order   = sa_list_order[1]
                    EOM.smallAngleApprox(sa_list, order=order)
            with Timer('Simplify', silent=silentTimer):
                EOM.simplify()

        # --- Separate EOM into mass matrix and forcing
        with Timer('Mass forcing term', silent=silentTimer):
            EOM.mass_forcing_form() # EOM.M and EOM.F

        # --- Linearize equation of motions 
        if linearize:
            with Timer('Linearization', silent=silentTimer):
                EOM.linearize(noAcc=True) # EOM.M0, EOM.K0, EOM.C0, EOM.B0

        # --- Points Motion
        #if len(self.Points)>0:
        #    self.computePointsMotion(noPosInJac=False)


        # --- Python export path
        if len(path)>0:
            folder = os.path.dirname(path)
            name= os.path.basename(path)
            if os.path.exists(folder):
                pass
            else:
                try:
                    os.makedirs(folder)
                    with open(os.path.join(folder,'__init__.py'),'w') as f:
                        pass
                except:
                    pass
        else:
            folder = './'
            name   = self.name

        # --- Export equations
        with Timer('Export to python', silent=silentTimer):
            outFileName=EOM.savePython(name=name, folder=folder, replaceDict=replaceDict)

        if len(self.PointsMotions)>0:
            with Timer('Export Point Motion', silent=silentTimer):
                with open(outFileName, 'a') as f:
                    for p,pm in zip(self.Points, self.PointsMotions):
                        pos,vel,acc,Ma,Ca,Ka = pm
                        s = PointAcc2Py(p, acc, self.q)
                        f.write(s)
                        s = PointAccLin2Py(p, Ma, Ca, Ka, self.q)
                        f.write(s)





        if pathtex is not None:
            folder = os.path.dirname(pathtex)
            name= os.path.basename(pathtex)
            try:
                os.makedirs(folder)
            except:
                pass
            with Timer('Export to latex', silent=silentTimer):
                EOM.saveTex(name=name, folder=folder, variables=['M','F','M0','K0','C0','B0'], fullPage=fullPage)


    def saveTex(self, name='', prefix='', suffix='', folder='./', extraSubs=None, header=True, extraHeader=None, variables=['MM','FF','M','C','K','B','MMsa','FFsa','Msa','Csa','Ksa','Bsa','body_details'], doSimplify=False, velSubs=[(0,0)], fullPage=True):
        """ 
        Save forcing and mass matrix to latex file
        """
        extraSubs = [] if extraSubs is None else extraSubs
        if len(name)==0:
            name=self.name
        name=prefix
        if len(suffix)>0:
            name=name+'_'+suffix.strip('_')

        filename=os.path.join(folder,name+'.tex')
        with Timer('Latex to {}'.format(filename),True,silent=True):
            with open(filename,'w') as f:
                if header:
                    f.write('Model: {}, \n'.format(self.name.replace('_','\_')))
                    f.write('Degrees of freedom: ${}$, \n'.format(cleantex(self.coordinates)))
                    try:
                        f.write('Small angles:       ${}$\\\\ \n'.format(cleantex(self.smallAnglesUsed)))
                    except:
                        pass
                    f.write('Free vars:       ${}$\\\\ \n'.format(cleantex(self.var)))
                if extraHeader:
                    f.write('\\clearpage\n{}\\\\\n'.format(extraHeader))
                if 'F' in variables:
                    FF = self.kane.forcing.subs(self.kdeqsSubs)
                    toTex(f, FF, label='Forcing', fullPage=fullPage, extraSubs=extraSubs, velSubs=velSubs, doSimplify=doSimplify)
                if 'M' in variables:
                    MM = self.kane.mass_matrix.subs(self.kdeqsSubs) # NOTE this is wrong
                    toTex(f, MM, label='Mass matrix', fullPage=fullPage, extraSubs=extraSubs, velSubs=velSubs, doSimplify=doSimplify)
                if 'M0' in variables:
                    toTex(f, self.M0, label='Linearized mass matrix', fullPage=fullPage, extraSubs=extraSubs, velSubs=velSubs, doSimplify=doSimplify)
                if 'C0' in variables:
                    toTex(f, self.C0, label='Linearized damping matrix', fullPage=fullPage, extraSubs=extraSubs, velSubs=velSubs, doSimplify=doSimplify)
                if 'K0' in variables:
                    toTex(f, self.K0, label='Linearized stiffness matrix', fullPage=fullPage, extraSubs=extraSubs, velSubs=velSubs, doSimplify=doSimplify)
                if 'B0' in variables:
                    toTex(f, self.B0, label='Linearized forcing matrix', fullPage=fullPage, extraSubs=extraSubs, velSubs=velSubs, doSimplify=doSimplify)
                if 'FFsa' in variables:
                    FF=subs_no_diff(self._sa_forcing,extraSubs)
                    toTex(f, FF, label='Forcing small angle', fullPage=fullPage, extraSubs=extraSubs, velSubs=velSubs, doSimplify=doSimplify)
                if 'MMsa' in variables:
                    toTex(f, self._sa_mass_matrix, label='Mass matrix small angle', fullPage=fullPage, extraSubs=extraSubs, velSubs=velSubs, doSimplify=doSimplify)
                if 'Msa' in variables:
                    toTex(f, self._sa_M, label='Linearized mass matrix small angle', fullPage=fullPage, extraSubs=extraSubs, velSubs=velSubs, doSimplify=doSimplify)
                if 'Csa' in variables:
                    toTex(f, self._sa_C, label='Linearized damping matrix small angle', fullPage=fullPage, extraSubs=extraSubs, velSubs=velSubs, doSimplify=doSimplify)
                if 'Ksa' in variables:
                    toTex(f, self._sa_K, label='Linearized stiffness matrix small angle', fullPage=fullPage, extraSubs=extraSubs, velSubs=velSubs, doSimplify=doSimplify)
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
    def toPython(self, extraSubs=None, variables=['MM','FF','MMsa','FFsa','M','C','K','B','Msa','Csa','Ksa','Bsa'], replaceDict=None, doSimplify=False, velSubs=[(0,0)]):
        """ 
        Save forcing, mass matrix and linear model to python package
        """
        extraSubs = [] if extraSubs is None else extraSubs

        # --- replace  dict
        if replaceDict is None: 
            replaceDict=OrderedDict()
        for b in self.bodies:
            if isinstance(b, YAMSFlexibleBody):
                b.replaceDict(replaceDict)
        s=''
        s+='"""\n'
        s+='{}\n'.format(self.__repr__())
        s+='"""\n'
        s+='import numpy as np\n'
        s+='from numpy import cos, sin, pi, sqrt\n'

        s += infoToPy(self.name, self.coordinates, self.var)

        if 'FF' in variables:
            forcing = self.kane.forcing.subs(self.kdeqsSubs)
            s+=forcingToPy(self.coordinates, forcing, replaceDict=replaceDict, extraSubs=extraSubs, velSubs=velSubs, doSimplify=doSimplify)
        if 'MM' in variables:
            MM = self.kane.mass_matrix.subs(self.kdeqsSubs)
            s+=MMToPy(self.coordinates, MM, replaceDict=replaceDict, extraSubs=extraSubs, velSubs=velSubs, doSimplify=doSimplify)
        if 'M0' in variables:
            #M0ToPy(f, self.q, self.M0, replaceDict=replaceDict, extraSubs=extraSubs, velSubs=velSubs, doSimplify=doSimplify)
            s+=M0ToPy(self.coordinates, self.M0, replaceDict=replaceDict, extraSubs=extraSubs, velSubs=velSubs, doSimplify=doSimplify)
        if 'C0' in variables:
            s+=C0ToPy(self.coordinates, self.C0, replaceDict=replaceDict, extraSubs=extraSubs, velSubs=velSubs, doSimplify=doSimplify)
        if 'K0' in variables:
            s+=K0ToPy(self.coordinates, self.K0, replaceDict=replaceDict, extraSubs=extraSubs, velSubs=velSubs, doSimplify=doSimplify)
        if 'B0' in variables:
            s+=B0ToPy(self.coordinates, self.B0, self.var, replaceDict=replaceDict, extraSubs=extraSubs, velSubs=velSubs, doSimplify=doSimplify)

        if 'FFsa' in variables and self._sa_forcing is not None:
            forcing = self._sa_forcing
            s+=forcingToPy(self.coordinates, forcing, replaceDict=replaceDict, extraSubs=extraSubs, velSubs=velSubs, doSimplify=doSimplify, fname='forcing_sa',extraComment='with small angle approximation')
        if 'MMsa' in variables and self._sa_mass_matrix is not None:
            MM = self._sa_mass_matrix
            s+=MMToPy(self.coordinates, MM, replaceDict=replaceDict, extraSubs=extraSubs, velSubs=velSubs, doSimplify=doSimplify, fname='mass_matrix_sa', extraComment='with small angle approximation')
        if 'Msa' in variables and self._sa_M is not None:
            s+=M0ToPy(self.coordinates, self._sa_M, replaceDict=replaceDict, extraSubs=extraSubs, velSubs=velSubs, doSimplify=doSimplify, fname='M_lin_sa', extraComment='with small angle approximation')
        if 'Csa' in variables and self._sa_K is not None:
            s+=C0ToPy(self.coordinates, self._sa_C, replaceDict=replaceDict, extraSubs=extraSubs, velSubs=velSubs, doSimplify=doSimplify, fname='C_lin_sa', extraComment='with small angle approximation')
        if 'Ksa' in variables and self._sa_K is not None:
            s+=K0ToPy(self.coordinates, self._sa_K, replaceDict=replaceDict, extraSubs=extraSubs, velSubs=velSubs, doSimplify=doSimplify, fname='K_lin_sa', extraComment='with small angle approximation')
        if 'Bsa' in variables and self._sa_B is not None:
            s+=B0ToPy(self.coordinates, self._sa_B, replaceDict=replaceDict, extraSubs=extraSubs, velSubs=velSubs, doSimplify=doSimplify, fname='B_lin_sa', extraComment='with small angle approximation')
        return s

        #print(replaceDict)
    def savePython(self, name='', prefix='', suffix='', folder='./', **kwargs):
        """ model save to python, see toPython for arguments"""
        if len(name)==0:
            name=self.name
        name=prefix+name
        if len(suffix)>0:
            name=name+'_'+suffix.strip('_')
        filename=os.path.join(folder,name+'.py')
        with Timer('Python to {}'.format(filename),True,silent=True):
            with open(filename,'w') as f:
                f.write(self.toPython(*args, **kwargs))


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

    def __init__(self, EOM, q, name, bodyReplaceDict=None):
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
        self.smallAnglesUsed=[]

        self.name=name
        self.bodyReplaceDict=bodyReplaceDict

        self.input_vars=findInputs(EOM, q)     
        #self.mass_forcing_form(self, extraSubs=None) # M and F
        #self.linearize(self, op_point, noAcc, noVel=False, extraSubs=None): # M0,K0,B0
 

    def __repr__(self):
        s='<{} object "{}" with attributes:>\n'.format(type(self).__name__,self.name)
        s+=' - name:    {}\n'.format(self.name)
        s+=' - q:       {}\n'.format(self.q)
        s+=' - input_vars: {}\n'.format(self.input_vars)
        s+=' - bodyReplaceDict: {}\n'.format(self.bodyReplaceDict)
        s+=' - smallAnglesUsed: {}\n'.format(self.smallAnglesUsed)
        s+='attributes: EOM\n'.format(self.smallAnglesUsed)
        s+='attributes: M,F          (call mass_forcing_form) \n'.format(self.smallAnglesUsed)
        s+='attributes: M0,K0,C0,B0  (call linearize) \n'.format(self.smallAnglesUsed)
        s+='methods that act on EOM in place (not M/F,M0): subs, simplify, trigsimp, expand\n'
        return s


    def subs(self, subs_list, inPlace=True):
        """ Apply substitutions to equations of motion """
        if inPlace:
            self.EOM = self.EOM.subs(subs_list)
        else:
            return self.EOM.subs(subs_list)

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

    def smallAngleApprox(self, angle_list, order=1, inPlace=True):
        """ 
        Apply small angle approximation to EOM 

        NOTE: inPlace!
        """
        with Timer('Small angle approx',True,silent=True):
            if inPlace:
                self.EOM = smallAngleApprox(self.EOM, angle_list, order=order)
                self.smallAnglesUsed+=angle_list
            else:
                return smallAngleApprox(self.EOM, angle_list, order=order)

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
        """
        Linearize EOM
        """
        op_point  = [] if op_point is None else op_point
        extraSubs = [] if extraSubs is None else extraSubs
        self.M0,self.C0,self.K0,self.B0, self.input_vars = linearizeQ(self.EOM, self.q, op_point=op_point, noAcc=noAcc, noVel=noVel, extraSubs=extraSubs)
        return self.M0, self.C0, self.K0, self.B0

    def saveTex(self, name='', prefix='', suffix='', folder='./', extraSubs=None, header=True, extraHeader=None, variables=['M','F','M0','C0','K0','B0'], doSimplify=False, velSubs=[(0,0)], fullPage=True):
        """ 
        Save EOM to a latex file
        """
        extraSubs = [] if extraSubs is None else extraSubs
        if len(name)==0:
            name=self.name
        name=prefix+name
        if len(suffix)>0:
            name=name+'_'+suffix.strip('_')

        filename=os.path.join(folder,name+'.tex')
        with Timer('Latex to {}'.format(filename),True,silent=True):
            with open(filename,'w') as f:
                if header:
                    f.write('Model: {}, \n'.format(self.name.replace('_','\_')))
                    f.write('Degrees of freedom: ${}$, \n'.format(cleantex(self.q)))
                    try:
                        f.write('Small angles:       ${}$\\\\ \n'.format(cleantex(self.smallAnglesUsed)))
                    except:
                        pass
                    f.write('Free vars:       ${}$\\\\ \n'.format(cleantex(self.input_vars)))

                if extraHeader:
                    f.write('\\clearpage\n{}\\\\\n'.format(extraHeader))

                if 'F' in variables:
                    toTex(f, self.F, label='Forcing', fullPage=fullPage, extraSubs=extraSubs, velSubs=velSubs, doSimplify=doSimplify)
                if 'M' in variables:
                    toTex(f, self.M, label='Mass matrix', fullPage=fullPage, extraSubs=extraSubs, velSubs=velSubs, doSimplify=doSimplify)
                if 'M0' in variables:
                    toTex(f, self.M0, label='Linearized mass matrix', fullPage=fullPage, extraSubs=extraSubs, velSubs=velSubs, doSimplify=doSimplify)
                if 'C0' in variables:
                    toTex(f, self.C0, label='Linearized damping matrix', fullPage=fullPage, extraSubs=extraSubs, velSubs=velSubs, doSimplify=doSimplify)
                if 'K0' in variables:
                    toTex(f, self.K0, label='Linearized stiffness matrix', fullPage=fullPage, extraSubs=extraSubs, velSubs=velSubs, doSimplify=doSimplify)
                if 'B0' in variables:
                    toTex(f, self.B0, label='Linearized forcing matrix', fullPage=fullPage, extraSubs=extraSubs, velSubs=velSubs, doSimplify=doSimplify)


    def toPython(self, extraSubs=None, replaceDict=None, doSimplify=False, velSubs=[(0,0)]):
        """ Export EOM to python string """
        extraSubs = [] if extraSubs is None else extraSubs
        # --- replace  dict
        if replaceDict is None: 
            replaceDict=OrderedDict()
        replaceDict.update(self.bodyReplaceDict)
        s=''
        s+='"""\n'
        s+='Equations of motion\n'
        s+='model name: {}\n'.format(self.name)
        s+='"""\n'
        s+='import numpy as np\n'
        s+='from numpy import cos, sin, pi, sqrt\n'

        s+=infoToPy(self.name, self.q, self.input_vars)

        if self.M is not None:
            s+=forcingToPy(self.q, self.F, replaceDict=replaceDict, extraSubs=extraSubs, velSubs=velSubs, doSimplify=doSimplify)
            s+=MMToPy(self.q, self.M, replaceDict=replaceDict, extraSubs=extraSubs, velSubs=velSubs, doSimplify=doSimplify)

        if self.M0 is not None:
            s+=M0ToPy(self.q, self.M0, replaceDict=replaceDict, extraSubs=extraSubs, velSubs=velSubs, doSimplify=doSimplify)
            s+=C0ToPy(self.q, self.C0, replaceDict=replaceDict, extraSubs=extraSubs, velSubs=velSubs, doSimplify=doSimplify)
            s+=K0ToPy(self.q, self.K0, replaceDict=replaceDict, extraSubs=extraSubs, velSubs=velSubs, doSimplify=doSimplify)
            s+=B0ToPy(self.q, self.B0, self.input_vars, replaceDict=replaceDict, extraSubs=extraSubs, velSubs=velSubs, doSimplify=doSimplify)
        return s
    
    def savePython(self, name='', prefix='', suffix='', folder='./', **kwargs):
        """ Write EOM to python. See toPython for extra arguments"""
        if len(name)==0:
            name=self.name
        name=prefix+name
        if len(suffix)>0:
            name=name+'_'+suffix.strip('_')
        filename=os.path.join(folder,name+'.py')
        with Timer('Python to {}'.format(filename),True,silent=True):
            with open(filename, 'w') as f:
                f.write(self.toPython(**kwargs))
        return filename


# --------------------------------------------------------------------------------}
# --- Main helper functions  
# --------------------------------------------------------------------------------{
def findInputs(EOM, q):
    """ Find variables within the EOM that are not """
    qd  = [qi.diff(dynamicsymbols._t) for qi in q]
    qdd = [qdi.diff(dynamicsymbols._t) for qdi in qd]
    dyn_symbols = find_dynamicsymbols(EOM)
    all_q       = q + qd + qdd
    inputs      = [s for s in dyn_symbols if s not in all_q]
    II = np.argsort([str(v) for v in inputs]) # sorting for consistency
    inputs = list(np.array(inputs)[II])
    return inputs

def linearizeQ(EOM, q, u=None, op_point=None, noAcc=True, noVel=False, extraSubs=None):
    """ Linearize the equations of motions using state "q" and derivatives
    The alternative is to use the kinematic equations using linearizer.linearize.
    """
    op_point  = [] if op_point is None else op_point
    extraSubs = [] if extraSubs is None else extraSubs
    qd  = [qi.diff(dynamicsymbols._t) for qi in q]
    qdd = [qdi.diff(dynamicsymbols._t) for qdi in qd]

    # NOTE: order important
    op_point0=[]
    if noAcc: 
        op_point0=[(qddi,0) for qddi in qdd]
    if noVel: 
        op_point0=[(qdi,0) for qdi in qd]
    op_point= op_point0+op_point # order might matter
    print('>>> TODO sort op point so that diff wrt time are first, or do the trick with symbols')
    # use if isinstance sympy.core.function.Derivative

    # --- Inputs are dynamic symbols that are not coordinates
    if u is None:
        u = findInputs(EOM, q)
    # KEEP ME Alternative
    #M, A, B = linearizer.linearize(op_point=op_point ) #o

    #print('>>>> op_point',op_point)
    #print('>>>> extraSubs',extraSubs)
    M =-EOM.jacobian(qdd).subs(extraSubs).subs(op_point)
    C =-EOM.jacobian(qd ).subs(extraSubs).subs(op_point)
    K =-EOM.jacobian(q  ).subs(extraSubs).subs(op_point)
    #M =-EOM.jacobian(qdd).subs(extraSubs)
    #C =-EOM.jacobian(qd ).subs(extraSubs)
    #K =-EOM.jacobian(q  ).subs(extraSubs)
    #M =subs_no_diff(M, op_point)
    #C =subs_no_diff(C, op_point)
    #K =subs_no_diff(K, op_point)
    if len(u)>0:
        B = EOM.jacobian(u).subs(extraSubs).subs(op_point)
    else:
        B=Matrix([])
    return M,C,K,B,u


def forcingToPy(q, forcing, replaceDict=None, extraSubs=None, velSubs=[(0,0)], doSimplify=False, fname='forcing', extraComment=''):
    extraSubs = [] if extraSubs is None else extraSubs
    forcing = subs_no_diff(forcing, extraSubs)
    forcing = forcing.subs(velSubs)
    if doSimplify:
        forcing=trigsimp(forcing)

    s=''
    s0, params, inputs, sdofs  = cleanPy(forcing, varname='FF', dofs = q, indent=4, replDict=replaceDict)
    s += 'def {}(t,q=None,qd=None,p=None,u=None,z=None):\n'.format(fname)
    s += '    """ Non linear mass forcing{} \n'.format(extraComment)
    s += '    q:  degrees of freedom, array-like: {}\n'.format(sdofs)
    s += '    qd: dof velocities, array-like\n'
    s += '    p:  parameters, dictionary with keys: {}\n'.format(params)
    s += '    u:  inputs, dictionary with keys: {}\n'.format(inputs)
    s += '           where each values is a function of time\n'
    s += '    """\n'
    s += '    if z is not None:\n'
    s += '        q  = z[0:int(len(z)/2)] \n'
    s += '        qd = z[int(len(z)/2): ] \n'
    s += s0
    s += '    return FF\n\n'
    return s


def MMToPy(q, MM, replaceDict=None, extraSubs=None, velSubs=[(0,0)], doSimplify=False, fname='mass_matrix', extraComment=''):
    extraSubs = [] if extraSubs is None else extraSubs
    MM = subs_no_diff(MM,extraSubs)
    MM = MM.subs(velSubs)
    if doSimplify:
        MM=trigsimp(MM)
        #MM.simplify()
    s=''
    s0, params, inputs, sdofs  = cleanPy(MM, varname='MM', dofs = q, indent=4, replDict=replaceDict)
    s += 'def {}(q=None,p=None,z=None):\n'.format(fname)
    s += '    """ Non linear mass matrix {}\n'.format(extraComment)
    s += '     q:  degrees of freedom, array-like: {}\n'.format(sdofs)
    s += '     p:  parameters, dictionary with keys: {}\n'.format(params)
    s += '    """\n'
    s += '    if z is not None:\n'
    s += '        q  = z[0:int(len(z)/2)] \n'
    s += s0
    s += '    return MM\n\n'
    return s


def M0ToPy(q, MM, replaceDict=None, extraSubs=None, velSubs=[(0,0)], doSimplify=False, fname='M_lin', extraComment=''):
    extraSubs = [] if extraSubs is None else extraSubs
    MM=subs_no_diff(MM,extraSubs)
    MM = MM.subs(velSubs)
    if doSimplify:
        MM.simplify()
    s=''
    s0, params, inputs, sdofs  = cleanPy(MM, varname='MM', dofs = q, indent=4, replDict=replaceDict, noTimeDep=True)
    s+='def {}(q=None,p=None,z=None):\n'.format(fname)
    s+='    """ Linear mass matrix {}\n'.format(extraComment)
    s+='    q:  degrees of freedom at operating point, array-like: {}\n'.format(sdofs)
    s+='    p:  parameters, dictionary with keys: {}\n'.format(params)
    s+='    """\n'
    s+='    if z is not None:\n'
    s+='        q  = z[0:int(len(z)/2)] \n'
    s+=s0
    s+='    return MM\n\n'
    return s

def C0ToPy(q, MM, replaceDict=None, extraSubs=None, velSubs=[(0,0)], doSimplify=False, fname='C_lin', extraComment=''):
    extraSubs = [] if extraSubs is None else extraSubs
    MM=subs_no_diff(MM,extraSubs)
    MM = MM.subs(velSubs)
    if doSimplify:
        MM.simplify()
    s=''
    s0, params, inputs, sdofs  = cleanPy(MM, varname='CC', dofs = q, indent=4, replDict=replaceDict, noTimeDep=True)
    s += 'def {}(q=None,qd=None,p=None,u=None,z=None):\n'.format(fname)
    s += '    """ Linear damping matrix {}\n'.format(extraComment)
    s += '    q:  degrees of freedom at operating point, array-like: {}\n'.format(sdofs)
    s += '    qd: dof velocities at operating point, array-like\n'
    s += '    p:  parameters, dictionary with keys: {}\n'.format(params)
    s += '    u:  inputs at operating point, dictionary with keys: {}\n'.format(inputs)
    s += '           where each values is a constant!\n'
    s += '    """\n'
    s += '    if z is not None:\n'
    s += '        q  = z[0:int(len(z)/2)] \n'
    s += '        qd = z[int(len(z)/2): ] \n'
    s += s0
    s += '    return CC\n\n'
    return s

def K0ToPy(q, MM, replaceDict=None, extraSubs=None, velSubs=[(0,0)], doSimplify=False, fname='K_lin', extraComment=''):
    extraSubs = [] if extraSubs is None else extraSubs
    MM=subs_no_diff(MM,extraSubs)
    MM = MM.subs(velSubs)
    if doSimplify:
        MM.simplify()
    s=''
    s0, params, inputs, sdofs  = cleanPy(MM, varname='KK', dofs = q, indent=4, replDict=replaceDict, noTimeDep=True)
    s +='def {}(q=None,qd=None,p=None,u=None,z=None):\n'.format(fname)
    s +='    """ Linear stiffness matrix {}\n'.format(extraComment)
    s +='    q:  degrees of freedom, array-like: {}\n'.format(sdofs)
    s +='    qd: dof velocities, array-like\n'
    s +='    p:  parameters, dictionary with keys: {}\n'.format(params)
    s +='    u:  inputs at operating point, dictionary with keys: {}\n'.format(inputs)
    s +='           where each values is a constant!\n'
    s +='    """\n'
    s +='    if z is not None:\n'
    s +='        q  = z[0:int(len(z)/2)] \n'
    s +='        qd = z[int(len(z)/2): ] \n'
    s +=s0
    s +='    return KK\n\n'
    return s


def B0ToPy(q, MM, input_vars, replaceDict=None, extraSubs=None, velSubs=[(0,0)], doSimplify=False, fname='B_lin', extraComment=''):
    extraSubs = [] if extraSubs is None else extraSubs
    MM=subs_no_diff(MM,extraSubs)
    MM = MM.subs(velSubs)
    if doSimplify:
        MM.simplify()
    s=''
    s0, params, inputs, sdofs  = cleanPy(MM, varname='BB', dofs = q, indent=4, replDict=replaceDict, noTimeDep=True)
    s += 'def {}(q=None,qd=None,p=None,u=None):\n'.format(fname)
    s += '    """ Linear mass matrix {}\n'.format(extraComment)
    s += '    q:  degrees of freedom at operating point, array-like: {}\n'.format(sdofs)
    s += '    qd: dof velocities at operating point, array-like\n'
    s += '    p:  parameters, dictionary with keys: {}\n'.format(params)
    s += '    u:  inputs at operating point, dictionary with keys: {}\n'.format(inputs)
    s += '           where each values is a constant!\n'
    s += '    The columns of B correspond to:   {}\\\\ \n'.format(input_vars)
    s += '    """\n'
    s += s0
    s += '    return BB\n\n'
    return s

def infoToPy(name, q, u):
    s =  'def info():\n'
    s += '    """ Return information about current model present in this package """\n'
    s += '    I=dict()\n'
    s += '    I[\'name\']=\'{}\'\n'.format(name)
    s += '    I[\'nq\']={}\n'.format(len(q))
    s += '    I[\'nu\']={}\n'.format(len(u))
    s += '    I[\'sq\']=[{}]\n'.format(','.join(['\''+repr(qi).replace('(t)','')+'\'' for qi in q]))
    s += '    I[\'su\']=[{}]\n'.format(','.join(['\''+repr(ui).replace('(t)','')+'\'' for ui in u]))
    s += '    return I\n\n'
    return s

def PointAcc2Py(P, acc, q, replaceDict=None, extraSubs=None, velSubs=[(0,0)], doSimplify=False, fname='forcing', extraComment=''):
    extraSubs = [] if extraSubs is None else extraSubs
    acc = subs_no_diff(acc, extraSubs)
    acc = acc.subs(velSubs)
    if doSimplify:
        acc=trigsimp(acc)
    s0, params, inputs, sdofs  = cleanPy(acc, varname='acc', dofs = q, indent=4, replDict=replaceDict)
    Pname=str(P)
    s=''
    s += 'def Acceleration{}(q=None,qd=None,qdd=None,p=None):\n'.format(Pname)
    s += '    """ Acceleration of point {} {} \n'.format(Pname,extraComment)
    s += '    q:  degrees of freedom, array-like: {}\n'.format(sdofs)
    s += '    qd: dof velocities, array-like\n'
    s += '    qdd:dof accelerations, array-like\n'
    s += '    p:  parameters, dictionary with keys: {}\n'.format(params)
    s += '    """\n'
    s += s0
    s += '    return acc\n\n'
    return s

def PointAccLin2Py(P, Ma, Ca, Ka, q, replaceDict=None, extraSubs=None, velSubs=[(0,0)], doSimplify=False, fname='forcing', extraComment=''):
    sMa, params, inputs, sdofs  = cleanPy(Ma, varname='Ma', dofs = q, indent=4, replDict=replaceDict)
    sCa, params, inputs, sdofs  = cleanPy(Ca, varname='Ca', dofs = q, indent=4, replDict=replaceDict)
    sKa, params, inputs, sdofs  = cleanPy(Ca, varname='Ka', dofs = q, indent=4, replDict=replaceDict)
    Pname=str(P)
    s=''
    s += 'def AccLin{}(q=None,qd=None,p=None):\n'.format(Pname)
    s += '    """ Acceleration of point {} {} \n'.format(Pname,extraComment)
    s += '    q:  degrees of freedom, array-like: {}\n'.format(sdofs)
    s += '    qd: dof velocities, array-like\n'
    s += '    p:  parameters, dictionary with keys: {}\n'.format(params)
    s += '    """\n'
    s += sMa
    s += sCa
    s += sKa
    s += '    return Ma, Ca, Ka\n\n'
    return s



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

