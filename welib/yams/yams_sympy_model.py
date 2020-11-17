
from .yams_sympy_tools import smallAngleApprox 
from .yams_kane import YAMSKanesMethod
# --------------------------------------------------------------------------------}
# ---  
# --------------------------------------------------------------------------------{
class YAMSModel(object):
    def __init__(self):
        self.ref         = None
        self.coordinates = None
        self.speeds      = None
        self.kdeqs       = None
        self.kdeqsSubs   = None
        self.bodies      = None
        self.loads       = None
        self.kane        = None
        self.forcing     = None
        self.mass_matrix = None

    def kaneEquations(self, Mform='symbolic'):
        from welib.tools.tictoc import Timer
        for sa in ['ref', 'coordinates', 'speeds','kdeqs','bodies','loads']:
            if getattr(self,sa) is None:
                raise Exception('Attribute {} needs to be set before calling `kane` method'.format(sa))

        with Timer('Kane step1',True):
            self.kane = YAMSKanesMethod(self.ref.frame, self.coordinates, self.speeds, self.kdeqs)

        # --- Expensive kane step
        with Timer('Kane step 2',True):
            #(use  Mform ='symbolic' or 'TaylorExpanded'), Mform='symbolic'
            self.fr, self.frstar  = self.kane.kanes_equations(self.bodies, self.loads)
        self.kane.fr     = self.fr
        self.kane.frstar = self.frstar

    def smallAngleApprox(self, angle_list, extraSubs=[]):
        """ 
        Apply small angle approximation
        """
        from welib.tools.tictoc import Timer
        # Forcing
        with Timer('Small angle approx. forcing',True):
            if self.forcing is None:
                self.forcing=self.kane.forcing
            self.forcing = self.forcing.subs(self.kdeqsSubs).subs(extraSubs)
            self.forcing = smallAngleApprox(self.forcing, angle_list).subs(extraSubs)
            self.forcing.simplify()
        # Mass matrix
        with Timer('Small angle approx. mass matrix',True):
            if self.mass_matrix is None:
                self.mass_matrix=self.kane.mass_matrix
            self.mass_matrix = self.mass_matrix.subs(self.kdeqsSubs).subs(extraSubs)
            self.mass_matrix = smallAngleApprox(self.mass_matrix, angle_list).subs(extraSubs)
            self.mass_matrix.simplify()
