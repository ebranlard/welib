from __future__ import print_function, division

from sympy.core.backend import zeros, Matrix, diff, eye
from sympy import solve_linear_system_LU
# try:
#     from sympy.core.compatibility import range
# except:
#     pass
try:
    from sympy.core.sorting import default_sort_key
except:
    from sympy.utilities import default_sort_key
from sympy.physics.vector import (ReferenceFrame, dynamicsymbols)
                             
from sympy.physics.vector import partial_velocity
from sympy.physics.mechanics.particle import Particle
from sympy.physics.mechanics.rigidbody import RigidBody
from sympy.physics.mechanics.functions import (msubs, find_dynamicsymbols,
                                               _f_list_parser)

from sympy.physics.mechanics import inertia
from sympy.physics.mechanics import Point
from sympy.physics.mechanics.linearize import Linearizer
from sympy.utilities.exceptions import SymPyDeprecationWarning
from sympy.utilities.iterables import iterable

from welib.yams.yams_sympy import YAMSFlexibleBody, YAMSRigidBody, coord2vec


__all__ = ['YAMSKanesMethod']
__all__+= ['kane_frstar','kane_fr','kane_fr_alt','kane_frstar_alt']


class YAMSKanesMethod(object):
    """Kane's method object.

    Explanation
    ===========

    This object is used to do the "book-keeping" as you go through and form
    equations of motion in the way Kane presents in:
    Kane, T., Levinson, D. Dynamics Theory and Applications. 1985 McGraw-Hill

    The attributes are for equations in the form [M] udot = forcing.

    Attributes
    ==========

    q, u : Matrix
        Matrices of the generalized coordinates and speeds
    bodylist : iterable
        Iterable of Point and RigidBody objects in the system.
    forcelist : iterable
        Iterable of (Point, vector) or (ReferenceFrame, vector) tuples
        describing the forces on the system.
    auxiliary : Matrix
        If applicable, the set of auxiliary Kane's
        equations used to solve for non-contributing
        forces.
    mass_matrix : Matrix
        The system's mass matrix
    forcing : Matrix
        The system's forcing vector
    mass_matrix_full : Matrix
        The "mass matrix" for the u's and q's
    forcing_full : Matrix
        The "forcing vector" for the u's and q's

    Examples
    ========

    This is a simple example for a one degree of freedom translational
    spring-mass-damper.

    In this example, we first need to do the kinematics.
    This involves creating generalized speeds and coordinates and their
    derivatives.
    Then we create a point and set its velocity in a frame.

        >>> from sympy import symbols
        >>> from sympy.physics.mechanics import dynamicsymbols, ReferenceFrame
        >>> from sympy.physics.mechanics import Point, Particle, KanesMethod
        >>> q, u = dynamicsymbols('q u')
        >>> qd, ud = dynamicsymbols('q u', 1)
        >>> m, c, k = symbols('m c k')
        >>> N = ReferenceFrame('N')
        >>> P = Point('P')
        >>> P.set_vel(N, u * N.x)

    Next we need to arrange/store information in the way that KanesMethod
    requires.  The kinematic differential equations need to be stored in a
    dict.  A list of forces/torques must be constructed, where each entry in
    the list is a (Point, Vector) or (ReferenceFrame, Vector) tuple, where the
    Vectors represent the Force or Torque.
    Next a particle needs to be created, and it needs to have a point and mass
    assigned to it.
    Finally, a list of all bodies and particles needs to be created.

        >>> kd = [qd - u]
        >>> FL = [(P, (-k * q - c * u) * N.x)]
        >>> pa = Particle('pa', P, m)
        >>> BL = [pa]

    Finally we can generate the equations of motion.
    First we create the KanesMethod object and supply an inertial frame,
    coordinates, generalized speeds, and the kinematic differential equations.
    Additional quantities such as configuration and motion constraints,
    dependent coordinates and speeds, and auxiliary speeds are also supplied
    here (see the online documentation).
    Next we form FR* and FR to complete: Fr + Fr* = 0.
    We have the equations of motion at this point.
    It makes sense to rearrange them though, so we calculate the mass matrix and
    the forcing terms, for E.o.M. in the form: [MM] udot = forcing, where MM is
    the mass matrix, udot is a vector of the time derivatives of the
    generalized speeds, and forcing is a vector representing "forcing" terms.

        >>> KM = KanesMethod(N, q_ind=[q], u_ind=[u], kd_eqs=kd)
        >>> (fr, frstar) = KM.kanes_equations(BL, FL)
        >>> MM = KM.mass_matrix
        >>> forcing = KM.forcing
        >>> rhs = MM.inv() * forcing
        >>> rhs
        Matrix([[(-c*u(t) - k*q(t))/m]])
        >>> KM.linearize(A_and_B=True)[0]
        Matrix([
        [   0,    1],
        [-k/m, -c/m]])

    Please look at the documentation pages for more information on how to
    perform linearization and how to deal with dependent coordinates & speeds,
    and how do deal with bringing non-contributing forces into evidence.

    """

    def __init__(self, frame, q_ind, u_ind, kd_eqs=None, q_dependent=None,
            configuration_constraints=None, u_dependent=None,
            velocity_constraints=None, acceleration_constraints=None,
            u_auxiliary=None, bodies=None, forcelist=None):

        """Please read the online documentation. """
        if not q_ind:
            q_ind = [dynamicsymbols('dummy_q')]
            kd_eqs = [dynamicsymbols('dummy_kd')]

        if not isinstance(frame, ReferenceFrame):
            raise TypeError('An inertial ReferenceFrame must be supplied')
        self._inertial = frame

        self._fr = None
        self._frstar = None

        self._forcelist = forcelist
        self._bodylist = bodies

        self._initialize_vectors(q_ind, q_dependent, u_ind, u_dependent,
                u_auxiliary)
        self._initialize_kindiffeq_matrices(kd_eqs)
        self._initialize_constraint_matrices(configuration_constraints,
                velocity_constraints, acceleration_constraints)

    def _initialize_vectors(self, q_ind, q_dep, u_ind, u_dep, u_aux):
        """Initialize the coordinate and speed vectors."""

        none_handler = lambda x: Matrix(x) if x else Matrix()

        # Initialize generalized coordinates
        q_dep = none_handler(q_dep)
        if not iterable(q_ind):
            raise TypeError('Generalized coordinates must be an iterable.')
        if not iterable(q_dep):
            raise TypeError('Dependent coordinates must be an iterable.')
        q_ind = Matrix(q_ind)
        self._qdep = q_dep
        self._q = Matrix([q_ind, q_dep])
        self._qdot = self.q.diff(dynamicsymbols._t)

        # Initialize generalized speeds
        u_dep = none_handler(u_dep)
        if not iterable(u_ind):
            raise TypeError('Generalized speeds must be an iterable.')
        if not iterable(u_dep):
            raise TypeError('Dependent speeds must be an iterable.')
        u_ind = Matrix(u_ind)
        self._udep = u_dep
        self._u = Matrix([u_ind, u_dep])
        self._udot = self.u.diff(dynamicsymbols._t)
        self._uaux = none_handler(u_aux)

    def _initialize_constraint_matrices(self, config, vel, acc):
        """Initializes constraint matrices."""

        # Define vector dimensions
        o = len(self.u)
        m = len(self._udep)
        p = o - m
        none_handler = lambda x: Matrix(x) if x else Matrix()

        # Initialize configuration constraints
        config = none_handler(config)
        if len(self._qdep) != len(config):
            raise ValueError('There must be an equal number of dependent '
                             'coordinates and configuration constraints.')
        self._f_h = none_handler(config)

        # Initialize velocity and acceleration constraints
        vel = none_handler(vel)
        acc = none_handler(acc)
        if len(vel) != m:
            raise ValueError('There must be an equal number of dependent '
                             'speeds and velocity constraints.')
        if acc and (len(acc) != m):
            raise ValueError('There must be an equal number of dependent '
                             'speeds and acceleration constraints.')
        if vel:
            u_zero = dict((i, 0) for i in self.u)
            udot_zero = dict((i, 0) for i in self._udot)

            # When calling kanes_equations, another class instance will be
            # created if auxiliary u's are present. In this case, the
            # computation of kinetic differential equation matrices will be
            # skipped as this was computed during the original KanesMethod
            # object, and the qd_u_map will not be available.
            if self._qdot_u_map is not None:
                vel = msubs(vel, self._qdot_u_map)

            self._f_nh = msubs(vel, u_zero)
            self._k_nh = (vel - self._f_nh).jacobian(self.u)
            # If no acceleration constraints given, calculate them.
            if not acc:
                self._f_dnh = (self._k_nh.diff(dynamicsymbols._t) * self.u +
                               self._f_nh.diff(dynamicsymbols._t))
                self._k_dnh = self._k_nh
            else:
                if self._qdot_u_map is not None:
                    acc = msubs(acc, self._qdot_u_map)
                self._f_dnh = msubs(acc, udot_zero)
                self._k_dnh = (acc - self._f_dnh).jacobian(self._udot)

            # Form of non-holonomic constraints is B*u + C = 0.
            # We partition B into independent and dependent columns:
            # Ars is then -B_dep.inv() * B_ind, and it relates dependent speeds
            # to independent speeds as: udep = Ars*uind, neglecting the C term.
            B_ind = self._k_nh[:, :p]
            B_dep = self._k_nh[:, p:o]
            self._Ars = -B_dep.LUsolve(B_ind)
        else:
            self._f_nh = Matrix()
            self._k_nh = Matrix()
            self._f_dnh = Matrix()
            self._k_dnh = Matrix()
            self._Ars = Matrix()

    def _initialize_kindiffeq_matrices(self, kdeqs):
        """Initialize the kinematic differential equation matrices.

        Parameters
        ==========
        kdeqs : sequence of sympy expressions
            Kinematic differential equations in the form of f(u,q',q,t) where
            f() = 0. The equations have to be linear in the generalized
            coordinates and generalized speeds.

        """

        if kdeqs:
            if len(self.q) != len(kdeqs):
                raise ValueError('There must be an equal number of kinematic '
                                 'differential equations and coordinates.')

            u = self.u
            qdot = self._qdot

            kdeqs = Matrix(kdeqs)
            # Dictionaries setting things to zero
            u_zero = dict((i, 0) for i in u)
            uaux_zero = dict((i, 0) for i in self._uaux)
            qdot_zero = dict((i, 0) for i in qdot)

            f_k = msubs(kdeqs, u_zero, qdot_zero)
            k_ku = (msubs(kdeqs, qdot_zero) - f_k).jacobian(u)
            k_kqdot = (msubs(kdeqs, u_zero) - f_k).jacobian(qdot)

            f_k = k_kqdot.LUsolve(f_k)
            k_ku = k_kqdot.LUsolve(k_ku)
            k_kqdot = eye(len(qdot))

            self._qdot_u_map = solve_linear_system_LU(
                    Matrix([k_kqdot.T, -(k_ku * u + f_k).T]).T, qdot)

            self._f_k = msubs(f_k, uaux_zero)
            self._k_ku = msubs(k_ku, uaux_zero)
            self._k_kqdot = k_kqdot
        else:
            self._qdot_u_map = None
            self._f_k = Matrix()
            self._k_ku = Matrix()
            self._k_kqdot = Matrix()

    def _form_fr(self, fl):
        """Form the generalized active force."""
        if fl is not None and (len(fl) == 0 or not iterable(fl)):
            raise ValueError('Force pairs must be supplied in an '
                'non-empty iterable or None.')
        def _my_f_list_parser(fl, ref_frame):
            """Parses the provided forcelist composed of items of the form (obj, force).
            Returns a tuple containing:
                vel_list: The velocity (ang_vel for Frames, vel for Points) in the provided reference frame.
                f_list: The forces.
            Used internally in the KanesMethod and LagrangesMethod classes.
            """
            genF=[]
            def flist_iter():
                for pair in fl:
                    obj, force = pair
                    if isinstance(obj, ReferenceFrame):
                        yield obj.ang_vel_in(ref_frame), force
                    elif isinstance(obj, Point):
                        yield obj.vel(ref_frame), force
                    elif isinstance(obj, YAMSFlexibleBody):
                        genF.append(force) # TODO
                    else:
                        raise TypeError('First entry in each forcelist pair must be a point or frame.')
            if not fl:
                vel_list, f_list = (), ()
            else:
                unzip = lambda l: list(zip(*l)) if l[0] else [(), ()]
                vel_list, f_list = unzip(list(flist_iter()))
            return vel_list, f_list, genF


        N = self._inertial
        # pull out relevant velocities for constructing partial velocities
        vel_list, f_list, genF = _my_f_list_parser(fl, N)
        vel_list = [msubs(i, self._qdot_u_map) for i in vel_list]
        f_list = [msubs(i, self._qdot_u_map) for i in f_list]

        # Fill Fr with dot product of partial velocities and forces
        o = len(self.u)
        b = len(f_list)
        FR = zeros(o, 1)
        partials = partial_velocity(vel_list, self.u, N)
        self._fr_products=[]
        for i in range(o):
            for j in range(b):
                self._fr_products.append((partials[j][i] , f_list[j]))
            FR[i] = sum(partials[j][i] & f_list[j] for j in range(b))

        # In case there are dependent speeds
        if self._udep:
            p = o - len(self._udep)
            FRtilde = FR[:p, 0]
            FRold = FR[p:o, 0]
            FRtilde += self._Ars.T * FRold
            FR = FRtilde

        self._forcelist = fl
        self._fr = FR
        return FR

    def _form_frstar(self, bl, Mform='TaylorExpanded', addGravity=True, g_vect=None):
        """Form the generalized inertia force.

        Mform: which form to use for flexible body mass matrix: TaylorExpanded or symbolic
        
        """

        if not iterable(bl):
            raise TypeError('Bodies must be supplied in an iterable.')

        t = dynamicsymbols._t
        N = self._inertial
        # Dicts setting things to zero
        udot_zero = dict((i, 0) for i in self._udot)
        uaux_zero = dict((i, 0) for i in self._uaux)
        uauxdot = [diff(i, t) for i in self._uaux]
        uauxdot_zero = dict((i, 0) for i in uauxdot)
        # Dictionary of q' and q'' to u and u'
        q_ddot_u_map = dict((k.diff(t), v.diff(t)) for (k, v) in
                self._qdot_u_map.items())
        q_ddot_u_map.update(self._qdot_u_map)

        # Compute fr_star in two components:
        # fr_star = -(MM*u' + nonMM)
        o = len(self.u)
        MM = zeros(o, o)
        nonMM = zeros(o, 1)
        zero_uaux = lambda expr: msubs(expr, uaux_zero)
        zero_udot_uaux = lambda expr: msubs(msubs(expr, udot_zero), uaux_zero)
        for i, body in enumerate(bl):
            bodyMM    = zeros(o, o)
            bodynonMM = zeros(o, 1)

            M = zero_uaux(body.mass)
            if isinstance(body, YAMSFlexibleBody):
                P = body.origin
                I = None
            else:
                P = body.masscenter
                I = zero_uaux(body.central_inertia)
            # --- Step 2: Positions and orientation
            #r = P.pos_from(O_E)
            #R = N.dcm(body.frame) # from body to inertial
            # --- Step 3/4: Velocities and accelerations
            omega = body.frame.ang_vel_in(N)
            vel = zero_uaux(P.vel(N))
            acc = zero_udot_uaux(P.acc(N))
            if not isinstance(body, Particle):
                omega = zero_uaux(body.frame.ang_vel_in(N))
                alpha = omega.diff(t, N)
            # NOTE: Keep me Alternative: vel from r and omega from identification:
            #vel_drdt = r.diff(t, N).simplify()
            #OmSkew = (R.diff(t) *  R.transpose()).simplify()
            #omega_ident = OmSkew[2,1] * N.x + OmSkew[0,2]*N.y + OmSkew[1,0] * N.z
            #acc_drdt2 = r.diff(t, N).diff(t,N).simplify()
            # --- Step 5 Partial velocities
            # Method 1: use "partial_velocity" function, which returns a vector
            vel_sub   = msubs(vel, self._qdot_u_map)
            Jv_vect   = partial_velocity([vel_sub], self.u, N)[0]
            if not isinstance(body, Particle):
                omega_sub = msubs(omega, self._qdot_u_map)
                Jo_vect = partial_velocity([omega_sub], self.u, N)[0]
            # NOTE: Keep me: Method 2: express everything in ref frame, and use "jacobians"
            #v  = vel.subs(Omega_Subs).to_matrix(N)
            #om = omega.subs(Omega_Subs).to_matrix(N)
            #Jv = v.jacobian(qspeeds)
            #Jo = om.jacobian(qspeeds)

            if isinstance(body,YAMSRigidBody) or isinstance(body, RigidBody):
                # --- Mass Matrix 
                for j in range(o):
                    tmp_vel = Jv_vect[j]      # Jv[:,j]
                    tmp_ang = I & Jo_vect[j]  # Jo[:,j]
                    for k in range(o):
                        # translational
                        bodyMM[j, k] += M * (tmp_vel & Jv_vect[k]) # M * Jv[:,j] dot Jv[:,k]
                        # rotational
                        bodyMM[j, k] +=     (tmp_ang & Jo_vect[k]) # I dot Jo[:,j] dot Jo[:,k]
                # --- Full inertial loads
                inertial_force = (M.diff(t) * vel + M * acc)
                inertial_torque = zero_uaux((I.dt(body.frame) & omega) + msubs(I & body.frame.ang_acc_in(N), udot_zero) + (omega ^ (I & omega)))

                # NOTE KEEP ME: Alternative formulation using "matrices" 
                #inertial_force = inertial_force.subs(Omega_Subs)
                #RIRt  = R*I.to_matrix(body.frame)*R.transpose()
                #inertial_torque_2 = - RIRt * alpha.to_matrix(N) \
                #        - coord2vec(om, N).cross( coord2vec(RIRt *om, N)).to_matrix(N)
                # Computing generatlized force Jv.f + Jo*M
                for j in range(o):
                    bodynonMM[j] += inertial_force  & Jv_vect[j]
                    bodynonMM[j] += inertial_torque & Jo_vect[j]
            elif isinstance(body,YAMSFlexibleBody):
                if addGravity:
                    acc=acc-g_vect
                MMloc = body.bodyMassMatrix(form=Mform)
                body.h_omega = body.bodyQuadraticForce(omega.to_matrix(body.frame), body.q, body.qdot)
                body.h_elast = body.bodyElasticForce(body.q, body.qdot)
                #body.h_g     = body.bodyGravitationalForce(g_vect, body.q, form='TaylorExpanded')
                inertial_force=0 # Fstar
                inertial_torque=0 # Tstar
                inertial_force_coord =MMloc[0:3,0:3] * acc.to_matrix(body.frame) 
                inertial_force_coord+=MMloc[0:3,3:6] * alpha.to_matrix(body.frame)
                inertial_force_coord+=MMloc[0:3,6:]  * Matrix(body.qddot)
                inertial_force_coord+=body.h_omega[0:3,0]
                inertial_force_coord+=body.h_elast[0:3,0]
                inertial_torque_coord =MMloc[3:6,0:3] * acc.to_matrix(body.frame) 
                inertial_torque_coord+=MMloc[3:6,3:6] * alpha.to_matrix(body.frame)
                inertial_torque_coord+=MMloc[3:6,6:]  * Matrix(body.qddot)
                inertial_torque_coord+=body.h_omega[3:6,0]
                inertial_torque_coord+=body.h_elast[3:6,0]
                inertial_elast_coord =MMloc[6:,0:3] * acc.to_matrix(body.frame) 
                inertial_elast_coord+=MMloc[6:,3:6] * alpha.to_matrix(body.frame)
                inertial_elast_coord+=MMloc[6:,6:]  * Matrix(body.qddot)
                inertial_elast_coord+=body.h_omega[6:,0]
                inertial_elast_coord+=body.h_elast[6:,0]
                body.inertial_elast=inertial_elast_coord

                inertial_force  = coord2vec(inertial_force_coord,body.frame) 
                inertial_torque = coord2vec(inertial_torque_coord,body.frame) 

                # Computing generatlized force Jv.f + Jo*M
                for j in range(o):
                    bodynonMM[j] += inertial_force  & Jv_vect[j]
                    bodynonMM[j] += inertial_torque & Jo_vect[j]
                    for k in range(len(body.q)):
                        if self.q[j] == body.q[k]:
                            bodynonMM[j] +=  inertial_elast_coord[k]
                bnMMSubs = msubs(bodynonMM, q_ddot_u_map)
                bodyMM = bnMMSubs.jacobian(self._udot)
            else:
                M = zero_uaux(body.mass)
                vel = zero_uaux(body.point.vel(N))
                acc = zero_udot_uaux(body.point.acc(N))
                vel_sub   = msubs(vel, qdot_u_map)
                Jv_vect   = partial_velocity([vel_sub], u, N)[0]
                inertial_force = (M.diff(t) * vel + M * acc)
                for j in range(o):
                    tmp_vel = Jv_vect[j]      # Jv[:,j]
                    for k in range(o):
                        bodyMM[j, k] += M * (tmp_vel & Jv_vect[k]) # M * Jv[:,j] dot Jv[:,k]
                    bodynonMM[j] += inertial_force & Jv_vect[j]
            # Perform important substitution and store body contributions
            body.MM      = msubs(bodyMM, q_ddot_u_map)
            body.nonMM   = bodynonMM

            # Cumulative MM and nonMM over all bodies
            MM   +=bodyMM
            nonMM+=bodynonMM
            # --- Storing for debug
            body._acc             = acc
            body._vel             = vel
            body._omega           = omega
            body._inertial_force  = inertial_force
            body._inertial_torque = inertial_torque
            body._Jv_vect=Jv_vect
            body._Jo_vect=Jo_vect

        # Compose fr_star out of MM and nonMM
        MM = zero_uaux(msubs(MM, q_ddot_u_map))
        nonMM = msubs(msubs(nonMM, q_ddot_u_map),
                udot_zero, uauxdot_zero, uaux_zero)
        fr_star = -(MM * msubs(Matrix(self._udot), uauxdot_zero) + nonMM)

        # If there are dependent speeds, we need to find fr_star_tilde
        if self._udep:
            p = o - len(self._udep)
            fr_star_ind = fr_star[:p, 0]
            fr_star_dep = fr_star[p:o, 0]
            fr_star = fr_star_ind + (self._Ars.T * fr_star_dep)
            # Apply the same to MM
            MMi = MM[:p, :]
            MMd = MM[p:o, :]
            MM = MMi + (self._Ars.T * MMd)

        self._bodylist = bl
        self._frstar = fr_star
        self._k_d = MM
        self._f_d = -msubs(self._fr + self._frstar, udot_zero)
        return fr_star

    def to_linearizer(self):
        """Returns an instance of the Linearizer class, initiated from the
        data in the KanesMethod class. This may be more desirable than using
        the linearize class method, as the Linearizer object will allow more
        efficient recalculation (i.e. about varying operating points)."""

        if (self._fr is None) or (self._frstar is None):
            raise ValueError('Need to compute Fr, Fr* first.')

        # Get required equation components. The Kane's method class breaks
        # these into pieces. Need to reassemble
        f_c = self._f_h
        if self._f_nh and self._k_nh:
            f_v = self._f_nh + self._k_nh*Matrix(self.u)
        else:
            f_v = Matrix()
        if self._f_dnh and self._k_dnh:
            f_a = self._f_dnh + self._k_dnh*Matrix(self._udot)
        else:
            f_a = Matrix()
        # Dicts to sub to zero, for splitting up expressions
        u_zero = dict((i, 0) for i in self.u)
        ud_zero = dict((i, 0) for i in self._udot)
        qd_zero = dict((i, 0) for i in self._qdot)
        qd_u_zero = dict((i, 0) for i in Matrix([self._qdot, self.u]))
        # Break the kinematic differential eqs apart into f_0 and f_1
        f_0 = msubs(self._f_k, u_zero) + self._k_kqdot*Matrix(self._qdot)
        f_1 = msubs(self._f_k, qd_zero) + self._k_ku*Matrix(self.u)
        # Break the dynamic differential eqs into f_2 and f_3
        f_2 = msubs(self._frstar, qd_u_zero)
        f_3 = msubs(self._frstar, ud_zero) + self._fr
        f_4 = zeros(len(f_2), 1)

        # Get the required vector components
        q = self.q
        u = self.u
        if self._qdep:
            q_i = q[:-len(self._qdep)]
        else:
            q_i = q
        q_d = self._qdep
        if self._udep:
            u_i = u[:-len(self._udep)]
        else:
            u_i = u
        u_d = self._udep

        # Form dictionary to set auxiliary speeds & their derivatives to 0.
        uaux = self._uaux
        uauxdot = uaux.diff(dynamicsymbols._t)
        uaux_zero = dict((i, 0) for i in Matrix([uaux, uauxdot]))

        # Checking for dynamic symbols outside the dynamic differential
        # equations; throws error if there is.
        sym_list = set(Matrix([q, self._qdot, u, self._udot, uaux, uauxdot]))
        if any(find_dynamicsymbols(i, sym_list) for i in [self._k_kqdot,
                self._k_ku, self._f_k, self._k_dnh, self._f_dnh, self._k_d]):
            raise ValueError('Cannot have dynamicsymbols outside dynamic \
                             forcing vector.')

        # Find all other dynamic symbols, forming the forcing vector r.
        # Sort r to make it canonical.
        r = list(find_dynamicsymbols(msubs(self._f_d, uaux_zero), sym_list))
        r.sort(key=default_sort_key)
        self._lin_r = r
        self._lin_q   = q
        self._lin_q_i = q_i
        self._lin_q_d = q_d
        self._lin_u   = u
        self._lin_u_i = u_i
        self._lin_u_d = u_d

        # Check for any derivatives of variables in r that are also found in r.
        for i in r:
            if diff(i, dynamicsymbols._t) in r:
                print('i',i)
                print('r',r)
                raise ValueError('Cannot have derivatives of specified \
                                 quantities when linearizing forcing terms.')
        return Linearizer(f_0, f_1, f_2, f_3, f_4, f_c, f_v, f_a, q, u, q_i,
                q_d, u_i, u_d, r)

    def linearize(self, **kwargs):
        """ Linearize the equations of motion about a symbolic operating point.

        Explanation
        ===========

        If kwarg A_and_B is False (default), returns M, A, B, r for the
        linearized form, M*[q', u']^T = A*[q_ind, u_ind]^T + B*r.

        If kwarg A_and_B is True, returns A, B, r for the linearized form
        dx = A*x + B*r, where x = [q_ind, u_ind]^T. Note that this is
        computationally intensive if there are many symbolic parameters. For
        this reason, it may be more desirable to use the default A_and_B=False,
        returning M, A, and B. Values may then be substituted in to these
        matrices, and the state space form found as
        A = P.T*M.inv()*A, B = P.T*M.inv()*B, where P = Linearizer.perm_mat.

        In both cases, r is found as all dynamicsymbols in the equations of
        motion that are not part of q, u, q', or u'. They are sorted in
        canonical form.

        The operating points may be also entered using the ``op_point`` kwarg.
        This takes a dictionary of {symbol: value}, or a an iterable of such
        dictionaries. The values may be numeric or symbolic. The more values
        you can specify beforehand, the faster this computation will run.

        For more documentation, please see the ``Linearizer`` class."""
        linearizer = self.to_linearizer()
        result = linearizer.linearize(**kwargs)
        return result + (linearizer.r,)

    def kanes_equations(self, bodies, loads=None, Mform='TaylorExpanded',addGravity=True, g_vect=None):
        """ Method to form Kane's equations, Fr + Fr* = 0.

        Explanation
        ===========

        Returns (Fr, Fr*). In the case where auxiliary generalized speeds are
        present (say, s auxiliary speeds, o generalized speeds, and m motion
        constraints) the length of the returned vectors will be o - m + s in
        length. The first o - m equations will be the constrained Kane's
        equations, then the s auxiliary Kane's equations. These auxiliary
        equations can be accessed with the auxiliary_eqs().

        Parameters
        ==========

        bodies : iterable
            An iterable of all RigidBody's and Particle's in the system.
            A system must have at least one body.
        loads : iterable
            Takes in an iterable of (Particle, Vector) or (ReferenceFrame, Vector)
            tuples which represent the force at a point or torque on a frame.
            Must be either a non-empty iterable of tuples or None which corresponds
            to a system with no constraints.
        """
        if (bodies is None and loads is not None) or isinstance(bodies[0], tuple):
            # This switches the order if they use the old way.
            bodies, loads = loads, bodies
            SymPyDeprecationWarning(value='The API for kanes_equations() has changed such '
                    'that the loads (forces and torques) are now the second argument '
                    'and is optional with None being the default.',
                    feature='The kanes_equation() argument order',
                    useinstead='switched argument order to update your code, For example: '
                    'kanes_equations(loads, bodies) > kanes_equations(bodies, loads).',
                    issue=10945, deprecated_since_version="1.1").warn()

        if not self._k_kqdot:
            raise AttributeError('Create an instance of KanesMethod with '
                    'kinematic differential equations to use this method.')
        fr = self._form_fr(loads)
        frstar = self._form_frstar(bodies, Mform=Mform, addGravity=addGravity, g_vect=g_vect)
        if self._uaux:
            if not self._udep:
                km = KanesMethod(self._inertial, self.q, self._uaux,
                             u_auxiliary=self._uaux)
            else:
                km = KanesMethod(self._inertial, self.q, self._uaux,
                        u_auxiliary=self._uaux, u_dependent=self._udep,
                        velocity_constraints=(self._k_nh * self.u +
                        self._f_nh))
            km._qdot_u_map = self._qdot_u_map
            self._km = km
            fraux = km._form_fr(loads)
            frstaraux = km._form_frstar(bodies)
            self._aux_eq = fraux + frstaraux
            self._fr = fr.col_join(fraux)
            self._frstar = frstar.col_join(frstaraux)
        return (self._fr, self._frstar)

    def _form_eoms(self):
        fr, frstar = self.kanes_equations(self.bodylist, self.forcelist)
        return fr + frstar

    def rhs(self, inv_method=None):
        """Returns the system's equations of motion in first order form. The
        output is the right hand side of::

           x' = |q'| =: f(q, u, r, p, t)
                |u'|

        The right hand side is what is needed by most numerical ODE
        integrators.

        Parameters
        ==========

        inv_method : str
            The specific sympy inverse matrix calculation method to use. For a
            list of valid methods, see
            :meth:`~sympy.matrices.matrices.MatrixBase.inv`

        """
        rhs = zeros(len(self.q) + len(self.u), 1)
        kdes = self.kindiffdict()
        for i, q_i in enumerate(self.q):
            rhs[i] = kdes[q_i.diff()]

        if inv_method is None:
            rhs[len(self.q):, 0] = self.mass_matrix.LUsolve(self.forcing)
        else:
            rhs[len(self.q):, 0] = (self.mass_matrix.inv(inv_method,
                                                         try_block_diag=True) *
                                    self.forcing)

        return rhs

    def kindiffdict(self):
        """Returns a dictionary mapping q' to u."""
        if not self._qdot_u_map:
            raise AttributeError('Create an instance of KanesMethod with '
                    'kinematic differential equations to use this method.')
        return self._qdot_u_map

    @property
    def auxiliary_eqs(self):
        """A matrix containing the auxiliary equations."""
        if not self._fr or not self._frstar:
            raise ValueError('Need to compute Fr, Fr* first.')
        if not self._uaux:
            raise ValueError('No auxiliary speeds have been declared.')
        return self._aux_eq

    @property
    def mass_matrix(self):
        """The mass matrix of the system."""
        if not self._fr or not self._frstar:
            raise ValueError('Need to compute Fr, Fr* first.')
        return Matrix([self._k_d, self._k_dnh])

    @property
    def mass_matrix_full(self):
        """The mass matrix of the system, augmented by the kinematic
        differential equations."""
        if not self._fr or not self._frstar:
            raise ValueError('Need to compute Fr, Fr* first.')
        o = len(self.u)
        n = len(self.q)
        return ((self._k_kqdot).row_join(zeros(n, o))).col_join((zeros(o,
                n)).row_join(self.mass_matrix))

    @property
    def forcing(self):
        """The forcing vector of the system."""
        if not self._fr or not self._frstar:
            raise ValueError('Need to compute Fr, Fr* first.')
        return -Matrix([self._f_d, self._f_dnh])

    @property
    def forcing_full(self):
        """The forcing vector of the system, augmented by the kinematic
        differential equations."""
        if not self._fr or not self._frstar:
            raise ValueError('Need to compute Fr, Fr* first.')
        f1 = self._k_ku * Matrix(self.u) + self._f_k
        return -Matrix([f1, self._f_d, self._f_dnh])

    @property
    def q(self):
        return self._q

    @property
    def u(self):
        return self._u

    @property
    def bodylist(self):
        return self._bodylist

    @property
    def forcelist(self):
        return self._forcelist

    @property
    def bodies(self):
        return self._bodylist

    @property
    def loads(self):
        return self._forcelist


# --------------------------------------------------------------------------------}
# --- Kane's method 
# --------------------------------------------------------------------------------{
def _initialize_kindiffeq_matrices(coordinates, speeds, kdeqs, uaux=Matrix()):
    """Initialize the kinematic differential equation matrices.
    See sympy.mechanics.kane

    kdeqs: kinematic differential equations
    """
    from sympy import solve_linear_system_LU

    if kdeqs:
        if len(coordinates) != len(kdeqs):
            raise ValueError('There must be an equal number of kinematic differential equations and coordinates.')
        coordinates = Matrix(coordinates)
        kdeqs = Matrix(kdeqs)
        u     = Matrix(speeds)
        qdot  = coordinates.diff(dynamicsymbols._t)
        # Dictionaries setting things to zero
        u_zero = dict((i, 0) for i in u)
        uaux_zero = dict((i, 0) for i in uaux)
        qdot_zero = dict((i, 0) for i in qdot)

        f_k = msubs(kdeqs, u_zero, qdot_zero)
        k_ku = (msubs(kdeqs, qdot_zero) - f_k).jacobian(u)
        k_kqdot = (msubs(kdeqs, u_zero) - f_k).jacobian(qdot)

        f_k = k_kqdot.LUsolve(f_k)
        k_ku = k_kqdot.LUsolve(k_ku)
        k_kqdot = eye(len(qdot))

        _qdot_u_map = solve_linear_system_LU( Matrix([k_kqdot.T, -(k_ku * u + f_k).T]).T, qdot)

        _f_k = msubs(f_k, uaux_zero)
        _k_ku = msubs(k_ku, uaux_zero)
        _k_kqdot = k_kqdot
    else:
        _qdot_u_map = None
        _f_k = Matrix()
        _k_ku = Matrix()
        _k_kqdot = Matrix()
    return _qdot_u_map, _f_k, _k_ku, _k_kqdot



def kane_frstar_alt(bodies, coordinates, speeds, kdeqs, inertial_frame, uaux=Matrix(), udep=None, Ars=None):
    """Form the generalized inertia force."""

    t = dynamicsymbols._t
    N = inertial_frame

    # Derived inputs
    q = Matrix(coordinates) # q
    u = Matrix(speeds) # u
    udot = u.diff(t)
    qdot_u_map,_,_,_k_kqdot = _initialize_kindiffeq_matrices(q, u, kdeqs, uaux=Matrix())
                                                  
    # Dicts setting things to zero
    udot_zero = dict((i, 0) for i in udot)
    uaux_zero = dict((i, 0) for i in uaux)
    uauxdot   = [diff(i, t) for i in uaux]
    uauxdot_zero = dict((i, 0) for i in uauxdot)
    # Dictionary of q' and q'' to u and u'
    q_ddot_u_map = dict((k.diff(t), v.diff(t)) for (k, v) in qdot_u_map.items())
    q_ddot_u_map.update(qdot_u_map)

    # Fill up the list of partials: format is a list with num elements
    # equal to number of entries in body list. Each of these elements is a
    # list - either of length 1 for the translational components of
    # particles or of length 2 for the translational and rotational
    # components of rigid bodies. The inner most list is the list of
    # partial velocities.
    def get_partial_velocity(body):
        if isinstance(body,YAMSRigidBody) or isinstance(body, RigidBody):
            vlist = [body.masscenter.vel(N), body.frame.ang_vel_in(N)]
        elif isinstance(body, Particle):
            vlist = [body.point.vel(N),]
        elif isinstance(body,YAMSFlexibleBody):
            print('>>>> FlexibleBody TODO, Jv Jo to partials')
            vlist=[body.masscenter.vel(N), body.frame.ang_vel_in(N)]
        else:
            raise TypeError('The body list may only contain either ' 'RigidBody or Particle as list elements.')
        v = [msubs(vel, qdot_u_map) for vel in vlist]
        return partial_velocity(v, u, N)

    partials = [get_partial_velocity(body) for body in bodies]

    # Compute fr_star in two components:
    # fr_star = -(MM*u' + nonMM)
    o = len(u)
    MM = zeros(o, o)
    nonMM = zeros(o, 1)
    zero_uaux      = lambda expr: msubs(expr, uaux_zero)
    zero_udot_uaux = lambda expr: msubs(msubs(expr, udot_zero), uaux_zero)
    for i, body in enumerate(bodies):
        bodyMM = zeros(o, o)
        bodynonMM = zeros(o, 1)
        if isinstance(body,YAMSRigidBody) or isinstance(body, RigidBody):
            # Rigid Body (see sympy.mechanics.kane)
            M     = zero_uaux(       body.mass                )
            I     = zero_uaux(       body.central_inertia     )
            vel   = zero_uaux(       body.masscenter.vel(N)   )
            omega = zero_uaux(       body.frame.ang_vel_in(N) )
            acc   = zero_udot_uaux(  body.masscenter.acc(N)   )
            # --- Mas Matrix
            for j in range(o):
                tmp_vel = zero_uaux(partials[i][0][j])
                tmp_ang = zero_uaux(I & partials[i][1][j])
                for k in range(o):
                    # translational
                    bodyMM[j, k] += M * (tmp_vel & partials[i][0][k])
                    # rotational
                    bodyMM[j, k] += (tmp_ang & partials[i][1][k])
            # --- Full inertial loads Matrix
            inertial_force  = (M.diff(t) * vel + M * acc)
            inertial_torque = zero_uaux((I.dt(body.frame) & omega) + msubs(I & body.frame.ang_acc_in(N), udot_zero) + (omega ^ (I & omega)))  # "&" = dot, "^"=cross
            for j in range(o):
                bodynonMM[j] += inertial_force & partials[i][0][j]
                bodynonMM[j] += inertial_torque & partials[i][1][j]

        elif isinstance(body,YAMSFlexibleBody):
            print('>>>> FlexibleBody TODO')
            M     = zero_uaux(body.mass)
            #I     = zero_uaux(body.central_inertia)
            vel   = zero_uaux(body.origin.vel(N))
            omega = zero_uaux(body.frame.ang_vel_in(N))
            acc   = zero_udot_uaux(body.origin.acc(N))
            inertial_force=0 # Fstar  !<<<< TODO
            inertial_torque=0 # Tstar  !<<<< TODO

        else:
            # Particles
            M = zero_uaux(body.mass)
            vel = zero_uaux(body.point.vel(N))
            acc = zero_udot_uaux(body.point.acc(N))
            inertial_force = (M.diff(t) * vel + M * acc)
            inertial_torque=0 # Tstar
            for j in range(o):
                temp = zero_uaux(partials[i][0][j])
                for k in range(o):
                    bodyMM[j, k] += M * (temp & partials[i][0][k])
                bodynonMM[j] += inertial_force & partials[i][0][j]

        # Perform important substitution and store body contributions
        body.MM_alt     = zero_uaux(msubs(bodyMM, q_ddot_u_map))
        body.nonMM_alt_bef = bodynonMM
        #body.nonMM_alt  = msubs(msubs(bodynonMM, q_ddot_u_map), udot_zero, uauxdot_zero, uaux_zero)
        # Cumulative MM and nonMM over all bodies
        MM    += bodyMM
        nonMM += bodynonMM
        # --- Storing for debug
        body.acc_alt             = acc
        body.vel_alt             = vel
        body.omega_alt           = omega
        body.inertial_force_alt  = inertial_force
        body.inertial_torque_alt = inertial_torque
        body.Jv_vect_alt=partials[i][0]
        body.Jo_vect_alt=partials[i][1]
    # End loop on bodies
    nonMM = msubs(msubs(nonMM, q_ddot_u_map), udot_zero, uauxdot_zero, uaux_zero)

    # Compose fr_star out of MM and nonMM
    fr_star = -(MM * msubs(Matrix(udot), uauxdot_zero) + nonMM)

    # If there are dependent speeds, we need to find fr_star_tilde
    if udep:
        p = o - len(udep)
        fr_star_ind = fr_star[:p, 0]
        fr_star_dep = fr_star[p:o, 0]
        fr_star = fr_star_ind + (Ars.T * fr_star_dep)
        # Apply the same to MM
        MMi = MM[:p, :]
        MMd = MM[p:o, :]
        MM = MMi + (Ars.T * MMd)

    MM_full= mass_matrix_full(MM,_k_kqdot)
    #self._bodylist = bodies
    #self._frstar = fr_star
    #self._k_d = MM
    #self._f_d = -msubs(self._fr + self._frstar, udot_zero)
    return fr_star, MM, MM_full

def kane_frstar(bodies, coordinates, speeds, kdeqs, origin, inertial_frame, 
        Omega_Subs=[(None,None)], Mform='TaylorExpanded', addGravity=False, g_vect=None):
    """ 
    coordinates "q"
    speeds      "u"
    kdeqs:   relates qdot and udot


    Computes Jv and Jo at a reference point:
      - For rigid bodies centere of mass
      - For flexible bodies origin
    
    """ 
    from sympy.physics.vector import partial_velocity

    nq = len(coordinates)
    MM    = zeros(nq, nq)
    nonMM = zeros(nq, 1)

    O_E = origin
    N = inertial_frame
    t = dynamicsymbols._t

    # Derived inputs
    q = Matrix(coordinates) # q
    u = Matrix(speeds) # u
    udot = u.diff(t)
    qspeeds = q.diff(t)
    qacc = qspeeds.diff(t)
    qdot_u_map,_,_,_k_kqdot = _initialize_kindiffeq_matrices(q, u, kdeqs, uaux=Matrix())

    # Dicts setting things to zero
    udot_zero = dict((i, 0) for i in udot)
    qdot_zero = dict((diff(qd,t), 0) for qd in qspeeds)
    # Dictionary of q' and q'' to u and u'
    q_ddot_u_map = dict((k.diff(t), v.diff(t)) for (k, v) in qdot_u_map.items())
    q_ddot_u_map.update(qdot_u_map)

    for i,body in enumerate(bodies):
        bodyMM    = zeros(nq, nq)
        bodynonMM = zeros(nq, 1)
        M     = body.mass
        #print(type(body),isinstance(body, YAMSFlexibleBody), isinstance(body, YAMSRigidBody), isinstance(body, YAMSBody), isinstance(body, SympyBody))
        if isinstance(body, YAMSFlexibleBody):
            P = body.origin
            I = None
        else:
            P = body.masscenter
            I = body.central_inertia

        # --- Step 2: Positions and orientation
        r = P.pos_from(O_E)
        R = N.dcm(body.frame) # from body to inertial

        # --- Step 3/4: Velocities and accelerations
        vel   = P.vel(N)
        omega = body.frame.ang_vel_in(N)
        acc   = P.acc(N)
        alpha = omega.diff(t, N)
        # NOTE: Keep me Alternative: vel from r and omega from identification:
        #vel_drdt = r.diff(t, N).simplify()
        #OmSkew = (R.diff(t) *  R.transpose()).simplify()
        #omega_ident = OmSkew[2,1] * N.x + OmSkew[0,2]*N.y + OmSkew[1,0] * N.z
        #acc_drdt2 = r.diff(t, N).diff(t,N).simplify()

        # --- Step 5 Partial velocities
        # Method 1: use "partial_velocity" function, which returns a vector
        vel_sub   = msubs(vel, qdot_u_map)
        Jv_vect   = partial_velocity([vel_sub], u, N)[0]
        omega_sub = msubs(omega, qdot_u_map)
        Jo_vect = partial_velocity([omega_sub], u, N)[0]
        # NOTE: Keep me: Method 2: express everything in ref frame, and use "jacobians"
        #v  = vel.subs(Omega_Subs).to_matrix(N)
        #om = omega.subs(Omega_Subs).to_matrix(N)
        #Jv = v.jacobian(qspeeds)
        #Jo = om.jacobian(qspeeds)

        # --- Step 6 Inertia forces
        if isinstance(body,YAMSRigidBody) or isinstance(body, SympyRigidBody):
            # --- Mass Matrix 
            for j in range(nq):
                tmp_vel = Jv_vect[j]      # Jv[:,j]
                tmp_ang = I & Jo_vect[j]  # I Jo[:,j]
                for k in range(nq):
                    # translational
                    bodyMM[j, k] += M * (tmp_vel & Jv_vect[k]) # M * Jv[:,j] dot Jv[:,k]
                    # rotational
                    bodyMM[j, k] +=     (tmp_ang & Jo_vect[k]) # I dot Jo[:,j] dot Jo[:,k]

            # --- Full inertial loads
            inertial_force  = (M.diff(t) * vel + M * acc) # "Fstar"
            inertial_torque = (I.dt(body.frame) & omega) + msubs(I & body.frame.ang_acc_in(N), udot_zero) + (omega ^ (I & omega))  # "&" = dot, "^"=cross

            # NOTE KEEP ME: Alternative formulation using "matrices" 
            #inertial_force = inertial_force.subs(Omega_Subs)
            #RIRt  = R*I.to_matrix(body.frame)*R.transpose()
            #inertial_torque_2 = - RIRt * alpha.to_matrix(N) \
            #        - coord2vec(om, N).cross( coord2vec(RIRt *om, N)).to_matrix(N)

            # Computing generatlized force Jv.f + Jo*M
            for j in range(nq):
                bodynonMM[j] += inertial_force  & Jv_vect[j]
                bodynonMM[j] += inertial_torque & Jo_vect[j]
            # Alternative:
            #frstar_t+ = Jv.transpose() * inertial_force.to_matrix(N)
            #frstar_o+ = Jo.transpose() * inertial_torque

        elif isinstance(body,YAMSFlexibleBody):
            if addGravity:
                print('>>> add Gravity')
                acc=acc-g_vect
            MMloc = body.bodyMassMatrix(form=Mform)
            body.h_omega = body.bodyQuadraticForce(omega.to_matrix(body.frame), body.q, body.qdot)
            body.h_elast = body.bodyElasticForce(body.q, body.qdot)
            #body.h_g     = body.bodyGravitationalForce(g_vect, body.q, form='TaylorExpanded')
            #  Inertial forces and torques at the body origin (wher Jv and Jo were computed)
            inertial_force=0 # Fstar
            inertial_torque=0 # Tstar
            inertial_force_coord =MMloc[0:3,0:3] * acc.to_matrix(body.frame) 
            inertial_force_coord+=MMloc[0:3,3:6] * alpha.to_matrix(body.frame)
            inertial_force_coord+=MMloc[0:3,6:]  * Matrix(body.qddot)
            inertial_force_coord+=body.h_omega[0:3,0]
            inertial_force_coord+=body.h_elast[0:3,0]
            inertial_torque_coord =MMloc[3:6,0:3] * acc.to_matrix(body.frame) 
            inertial_torque_coord+=MMloc[3:6,3:6] * alpha.to_matrix(body.frame)
            inertial_torque_coord+=MMloc[3:6,6:]  * Matrix(body.qddot)
            inertial_torque_coord+=body.h_omega[3:6,0]
            inertial_torque_coord+=body.h_elast[3:6,0]
            inertial_elast_coord =MMloc[6:,0:3] * acc.to_matrix(body.frame) 
            inertial_elast_coord+=MMloc[6:,3:6] * alpha.to_matrix(body.frame)
            inertial_elast_coord+=MMloc[6:,6:]  * Matrix(body.qddot)
            inertial_elast_coord+=body.h_omega[6:,0]
            inertial_elast_coord+=body.h_elast[6:,0]
            body.inertial_elast=inertial_elast_coord

            inertial_force  = coord2vec(inertial_force_coord,body.frame) 
            inertial_torque = coord2vec(inertial_torque_coord,body.frame) 

            # Computing generatlized force Jv.f + Jo*M
            for j in range(nq):
                bodynonMM[j] += inertial_force  & Jv_vect[j]
                bodynonMM[j] += inertial_torque & Jo_vect[j]
                for k in range(len(body.q)):
                    if q[j] == body.q[k]:
                        bodynonMM[j] +=  inertial_elast_coord[k]
            bnMMSubs = msubs(bodynonMM, q_ddot_u_map)
            bodyMM = bnMMSubs.jacobian(udot) # TODO why not using MMloc????? <<<<<<

        else:
            raise Exception('Unsupported body type: {}'.format(type(body)))

        # Perform important substitution and store body contributions
        body.MM      = msubs(bodyMM, q_ddot_u_map)
        body.nonMM   = bodynonMM

        # Cumulative MM and nonMM over all bodies
        MM   +=bodyMM
        nonMM+=bodynonMM
        # --- Storing for debug
        body._acc             = acc
        body._vel             = vel
        body._omega           = omega
        body._inertial_force  = inertial_force
        body._inertial_torque = inertial_torque
        body._Jv_vect=Jv_vect
        body._Jo_vect=Jo_vect
    # End loop on bodies
    # NOTE: substitution needs to be done at the end
    nonMM  = msubs(msubs(nonMM, q_ddot_u_map), udot_zero) #, uauxdot_zero, uaux_zero)

    # Compose fr_star out of MM and nonMM
    fr_star = -(MM *Matrix(udot) + nonMM)

    MM_full= mass_matrix_full(MM,_k_kqdot)

    return fr_star, MM, MM_full

# --------------------------------------------------------------------------------}
# --- Kane fr 
# --------------------------------------------------------------------------------{
def kane_fr_alt(loads, coordinates, speeds, kdeqs, inertial_frame, uaux=Matrix(), udep=None):
    """
      - For each force: compute the velocity at the point of application, v_P, and then do
        fr =  [d v/ dqdot]^t F
      - For each moment: compute the angular velocity of the frame (in E)
        fr =  [d om/ dqdot]^t M
    """
    from sympy.physics.vector import partial_velocity
 

    def _f_list_parser(fl, ref_frame):
        """Parses the provided forcelist composed of items of the form (obj, force).
        Returns a tuple containing:
            vel_list: The velocity (ang_vel for Frames, vel for Points) in the provided reference frame.
            f_list: The forces.
        Used internally in the KanesMethod and LagrangesMethod classes.
        """
        def flist_iter():
            for pair in fl:
                obj, force = pair
                if isinstance(obj, ReferenceFrame):
                    yield obj.ang_vel_in(ref_frame), force
                elif isinstance(obj, Point):
                    yield obj.vel(ref_frame), force
                else:
                    raise TypeError('First entry in each forcelist pair must be a point or frame.')
        if not fl:
            vel_list, f_list = (), ()
        else:
            unzip = lambda l: list(zip(*l)) if l[0] else [(), ()]
            vel_list, f_list = unzip(list(flist_iter()))
        return vel_list, f_list

    """
    Form the generalized active force.
        See _form_fr in sympy.mechanics
    """
    N = inertial_frame
    # Derived inputs
    speeds = Matrix(speeds) # u
    qdot_u_map,_,_,_k_kqdot = _initialize_kindiffeq_matrices(coordinates, speeds, kdeqs, uaux=Matrix())
    
    # pull out relevant velocities for constructing partial velocities
    vel_list, f_list = _f_list_parser(loads, N)
    vel_list = [msubs(i, qdot_u_map) for i in vel_list]
    f_list   = [msubs(i, qdot_u_map) for i in f_list]

    # Fill Fr with dot product of partial velocities and forces
    o = len(speeds)
    b = len(f_list)
    FR = zeros(o, 1)
    partials = partial_velocity(vel_list, speeds, N)
    for i in range(o):
        FR[i] = sum(partials[j][i] & f_list[j] for j in range(b))

    # In case there are dependent speeds
    if udep:
        p = o - len(udep)
        FRtilde = FR[:p, 0]
        FRold = FR[p:o, 0]
        FRtilde += Ars.T * FRold
        FR = FRtilde
    #self._forcelist = loads
    #self._fr = FR
    return FR


def kane_fr(body_loads, speeds, inertial_frame):
    """
    Compute Kane's "fr" terms, using a list of bodies and external loads

    For each body:  fr = Jv * F@refP  + Jo * M
    where Jv and Jo are the jacoban of the linear velocity of the point (for a force) and angular velocity of the frame (for a moment). The point "refP" should be the one used in the calculation of Jv (typically the center of mass for rigid body)

    Right now, Jv, and Jo are computed when calling kane_frstar...

    Alternatively (see kane_fr_alt and _form_fr):
      - For each force: compute the velocity at the point of application, v_P, and then do
        fr =  [d v/ dqdot]^t F
      - For each moment: compute the angular velocity of the frame (in E)
        fr =  [d om/ dqdot]^t M

    INPUTS:
        body_loads: a list of tuples of the form  (body, (point_or_frame, force_or_moment ) )
            The tuples (point_or_frame, force_or_moment) are the ones needed when calling sympy's kane
            For instance:
               body_loads = [
                    (nac, (N        , Thrust*N.x)),
                    (nac, (nac.frame, Qaero*N.x ))
                    ]
        
    """
    import sympy.physics.vector as vect

    fr_t = zeros(len(speeds), 1)
    fr_o = zeros(len(speeds), 1)

    N = inertial_frame

    for bl in body_loads:
        body, (point_or_frame, force_or_moment) = bl
        if not hasattr(body,'_Jv_vect') or not hasattr(body, '_Jo_vect'):
            raise Exception('Jacobians matrices need to be computed for body {}. (Call frstart first)'.format(body.name))
        if isinstance(point_or_frame, ReferenceFrame):
            # Moment and frame
            Moment = force_or_moment
            for j in range(len(speeds)):
                fr_o[j] += body._Jo_vect[j] & Moment # Jo^t * M
            pass
        else:
            # Force and point
            Force = force_or_moment
            point = point_or_frame
            r = point.pos_from(body.masscenter)
            for j in range(len(speeds)):
                fr_t[j] += body._Jv_vect[j] & Force  # Jv^t * F
            # Need to add moment if r/=0
            for j in range(len(speeds)):
                fr_o[j] += body._Jo_vect[j] & (vect.cross( r, Force)) # Jo^t * M
    return fr_t+fr_o



# --------------------------------------------------------------------------------}
# --- Mass matrix and forcing 
# --------------------------------------------------------------------------------{
def mass_matrix_full(MM,k_kqdot):
    """The mass matrix of the system, augmented by the kinematic
    differential equations.
    x'= [q'; u']
    """
    n = MM.shape[0]
    o = k_kqdot.shape[0]
    return ((k_kqdot).row_join(zeros(n, o))).col_join((zeros(o, n)).row_join(MM))

def forcing(fr, frstar, speeds):
    """The forcing vector of the system."""
    u = Matrix(speeds) # u
    udot = u.diff(dynamicsymbols._t)
    udot_zero = dict((i, 0) for i in udot)
    _f_d = -msubs(fr + frstar, udot_zero)
    return -Matrix([_f_d])

def forcing_full(fr, frstar, coordinates, speeds, kdeqs):
    """The forcing vector of the system, augmented by the kinematic
    differential equations.
    x'= [q'; u']
    """
    _f_d = -forcing(fr, frstar, speeds)
    _, _f_k, _k_ku, _= _initialize_kindiffeq_matrices(coordinates, speeds, kdeqs)
    f1 = _k_ku * Matrix(speeds) +_f_k
    return -Matrix([f1, _f_d])
# 
