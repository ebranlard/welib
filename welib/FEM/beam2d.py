
def beam2d_KeMe(EI,L,Mass,T=0, theta=None, MMFormulation='consistent' ):
    """ 
    Stiffness and mass matrices for Hermitian beam element with 2DOF per node.
      - Bending in one transverse direction
    Euler-Bernoulli beam model. 
        
    The beam coordinate system is such that the cross section is assumed to be in the y-z plane
    (along x)
        
    Nodal DOF:    (u theta)
    Element DOFs: (u1 t1 u2 t2) 
        
    INPUTS
        EI  : Young Modulus times Planar second moment of area,local y-axis. Iy=\iint z^2 dy dz [m4]
        L   :    Element length
        Mass :    Element mass = rho * A * L [kg]

    OPTIONAL INPUTS
        theta: orientation angle of element in global coordinate system
               if provided, element matrix is provided in global coord
        T   : Axial load to be used for the compuation of the geometrical stiffness

    OUTPUTS
        Ke: Element stiffness matrix (4x4)
        Me: Element mass matrix      (4x4)
        Kg: Element geometrical stiffness matrix (4,4)
        
    AUTHOR: E. Branlard
    """
    # NOTE: matrices determined using sympy, see scripts in current folder 

    # --- Stiffness matrices
    # ke=EI/(L^3)*[12      6*L   -12      6*L ; ...
    #              6*L   4*L^2   -6*L   2*L^2 ; ...
    #             -12     -6*L    12     -6*L ; ...
    #              6*L   2*L^2   -6*L   4*L^2];
    Ke = np.array( [
            [ 12*EI/L**3  , 6*EI/L**2   , -12*EI/L**3 , 6*EI/L**2]  , 
            [ 6*EI/L**2   , 4*EI/L      , -6*EI/L**2  , 2*EI/L]     , 
            [ -12*EI/L**3 , -6*EI/L**2  , 12*EI/L**3  , -6*EI/L**2] , 
            [ 6*EI/L**2   , 2*EI/L      , -6*EI/L**2  , 4*EI/L]
            ])
    Kg = np.array([
            [ 6*T/(5*L)  , T/10     , -6*T/(5*L) , T/10]     , 
            [ T/10       , 2*L*T/15 , -T/10      , -L*T/30]  , 
            [ -6*T/(5*L) , -T/10    , 6*T/(5*L)  , -T/10]    , 
            [ T/10       , -L*T/30  , -T/10      , 2*L*T/15]
        ])

    # Mass matrix
    if MMFormulation=='consistent':
        # Consistent Formulation
        # me=Mass/420*[156    22*L   54    -13*L  ; ...
        #              22*L  4*L^2  13*L  -3*L^2  ; ...
        #              54    13*L   156      -22*L; ...
        #             -13*L -3*L^2 -22*L   4*L^2] ; 
        me = Mass / 420 * np.array([
            [156      , 22 * L       , 54       , - 13 * L]     , 
            [22 * L   , 4 * L ** 2   , 13 * L   , - 3 * L ** 2] , 
            [54       , 13 * L       , 156      , - 22 * L]     , 
            [- 13 * L , - 3 * L ** 2 , - 22 * L , 4 * L ** 2]
            ])
    elif MMFormulation=='lumped':
        # Lumped formulation
        me = np.diag([Mass/2, 0, Mass/2, 0]))
        #  TODO?
        #     alpha = 17.5;
        #   me = rho*A*L / 2 * ...  # lumped
        #   [ 1   0              0  0
        #     0   alpha*L^2/210  0  0
        #     0   0              1  0
        #     0   0              0  alpha*L^2/210 ];
    elif MMFormulation=='diagonal':
        # Diagonal formulation
        me = Mass * np.diag([1/2, L**2/78, 1/2, L**2/78])
    else:
        raise Exception('Unknown  mass matrix formulation {}'.format(MMFormulation))


    # --- Conversion to global system if requested
    if theta is not None:
        # TODO
        #R = np.array([
        #    [np.cos(theta)   , np.sin(theta) , 0 , 0               , 0             , 0]   , 
        #    [- np.sin(theta) , np.cos(theta) , 0 , 0               , 0             , 0]   , 
        #    [0               , 0             , 1 , 0               , 0             , 0]   , 
        #    [0               , 0             , 0 , np.cos(theta)   , np.sin(theta) , 0]   , 
        #    [0               , 0             , 0 , - np.sin(theta) , np.cos(theta) , 0]   , 
        #    [0               , 0             , 0 , 0               , 0             , 1]])
        #Me = np.transpose(R) * Me * R
        #Ke = np.transpose(R) * Ke * R
        #Kg = np.transpose(R) * Kg * R

    return Ke, Me, Kg
