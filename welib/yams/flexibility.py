import numpy as np
import scipy.integrate as sciint
'''
Flexible beam tools:
    - computation of generalized mass and stiffness matrix

Reference:
     [1]: Flexible multibody dynamics using joint coordinates and the Rayleigh-Ritz approximation: the general framework behind and beyond Flex
'''

def fcumtrapzlr(s_span, p):
    r""" Cumulative trapezoidal integration, flipped left-right 
    Useful to return the following:
         P(x) = \int_x^R p(r) dr
    """
    P = - sciint.cumtrapz( p[-1::-1], s_span[-1::-1],)[-1::-1] 
    P = np.concatenate((P,[0]))
    return P


def polyshape(x, coeff, exp=None, x_max=None, doscale=True):
    """ 
    Computes a shape function described as a polynomial expression y = a_i x^e_i
        where the a_i are given by `coeff`
              the e_i are given by `exp`
    The shape function is normalized such as to have a unitary tip deflection

    INPUTS:
      - x : spanwise dimension, from 0 to L, not dimensionless!
            The points 0 and L may not be present in x but in that case prescribing 
            xmax is very important)
      - coeff : polynomial coefficients
      - exp: exponents of the polynomial. Should be length of coeff. 
            If None, exp = [2,3,4,5,6] as used in OpenFAST
      - xmax:value used to non-dimensionlize x:  x_bar = x/x_max. default: max(x)
      - doscale: used the tip value to scale the shapes and ensure a unit value at the tip
               In this case the "tip" is defined as xbar=1 
    Returns:
      - U, dU, ddU the shape, slope and curvature
    """
    if exp is None:
        exp = np.arange(2,7)
    mode   = np.zeros(x.shape)
    dmode  = np.zeros(x.shape)
    ddmode = np.zeros(x.shape)
    # Polynomials assume x to be dimensionless
    # TODO TODO TODO substract x[0]
    if x_max is None: 
        x_max= (x[-1]-0)  # Might not always be desired if "x" are mid-points
    x_bar=(x-0)/x_max
    mode_max =0   #  value of shape function at max x
    for i in range(0,len(coeff)):
        mode     += coeff[i]*x_bar**exp[i]
        mode_max += coeff[i]*1**exp[i]
        if exp[i]-1>=0:
            dmode += coeff[i]*exp[i]* x_bar**(exp[i]-1)
        if exp[i]-2>=0:
            ddmode += coeff[i]*exp[i]*(exp[i]-1) * x_bar**(exp[i]-2)
    # Scaling by the tip deflection, and include x_max for derivatives since derivates were computed w.r.t. x_bar not x
    #scale= mode[-1]
    if doscale:
        scale= mode_max
    else:
        scale=1
    return mode/scale, dmode/(scale*x_max), ddmode/(scale*x_max*x_max)


def integrationWeights(s_span,m):
    r""" Returns integration weights convenient to integrate functions along the span of the beam
    The equations are written such that s_span(1) is not necessary 0
    
    - Span Integration weights  IW and IW_x 
      Assuming a fonction f that varies linearly
         IW   is such that \int   f(x) dx = \Sum IW  [i] F[i]
         IW_x is such that \int x.f(x) dx = \Sum IW_x[i] F[i]

    - Mass integration weights iW_m iW_xm 
         iW_m  is such that \int   f(x).m(x) dx = \Sum iW_m [i] F[i]
         iW_xm is such that \int x.f(x).m(x) dx = \Sum iW_xm[i] F[i]
    """
    IW  =np.zeros(s_span.shape); 
    IW_x=np.zeros(s_span.shape); 
    for i in np.arange(len(s_span)-1):
        L         = s_span[i+1] - s_span[i]
        IW  [i]   = IW[i] + L/2
        IW  [i+1] = L/2
        IW_x[i]   = IW_x[i] + (s_span[i]/2 + L/6)*L
        IW_x[i+1] = (s_span[i]/2 + L/3)*L
    IW_m  = IW*m   
    IW_xm = IW_x*m 
    return IW,IW_x,IW_m,IW_xm


def GKBeamStiffnening(s_span, dU, gravity, m, Mtop, Omega=0, bSelfWeight=True, bMtop=True, bRot=True, main_axis='x'):
    """ 
       Computes geometrical stiffnening for a beam

    OUTPUTS:
     - KKg: Geometrical stiffeness matrix. 
    """
    if gravity is None and (bMtop or bSelfWeight):
        raise Exception('`gravity` is none, but Mtop or SelfWeight is true. Please provide `gravity`')
    nSpan = len(s_span)
    nf    = len(dU)
    KKg = np.zeros((6+nf,6+nf))
    # --- Axial force 
    Pacc    = np.zeros(nSpan) 
    # TopMass contribution to Pacc
    if bMtop:
        Pacc_MT = -Mtop * gravity*np.ones(nSpan)
        Pacc=Pacc+Pacc_MT
    if bSelfWeight:
        Pacc_SW  = fcumtrapzlr(s_span, -m * gravity)
        Pacc=Pacc+Pacc_SW
    if bRot:
        Pacc_Rot = fcumtrapzlr(s_span,  m * Omega**2 * s_span)
        Pacc=Pacc+Pacc_Rot
    # Method 2
    KKCorr = np.zeros((nf,nf))
    for i in range(0,nf):
        for j in range(0,nf):
            #xx=trapz(s_span, Pacc .* PhiV{i}(1,:).* o.PhiV{j}(1,:));
            if main_axis=='x':
                yy=np.trapz(Pacc * dU[i][1,:] * dU[j][1,:] , s_span )
                zz=np.trapz(Pacc * dU[i][2,:] * dU[j][2,:] , s_span )
                KKCorr[i,j]=yy+zz
            elif main_axis=='z':
                xx=np.trapz(Pacc * dU[i][0,:] * dU[j][0,:] , s_span )
                yy=np.trapz(Pacc * dU[i][1,:] * dU[j][1,:] , s_span )
                KKCorr[i,j]=yy+xx
            else:
                raise Exception('Axis not suported')
    #print('KKCorr\n',KKCorr)
    KKg[6:,6:] = KKCorr
    return KKg


def GKBeam(s_span, EI, ddU, bOrth=False):
    """ 
       Computes generalized stiffness matrix for a beam
       Eq.(20) from [1]
       TODO torsion

    OPTIONAL INPUTS:
     - bOrth : if true, enforce orthogonality of modes

    OUTPUTS:
     - KK0: Stiffness matrix without geometrical stiffening
     - The total stiffness matrix should then be KK0+KKg
    """
    nf = len(ddU)
    KK0 = np.zeros((6+nf,6+nf))
    Kgg = np.zeros((nf,nf))
    for i in range(0,nf):
        for j in range(0,nf):
            Kgg[i,j] = np.trapz(EI[0,:]*ddU[i][0,:]*ddU[j][0,:] + EI[1,:]*ddU[i][1,:]*ddU[j][1,:] + EI[2,:]*ddU[i][2,:]*ddU[j][2,:],s_span)
    if bOrth:
        Kgg=Kgg*np.eye(nf)
    #print('Kgg\n',Kgg)
    KK0[6:,6:] = Kgg
    return KK0
    
def GMBeam(s_G, s_span, m, U=None, V=None, jxxG=None, bOrth=False, bAxialCorr=False, IW=None, IW_xm=None, main_axis='x', bUseIW=True, V_tot=None, Peq_tot=None, split_outputs=False, rot_terms=False, method='trapz', U_untwisted=None, M1=False):
    r"""
    Computes generalized mass matrix for a beam.
    Eq.(2) from [1]

    Performing full integration of mass matrix without shape integral functions
    NOTE: Beam assumed to be along x for now (only because of Jxx)
    
    INPUTS
     - s_G    : [m] 3 x nSpan , location of cross sections COG
     - s_span : [m] span integration variable (e.g. s_G(1,:))
     - m      : [kg/m] cross section mass along the beam

    OPTIONAL INPUTS:
     - jxxG   : [kg.m] second moment of inertia of cross section # TODO
     - bOrth : if true, enforce orthogonality of modes
     - U , if omitted, then rigid body (6x6) mass matrix is returned
     - split_outputs: if false (default) return MM, else returns Mxx, Mtt, Mxt, Mtg, Mxg, Mgg
     - rot_terms : if True, outputs the rotational terms as well
     - method: 'trapz', 'Flex', 'OpenFAST' (see below)
     - U_untwsited: untwisted shape functions, used with OpenFAST method only

    OpenFAST method:
      - s_span needs to be [0, np.arange(L/n/2, L, L./n), L]
      - The values for indices 1:-1 represent "elements" (at mid points). 
        First element extends from 0 to L/n, with mid point at L/n/2 
      - Integrals are obtained using summations
      - If provided, the untwisted spape functions "U_untwisted" is used for the mass matrix
        Mgg. The twisted shape functions (likely in U) are used for the coupling terms Ct Cr
    """

    # --- Sanity check on method
    if method=='OpenFAST':
        n = len(s_span)-2
        L = s_span[-1]-s_span[0]
        ds = L/n
        dr = np.diff(s_span[1:-1])
        vds = np.unique(np.around(dr,5))
        dr = np.concatenate( (dr, [ds]))
        melem = m[1:-1]*dr
        if len(vds)>1:
            raise Exception('When using `OpenFAST` method, the user should provide inputs at 0, L, and mid nodes similar to OpenFAST ElastoDYn nodes')
        if vds[0]!=np.around(ds,5):
            raise Exception('When using `OpenFAST` method, the mid nodes should have a spacing equal to L/n')
        # We will only use the inner nodes ("elements")
        U   = U[:,:,1:-1]
        s_G = s_G[:,1:-1]
        #m      = m[1:-1]  # NOTE: temporary, m shouldn't me used with this method
        s_span = s_span[1:-1]  # NOTE: temporary, m shouldn't me used with this method
        m = melem # Important Hack we replace m by melem
        if jxxG is not None:
            jxxG = jxxG[1:-1]*dr # Important Hack
        if U_untwisted is not None:
            U_untwisted = U_untwisted[:,:,1:-1]
        if V is not None:
            V = V[:,:,1:-1]
        bUseIW=False


    elif method=='trapz':
        pass
    elif method=='Flex':
        bUseIW=True
    else:
        raise NotImplementedError()

    if method=='OpenFAST':
        # OpenFAST integration is simple mid-rule summation
        # NOTE: yy is hacked to include "dr" in it already
        def trapzs(yy):
            return np.sum(yy) 
    else:
        # Speed up integration along the span, using integration weight
        def trapzs(yy,**args):
            return np.sum(yy*IW) # NOTE: this is equivalent to trapezoidal integration
    if IW is None or IW_xm is None:
        IW,_,_,IW_xm=integrationWeights(s_span,m)

    if U is not None:
        nf = len(U)
    else:
        nf=0

    # --- Torsion-related variables - Zeros by default
    if jxxG is not None:
        Jxx = trapzs(jxxG) # Imomx  OR Imomz is along z
    else:
        Jxx = 0
    GMJxx=np.zeros(nf);
    I_Jxx=np.zeros(nf);
    if V is not None:
        if main_axis=='x':
            for j in range(nf):
                VJ       = jxxG*V[j][0,:]
                GMJxx[j] = trapzs(V[j][0,:]*VJ)
                I_Jxx[j] = trapzs(VJ)
        elif main_axis=='z':
            # TODO verify me
            for j in range(nf):
                VJ       = jxxG*V[j][2,:]
                GMJxx[j] = trapzs(V[j][2,:]*VJ)
                I_Jxx[j] = trapzs(VJ)

    # --- Mxx
    M = trapzs(m)
    Mxx = np.identity(3)*M
    #print('Mxx\n',Mxx)

    # --- Mxt = -\int [~s] dm    =  -Skew(sigma+Psi g)    Or: +/- Skew(mdCM)
    C_x = trapzs(s_G[0,:]*m)
    C_y = trapzs(s_G[1,:]*m)
    C_z = trapzs(s_G[2,:]*m)
    Mxt = np.array([[0, C_z, -C_y],[-C_z, 0, C_x],[C_y, -C_x, 0]])

    if bAxialCorr:
        # TODO TODO TODO m15 and m16 may need to be additive!
        # --- Variables for axial correction
        # FT=fcumtrapzlr(s_span,m);
        FT = - sciint.cumtrapz( m[-1::-1], s_span[-1::-1],)[-1::-1] 
        FT = np.concatenate((FT,[0]))
        if V_tot is None: 
            raise Exception('Please provide Vtot for axial correction'); end
        if main_axis=='x':
            Mxt[0,1]=+trapzs(V_tot[2,:]*FT) # m15
            Mxt[0,2]=-trapzs(V_tot[1,:]*FT) # m16
        else:
            # TODO TODO TODO VERIFY ME
            Mxt[2,0]=+trapzs(V_tot[1,:]*FT) # m15
            Mxt[2,1]=-trapzs(V_tot[0,:]*FT) # m16
    #print('Mxt\n',Mxt)

    # --- Mxg = \int Phi dm     Or:  Psi , Ct^T
    Mxg      = np.zeros((3,nf))
    for j in range(nf):
        Mxg[0,j] = trapzs(U[j][0,:]*m)
        Mxg[1,j] = trapzs(U[j][1,:]*m)
        Mxg[2,j] = trapzs(U[j][2,:]*m)
    if bAxialCorr:
        # NOTE: not implement with method OpenFAST for now
        # TODO TODO TODO correction may need to be additive
        if (V_tot is not None) and (Peq_tot is not None):
            raise Exception('Provide either V_tot or Peq_tot')
        if V_tot is not None:
            if main_axis=='x':
                for j in range(nf):
                    Mxg[0,j]= trapzs(-V[j][1,:]*V_tot[1,:]*FT - V[j][2,:]*V_tot[2,:]*FT); 
            else:
                for j in range(nf):
                    Mxg[2,j]= trapzs(-V[j][0,:]*V_tot[0,:]*FT - V[j][1,:]*V_tot[1,:]*FT); 
        elif Peq_tot is not None:
            if main_axis=='x':
                for j in range(nf):
                    Mxg[0,j] = trapzs(U[j][1,:]*Peq_tot[1,:] + U[j][2,:]*Peq_tot[2,:] );
            else:
                for j in range(nf):
                    Mxg[2,j] = trapzs(U[j][0,:]*Peq_tot[0,:] + U[j][1,:]*Peq_tot[1,:] );
        else:
            raise Exception('Please provide Vtot of Peq_tot for axial correction');
    #print('Mxg\n',Mxg)
        
    # --- Mtt = - \int [~s][~s] dm  - Or: J, Mrr
    if bUseIW:
        if main_axis=='x':
            s00= np.sum(IW_xm * s_G[0,:]);
            s01= np.sum(IW_xm * s_G[1,:]);
            s02= np.sum(IW_xm * s_G[2,:]);
            s11 = trapzs(s_G[1,:]*s_G[1,:]*m)
            s12 = trapzs(s_G[1,:]*s_G[2,:]*m)
            s22 = trapzs(s_G[2,:]*s_G[2,:]*m)
        elif main_axis=='z':
            s02= np.sum(IW_xm * s_G[0,:]);
            s12= np.sum(IW_xm * s_G[1,:]);
            s22= np.sum(IW_xm * s_G[2,:]);
            s11 = trapzs(s_G[1,:]*s_G[1,:]*m)
            s00 = trapzs(s_G[0,:]*s_G[0,:]*m)
            s01 = trapzs(s_G[0,:]*s_G[1,:]*m)
    else:
        # OpenFAST & trapz (unified via trapzs & m=melem)
        s00 = trapzs(s_G[0,:]*s_G[0,:]*m)
        s01 = trapzs(s_G[0,:]*s_G[1,:]*m)
        s02 = trapzs(s_G[0,:]*s_G[2,:]*m)
        s11 = trapzs(s_G[1,:]*s_G[1,:]*m)
        s12 = trapzs(s_G[1,:]*s_G[2,:]*m)
        s22 = trapzs(s_G[2,:]*s_G[2,:]*m)
    Mtt = np.zeros((3,3))
    Mtt[0,0] = s11 + s22    ;     Mtt[0,1] = -s01;       Mtt[0,2] = -s02
    Mtt[1,0] = -s01;              Mtt[1,1] = s00 + s22;  Mtt[1,2] = -s12
    Mtt[2,0] = -s02;              Mtt[2,1] = -s12;       Mtt[2,2] = s00+s11
    if main_axis=='x':
        Mtt[0,0] += Jxx
    else:
        Mtt[2,2] += Jxx
    #print('Mtt\n',Mtt)

    # --- Mtg  = \int [~s] Phi dm -  Or: Mrg, Cr^T
    # [~s]=[ 0 -z  y]
    #      [ z  0 -x]
    #      [-y  x  0]
    Mtg      = np.zeros((3,nf))
    if bUseIW:
        if main_axis=='x':
            for j in range(nf):
                Mtg[0,j] = trapzs(  (-s_G[2,:]*U[j][1,:] + s_G[1,:]*U[j][2,:])*m)
                Mtg[1,j] = trapzs(  (+s_G[2,:]*U[j][0,:]*m)) - sum(IW_xm*U[j][2,:]);
                Mtg[2,j] = trapzs(  (-s_G[1,:]*U[j][0,:]*m)) + sum(IW_xm*U[j][1,:]);
        elif main_axis=='z':
            for j in range(nf):
                Mtg[0,j] = -sum(IW_xm*U[j][1,:])+trapzs((+ s_G[1,:]*U[j][2,:])*m) 
                Mtg[1,j] =  sum(IW_xm*U[j][0,:])+trapzs((- s_G[0,:]*U[j][2,:])*m)
                Mtg[2,j] = trapzs((-s_G[1,:]*U[j][0,:]   + s_G[0,:]*U[j][1,:])*m)
    else:
        # OpenFAST & trapz (unified via trapzs & m=melem)
        for j in range(nf):
            Mtg[0,j] = trapzs((-s_G[2,:]*U[j][1,:] + s_G[1,:]*U[j][2,:])*m)
            Mtg[1,j] = trapzs(( s_G[2,:]*U[j][0,:] - s_G[0,:]*U[j][2,:])*m)
            Mtg[2,j] = trapzs((-s_G[1,:]*U[j][0,:] + s_G[0,:]*U[j][1,:])*m)
    if main_axis=='x':
        Mtg[0,:] +=I_Jxx[:]
    elif main_axis=='z':
        Mtg[2,:] +=I_Jxx[:]

    #print('Mtg\n',Mtg)
        
    # --- Mgg  = \int Phi^t Phi dm  =  Sum Upsilon_kl(i,i)  Or: Me
    Mgg = np.zeros((nf,nf))
    if method=='OpenFAST' and U_untwisted is not None:
        U0=U_untwisted[:,:,:]
    else:
        U0=U[:,:,:]
    for i in range(nf):
        for j in range(nf): # NOTE: we could remove cross couplings here
            Mgg[i,j] = trapzs((U0[i][0,:]*U0[j][0,:] + U0[i][1,:]*U0[j][1,:] + U0[i][2,:]*U0[j][2,:])*m)

    # Adding torsion contribution if any
    Mgg=Mgg+np.diag(GMJxx)
    if bOrth:
        Mgg=Mgg*np.eye(nf)
    #print('Mgg\n',Mgg)

    # --- Build complete mass matrix
    MM = np.zeros((6+nf,6+nf))
    MM[:3,:3]   = Mxx; MM[:3,3:6] = Mxt; MM[:3,6:] = Mxg
    MM[3:6,3:6] = Mtt; MM[3:6,6:] = Mtg
    MM[6:,6:]   = Mgg

    i_lower     = np.tril_indices(len(MM), -1)
    MM[i_lower] = MM.T[i_lower]



    # --- Additional shape functions
    if rot_terms:
        # --- Gr_j = - 2*\int [~s] [~Phi_j]
        Gr = np.zeros((nf,3,3))

        # --- Oe_j = \int [~Phi_j] [~s] = { \int [~s] [~Phi_j] }^t = -1/2 *(Gr_j)^t
        Oe = np.zeros((nf,3,3))
        Oe6= np.zeros((nf,6))

        for j in range(nf):
            sxx = trapzs(s_G[0,:]*U[j][0,:]*m)
            sxy = trapzs(s_G[0,:]*U[j][1,:]*m)
            sxz = trapzs(s_G[0,:]*U[j][2,:]*m)
            syx = trapzs(s_G[1,:]*U[j][0,:]*m)
            syy = trapzs(s_G[1,:]*U[j][1,:]*m)
            syz = trapzs(s_G[1,:]*U[j][2,:]*m)
            szx = trapzs(s_G[2,:]*U[j][0,:]*m)
            szy = trapzs(s_G[2,:]*U[j][1,:]*m)
            szz = trapzs(s_G[2,:]*U[j][2,:]*m)
            Gr[j][0,:] = 2*np.array([ syy+szz, -syx  , -szx     ])
            Gr[j][1,:] = 2*np.array([ -sxy   ,sxx+szz, -szy     ])
            Gr[j][2,:] = 2*np.array([ -sxz   , -syz  , sxx+syy  ])

            Oe[j] = -0.5*Gr[j].T
            Oe6[j][0] = Oe[j][0,0]
            Oe6[j][1] = Oe[j][1,1]
            Oe6[j][2] = Oe[j][2,2]
            Oe6[j][3] = Oe[j][0,1] + Oe[j][1,0]
            Oe6[j][4] = Oe[j][1,2] + Oe[j][2,1]
            Oe6[j][5] = Oe[j][0,2] + Oe[j][2,0]

        # --- Ge_j = - 2*\int [Phi]^t  [~Phi_j]
        # [Phi]: 3xnf
        Ge = np.zeros((nf,nf,3))
        for j in range(nf):
            for k in range(nf):
                Ge[j][k,0] = -2*( trapzs(U[k][1,:]*U[j][2,:]*m) - trapzs(U[k][2,:]*U[j][1,:]*m))
                Ge[j][k,1] = -2*(-trapzs(U[k][0,:]*U[j][2,:]*m) + trapzs(U[k][2,:]*U[j][0,:]*m))
                Ge[j][k,2] = -2*( trapzs(U[k][0,:]*U[j][1,:]*m) - trapzs(U[k][1,:]*U[j][0,:]*m))	

    # --- M1 terms
    # Computing the M1 terms, this assumes that "s_G" is the undisplaced position!
    # The displaced position for each dof l is then s_G+ U[l]q_l with q_l=1
    if M1:
        OeM1 = np.zeros((nf,3,3,nf))
        Oe6M1= np.zeros((nf,6,nf))
        GrM1 = np.zeros((nf,3,3)) # we do not really store all of them
        for j in range(nf):
            for l in range(nf):
                sxx = trapzs((s_G[0,:]+U[l][0,:])*U[j][0,:]*m)
                sxy = trapzs((s_G[0,:]+U[l][0,:])*U[j][1,:]*m)
                sxz = trapzs((s_G[0,:]+U[l][0,:])*U[j][2,:]*m)
                syx = trapzs((s_G[1,:]+U[l][1,:])*U[j][0,:]*m)
                syy = trapzs((s_G[1,:]+U[l][1,:])*U[j][1,:]*m)
                syz = trapzs((s_G[1,:]+U[l][1,:])*U[j][2,:]*m)
                szx = trapzs((s_G[2,:]+U[l][2,:])*U[j][0,:]*m)
                szy = trapzs((s_G[2,:]+U[l][2,:])*U[j][1,:]*m)
                szz = trapzs((s_G[2,:]+U[l][2,:])*U[j][2,:]*m)
                GrM1[j][0,:] = 2*np.array([ syy+szz, -syx  , -szx     ])
                GrM1[j][1,:] = 2*np.array([ -sxy   ,sxx+szz, -szy     ])
                GrM1[j][2,:] = 2*np.array([ -sxz   , -syz  , sxx+syy  ])

                OeM1[j][l] = -0.5*GrM1[j].T
                Oe6M1[j][0][l] = OeM1[j][l][0,0]
                Oe6M1[j][1][l] = OeM1[j][l][1,1]
                Oe6M1[j][2][l] = OeM1[j][l][2,2]
                Oe6M1[j][3][l] = OeM1[j][l][0,1] + Oe[j][1,0]
                Oe6M1[j][4][l] = OeM1[j][l][1,2] + Oe[j][2,1]
                Oe6M1[j][5][l] = OeM1[j][l][0,2] + Oe[j][2,0]
# --- TODO
#         MxgM1      = np.zeros((3,nf,nf))
#         for j in range(nf):
#             for l in range(nf):
#                 MxgM1[0,j,l] = trapzs(U[j][0,:]*m)
#                 MxgM1[1,j,l] = trapzs(U[j][1,:]*m)
#                 MxgM1[2,j,l] = trapzs(U[j][2,:]*m)
    if split_outputs:
        if rot_terms:
            if M1:
                return Mxx, Mtt, Mxt, Mtg, Mxg, Mgg, Gr, Ge, Oe, Oe6, Oe6M1
            else:
                return Mxx, Mtt, Mxt, Mtg, Mxg, Mgg, Gr, Ge, Oe, Oe6
        else:
            return Mxx, Mtt, Mxt, Mtg, Mxg, Mgg
    else:
        if rot_terms:
            if M1:
                return MM, Gr, Ge, Oe, Oe6, Oe6M1
            else:
                return MM, Gr, Ge, Oe, Oe6
        else:
            return MM



def beamSectionLoads1D(z, p, Ftop=0, Mtop=0, s=1, F_lumped=None, method='plin'):
    r""" 
    Integrate section loads along a beam based on inline and lumped loads and top load.

    S(z) = int_z^L p(z') dz',   dS/dz = - p(z)
    M(z) =-int_z^L S(z') dz',   dM/dz = - S(z)

    Lumped forces can be inserted in the integral p= F_i \delta(x_i) with delta a Dirac

    Loads are integrated from "top" (z=L) to "bottom" (z=0)

    - z: array, linear station along the beam (typically from 0 to L)
    - p: array, inline load, assumed to go from z=0 to z=L
    - Ftop: Mtop: force and moment at top of Beam
    - s: sign +/-1
    - F: array, lumped forces along the span of the beam
    
    """
    n=len(z)
    Fsec=np.zeros(n)
    Msec=np.zeros(n)
    if F_lumped is None:
        F_lumped=np.zeros(n)

    if method=='plin':
        """ 
        Analytical results assuming linear variation of p:
            p(z) = (p_i-p_{i-1}) \tilde{z} + p_{i-1}
        """
        Fsec[-1]=Ftop 
        Msec[-1]=Mtop 
        for i in np.arange(n-2,-1,-1): # TODO vectorize me
            i1      = i+1
            dz      = z[i1]-z[i]
            Fsec[i] = Fsec[i1] +                 (p[i1]   +p[i] )/2.*dz   + F_lumped[i]
            Msec[i] = Msec[i1] + s*Fsec[i1]*dz+s*(p[i1]/3.+p[i]/6.)*dz*dz

    elif method=='manual':

        for i in np.arange(len(z)):
            zcur=z[i]
            Iabove = z>=zcur
            zabove = z[Iabove]
            pabove = p[Iabove]
            Fabove = F_lumped[Iabove]
            Fsec[i]=np.trapz(pabove, zabove) + Ftop + np.sum(Fabove)
            Msec[i]=np.trapz(pabove*(zabove-zcur), zabove)+Ftop*(z[-1]-zcur) + np.sum(Fabove*(zabove-zcur))
        Msec+=Mtop

    elif method=='cumtrapz':
        # NOTE: might not work as well when lumped forces are present
        zn = z[-1::-1]  # flip z so it goes from top to bottom for cumtrapz
        Fsec[:]  = Ftop
        Msec[:]  = Mtop
        Fsec      += np.cumsum(F_lumped[-1::-1])[-1::-1]
        Fsec[:-1] +=  - sciint.cumtrapz(p[-1::-1], zn)[-1::-1]    
        Msec[:-1] +=  - sciint.cumtrapz(Fsec[-1::-1], zn)[-1::-1] 

    else:
        raise NotImplementedError()

    return Fsec, Msec



def beamSectionLoads(x, xd, xdd, p_ext, F_top, M_top, s_span, PhiU, PhiV, m, 
        M_lumped=None, m_hydro=None, a_ext=None, F_ext_lumped=None, corrections=1):
    """ 
    
    """
    # Main dimensions
    shapeDisp = PhiU[0].shape
    nf        = len(PhiU)
    nSpan     = shapeDisp[1]

    # Default values
    if m_hydro is None:
        m_hydro  = np.zeros(nSpan)           # added mass, only on wet surface of structure
    if M_lumped is None:
        M_lumped = np.zeros(nSpan)
    if F_ext_lumped is None:
        F_ext_lumped=np.zeros(nSpan)

    # Linear Translation, Velocity, Acceleration
    U        = np.zeros(shapeDisp)
    V        = np.zeros(shapeDisp)
    #v_struct = np.zeros(shapeDisp)
    a_struct = np.zeros(shapeDisp)
    for j in np.arange(nf):
        U         += x  [j] * PhiU[j] # Deflections
        V         += x  [j] * PhiV[j] # Slopes
        #v_struct  += xd [j] * PhiU[j]
        a_struct  += xdd[j] * PhiU[j]
    if a_ext is not None:
        # Typically gravity
        # TODO Body root acceleration!
        a_struct[0,:] -= a_ext[0]
        a_struct[1,:] -= a_ext[1]
        a_struct[2,:] -= a_ext[2]

    #print('a2',a_struct[0,:])
    # --- Inertial loads
    try:
        m_struct = m
        m_tot = m_struct + m_hydro
        p_inertia        = m_tot    * a_struct # TODO is it really m_tot
        F_inertia_lumped = M_lumped * a_struct
    except:
        import pdb; pdb.set_trace()

    # --- Total loads from external forces and inertia (inline and lumped)
    p_all        = p_ext        - p_inertia 
    F_lumped_all = F_ext_lumped - F_inertia_lumped
#     print('p2',p_all[0,:])

    # --- Section Loads
    z  = s_span-s_span[0]
    zn = z[-1::-1]  # flip z so it goes from top to bottom for cumtrapz
    F_sec=np.zeros((3,len(z)))
    M_sec=np.zeros((3,len(z)))
    # Bending momemts 
    #print('p_inertia[0,:]',p_inertia[0,:])
    #print('p_all[0,:]',p_all[0,:])
    F_sec[0,:], M_sec[1,:] = beamSectionLoads1D(z, p_all[0,:], F_top[0], M_top[1], s=1,  F_lumped = F_lumped_all[0,:])
    F_sec[1,:], M_sec[0,:] = beamSectionLoads1D(z, p_all[1,:], F_top[1], M_top[0], s=-1, F_lumped = F_lumped_all[1,:])
    # Axial force
    F_sec[2,:-1] =- sciint.cumtrapz(p_all[2, -1::-1], zn)[-1::-1] # NOTE: mostly m*acc, can use FXG
    F_sec[2,:] += F_top[2] 
    # Torsion moment
    M_sec[2,:] += M_top[2]  # TODO integrate external torsions - torsional inertias and contributions from sectionn loads due to lever arm of deflection

    # Additional forces and moments from top loads due to deflections (ExtraLeverArm)
    if corrections>=1:
        F_sec[0,:] +=  F_top[2] * V[1,:] # Fx = Fz v_y 
        F_sec[1,:] +=  F_top[2] * V[0,:] # Fy = Fz v_x  # TODO check sign
        dx = U[0,-1] - U[0,:]
        dy = U[1,-1] - U[1,:]
        M_sec[1,:] += -F_top[2] * dx # My =-Fz dx 
        M_sec[0,:] += +F_top[2] * dy # Mx = Fz dy 
        M_sec[2,:] +=  F_top[1]*dx - F_top[0]*dy # Mz = Fy dx - Fx dy 

    # Torsion correction
    if corrections>=2:
        M_sec[2,1:] +=- V[1,1:]*M_sec[0,1:]-V[0,1:]*M_sec[1,1:] # Mx = - Vy Mx - Vx My # TODO check sign

    # KEEP ME: M_y approximation
    #M_sec[1,0] = F_top[1]*z[-1] + M_top[1] # approximation

    return F_sec, M_sec




def GeneralizedMCK_PolyBeam(s_span, m, EIFlp, EIEdg, coeffs, exp, damp_zeta, jxxG=None, gravity=None, Mtop=0, Omega=0, nSpan=None, bAxialCorr=False, bStiffening=True, main_axis='z', shapes=[0,1,2,3], algo=''):
    """ 
    Compute generalized mass, stiffness and damping matrix, and shape integrals for a beam defined using polynomial coefficients
    The shape functions of the beam are defined as:
       (U_j = a_ij s^i,  i=1..nExponents),  j=1..nShapes

    INPUTS:
       s_span    : spanwise coordinate (blade radius, or tower height)
       m         : mass per length [kg/m], array-like of length n
       EIFlp     : Flap EI, (FA EI for tower), array-like of length n
       EIEdg     : Edge EI  (SS EI for tower), array-like of length n
       coeffs    : coefficient used to define beam shape functions, array of shape nCoeffs x nShapes
       exp       : exponents used to defined beam shape functions, array of shape nExponent
       damp_zeta : damping ratio for the shape functions, array-like of length nShapes

    OPTIONAL INPUTS:
       jxxG      : cross section moment of inertia about G
       gravity   : acceleration of gravity
       Mtop      : mass on top of beam if any
       Omega     : rotational speed of the beam, if any [rad/s]
       nSpan     : number of spanwise station desired (will be interpolated)
    """
    def skew(x):
        x=np.asarray(x).ravel()
        """ Returns the skew symmetric matrix M, such that: cross(x,v) = M v """
        return np.array([[0, -x[2], x[1]],[x[2],0,-x[0]],[-x[1],x[0],0]])

    # --- Reading properties, coefficients
    nShapes = coeffs.shape[1]

    # --- Interpolating structural properties
    if jxxG is None:
        jxxG=m*0

    s_min = np.min(s_span)
    s_max = np.max(s_span)
    if algo=='ElastoDyn':
        # Settting up nodes like ElastoDyn
        s_span0 = s_span

        if nSpan is None:
            nSpan=len(s_span)
        length = s_span0[-1]-s_span0[0]
        fract  = np.arange(1./nSpan/2., 1, 1./nSpan)
        s_span = fract*length
        m       = np.interp(s_span, s_span0, m)
        EIFlp   = np.interp(s_span, s_span0, EIFlp)
        EIEdg   = np.interp(s_span, s_span0, EIEdg)
        jxxG    = np.interp(s_span, s_span0, jxxG)

    else:
        if nSpan is not None:
            if s_min>0:
                print('TODO flexibility.py, support s_span0/=0')
            s_span0 = s_span
            s_span  = np.linspace(0,np.max(s_span),nSpan)
            m       = np.interp(s_span, s_span0, m)
            EIFlp   = np.interp(s_span, s_span0, EIFlp)
            EIEdg   = np.interp(s_span, s_span0, EIEdg)
            jxxG    = np.interp(s_span, s_span0, jxxG)
    nSpan=len(s_span)

    # --- Definition of main directions
    ShapeDir=np.zeros(nShapes).astype(int)
    EI =np.zeros((3,nSpan))
    shapes = shapes[:nShapes]

    if main_axis=='x':
        iMain = 0  # longitudinal axis along x
         # First two shapes are along z (flapwise/Fore-Aft)
         # Third shape along along y (edgewise/Side-Side) Sign...
        shape2Dir = np.array([2,2,1,1])
        ShapeDir = shape2Dir[shapes]
        EI[2,:] = EIFlp
        EI[1,:] = EIEdg
    elif main_axis=='z':
        iMain = 2  # longitudinal axis along z
        # First two shapes are along x (flapwise/Fore-Aft)
        # Third shape along along y (edgewise/Side-Side) Sign...
        shape2Dir = np.array([0,0,1,1])
        ShapeDir = shape2Dir[shapes]
        EI[0,:] = EIFlp # EIFlp is EIy But formulae are reversed in GKBeam
        EI[1,:] = EIEdg
    else:
        raise NotImplementedError()

    # --- Undeflected shape
    s_P0          = np.zeros((3,nSpan))
    s_P0[iMain,:] = s_span 
    s_G0 = s_P0
    # TODO blade COG

    # --- Shape functions
    PhiU = np.zeros((nShapes,3,nSpan)) # Shape
    PhiV = np.zeros((nShapes,3,nSpan)) # Slope
    PhiK = np.zeros((nShapes,3,nSpan)) # Curvature
    for j in np.arange(nShapes):
        iAxis = ShapeDir[j]
        PhiU[j][iAxis,:], PhiV[j][iAxis,:], PhiK[j][iAxis,:] = polyshape(s_span,coeffs[:,j],exp, x_max=s_max)


    # --- Generalized stiffness
    KK0 = GKBeam(s_span, EI, PhiK, bOrth=False)
    if bStiffening:
        KKg     = GKBeamStiffnening(s_span, PhiV, gravity, m, Mtop, Omega, main_axis=main_axis)
        KKg_self= GKBeamStiffnening(s_span, PhiV, gravity, m, Mtop, Omega, main_axis=main_axis, bSelfWeight=True , bMtop=False, bRot=False)
        KKg_Mtop= GKBeamStiffnening(s_span, PhiV, gravity, m, Mtop, Omega, main_axis=main_axis, bSelfWeight=False, bMtop=True,  bRot=False)
        KKg_rot = GKBeamStiffnening(s_span, PhiV, gravity, m, Mtop, Omega, main_axis=main_axis, bSelfWeight=False, bMtop=False, bRot=True)
    else:
        KKg      = KK0*0
        KKg_self = KK0*0
        KKg_Mtop = KK0*0
        KKg_rot = KK0*0

    KK=KK0+KKg

    # --- Generalized mass
    MM, Gr, Ge, Oe, Oe6 = GMBeam(s_G0, s_span, m, PhiU, jxxG=jxxG, bUseIW=True, main_axis=main_axis, bAxialCorr=bAxialCorr, bOrth=False, rot_terms=True)

    # Beam COG
    s_COG = np.trapz(m*s_G0,s_span)/MM[0,0]

    # J at COG
    J_O  = MM[3:6, 3:6]
    r_PG = s_COG
    J_G  = J_O + MM[0,0] * np.dot(skew(r_PG), skew(r_PG))

    # --- Generalized damping
    DD = np.zeros((6+nShapes,6+nShapes))
    for j in range(nShapes):
        gm             = MM[6+j,6+j]
        gk             = KK[6+j,6+j]
        om            = np.sqrt(gk/gm)
        xi            = damp_zeta[j]*2*np.pi
        c             = xi * gm * om / np.pi
        DD[6+j,6+j] = c

    # --- alpha couplings 
    alpha = np.zeros((3,nShapes))
    for j in np.arange(nShapes):
        if main_axis=='x':
            alpha[0,j]=0                      # torsion
            alpha[1,j]=-PhiV[j][2,-1]
            alpha[2,j]= PhiV[j][1,-1]
        elif main_axis=='z':
            alpha[0,j]=-PhiV[j][1,-1]
            alpha[1,j]= PhiV[j][0,-1]
            alpha[2,j]= 0                     # torsion


    # --- Return dict
    return {'MM':MM, 'KK':KK, 'DD':DD, 'KK0':KK0, 'KKg':KKg, 'KKg_self':KKg_self,'KKg_Mtop':KKg_Mtop, 'KKg_rot':KKg_rot,
            'Oe':Oe, 'Oe6':Oe6, 'Gr':Gr, 'Ge':Ge,
            'PhiU':PhiU, 'PhiV':PhiV, 'PhiK':PhiK,
            's_P0':s_P0, 's_G':s_G0, 's_span':s_span, 's_min':s_min, 's_max':s_max,
            'm':m, 'EI':EI, 'jxxG':jxxG, 'ShapeDir':ShapeDir,
            's_OG':s_COG, 'J_G':J_G, 'alpha':alpha
            }



if __name__=='__main__':
    pass

