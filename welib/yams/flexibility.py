import numpy as np
import scipy.integrate as sciint
'''
Flexible beam tools:
    - computation of generalized mass and stiffness matrix

Reference:
     [1]: Flexible multibody dynamics using joint coordinates and the Rayleigh-Ritz approximation: the general framework behind and beyond Flex
     [2]: Standard Input Data of Flexible Members in Multibody Systems. Wallrapp 1993
'''

def skew(x):
    x=np.asarray(x).ravel()
    """ Returns the skew symmetric matrix M, such that: cross(x,v) = M v """
    return np.array([[0, -x[2], x[1]],[x[2],0,-x[0]],[-x[1],x[0],0]])

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

def checkRegularNode(s_span):
    """ 
    Check if a regular discretization is used, typically used for OpenFAST method
    """
    n = len(s_span)-2
    L = s_span[-1]-s_span[0]
    ds = L/n
    dr = np.diff(s_span[1:-1])
    def round_sigdigit(ref, v, ndigit=2):
        nround=-int(np.floor(np.log10(abs(ref))))+ndigit
        return np.around(v, nround)
    #vds = np.unique(np.around(dr,5))
    vds = np.unique(round_sigdigit(ds, dr)) # rounding dr to a precision given by ds
    dr = np.concatenate( (dr, [ds]))
    if len(vds)>1:
        raise Exception('When using `OpenFAST` method, the user should provide inputs at 0, L, and mid nodes similar to OpenFAST ElastoDYn nodes')
    if vds[0]!=round_sigdigit(ds, ds):
        raise Exception('When using `OpenFAST` method, the mid nodes should have a spacing equal to L/n')
    return dr


def GKBeamStiffnening(s_span, dU, gravity, m, Mtop=0, Omega=0, bSelfWeight=True, bMtop=True, bRot=True, main_axis='x', method='trapz'):
    """ 
    TODO swap gravity and m
    TODO Use Mtop Omega as flag, would still need a flag for selfweight

       Computes geometrical stiffnening for a beam

    OUTPUTS:
     - KKg: Geometrical stiffeness matrix. 
    """
    if method=='OpenFAST':
        dr = checkRegularNode(s_span)
        # We will only use the inner nodes ("elements")
        dU     = dU[:,:,1:-1]
        s_span = s_span[1:-1] # NOTE: temporary, m shouldn't me used with this method
        m = m[1:-1] *dr   # Important HACK 
        # OpenFAST integration is simple mid-rule summation
        # NOTE: yy is hacked to include "dr" in it already
        def trapzs(yy):
            return np.sum(yy) 
    elif method=='trapz' or method=='Flex':
        def trapzs(yy,**args):
            return np.trapz(yy, s_span)
    else:
        raise NotImplementedError()


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
                yy=trapzs(Pacc * dU[i][1,:] * dU[j][1,:])
                zz=trapzs(Pacc * dU[i][2,:] * dU[j][2,:])
                KKCorr[i,j]=yy+zz
            elif main_axis=='z':
                xx=trapzs(Pacc * dU[i][0,:] * dU[j][0,:])
                yy=trapzs(Pacc * dU[i][1,:] * dU[j][1,:])
                KKCorr[i,j]=yy+xx
            else:
                raise Exception('Axis not suported')
    #print('KKCorr\n',KKCorr)
    KKg[6:,6:] = KKCorr
    return KKg


def GKBeamStiffneningSplit(s_G, s_span, dU, m, main_axis='x', method='trapz'):
    """ 
    TODO swap gravity and m
    TODO Use Mtop Omega as flag, would still need a flag for selfweight

       Computes geometrical stiffnening for a beam

    OUTPUTS:
     - KKg: Geometrical stiffeness matrix. 


    xx=trapz(s_span, Pacc .* PhiV{i}(1,:).* o.PhiV{j}(1,:));
     where Pacc is a normal load, typically accumulated over the blade length


    """

    DR  = np.diff(s_G) # TODO
    m0      = m.copy()
    s_span0 = s_span.copy()
    dU0     = dU.copy()
    if method=='OpenFAST':
        dr = checkRegularNode(s_span)
        DR=DR[:,1:] 
        DR[2,: ] =dr
        # We will only use the inner nodes ("elements")
        dU     = dU[:,:,1:-1]
        s_span = s_span[1:-1]
        s_G    = s_G[:, 1:-1]
        m      = m[1:-1]
        # OpenFAST integration is simple mid-rule summation
        def trapzs(yy):
            return np.sum(yy * dr)  # NOTE: adding dr

        def fcumtrapzlrs(s, p):
            p_cum= np.cumsum(p[-1::-1])[-1::-1] - 0.5*p
            return p_cum*dr


    elif method=='trapz' or method=='Flex':
        dr = np.ones(s_span.shape)

        def trapzs(yy,**args):
            return np.trapz(yy, s_span)

        fcumtrapzlrs = fcumtrapzlr

    else:
        raise NotImplementedError()

    nSpan = len(s_span)
    nf    = len(dU)

    # --- 
    p=dict()

    if main_axis=='x':
        i1=1
        i2=2
        i3=0
    elif main_axis=='z':
        i1=0
        i2=1
        i3=2
    else:
        raise Exception('Axis not suported')



    # --- Translational acceleration ("self weight" with unit gravity)
    # TODO not too sure what to do in other directions
    Pacc = fcumtrapzlrs(s_span, m)
    K0t = np.zeros((3,nf,nf)) # xyz, nf x nf
    for i in range(0,nf):
        for j in range(0,nf):
            K0t[i3, i,j]=trapzs(Pacc * (dU[i][i1,:] * dU[j][i1,:] + dU[i][i2,:] * dU[j][i2,:]))

    # --- Point loads 
    K0F = np.zeros((nSpan,3,nf,nf)) # nS, xyz, nf x nf
    Pacc = np.zeros(nSpan)
    for iS in range(0,nSpan):
        Pacc[:iS+1]=1
        for i in range(0,nf):
            for j in range(0,nf):
                K0F[iS, i3, i,j]=trapzs(Pacc * (dU[i][i1,:] * dU[j][i1,:] + dU[i][i2,:] * dU[j][i2,:]))

    # --- Rotational velocity K0omega
    I_acc = np.zeros((3,nSpan))
    if method=='OpenFAST':
        for i_elem in list(range(nSpan))[-1::-1]:
            I_acc[:,i_elem]=                           0.5*m[i_elem]  *(  s_G[:,i_elem  ] + 0.5*DR[:,i_elem  ])
            if i_elem<nSpan-1:
                I_acc[:, i_elem]+= I_acc[:,i_elem+1] + 0.5*m[i_elem+1]*(  s_G[:,i_elem+1] - 0.5*DR[:,i_elem+1])
        I_acc[0,:]*=dr
        I_acc[1,:]*=dr
        I_acc[2,:]*=dr
    else:
        I_acc[0,:] = fcumtrapzlrs(s_span,  m * s_G[0,:])
        I_acc[1,:] = fcumtrapzlrs(s_span,  m * s_G[1,:])
        I_acc[2,:] = fcumtrapzlrs(s_span,  m * s_G[2,:])
    # Using sum 
    #Kom = np.zeros((3,3, nf, nf))
    #for a in [0,1,2]:
    #    ea = np.zeros(3); ea[a]= 1;
    #    for b in [0,1,2]:
    #        eb = np.zeros(3); eb[b]= 1;
    #        for i_elem in range(nSpan):
    #            AA = -skew(ea).dot(skew(eb)).dot(I_acc[:,i_elem]).dot(DR[:,i_elem])
    #            for j in range(nf):
    #                for k in range(nf):
    #                    BB =  (sum( dU[j][:,i_elem]* dU[k][:,i_elem] )  )
    #                    Kom[a,b,j,k] += AA*BB
    # Using trapz 
    Pacc=np.zeros((3,3,nSpan))
    K0om = np.zeros((3,3, nf, nf))
    for a in [0,1,2]:
        ea = np.zeros(3); ea[a]= 1;
        for b in [0,1,2]:
            eb = np.zeros(3); eb[b]= 1;
            # Compute accumulation force
            for i_elem in range(nSpan):
                Pacc[a,b,i_elem] = -skew(ea).dot(skew(eb)).dot(I_acc[:,i_elem]).dot(DR[:,i_elem]/dr[i_elem])
            # Integrate
            for j in range(nf):
                for k in range(nf):
                    K0om[a,b,j,k] = trapzs(Pacc[a,b] * (dU[j][i1,:] * dU[k][i1,:] + dU[j][i2,:] * dU[k][i2,:]))

    p['K0t']     = K0t
    p['K0F']     = K0F
    p['K0omega'] = K0om
    return p 




def GKBeam(s_span, EI, ddU, bOrth=False, method='trapz'):
    """ 
       Computes generalized stiffness matrix for a beam
       Eq.(20) from [1]

    INPUTS:
     - s_span : array(n  )    span integration variable   [m]    
     - EI     : array(3 x n)  stiffnesses in each beam directions [Nm^2] 
                TODO torsion

    OPTIONAL INPUTS:
     - bOrth : if true, enforce orthogonality of modes
     - method: 'trapz', 'Flex', 'OpenFAST' (see GMBeam)

    OUTPUTS:
     - KK0: Stiffness matrix without geometrical stiffening
     - The total stiffness matrix should then be KK0+KKg
    """
    if method=='OpenFAST':
        dr = checkRegularNode(s_span)
        # We will only use the inner nodes ("elements")
        ddU   = ddU[:,:,1:-1]
        s_span = s_span[1:-1]  # NOTE: temporary, m shouldn't me used with this method
        EI = EI[:,1:-1] *dr   # Important HACK 
        # OpenFAST integration is simple mid-rule summation
        # NOTE: yy is hacked to include "dr" in it already
        def trapzs(yy):
            return np.sum(yy) 
    elif method=='trapz' or method=='Flex':
        def trapzs(yy,**args):
            return np.trapz(yy, s_span)
            #return np.sum(yy*IW) # NOTE: this is equivalent to trapezoidal integration
    else:
        raise NotImplementedError()

    nf = len(ddU)
    KK0 = np.zeros((6+nf,6+nf))
    Kgg = np.zeros((nf,nf))
    for i in range(0,nf):
        for j in range(0,nf):
            Kgg[i,j] = trapzs(EI[0,:]*ddU[i][0,:]*ddU[j][0,:] + EI[1,:]*ddU[i][1,:]*ddU[j][1,:] + EI[2,:]*ddU[i][2,:]*ddU[j][2,:])
    if bOrth:
        Kgg=Kgg*np.eye(nf)
    #print('Kgg\n',Kgg)
    KK0[6:,6:] = Kgg
    return KK0
    
def GMBeam(s_G, s_span, m, U=None, V=None, jxxG=None, bOrth=False, bAxialCorr=False, IW=None, IW_xm=None, main_axis='x', V_tot=None, Peq_tot=None, rot_terms=False, method='trapz', U_untwisted=None, M1=False):
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
    OUTPUTS:
      - MM: generalized mass matrix
      - IT: dictionary containing inertial terms
    """

    # --- Sanity check on method
    if method=='OpenFAST':
        dr = checkRegularNode(s_span)
        # We will only use the inner nodes ("elements")
        U   = U[:,:,1:-1]
        s_G = s_G[:,1:-1]
        #m      = m[1:-1]  # NOTE: temporary, m shouldn't me used with this method
        s_span = s_span[1:-1]  # NOTE: temporary, m shouldn't me used with this method
        m = m[1:-1]*dr # Important Hack we replace m by melem
        if jxxG is not None:
            jxxG = jxxG[1:-1]*dr # Important Hack
        if U_untwisted is not None:
            U_untwisted = U_untwisted[:,:,1:-1]
        if V is not None:
            V = V[:,:,1:-1]
    elif method=='trapz':
        pass
    elif method=='Flex':
        pass
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
    Mxt = -skew((C_x, C_y, C_z))
    #np.array([[0, C_z, -C_y],[-C_z, 0, C_x],[C_y, -C_x, 0]])

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
    if method=='Flex':
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
    if method=='Flex':
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


    IT=dict()
    IT['Mxx'] = Mxx
    IT['Mtt'] = Mtt
    IT['Mxt'] = Mxt
    IT['Mtg'] = Mtg
    IT['Mxg'] = Mxg
    IT['Mgg'] = Mgg

    # --- Additional shape functions
    if rot_terms:
        # --- Gr_j = - 2*\int [~s] [~Phi_j]
        #        [ 0  -z   y ]
        # [s~] = [ z   0  -x ]
        #        [-y   x   0 ]
        # 
        #           [ syy+szz, -syx  , -szx     ]
        #  Gr_j =  2[ -sxy   ,sxx+szz, -szy     ]
        #           [ -sxz   , -syz  , sxx+syy  ]
        #
        #           [-(syy+szz),    sxy    , sxz      ]
        #  Oe_j =   [   syx    ,-(sxx+szz), syz     ]
        #           [   szx    ,    szy   ,-(sxx+syy) ]
        #           
        #  Oe6_j=  [-(syy+szz), -(sxx+szz), -(sxx+syy), sxy+syx, syz+szy, sxz+szx] 
        #
        # NOTE: for a straight blade along z:
        #     s_Gx=0 (so sxx,sxy,sxz=0), 
        #     s_Gy=0 (so syx,syy,syz=0) 
        #     and no axial deflection Uz=0  (sxz,syz szz=0)
        #        
        #          [ 0   0  -szx ]
        #  Gr_j = 2[ 0   0  -szy ]
        #          [ 0   0   0   ]
        # 
        #          [ 0   0     0 ]
        #  Oe_j =  [ 0   0     0 ]
        #          [ szx szy   0 ]
        # 
        #  Oe6_j= [0, 0, 0, 0, szy, szx]
        # 
        # NOTE: for M1 we use s=Uj, then for a straight blade along z:
        #           we mostly have Uz=0 (now meaning: szx, szy, szz=0 and sxz,syz,szz=0)
        # 
        #  Oe6_j=  [-syy, -sxx, -(sxx+syy), sxy+syx, 0, 0] 
        # 
        # 
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

        # --- Ge_j = - 2*\int [Phi]^t [~Phi_j] dm
        #     Ge_j =   2*\int [Phi]^t [~Phi_j]^t dm  = [2C5_jk']
        #     Ge_j =   2 Kr
        # [Phi]: 3xnf
        # TODO revisit this
        Ge = np.zeros((nf,nf,3))
        for j in range(nf):
            for k in range(nf):
                Ge[j][k,0] =-2*( trapzs(U[k][1,:]*U[j][2,:]*m) - trapzs(U[k][2,:]*U[j][1,:]*m))
                Ge[j][k,1] =-2*(-trapzs(U[k][0,:]*U[j][2,:]*m) + trapzs(U[k][2,:]*U[j][0,:]*m))
                Ge[j][k,2] =-2*( trapzs(U[k][0,:]*U[j][1,:]*m) - trapzs(U[k][1,:]*U[j][0,:]*m))	

        ## --- C5_jk = \int [~Phi_j^0][Phi_k^0]
        #IT['C5'] = np.zeros((nf,3,nf))
        #for j in range(nf):
        #    for k in range(nf):
        #        IT['C5'][j][0,k]= trapzs( (-U[j][2,:]*U[k][1,:] + U[j][1,:]*U[k][2,:])*m)
        #        IT['C5'][j][1,k]= trapzs( ( U[j][2,:]*U[k][0,:] - U[j][0,:]*U[k][2,:])*m)
        #        IT['C5'][j][2,k]= trapzs( (-U[j][1,:]*U[k][0,:] + U[j][0,:]*U[k][1,:])*m)


        IT['Ge']  = Ge
        IT['Gr']  = Gr
        IT['Oe']  = Oe
        IT['Oe6'] = Oe6

    # --- M1 terms
    # Computing the M1 terms, this assumes that "s_G" is the undisplaced position!
    # The displaced position for each dof l is then s_G+ U[l]q_l with q_l=1
    if M1:
        # --- M1 Mxt = -\int [~s] dm  = -Skew(sigma+Psi g)  Or: +/- Skew(mdCM)
        # mdCM1 = C1n = int phi_n dm 
        mdCM_M1 = np.zeros((3,nf))
        Mxt_M1 = np.zeros((3,3,nf))
        for j in range(nf):
            mdCM_M1[0,j] = trapzs(U[j][0,:]*m)
            mdCM_M1[1,j] = trapzs(U[j][1,:]*m)
            mdCM_M1[2,j] = trapzs(U[j][2,:]*m)
            Mxt_M1[:,:,j] = - skew(mdCM_M1[:,j])

        # --- M1  Mtt = - \int [~s][~s] dm  - Or: J, Mrr
        Mtt_M1 = np.zeros((3,3,nf))
        for j in range(nf):
            s00 = trapzs(U[j][0,:]*s_G[0,:]*m)
            s01 = trapzs(U[j][0,:]*s_G[1,:]*m)
            s02 = trapzs(U[j][0,:]*s_G[2,:]*m)
            s11 = trapzs(U[j][1,:]*s_G[1,:]*m)
            s12 = trapzs(U[j][1,:]*s_G[2,:]*m)
            s22 = trapzs(U[j][2,:]*s_G[2,:]*m)
            Mtt_M1[0,0,j] = s11 + s22    ;     Mtt_M1[0,1,j] = -s01;       Mtt_M1[0,2,j] = -s02
            Mtt_M1[1,0,j] = -s01;              Mtt_M1[1,1,j] = s00 + s22;  Mtt_M1[1,2,j] = -s12
            Mtt_M1[2,0,j] = -s02;              Mtt_M1[2,1,j] = -s12;       Mtt_M1[2,2,j] = s00+s11
            #if main_axis=='x':
            #    Mtt[0,0] += Jxx
            #else:
            #    Mtt[2,2] += Jxx

        # --- M1 Mxg = \int Phi dm  Or:  Psi , Ct^T
        # Ct1n= N1nk^T 
        # NOTE: N1nk requires phi1 which we don't have...
        # TODO This should be geometrical stiffening, K0t
        Mxg_M1 = np.zeros((3,nf,nf))
        #for j in range(nf):
        #    for l in range(nf):
        #        Mxg_M1[0,j,l] = trapzs(U[j][0,:]*m)
        #        Mxg_M1[1,j,l] = trapzs(U[j][1,:]*m)
        #        Mxg_M1[2,j,l] = trapzs(U[j][2,:]*m)

        # --- M1 Mtg  = \int [~s] Phi dm -  Or: Mrg, Cr^T
        # See Kr 
        Mtg_M1 = np.zeros((3,nf,nf))
        for j in range(nf):
            for l in range(nf):
                Mtg_M1[0,j,l] = trapzs((-U[l][2,:]*U[j][1,:] + U[l][1,:]*U[j][2,:])*m)
                Mtg_M1[1,j,l] = trapzs(( U[l][2,:]*U[j][0,:] - U[l][0,:]*U[j][2,:])*m)
                Mtg_M1[2,j,l] = trapzs((-U[l][1,:]*U[j][0,:] + U[l][0,:]*U[j][1,:])*m)


        # --- Gr/Oe M1
        Oe_M1 = np.zeros((nf,3,3,nf))
        Oe6_M1= np.zeros((nf,6,nf))
        Gr_M1 = np.zeros((nf,3,3,nf)) 
        SS_M1=np.zeros((nf,nf,3,3))
        for j in range(nf):
            for l in range(nf):
                # NOTE: we remove s_G as this is second order effect, and we would need Phi1
                # Look at Wallrap 1993 [2]
                sxx = trapzs((         U[l][0,:])*U[j][0,:]*m)
                sxy = trapzs((         U[l][0,:])*U[j][1,:]*m)
                sxz = trapzs((         U[l][0,:])*U[j][2,:]*m)
                syx = trapzs((         U[l][1,:])*U[j][0,:]*m)
                syy = trapzs((         U[l][1,:])*U[j][1,:]*m)
                syz = trapzs((         U[l][1,:])*U[j][2,:]*m)
                szx = trapzs((         U[l][2,:])*U[j][0,:]*m)
                szy = trapzs((         U[l][2,:])*U[j][1,:]*m)
                szz = trapzs((         U[l][2,:])*U[j][2,:]*m)
                s=np.array([[ sxx, sxy, sxz], [ syx, syy, syz], [ szx, szy, szz]])
                SS_M1[j,l,:,:] = s

                Gr_M1[j,0,:,l] = 2*np.array([ syy+szz, -syx  , -szx     ])
                Gr_M1[j,1,:,l] = 2*np.array([ -sxy   ,sxx+szz, -szy     ])
                Gr_M1[j,2,:,l] = 2*np.array([ -sxz   , -syz  , sxx+syy  ])

                Oe_M1[j,:,:,l] = -0.5*Gr_M1[j,:,:,l].T
                Oe6_M1[j,0,l] = Oe_M1[j,0,0,l]
                Oe6_M1[j,1,l] = Oe_M1[j,1,1,l]
                Oe6_M1[j,2,l] = Oe_M1[j,2,2,l]
                Oe6_M1[j,3,l] = Oe_M1[j,0,1,l] + Oe_M1[j,1,0,l]
                Oe6_M1[j,4,l] = Oe_M1[j,1,2,l] + Oe_M1[j,2,1,l]
                Oe6_M1[j,5,l] = Oe_M1[j,0,2,l] + Oe_M1[j,2,0,l]
                # Oe6M1[j,:,l] = [-(s[1,1]+s[2,2]), -(s[0,0]+s[2,2]), -(s[0,0]+s[1,1]), s[0,1]+s[1,0], s[1,2]+s[2,1], s[0,2]+s[2,0] ]
        

        IT['mdCM_M1'] = mdCM_M1
        IT['Mxt_M1']  = Mxt_M1
        IT['Mtt_M1']  = Mtt_M1
        IT['Mxg_M1']  = Mxg_M1
        IT['Mtg_M1']  = Mtg_M1
        IT['Gr_M1']   = Gr_M1
        IT['Oe6_M1']  = Oe6_M1
        IT['SS_M1']   = SS_M1

    return MM, IT


def shapeIntegrals(s_G, s_span, m, U, dU, ddU, method='trapz', EI=None):
    """
    Compute shape integrals according to Wallrapp 1993 [2]


    INPUTS
     - s_G    : [m] 3 x nSpan , location of cross sections COG
     - s_span : [m] span integration variable (e.g. s_G(1,:))
     - m      : [kg/m] cross section mass along the beam
     - U      : nf x 3 x nSpan, shape functions
     - dU     : nf x 3 x nSpan, shape functions slopes
     - ddU    : nf x 3 x nSpan, shape functions curvatures


    Wallrapp notations:
         r_P = r_O + d
         d   = c + u   "c" is "s_G0"
         u     = Phi_u q , Phi_u =[Phi_ul]  Phi_ul = Phi_l^0 + 0.5 Phi^1_ln qn
         udot  = Phi q ,   Phi   =[Phi_l]   Phi_l  = Phi_l^0 +     Phi^1_ln qn

           [ 0  -z   y ]
    [s~] = [ z   0  -x ]
           [-y   x   0 ]

    """
    DR  = np.diff(s_G) # TODO
    m0      = m.copy()
    s_span0 = s_span.copy()
    dU0     = dU.copy()

    if method=='OpenFAST':
        DR=DR[:,1:] # TODO
        dr = checkRegularNode(s_span)
        DR[2,: ] =dr
        # We will only use the inner nodes ("elements")
        U      =   U[:,:,1:-1]
        dU     =  dU[:,:,1:-1]
        ddU    = ddU[:,:,1:-1]
        s_G    = s_G[:,1:-1]
        s_span = s_span[1:-1]
        m = m[1:-1]*dr # Important Hack we replace m by melem
        if EI is not None:
            EI=EI[:,1:-1]*dr # Important Hack
        # OpenFAST integration is simple mid-rule summation
        # NOTE: yy is hacked to include "dr" in it already
        def trapzs(yy):
            return np.sum(yy) 
    else:
        # Speed up integration along the span, using integration weight
        def trapzs(yy,**args):
            #return np.sum(yy*IW) # NOTE: this is equivalent to trapezoidal integration
            return np.trapz(yy, s_span)


    p     = dict()
    nf    = len(U)
    nSpan = len(s_span)

    # --- Mass
    p['mass'] = trapzs(m)

    # --- Center of mass location
    p['mc0'] = np.array([trapzs(s_G[0,:]*m), trapzs(s_G[1,:]*m), trapzs(s_G[2,:]*m)])
        
    # --- J0
    s00 = trapzs(s_G[0,:]*s_G[0,:]*m)
    s01 = trapzs(s_G[0,:]*s_G[1,:]*m)
    s02 = trapzs(s_G[0,:]*s_G[2,:]*m)
    s11 = trapzs(s_G[1,:]*s_G[1,:]*m)
    s12 = trapzs(s_G[1,:]*s_G[2,:]*m)
    s22 = trapzs(s_G[2,:]*s_G[2,:]*m)
    p['J0'] = np.zeros((3,3))
    p['J0'][0,0] = s11 + s22    ;     p['J0'][0,1] = -s01;       p['J0'][0,2] = -s02
    p['J0'][1,0] = -s01;              p['J0'][1,1] = s00 + s22;  p['J0'][1,2] = -s12
    p['J0'][2,0] = -s02;              p['J0'][2,1] = -s12;       p['J0'][2,2] = s00+s11

    # --- C1_k = \int Phi^0_k dm
    # 3 x nf
    p['C1'] = np.zeros((3, nf))
    for j in range(nf):
        p['C1'][0,j]= trapzs(m*U[j][0,:])
        p['C1'][1,j]= trapzs(m*U[j][1,:])
        p['C1'][2,j]= trapzs(m*U[j][2,:])

    # --- C2_k = \int [~c] [Phi^0_k] dm
    # 3 x nf
    p['C2'] = np.zeros((3, nf))
    for j in range(nf):
        p['C2'][0,j]= trapzs( (-s_G[2,:]*U[j][1,:] + s_G[1,:]*U[j][2,:])*m)
        p['C2'][1,j]= trapzs( ( s_G[2,:]*U[j][0,:] - s_G[0,:]*U[j][2,:])*m)
        p['C2'][2,j]= trapzs( (-s_G[1,:]*U[j][0,:] + s_G[0,:]*U[j][1,:])*m)

    # --- C3_kl = \int [Phi^0_k]^T [phi^0_l] dm
    # nf x nf x 3 x 3
    # 3 x 3 x nf x nf
    p['C3'] = np.zeros((3,3,nf,nf))
    for i in range(nf):
        for j in range(nf):
            for a in [0,1,2]:
                for b in [0,1,2]:
                    #p['C3'][i,j,a,b] = trapzs(U[i][a,:]*U[j][b,:]) 
                    p['C3'][a,b,i,j] = trapzs(U[i][a,:]*U[j][b,:]*m) 

    # --- C4_k = \int [~c] [~Phi^0_k] dm
    # nf x 3 x 3 
    p['C4'] = np.zeros((nf,3,3))
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
        p['C4'][j][0,:] = -np.array([ syy+szz, -syx  , -szx     ])
        p['C4'][j][1,:] = -np.array([ -sxy   ,sxx+szz, -szy     ])
        p['C4'][j][2,:] = -np.array([ -sxz   , -syz  , sxx+syy  ])

    # --- C5_kl = \int [~Phi_k^0][Phi_l^0]
    p['C5'] = np.zeros((nf,3,nf))
    for j in range(nf):
        for k in range(nf):
            p['C5'][j][0,k]= trapzs( (-U[j][2,:]*U[k][1,:] + U[j][1,:]*U[k][2,:])*m)
            p['C5'][j][1,k]= trapzs( ( U[j][2,:]*U[k][0,:] - U[j][0,:]*U[k][2,:])*m)
            p['C5'][j][2,k]= trapzs( (-U[j][1,:]*U[k][0,:] + U[j][0,:]*U[k][1,:])*m)

    # --- C6_kl = \int [~Phi_k^0][~Phi_l^0]
    p['C6'] = np.zeros((nf,nf,3,3))
    for j in range(nf):
        for k in range(nf):
            sxx = trapzs(U[j][0,:]*U[k][0,:]*m)
            sxy = trapzs(U[j][0,:]*U[k][1,:]*m)
            sxz = trapzs(U[j][0,:]*U[k][2,:]*m)
            syx = trapzs(U[j][1,:]*U[k][0,:]*m)
            syy = trapzs(U[j][1,:]*U[k][1,:]*m)
            syz = trapzs(U[j][1,:]*U[k][2,:]*m)
            szx = trapzs(U[j][2,:]*U[k][0,:]*m)
            szy = trapzs(U[j][2,:]*U[k][1,:]*m)
            szz = trapzs(U[j][2,:]*U[k][2,:]*m)
            p['C6'][j,k][0,:] = -np.array([ syy+szz, -syx  , -szx     ])
            p['C6'][j,k][1,:] = -np.array([ -sxy   ,sxx+szz, -szy     ])
            p['C6'][j,k][2,:] = -np.array([ -sxz   , -syz  , sxx+syy  ])

    # Komega
    p['Komega'] = np.zeros((3,3,nf,nf))
    for a in [0,1,2]:
        for b in [0,1,2]:
            if a!=b:
                p['Komega'][a, b,:,:]= p['C3'][a, b,:,:]
            else:
                c= np.mod(a+1, 3)
                d= np.mod(c+1, 3)
                p['Komega'][a, b, :, :]= -p['C3'][c, c, :,:] - p['C3'][d, d,:,:]
    p['Komega_b'] = np.zeros((3,3,nf,nf))
    for a in [0,1,2]:
        for b in [0,1,2]:
            for j in range(nf):
                for k in range(nf):
                    p['Komega_b'][a,b,j,k] = p['C6'][j,k,b,a] # Somehow transpose
    # Kr
    p['Kr'] = np.zeros((3,nf,nf))
    for a in [0,1,2]:
        c = np.mod(a+1, 3)
        d = np.mod(c+1, 3)
        p['Kr'][a]= -p['C3'][c, d] + p['C3'][c, d].T


    # --- Geometric stiffnesses
    # stiffening due to point load
    if method=='OpenFAST':
        p['K0F'] = np.zeros((nSpan, 3, nf, nf))
        for i_elem in range(nSpan):
            if i_elem>0:
                p['K0F'][i_elem,:,:,:]= p['K0F'][i_elem-1,:,:,:]
            else:
                p['K0F'][i_elem,:,:,:]= np.zeros((3, nf, nf))
            for a in [0,1,2]:
                for j in range(nf):
                    for k in range(nf):
                        p['K0F'][i_elem,a, j, k] += DR[a,i_elem]*(sum( dU[j][:,i_elem]* dU[k][:,i_elem] )  )
    else:
        print('TODO, shape integral K0F with trapz')

    # stiffening due to acceleration
    m_cum= np.cumsum(m[-1::-1])[-1::-1] - 0.5*m
    p['K0t'] = np.zeros((3,nf,nf))
    for i_elem in range(nSpan):
        for a in [0,1,2]:
            for j in range(nf):
                for k in range(nf):
                    p['K0t'][a,j,k]+=  m_cum[i_elem]*DR[a,i_elem]*(sum( dU[j][:,i_elem]* dU[k][:,i_elem] )  )

    # stiffening due to rotational acceleration
    p['K0r'] = np.zeros((3, 1, nf, nf))

    # stiffening due to rotational velocity (centrifugal stiffening)
    # force due to rotation results from cumulated moment of inertia I_cum
    # for i_elem= n_elem:-1:1
    #Kgrot= GKBeamStiffnening(s_span0, dU0, 0, m0, Mtop=0, Omega=1, bSelfWeight=False, bMtop=False, bRot=True, main_axis='z', method='OpenFAST')
    #p['Kgrot'] = Kgrot
    I_cum = np.zeros((nSpan,3))
    for i_elem in list(range(nSpan))[-1::-1]:
        I_cum[i_elem]= 0.5*m[i_elem] * (s_G[:,i_elem] + 0.5*DR[:,i_elem])
        if i_elem<nSpan-1:
            I_cum[i_elem,:]+= I_cum[i_elem+1] + 0.5*m[i_elem+1]*(  s_G[:,i_elem+1] - 0.5*DR[:,i_elem+1])

    # 
    p['K0omega'] = np.zeros((3,3, nf, nf))
    for i_elem in range(nSpan):
        for a in [0,1,2]:
            for b in [0,1,2]:
                w1 = np.zeros(3)
                w2 = np.zeros(3)
                w1[a]= 1;
                w2[b]= 1;
                AA = -skew(w1).dot(skew(w2).dot(I_cum[i_elem, :])).dot(DR[:,i_elem])
                for j in range(nf):
                    for k in range(nf):
                        BB =  (sum( dU[j][:,i_elem]* dU[k][:,i_elem] )  )
                        p['K0omega'][a,b,j,k] += AA*BB


    # --------------------------------------------------------------------------------}
    # --- Not strictly shape integrals..
    # --------------------------------------------------------------------------------{
    p['Me'] = p['C3'][0, 0] + p['C3'][1, 1] + p['C3'][2, 2]
    # stiffness and damping
    if EI is not None:
        p['Ke']= np.zeros((nf,nf))
        for i in range(0,nf):
            for j in range(0,nf):
                p['Ke'][i,j] = trapzs(EI[0,:]*ddU[i][0,:]*ddU[j][0,:] + EI[1,:]*ddU[i][1,:]*ddU[j][1,:] + EI[2,:]*ddU[i][2,:]*ddU[j][2,:])
    else:
        p['Ke'] = None
    #for i_elem in range(nSpan):
    #    pass
        #stiff= diag(IE{i_elem})*norm(DR{i_elem});
        #Ke= Ke + ddPhi{i_elem}' * stiff * ddPhi{i_elem};
    #EF= sqrt(diag(Ke)./diag(Me))/2/pi;
    #De= Ke*diag(damping(:)./(pi*EF));
    # --- Store inputs
    p['s_G'] = s_G
    p['U']   = U
    p['dU']  = dU
    p['ddU'] = ddU

    return p



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




def GeneralizedMCK_PolyBeam(s_span, m, EIFlp, EIEdg, coeffs, exp, damp_zeta, jxxG=None, gravity=None, Mtop=0, Omega=0, nSpan=None, 
        bAxialCorr=False, bStiffening=True, main_axis='z', shapes=[0,1,2,3], algo='', s_start=0):
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

    # --- Reading properties, coefficients
    nShapes = coeffs.shape[1]

    # --- Interpolating structural properties
    if jxxG is None:
        jxxG=m*0

    s_min = np.min(s_span)
    s_max = np.max(s_span)
    if algo=='OpenFAST':
        # Settting up nodes like ElastoDyn
        method= 'OpenFAST'
        s_span0 = s_span

        if nSpan is None:
            nSpan=len(s_span)
        length = s_span0[-1]-s_span0[0]
        fract  = np.arange(1./nSpan/2., 1, 1./nSpan)
        s_span = np.concatenate([[0],fract,[1]])*length
        m       = np.interp(s_span, s_span0, m)
        EIFlp   = np.interp(s_span, s_span0, EIFlp)
        EIEdg   = np.interp(s_span, s_span0, EIEdg)
        jxxG    = np.interp(s_span, s_span0, jxxG)
        s_span  = s_span + s_start
        s_min = np.min(s_span)
        s_max = np.max(s_span)

    else:
        method= 'Flex'
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
    KK0 = GKBeam(s_span, EI, PhiK, bOrth=False, method=method)
    if bStiffening:
        KKg     = GKBeamStiffnening(s_span, PhiV, gravity, m, Mtop, Omega, main_axis=main_axis, method=method)
        KKg_self= GKBeamStiffnening(s_span, PhiV, gravity, m, Mtop, Omega, main_axis=main_axis, bSelfWeight=True , bMtop=False, bRot=False, method=method)
        KKg_Mtop= GKBeamStiffnening(s_span, PhiV, gravity, m, Mtop, Omega, main_axis=main_axis, bSelfWeight=False, bMtop=True,  bRot=False, method=method)
        KKg_rot = GKBeamStiffnening(s_span, PhiV, gravity, m, Mtop, Omega, main_axis=main_axis, bSelfWeight=False, bMtop=False, bRot=True , method=method)
    else:
        KKg      = KK0*0
        KKg_self = KK0*0
        KKg_Mtop = KK0*0
        KKg_rot = KK0*0

    KK=KK0+KKg

    # --- Generalized mass
    MM, IT = GMBeam(s_G0, s_span, m, PhiU, jxxG=jxxG, main_axis=main_axis, bAxialCorr=bAxialCorr, bOrth=False, rot_terms=True, method=method)
    Gr, Ge, Oe, Oe6 = IT['Gr'], IT['Ge'], IT['Oe'], IT['Oe6']

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
        if gk<0:
            print('[WARN] Flexibility: Shape function has negative stiffness (likely due to geometrical stiffening)')
            om = 0
        else:
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

