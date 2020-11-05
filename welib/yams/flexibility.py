import numpy as np
import scipy.integrate as sciint
'''
Flexible beam tools:
    - computation of generalized mass and stiffness matrix

Reference:
     [1]: Flexible multibody dynamics using joint coordinates and the Rayleigh-Ritz approximation: the general framework behind and beyond Flex
'''

def fcumtrapzlr(s_span, p):
    """ Cumulative trapezoidal integration, flipped left-right 
    Useful to return the following:
         P(x) = \int_x^R p(r) dr
    """
    P = - sciint.cumtrapz( p[-1::-1], s_span[-1::-1],)[-1::-1] 
    P = np.concatenate((P,[0]))
    return P


def polymode(x,coeff,exp):
    """ 
    Computes a shape function described as a polynomial expression y = a_i x^e_i
        where the a_i are given by `coeff`
              the e_i are given by `exp`
    The shape function is normalized such as to have a unitary tip deflection

    INPUTS:
        x : spanwise dimension, from 0 to L, not dimensionless!

    Returns:
        U, dU, ddU the shape, slope and curvature
    """
    mode   = np.zeros(x.shape)
    dmode  = np.zeros(x.shape)
    ddmode = np.zeros(x.shape)
    # Polynomials assume x to be dimensionless
    x_max= x[-1] 
    x_bar=x/x[-1] 
    for i in range(0,len(coeff)):
        mode += coeff[i]*x_bar**exp[i]
        if exp[i]-1>=0:
            dmode += coeff[i]*exp[i]* x_bar**(exp[i]-1)
        if exp[i]-2>=0:
            ddmode += coeff[i]*exp[i]*(exp[i]-1) * x_bar**(exp[i]-2)
    # Scaling by the tip deflection, and include x_max for derivatives since derivates were computed w.r.t. x_bar not x
    scale= mode[-1]
    return mode/scale, dmode/(scale*x_max), ddmode/(scale*x_max*x_max)


def integrationWeights(s_span,m):
    """ Returns integration weights convenient to integrate functions along the span of the beam
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


def GKBeamStiffnening(s_span, dU, gravity, m, Mtop, bSelfWeight=True, bMtop=True, main_axis='x'):
    """ 
       Computes geometrical stiffnening for a beam

    OUTPUTS:
     - KKg: Geometrical stiffeness matrix. 
    """
    nSpan = len(s_span)
    nf    = len(dU)
    KKg = np.zeros((6+nf,6+nf))
    # --- Axial force 
    Pacc_SW = fcumtrapzlr(s_span, -m * gravity)
    Pacc_MT = -Mtop * gravity*np.ones(nSpan)
    Pacc    = np.zeros(nSpan) 
    # TopMass contribution to Pacc
    if bMtop:
        Pacc=Pacc+Pacc_MT
    if bSelfWeight:
        Pacc=Pacc+Pacc_SW
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
    
def GMBeam(s_G, s_span, m, U=None, V=None, jxxG=None, bOrth=False, bAxialCorr=False, IW=None, IW_xm=None, main_axis='x', bUseIW=True, V_tot=None, Peq_tot=None):
    """
    Computes generalized mass matrix for a beam.
    Eq.(2) from [1]

    Performing full integration of mass matrix without shape integral functions
    NOTE: Beam assumed to be along x for now (only because of Jxx)
    
    INPUTS
     - s_G    : [m] 3 x nSpan , location of cross sections COG
     - s_span : [m] span integration variable (e.g. s_G(1,:))
     - m      : [kg/m] cross section mass along the beam
     - jxxG   : [kg.m] second moment of inertia of cross section # TODO


    OPTIONAL INPUTS:
     - bOrth : if true, enforce orthogonality of modes
     - JxxG, if omitted, assumed to be 0 # TODO
     - U , if omitted, then rigid body (6x6) mass matrix is returned
    
    """
    # Speed up integration along the span, using integration weight
    def trapzs(yy,**args):
        return np.sum(yy*IW)
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

    # --- Mxt = -\int [~s] dm    =  -Skew(sigma+Psi g)
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

    # --- Mxg = \int Phi dm  =   Psi
    Mxg      = np.zeros((3,nf))
    for j in range(nf):
        Mxg[0,j] = trapzs(U[j][0,:]*m)
        Mxg[1,j] = trapzs(U[j][1,:]*m)
        Mxg[2,j] = trapzs(U[j][2,:]*m)
    if bAxialCorr:
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
        
    # --- Mtt = - \int [~s][~s] dm
    if main_axis=='x':
        if bUseIW:
            s00= np.sum(IW_xm * s_G[0,:]);
            s01= np.sum(IW_xm * s_G[1,:]);
            s02= np.sum(IW_xm * s_G[2,:]);
        else:
            s00 = trapzs(s_G[0,:]*s_G[0,:]*m)
            s01 = trapzs(s_G[0,:]*s_G[1,:]*m)
            s02 = trapzs(s_G[0,:]*s_G[2,:]*m)

        s11 = trapzs(s_G[1,:]*s_G[1,:]*m)
        s12 = trapzs(s_G[1,:]*s_G[2,:]*m)
        s22 = trapzs(s_G[2,:]*s_G[2,:]*m)
    elif main_axis=='z':
        if bUseIW:
            s02= np.sum(IW_xm * s_G[0,:]);
            s12= np.sum(IW_xm * s_G[1,:]);
            s22= np.sum(IW_xm * s_G[2,:]);
        else:
            s02 = trapzs(s_G[2,:]*s_G[0,:]*m)
            s12 = trapzs(s_G[2,:]*s_G[1,:]*m)
            s22 = trapzs(s_G[2,:]*s_G[2,:]*m)

        s11 = trapzs(s_G[1,:]*s_G[1,:]*m)
        s00 = trapzs(s_G[0,:]*s_G[0,:]*m)
        s01 = trapzs(s_G[0,:]*s_G[1,:]*m)

    Mtt = np.zeros((3,3))
    Mtt[0,0] = s11 + s22    ;     Mtt[0,1] = -s01;       Mtt[0,2] = -s02
    Mtt[1,0] = -s01;              Mtt[1,1] = s00 + s22;  Mtt[1,2] = -s12
    Mtt[2,0] = -s02;              Mtt[2,1] = -s12;       Mtt[2,2] = s00+s11
    if main_axis=='x':
        Mtt[0,0] += Jxx
    else:
        Mtt[2,2] += Jxx
    #print('Mtt\n',Mtt)

    # --- Mtg  = \int [~s] Phi dm  
    Mtg      = np.zeros((3,nf))
    if main_axis=='x':
        if bUseIW:
            for j in range(nf):
                Mtg[0,j] = trapzs((-s_G[2,:]*U[j][1,:] + s_G[1,:]*U[j][2,:])*m) + I_Jxx[j]
                Mtg[1,j] = trapzs(  (+s_G[2,:]*U[j][0,:]*m)) - sum(IW_xm*U[j][2,:]);
                Mtg[2,j] = trapzs(  (-s_G[1,:]*U[j][0,:]*m)) + sum(IW_xm*U[j][1,:]);
        else:
            for j in range(nf):
                Mtg[0,j] = trapzs((-s_G[2,:]*U[j][1,:] + s_G[1,:]*U[j][2,:])*m) + I_Jxx[j]
                Mtg[1,j] = trapzs(( s_G[2,:]*U[j][0,:] - s_G[0,:]*U[j][2,:])*m)
                Mtg[2,j] = trapzs((-s_G[1,:]*U[j][0,:] + s_G[0,:]*U[j][1,:])*m)
    elif main_axis=='z':
        if bUseIW:
            for j in range(nf):
                Mtg[0,j] = trapzs((-sum(IW_xm*U[j][1,:]) + s_G[1,:]*U[j][2,:])*m) 
                Mtg[1,j] = trapzs(( sum(IW_xm*U[j][0,:]) - s_G[0,:]*U[j][2,:])*m)
                Mtg[2,j] = trapzs((-s_G[1,:]*U[j][0,:] + s_G[0,:]*U[j][1,:])*m)+ I_Jxx[j]
        else:
            for j in range(nf):
                Mtg[0,j] = trapzs((-s_G[2,:]*U[j][1,:] + s_G[1,:]*U[j][2,:])*m) 
                Mtg[1,j] = trapzs(( s_G[2,:]*U[j][0,:] - s_G[0,:]*U[j][2,:])*m)
                Mtg[2,j] = trapzs((-s_G[1,:]*U[j][0,:] + s_G[0,:]*U[j][1,:])*m)+ I_Jxx[j]

    #print('Mtg\n',Mtg)
        
    # --- Mgg  = \int Phi^t Phi dm  =  Sum Upsilon_kl(i,i)
    Mgg = np.zeros((nf,nf))
    for i in range(nf):
        for j in range(nf):
            Mgg[i,j] = trapzs((U[i][0,:]*U[j][0,:] + U[i][1,:]*U[j][1,:] + U[i][2,:]*U[j][2,:])*m)

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
    return MM



if __name__=='__main__':
    pass

