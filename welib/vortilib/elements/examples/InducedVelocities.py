import numpy as np
from scipy.special import ellipk, ellipe

# --------------------------------------------------------------------------------}
# --- Induced velocity functions 
# --------------------------------------------------------------------------------{
def vc_tang_u_doublet(Xcp,Ycp,Zcp,gamma_t=-1,R=1,r_bar_Cut=6,polar_out=True):
    Xcp=np.asarray(Xcp)
    shape_in=Xcp.shape
    Xcp=Xcp.ravel()
    Ycp=Ycp.ravel()
    Zcp=np.asarray(Zcp).ravel()
    Rcp  = np.sqrt(Xcp ** 2 + Ycp **2)
    Rsph = np.sqrt(Rcp ** 2 + Zcp **2)
    ur   = np.full(Xcp.shape,np.nan)
    uz   = np.full(Xcp.shape,np.nan)
    bCut = Rsph>r_bar_Cut*R
    dmz_dz = gamma_t * R**2 * np.pi # doublet intensity per length

    ur[~bCut],uz[~bCut] = vc_tang_u           (Xcp[~bCut],Ycp[~bCut],Zcp[~bCut],gamma_t=gamma_t,R=R,polar_out=True)
    ur[ bCut],uz[ bCut] = doublet_line_polar_u(Rcp[ bCut],Zcp[ bCut],dmz_dz)

    ur = ur.reshape(shape_in)   
    uz = uz.reshape(shape_in)   

    if polar_out:
        return ur,uz
    else:
        psi = np.arctan2(Ycp,Xcp)     ;
        ux=ur*np.cos(psi)
        uy=ur*np.sin(psi)
        return ux,uy,uz

    return ur, uz

def doublet_line_polar_u(rcp,zcp,dmz_dz):
    if np.any(rcp<0):
        raise Exception('Script meant for positive r')
    r=np.asarray(rcp)
    z=np.asarray(zcp)

    bZ0    = np.abs(z)<1e-16
    bR0    = np.abs(r)<1e-16
    bZ0R0  = np.logical_and(bZ0,bR0)
    bZ0Rp  = np.logical_and(bZ0, np.abs(r)>1e-16)
    bR0Zp  = np.logical_and(bR0, z>1e-16)
    bR0Zm  = np.logical_and(bR0, z<-1e-16)
    bOK    = np.logical_and(~bZ0,~bR0)

    uz=np.zeros(r.shape)
    ur=np.zeros(r.shape)

    norm2 = r**2+z**2
    uz[bOK]  = dmz_dz/(4*np.pi) * 1/r[bOK]**2 * ( z[bOK]**3/(norm2[bOK])**(3/2) - z[bOK]/(norm2[bOK])**(1/2) )
    uz[bZ0Rp] = 0
    uz[bR0Zm] = dmz_dz/(4*np.pi) * 1/norm2[bR0Zm]
    uz[bR0Zp] = dmz_dz/(4*np.pi) * 1/norm2[bR0Zp]
    ur[bOK]   =-dmz_dz/(4*np.pi) * r[bOK]          *  1/(norm2[bOK]  )**(3/2)
    ur[bZ0Rp] =-dmz_dz/(4*np.pi) * r[bZ0Rp]        *  1/(norm2[bZ0Rp])**(3/2)
    ur[bR0Zm] = 0
    ur[bR0Zp] = 0

    ur[bZ0R0] = 0
    uz[bZ0R0] = 0


    return ur, uz

def vc_tang_u(Xcp,Ycp,Zcp,gamma_t=-1,R=1,polar_out=True,epsilon=0):
    """ Induced velocity from a semi infinite cylinder extending along the z axis, starting at z=0
    INPUTS:
      Xcp,Ycp,Zcp: vector or matrix of control points coordinates
      gamma_t: tangential vorticity sheet strength of the cylinder
      R: cylinder radius
      epsilon : Regularization parameter, e.g. epsilon=0.0001*R
    Reference: [1,2],  in particular, equations (7-8) from [1]"""
    EPSILON_AXIS=1e-7; # relative threshold for using axis formula

    # --- Main corpus
    Xcp=np.asarray(Xcp)
    Ycp=np.asarray(Ycp)
    Zcp=np.asarray(Zcp)
    if Xcp.shape==(0,):
        if polar_out:
            return np.array([]), np.array([])
        else:
            return np.array([]),np.array([]),np.array([]),

    r = np.sqrt(Xcp ** 2 + Ycp ** 2)
    z = Zcp
    ur = np.full(r.shape,np.nan)
    uz = np.full(r.shape,np.nan)
    # Enforcing axis formula : v_z=-Gamma/(2) *( 1 + z / sqrt(R^2+z^2))
    Iz = r < (EPSILON_AXIS * R)
    ur[Iz] = 0
    uz[Iz] = gamma_t/2 * (1 + z[Iz] / np.sqrt(z[Iz]** 2 + R**2))

    # --- From this point on, variables have the size of ~Iz..
    bnIz = np.logical_not(Iz)
    r = r[bnIz]
    z = z[bnIz]

    # Eliptic integrals
    if epsilon==0:
        k_2  = 4 * r * R / ((R + r)**2 + z**2)
        k0_2 = 4 * r * R/  ((R + r)**2       )
    else:
        epsilon2= r*0+epsilon**2
        epsilon2[z<0]=0 # No regularization when z<0 # TODO
        k_2  = 4 * r * R / ((R + r)**2 + z**2 + epsilon2)
        k0_2 = 4 * r * R/  ((R + r)**2        + epsilon2)
    k = np.sqrt(k_2)
    EE = ellipe(k_2)
    KK = ellipk(k_2)
    #     PI = ellippi(k0_2,k_2)
    PI = ellipticPiCarlson(k0_2,k_2)
    # --- Special values
    PI[PI==np.inf]==0
    PI[r==R]=0 ; # when r==R, PI=0 TODO, check
    KK[KK==np.inf]=0 ; # when r==R, K=0  TODO, check
    # ---
    ur[bnIz] = -gamma_t/(2*np.pi) * np.multiply(np.sqrt(R/r) , np.multiply((2-k_2)/k,KK) - np.multiply(2.0/k, EE))
    # Term 1 has a singularity at r=R, # T1 = (R-r + np.abs(R-r))/(2*np.abs(R-r))
    T1=np.zeros(r.shape) 
    T1[r==R] = 1/2
    T1[r<R]  = 1
    if (epsilon!=0):
        # TODO, more work needed on regularization
        epsilon2= r*0+epsilon**2
        b=z>=0
        T1[b]=1/2*(1 + (R-r[b])*np.sqrt(1+epsilon2[b]/(R+r[b])**2)/np.sqrt((R-r[b])**2 +epsilon2[b]))
    uz[bnIz] = gamma_t/2*( T1 + np.multiply(np.multiply(z,k) / (2 * np.pi * np.sqrt(r * R)),(KK + np.multiply((R - r)/(R + r),PI))))
    
    if polar_out:
        return ur,uz
    else:
        psi = np.arctan2(Ycp,Xcp)     ;
        ux=ur*np.cos(psi)
        uy=ur*np.sin(psi)
        return ux,uy,uz


def ring_u_polar(r,z,Gamma=-1,r0=1,z0=0,epsilon=0, reg_method='Saffman'):
    """ 
    Compute the induced velocity from a vortex ring of radius r0, located at z0, 
    """
    def ring_u_polar_singular(r,z,Gamma=-1,r0=1):
        # Formulation from Yoon 2004
        a = np.sqrt((r+r0)**2 + (z)**2)
        m = 4 * r * r0 / (a ** 2)
        A = (z)**2 + r**2 + r0**2
        B = - 2*r*r0
        I1 = 4.0 / a     * ellipk(m)
        I2 = 4.0 / a**3 * ellipe(m) / (1 - m)
        ur = Gamma/(4*np.pi)*r0 * z/B *(I1 - A*I2)
        uz = Gamma/(4*np.pi)*r0*( (r0 + r*A/B)*I2 - r/B*I1)
        return ur,uz
    EPSILON = 1e-07 # small value used for axis formula and singularity (when epsilon=0)
    # --- Main corpus
    r=np.asarray(r)
    z=np.asarray(z)-z0 # Note: from now on, ring is centered on 0
    ur = np.full(r.shape,np.nan)
    uz = np.full(r.shape,np.nan)

    # Enforcing  Axis formula : v_z=-Gamma/(2r0) *1 / (1+(z/r0)^2)^(3/2)  
    Iz = r < (EPSILON * r0)
    ur[Iz] = 0
    uz[Iz] = Gamma/(2*r0)*(1.0/((1 +(z[Iz]/r0)**2)**(3.0/2.0)))

    # Value on the neighborhood of the ring itself..
    if epsilon==0:
        Ir = np.logical_and(np.abs(r-r0)<(EPSILON*r0), np.abs(z)<EPSILON)
        ur[Ir]=0
        uz[Ir]=Gamma/(4*r0) # NOTE: this is arbitrary
    else:
        if reg_method=='Saffman':
            # Regularized velocity near the ring (see Saffman)
            rho = np.sqrt( (r-r0)**2 + z**2 )
            Ir = rho < epsilon/2 
            ur[Ir]=0
            uz[Ir]=Gamma/(4*np.pi*r0)*(np.log(16*r0/epsilon)-1/4) # Eq 35.36 from [1]
        else:
            raise NotImplementedError()


    # --- From this point on, variables have the size of ~Iz..
    bnIz = np.logical_and(np.logical_not(Iz), np.logical_not(Ir))
    r = r[bnIz]
    z = z[bnIz]
    ur[bnIz], uz[bnIz] = ring_u_polar_singular(r,z,Gamma,r0)

    return ur,uz


def ellipticPiCarlson(n,m):
    # Elliptic integral of the third kind using the method of Carlson
    #    PI(n,m)=int(1/((1-n*sin(t)^2)*sqrt(1-m*sin(t)^2)),t=0,pi/2)
    # AUTHOR: N. Troldborg
    # REF: B.C. Carlson (1979) "Computing Elliptic Integrals by Duplication"
    # --- Performance parameters
    RES   = 1.e-12  # 
    ITMAX = 20      # 20 iterations usually sufficient
    # --- Subfunctions RF, RJ, RC
    def ellipticRF(y = None): 
        xo = np.zeros(y.shape)
        yo = y
        zo = np.ones(y.shape)
        nIt = 0
        res = 1
        RFo = np.zeros(yo.shape)
        while res > RES and nIt < ITMAX:
            lambda_ = (np.multiply(xo,yo))**0.5 + (np.multiply(xo,zo))**0.5 + (np.multiply(yo,zo))**0.5
            mu = (xo + yo + zo) / 3.
            xn = (xo + lambda_) / 4.
            yn = (yo + lambda_) / 4.
            zn = (zo + lambda_) / 4.
            X = 1 - xo / mu
            Y = 1 - yo / mu
            Z = 1 - zo / mu
            X2 = X**2
            X3 = np.multiply(X2,X)
            Y2 = Y**2
            Y3 = np.multiply(Y2,Y)
            Z2 = Z**2
            Z3 = np.multiply(Z2,Z)
            s1 = (X2 + Y2 + Z2) / 4
            s2 = (X3 + Y3 + Z3) / 6
            s12 = s1 ** 2
            s13 = np.multiply(s12,s1)
            r = 5./26 * s13 + 3./ 26 * s2 ** 2
            RF = np.multiply(mu ** - 0.5,(1 + s1 / 5 + s2 / 7 + s12 / 6 + np.multiply(3 / 11 * s1,s2) + r))
            res = np.amax(np.abs(RF - RFo))
            RFo = RF
            xo = xn
            yo = yn
            zo = zn
            nIt = nIt + 1
        return RF
    def ellipticRJ(y = None,rho = None): 
        b1 = rho > 0
        xt = np.zeros(y.shape)
        yt = y
        zt = np.ones(y.shape)
        rhot = rho
        # --- Dealing first with positive values
        xo = xt[b1]
        yo = yt[b1]
        zo = zt[b1]
        rhoo = rhot[b1]
        RJ = np.full(y.shape,np.inf)
        if np.any(b1):
            nIt = 0
            res = 1
            RJo  = np.zeros(xo.shape)
            rhs1 = np.zeros(xo.shape)
            while res > RES and nIt < ITMAX:

                lambda_ = (np.multiply(xo,yo))**0.5 + (np.multiply(xo,zo))**0.5 + (np.multiply(yo,zo))**0.5
                mu = (xo + yo + zo + 2 * rhoo) / 5.
                xn = (xo + lambda_) / 4.
                yn = (yo + lambda_) / 4.
                zn = (zo + lambda_) / 4.
                rhon = (rhoo + lambda_) / 4.
                X = 1 - xo / mu
                Y = 1 - yo / mu
                Z = 1 - zo / mu
                RHO = 1 - rhoo / mu
                X2 = X ** 2
                X3 = np.multiply(X2,X)
                X4 = np.multiply(X3,X)
                X5 = np.multiply(X4,X)
                Y2 = Y ** 2
                Y3 = np.multiply(Y2,Y)
                Y4 = np.multiply(Y3,Y)
                Y5 = np.multiply(Y4,Y)
                Z2 = Z ** 2
                Z3 = np.multiply(Z2,Z)
                Z4 = np.multiply(Z3,Z)
                Z5 = np.multiply(Z4,Z)
                RHO2 = RHO ** 2
                RHO3 = np.multiply(RHO2,RHO)
                RHO4 = np.multiply(RHO3,RHO)
                RHO5 = np.multiply(RHO4,RHO)
                s1 = (X2 + Y2 + Z2 + 2 * RHO2) / 4.
                s2 = (X3 + Y3 + Z3 + 2 * RHO3) / 6.
                s3 = (X4 + Y4 + Z4 + 2 * RHO4) / 8.
                s4 = (X5 + Y5 + Z5 + 2 * RHO5) / 10.
                s12 = np.multiply(s1,s1)
                s13 = np.multiply(s12,s1)
                r = - 1./10 * s13 + 3./10*s2** 2 + np.multiply(3./5 * s1,s3)
                alfa = (np.multiply(rhoo,(xo ** 0.5 + yo ** 0.5 + zo ** 0.5)) + (np.multiply(np.multiply(xo,yo),zo)) ** 0.5) ** 2
                bet = np.multiply(rhoo,(rhoo + lambda_) ** 2)
                rhs1 = rhs1 + 3 * 4 ** - nIt * ellipticRC(alfa,bet)
                rhs2 = np.multiply(4**-(nIt+1) * mu**-1.5, (1+3/7*s1 + s2/3 + 3/22*s12 + 3/11*s3 + 3/13*(np.multiply(s1,s2) + s4) + r))
                RJLoc = rhs1 + rhs2
                res = np.amax(np.abs(RJLoc - RJo))
                RJo = RJLoc
                xo = xn
                yo = yn
                zo = zn
                rhoo = rhon
                nIt = nIt + 1
            RJ[b1] = RJLoc
        return RJ
    def ellipticRC(x = None,y = None): 
        # Computes Carlson's Degenerate Elliptic Integral
        # RC(x,y)=1/2int_0^infty (t+x)^-0.5*(t+y)^-1dt
        # Carlson, B.C. (1979) "Computing Elliptic Integrals by Duplication"
        # AUTHOR: N. Troldborg
        nIt = 0
        res = 1
        xo = x
        yo = y
        RCo = np.zeros(x.shape)
        while res > RES and nIt < ITMAX:
            lambda_ = 2 * (np.multiply(xo,yo)) ** 0.5 + yo
            xn = (xo + lambda_) / 4
            yn = (yo + lambda_) / 4
            mu = (xo + 2 * yo) / 3
            s = (yo - xo) / (3 * mu)
            s2 = s ** 2
            s3 = np.multiply(s2,s)
            s4 = np.multiply(s3,s)
            s5 = np.multiply(s4,s)
            s6 = np.multiply(s5,s)
            RC = np.multiply(mu**-0.5 ,(1+3/10*s2 + s3/7 + 3/8*s4 + 9/22*s5 + 159/208*s6))
            res = np.amax(np.abs(RC - RCo))
            RCo = RC
            xo = xn
            yo = yn
            nIt = nIt + 1
        return RC

    # --- Main corpus
    if type(m) is not np.ndarray:
        m=np.array(m)
        n=np.array(n)
    if m.shape==(0,):
        return np.array([])
    RF = ellipticRF(1-m)
    RJ = ellipticRJ(1-m,1-n)
    PI = RF + np.multiply(1 / 3 * n,RJ)
    return PI
