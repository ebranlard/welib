import numpy as np
import unittest
    

# --------------------------------------------------------------------------------}
# --- Point doublet 
# --------------------------------------------------------------------------------{
def doublet_u(Xcp,Ycp,Zcp,m,x0=[0,0,0]): 
    """
    Velocity field induced by one doublet located at x0 with vector intensity m
         u = 1/(4pi)       *( 3* r * m.r/|r|**5 - m/|r|**3)
         u = 1/(4pi |r|**3)*( 3* r * m.r/|r|**2 - m)
    Xcp, Ycp, Zcp: Cartesian coordinates of control points
    """
    
    Xcp=np.asarray(Xcp)
    Ycp=np.asarray(Ycp)
    Zcp=np.asarray(Zcp)
    
    DX = Xcp - x0[0]
    DY = Ycp - x0[1]
    DZ = Zcp - x0[2]

    r_norm3 = (DX**2 + DY**2 + DZ**2 )**(3/2) # |r|**3
    r_norm3=np.asarray(r_norm3)

    # --- Avoiding singularity by introducing a temporarily fake norm
    bSing= r_norm3<1e-8
    r_norm3[bSing]=1 # Temporary hack, replaced at the end

    r_norm5 = (r_norm3 )**(5/3)               # |r|**5
    m_dot_r = DX*m[0] + DY*m[1] + DZ*m[2]     # m . r
    Fact1 = 1./(4*np.pi)*3*m_dot_r/r_norm5
    Fact2 = 1./(4*np.pi)/r_norm3

    ui = np.asarray(Fact1*DX - Fact2*m[0])
    vi = np.asarray(Fact1*DY - Fact2*m[1])
    wi = np.asarray(Fact1*DZ - Fact2*m[2])

    # --- Singularity
    ui[bSing]=0
    vi[bSing]=0
    wi[bSing]=0
    
    return ui,vi,wi

def doublet_u_polar(rcp,zcp,m_z,z0=0): 
    """
    Velocity field induced by one doublet located at z0 on the z axis with vector intensity (0,0,mz)
         ur = 3 m_z /(4pi) r (z-z0) / |r|**5
         uz =   m_z /(4pi) (3 (z-z0)**2 / |r|**5 -1/|r|**3)
    Control points defined by polar coordinates `rcp` and `zcp`.
    """
    if np.any(rcp<0):
        raise Exception('Script meant for positive r')
    rcp=np.asarray(rcp)
    zcp=np.asarray(zcp)

    DZ = zcp-z0

    r_norm3 = (rcp**2 + DZ**2 )**(3/2) # |r|**3
    r_norm3=np.asarray(r_norm3)

    # --- Avoiding singularity by introducing a temporarily fake norm
    bSing= r_norm3<1e-8
    r_norm3[bSing]=1 # Temporary hack, replaced at the end

    r_norm5 = (r_norm3 )**(5/3)               # |r|**5

    ur = np.asarray( 3*rcp*m_z / (4*np.pi) * DZ / r_norm5                      )
    uz = np.asarray(       m_z / (4*np.pi) * (3 * DZ**2 / r_norm5 - 1/r_norm3) )

    # --- Singularity
    ur[bSing]= 0
    uz[bSing]= 0
    
    return ur,uz


# --------------------------------------------------------------------------------}
# --- Line doublet
# --------------------------------------------------------------------------------{
def doublet_line_polar_u_num(rcp,zcp,dmz_dz,z0=0,zmax=1000,nQuad=100):
    """
    Velocity field induced by a doublet line (on the z axis) of intensity `dmz_dz`.
    Control points defined by polar coordinates `rcp` and `zcp`.
    The line goes from `z0` to `zmax`.
    Numerical integration is used with `nQuad` quadrature points.

    """
    rcp=np.asarray(rcp)
    zcp=np.asarray(zcp)

    zq = np.linspace(z0,zmax,nQuad)

    # --- Summation
    dz = zq[1]-zq[0]
    mz = dmz_dz * dz
    ur =  np.zeros(rcp.shape)
    uz =  np.zeros(rcp.shape)
    for z0 in zq: 
        dur, duz = doublet_u_polar(rcp,zcp,mz,z0=z0)
        ur += dur
        uz += duz

    return ur, uz

def doublet_line_polar_u(rcp,zcp,dmz_dz, bSelfInd=False):
    """
    Velocity field induced by a semi-infinite doublet line (on the z axis) of intensity `dmz_dz`
    Control points defined by polar coordinates `rcp` and `zcp`.

     \int  1/(r^2 + (z-x)^2 )^(3/2) dx
     \int  1/(r^2 + (z-x)^2 )^(5/2) dx
    """
    if np.any(rcp<0):
        raise Exception('Script meant for positive r')
    r=np.asarray(rcp)
    z=np.asarray(zcp)

    # Vectorial "if" statements to isolate singular regions of the domain
    bZ0    = np.abs(z)<1e-8
    bR0    = np.abs(r)<1e-8
    bZ0R0  = np.logical_and(bZ0,bR0)
    bZ0Rp  = np.logical_and(bZ0, np.abs(r)>1e-8)
    bR0Zp  = np.logical_and(bR0, z>1e-8)
    bR0Zm  = np.logical_and(bR0, z<-1e-8)
    bOK    = np.logical_and(~bZ0,~bR0)

    uz=np.zeros(r.shape)
    ur=np.zeros(r.shape)

    norm2 = r**2+z**2
    uz[bOK]  = dmz_dz/(4*np.pi) * 1/r[bOK]**2 * ( z[bOK]**3/(norm2[bOK])**(3/2) - z[bOK]/(norm2[bOK])**(1/2) )
    uz[bZ0Rp] = 0
    uz[bR0Zm] = dmz_dz/(4*np.pi) * 1/norm2[bR0Zm]
    #uz[bR0Zp] = dmz_dz/(4*np.pi) * 1/norm2[bR0Zp] #<<< No singularity there, but we force it to 0
    ur[bOK]   =-dmz_dz/(4*np.pi) * r[bOK]          *  1/(norm2[bOK]  )**(3/2)
    ur[bZ0Rp] =-dmz_dz/(4*np.pi) * r[bZ0Rp]        *  1/(norm2[bZ0Rp])**(3/2)
    ur[bR0Zm] = 0
    ur[bR0Zp] = 0

    ur[bZ0R0] = 0
    uz[bZ0R0] = 0


    return ur, uz

def doublet_line_u(Xcp,Ycp,Zcp,dmz_dz,m=0,L=np.inf):
    """ 
    Velocity field induced by a double line, of intensity [0,0,dmz/dz]
    The line is defined by the coordinates (x=mz,y=0,z)
    By default a semi-infinite line is assumed (L=np.inf)
    INPUTS:
       Xcp,Ycp,Zcp: vector or matrix of control points Cartesian Coordinates
       m          : tangent of the skew angle (0 for a straight line along z)
       L          : length of the line (on the z axis!), True length is L*\sqrt{1+m^2}
    """
#     if m==0:
#         vr, vpsi = np.sqrt(Xcp**2+Ycp**2), np.arctan2(Ycp,Xcp) # polar coords
#         ur,uz=doublet_line_polar_u(vr,Zcp,dmz_dz)
#         ux=ur*np.cos(vpsi)
#         uy=ur*np.sin(vpsi)
#     else:
    if True:
        ux =  np.zeros(Xcp.shape)
        uy =  np.zeros(Xcp.shape)
        uz =  np.zeros(Xcp.shape)
        #
        a = np.sqrt(Xcp**2+Ycp**2+Zcp**2)
        b = m*Xcp+Zcp
        c = np.sqrt(m**2+1)

        # On the doublet line, we force the velocity to be zero even if a solution exist
        bOnLine = np.logical_and( np.abs(Ycp)<1e-8, np.abs(Xcp-m*Zcp)<1.e-8)
        bOnLine = np.logical_and(bOnLine, Zcp>=0)
        bOnLine = np.logical_and(bOnLine, Zcp<=L)
        ux[bOnLine]=0
        uy[bOnLine]=0
        uz[bOnLine]=0
        # Singularity condition of the integral:
#         bSing = np.abs(a*c-b)<1e-8
#         ux[bSing]=0
#         uy[bSing]=0
#         uz[bSing]=0
        # Standard formula
#         bOK   = np.logical_not(np.logical_or(bSing,bOnLine))
#         bOK   = np.logical_not(bSing)
        bOK   = np.logical_not(bOnLine)
        a=a[bOK]
        b=b[bOK]
        Xcp=Xcp[bOK]
        Ycp=Ycp[bOK]
        Zcp=Zcp[bOK]

        if L==np.inf:
            i30  = I30 (a,b,c) # NOTE: not scaling
            ti50 = tI50(a,b,c) # = 3 x I50
            ti51 = tI51(a,b,c) # = 3 x I50
            ti52 = tI52(a,b,c) # = 3 x I50
        else:
            i30, ti50, ti51, ti52 = _integrals_DL_finite(a,b,c,L) #= I30, 3*I50, 3*I51, 3*I52

        uz[bOK] = dmz_dz/(4*np.pi)  * (       Zcp**2 * ti50 +   ti52 - 2*Zcp * ti51 - i30 )
        ux[bOK] = dmz_dz/(4*np.pi)  * ( Xcp * Zcp    * ti50 + m*ti52 - (m*Zcp+Xcp) * ti51 )
        uy[bOK] = dmz_dz/(4*np.pi)  * ( Ycp * Zcp    * ti50         -        Ycp   * ti51 )

    return ux,uy,uz 


def doublet_line_u_num(Xcp,Ycp,Zcp,dmz_dz,m=0,zmax=1000,nQuad=100):
    """ 
    Velocity field induced by a double line, of intensity [0,0,dmz/dz]
    The line is defined by the coordinates (x=mz,y=0,z)

    INPUTS:
       Xcp,Ycp,Zcp: vector or matrix of control points Cartesian Coordinates
       m          : tangent of the skew angle (0 for a straight line along z)
    """
    ## --- Summation
    zq = np.linspace(0,zmax,nQuad)
    dz = zq[1]-zq[0]
    mz = dmz_dz * dz
    ux =  np.zeros(Xcp.shape)
    uy =  np.zeros(Xcp.shape)
    uz =  np.zeros(Xcp.shape)
    for z0 in zq: 
        dux, duy, duz = doublet_u(Xcp,Ycp,Zcp,m=[0,0,mz],x0=[m*z0,0,z0]) 
        ux += dux
        uy += duy
        uz += duz
    return ux,uy,uz 

# --------------------------------------------------------------------------------}
# --- Integrals useful for vortex doublet
# --------------------------------------------------------------------------------{
# --- Main Results for semi-infinite line
# a = sqrt(x^2 + y^2 + z^2)
# b = m x + z
# c = sqrt(m^2+1)
#  I30       = 1/(a (a c-b))   
#  I50       = (2 a c - b)/(3 a^3 (a c - b)^2)
#  I51       = 1/ (3 a (ac - b)^2)
#  I52       = 1/(3 c (a c - b)^2)

# --- INTEGRALS WRITTEN ALONG X AND USING PROPER a b c:   a^2 = x^2 + y^2 + z^2, b=my+x , c^2 = (1+m^2)
# II30       = \int  1    / ( a^2  - 2 b x'  + c^2 x'^2 )^(3/2) dx'
# II50       = \int  1    / ( a^2  - 2 b x'  + c^2 x'^2 )^(5/2) dx'
# II51       = \int  x'   / ( a^2  - 2 b x'  + c^2 x'^2 )^(5/2) dx'
# II52       = \int  x'^2 / ( a^2  - 2 b x'  + c^2 x'^2 )^(5/2) dx'

# --- I30
# indef         = (b - c^2 x)/((b^2 - a^2 c^2) sqrt(a^2 + x (-2 b + c^2 x)))
# indef_(b=-ac) = -1/(2 c (a + c x)^2) =  s/(2 c (a -sc x)^2) s==-1
# indef_(b=ac)  =  1/(2 c (a - c x)^2) =  s/(2 c (a -sc x)^2) s==1
#
# --- I50
# Indef  =  ((-b + c^2 x) (-b^2 + 2 c^4 x^2 + c^2 (3 a^2 - 4 b x)))/(3 (-b^2 + a^2 c^2)^2 (a^2 - x (2 b - c^2 x))^(3/2))
# indef_(b=-ac) = -1/(4 c (a + c x)^4) =  s/(4 c (a -sc x)^4) s==-1  b=-ac
# indef_(b= ac) =  1/(4 c (a - c x)^4) =  s/(4 c (a -sc x)^4) s==1


# --- I51
# Indef  =  (-(a^4 c^2) - a^2 b (b - 3 c^2 x) + b x (3 b^2 - 6 b c^2 x + 2 c^4 x^2))/(3 (b^2 - a^2 c^2)^2 (a^2 + x (-2 b + c^2 x))^(3/2)) 
# indef_(b=-ac) =  -(a + 4 c x)/(12 c^2 (a + c x)^4) =  -(a -s 4 c x)/(12 c^2 (a -sc x)^4) s==-1
# indef_(b= ac) =  -(a - 4 c x)/(12 c^2 (a - c x)^4) =  -(a -s 4 c x)/(12 c^2 (a -sc x)^4) s==1

# --- I52
# Indef         = (-2 a^4 b - b^2 x^2 (3 b - c^2 x) + a^2 x (6 b^2 - 3 b c^2 x + c^4 x^2))/(3 (-b^2 + a^2 c^2)^2 (a^2 - x (2 b - c^2 x))^(3/2))
# indef_(b=-ac) = -(a^2 + 4 a c x + 6 c^2 x^2)/(12 c^3 (a + c x)^4)= s(a^2 -s 4 a c x + 6 c^2 x^2)/(12 c^3 (a - s c x)^4) s==-1
# indef_(b= ac) =  (a^2 - 4 a c x + 6 c^2 x^2)/(12 c^3 (a - c x)^4)= s(a^2 -s 4 a c x + 6 c^2 x^2)/(12 c^3 (a - s c x)^4) s==1

def I30(a,b,c):
    # = I30
    # I30 = \int_0^\infty  1    / ( a^2  - 2 b x'  + c^2 x'^2 )^(3/2) dx'
#     if a==b/c:
#         return  -1/(2*c[b0]*a[b0]**2)
    return 1/(a*(a*c-b))

def tI50(a,b,c):
    # return 3 x I50
    # I50 = \int_0^\infty  1    / ( a^2  - 2 b x'  + c^2 x'^2 )^(5/2) dx'
#     if a==b/c:
#         return -3/(4*c*a**4)
    return (2*a*c-b)/(a**3*(a*c-b)**2)

def tI51(a,b,c):
    # return 3 x I51
    # I51 = \int_0^\infty  x'   / ( a^2  - 2 b x'  + c^2 x'^2 )^(5/2) dx'
#     if a==b/c:
#         return 3/(12*c**2*a**3)
    return 1/(a*(a*c-b)**2)

def tI52(a,b,c):
    # return 3 x I52
    # I52 = \int_0^\infty  x'^2 / ( a^2  - 2 b x'  + c^2 x'^2 )^(5/2) dx'
#     if a==b/c:
#         return -3/(12*c**3*a**2)
    return 1/(c*(a*c-b)**2)

def I30_num(x,y,z,m,nQuad=50000,zmax=5000):
    zp = np.linspace(0,zmax,nQuad)
    r2 = (x-m*zp)**2 + y**2 + (z-zp)**2
    dI = 1/ r2**(3/2)
    return np.trapz(dI, zp)

def I50_num(x,y,z,m,nQuad=50000,zmax=5000):
    zp = np.linspace(0,zmax,nQuad)
    r2 = (x-m*zp)**2 + y**2 + (z-zp)**2
    dI = 1/ r2**(5/2)
    return np.trapz(dI, zp)

def I51_num(x,y,z,m,nQuad=50000,zmax=5000):
    zp = np.linspace(0,zmax,nQuad)
    r2 = (x-m*zp)**2 + y**2 + (z-zp)**2
    dI = zp/ r2**(5/2)
    return np.trapz(dI, zp)

def I52_num(x,y,z,m,nQuad=50000,zmax=5000):
    zp = np.linspace(0,zmax,nQuad)
    r2 = (x-m*zp)**2 + y**2 + (z-zp)**2
    dI = zp**2/ r2**(5/2)
    return np.trapz(dI, zp)


# --------------------------------------------------------------------------------}
# --- FINITE LENGTH DOUBLET LINE
# --------------------------------------------------------------------------------{
def _integrals_DL_finite(a,b,c,L):
    """
    Compute integrals values for doublet line of finite length
    returns: I30,  3*I50, 3*I51, 3*I52
    """
    # I30 = \int_0^L  1    / ( a^2  - 2 b x'  + c^2 x'^2 )^(3/2) dx'
    def i30indef_sing(a,b,c,x):
        s=np.sign(b)
        return s/(2*c*(a -s*c*x)**2)
    def i30indef_norm(a,b,c,x):
        return (b - c**2*x)/((b**2 - a**2*c**2)*(a**2 - 2*b*x + c**2 * x**2)**(1/2))
    # I50 = \int_0^L  1    / ( a^2  - 2 b x'  + c^2 x'^2 )^(5/2) dx'
    def i50indef_sing(a,b,c,x):
        s=np.sign(b)
        return 3*s/(4*c*(a -s*c*x)**4)
    def i50indef_norm(a,b,c,x):
        return ( b*(b**2 - 3*a**2*c**2) + 2*c**2*(b**2 + a**2*c**2)*x - 6*b*c**4*x**2 + 2*c**6*x**3 )/ ((b**2 - a**2*c**2)**2*(a**2 - 2*b*x + c**2*x**2)**(3/2) ) 
    # I51 = \int_0^L  x'   / ( a^2  - 2 b x'  + c^2 x'^2 )^(5/2) dx'
    def i51indef_sing(a,b,c,x):
        s=np.sign(b)
        return -3*(a -s*4*c*x)/(12*c**2*(a -s*c*x)**4)
    def i51indef_norm(a,b,c,x):
        return ( -a**2*(b**2 + a**2*c**2) + 3*b*(b**2 + a**2*c**2)*x - 6*b**2*c**2*x**2 + 2*b*c**4*x**3 )/ ( (b**2 - a**2*c**2)**2 *(a**2 -2*b*x + c**2*x**2)**(3/2) ) 
    # I52 = \int  x'^2 / ( a^2  - 2 b x'  + c^2 x'^2 )^(5/2) dx'
    def i52indef_sing(a,b,c,x):
        s=np.sign(b)
        return 3*s*(a**2 -s*4*a*c*x + 6*c**2*x**2)/(12*c**3*(a - s*c*x)**4)
    def i52indef_norm(a,b,c,x):
        return ( -2*a**4*b + 6*a**2*b**2*x - 3*b*(b**2 +  a**2*c**2)*x**2 + c**2*(b**2  + a**2*c**2)*x**3 )/ ( (b**2 - a**2*c**2)**2 *(a**2 - 2*b*x + c**2*x**2)**(3/2) )

    a=np.asarray(a); b=np.asarray(b);
    S = np.abs((b**2 - a**2*c**2))<1e-8 # a== b/c
    N = np.logical_not(S)
    aS=a[S]; bS=b[S];
    aN=a[N]; bN=b[N];

    I30 = np.zeros(a.shape); I50 = np.zeros(a.shape); I51 = np.zeros(a.shape); I52 = np.zeros(a.shape)

    I30[S] = i30indef_sing(aS, bS, c, L)-i30indef_sing(aS, bS, c, 0)
    I30[N] = i30indef_norm(aN, bN, c, L)-i30indef_norm(aN, bN, c, 0)
    I50[S] = i50indef_sing(aS, bS, c, L)-i50indef_sing(aS, bS, c, 0)
    I50[N] = i50indef_norm(aN, bN, c, L)-i50indef_norm(aN, bN, c, 0)
    I51[S] = i51indef_sing(aS, bS, c, L)-i51indef_sing(aS, bS, c, 0)
    I51[N] = i51indef_norm(aN, bN, c, L)-i51indef_norm(aN, bN, c, 0)
    I52[S] = i52indef_sing(aS, bS, c, L)-i52indef_sing(aS, bS, c, 0)
    I52[N] = i52indef_norm(aN, bN, c, L)-i52indef_norm(aN, bN, c, 0)
    return I30, I50, I51, I52

class TestDoublet(unittest.TestCase):

    def test_Integrals(self):
        # --- random point, z>0
        x = 3; y = 1; z = 1.0; m = 0.3
        a = np.sqrt(x**2+y**2+z**2); b = m*x+z ; c = np.sqrt(m**2+1)
        np.testing.assert_almost_equal(I30(a,b,c), I30_num(x,y,z,m), 5)
        np.testing.assert_almost_equal(tI50(a,b,c)/3, I50_num(x,y,z,m), 5)
        np.testing.assert_almost_equal(tI51(a,b,c)/3, I51_num(x,y,z,m), 5)
        np.testing.assert_almost_equal(tI52(a,b,c)/3, I52_num(x,y,z,m), 5)
        # --- random point, z<0
        x = 1; y = 2; z = -0.5; m = 0.3
        a = np.sqrt(x**2+y**2+z**2); b = m*x+z ; c = np.sqrt(m**2+1)
        np.testing.assert_almost_equal(I30(a,b,c), I30_num(x,y,z,m), 5)
        np.testing.assert_almost_equal(tI50(a,b,c)/3, I50_num(x,y,z,m), 5)
        np.testing.assert_almost_equal(tI51(a,b,c)/3, I51_num(x,y,z,m), 5)
        np.testing.assert_almost_equal(tI52(a,b,c)/3, I52_num(x,y,z,m), 5)
        # --- On Axis, z<0 (a==-b/c), no need for singularity
        y = 0; z = -1 ; m = 0.3; x = z*m
        a = np.sqrt(x**2+y**2+z**2); b = m*x+z ; c = np.sqrt(m**2+1)
        np.testing.assert_almost_equal(I30(a,b,c), I30_num(x,y,z,m), 2)
        np.testing.assert_almost_equal(tI50(a,b,c)/3, I50_num(x,y,z,m), 2)
        np.testing.assert_almost_equal(tI51(a,b,c)/3, I51_num(x,y,z,m), 2)
        np.testing.assert_almost_equal(tI52(a,b,c)/3, I52_num(x,y,z,m), 2)
        # --- On Axis, z>0 (a==b/c) Can't match singularity with numerical integration
        y = 0; z = 1 ; m = 0.3; x = z*m
        a = np.sqrt(x**2+y**2+z**2); b = m*x+z ; c = np.sqrt(m**2+1)
        #np.testing.assert_almost_equal(I30(a,b,c), I30_num(x,y,z,m), 2)
        #np.testing.assert_almost_equal(tI50(a,b,c)/3, I50_num(x,y,z,m), 2)
        #np.testing.assert_almost_equal(tI51(a,b,c)/3, I51_num(x,y,z,m), 2)
        #np.testing.assert_almost_equal(tI52(a,b,c)/3, I52_num(x,y,z,m), 2)

    def test_Integrals_Finite(self):
        L = 1e+4
        # --- random point, z>0
        x = 3; y = 1; z = 1.0; m = 0.3
        a = np.sqrt(x**2+y**2+z**2); b = m*x+z ; c = np.sqrt(m**2+1)
        I30L, tI50L, tI51L, tI52L = _integrals_DL_finite(a,b,c,L)
        np.testing.assert_almost_equal(I30(a,b,c) , I30L  , 7)
        np.testing.assert_almost_equal(tI50(a,b,c), tI50L, 7)
        np.testing.assert_almost_equal(tI51(a,b,c), tI51L, 7)
        np.testing.assert_almost_equal(tI52(a,b,c), tI52L, 7)
        # --- random point, z<0
        x = 1; y = 2; z = -0.5; m = 0.3
        a = np.sqrt(x**2+y**2+z**2); b = m*x+z ; c = np.sqrt(m**2+1)
        I30L, tI50L, tI51L, tI52L = _integrals_DL_finite(a,b,c,L)
        np.testing.assert_almost_equal(I30(a,b,c) , I30L  , 7)
        np.testing.assert_almost_equal(tI50(a,b,c), tI50L, 7)
        np.testing.assert_almost_equal(tI51(a,b,c), tI51L, 7)
        np.testing.assert_almost_equal(tI52(a,b,c), tI52L, 7)
        # --- On Axis, z<0 (a==-b/c)
        y = 0; z = -1 ; m = 0.3; x = z*m
        a = np.sqrt(x**2+y**2+z**2); b = m*x+z ; c = np.sqrt(m**2+1)
        I30L, tI50L, tI51L, tI52L = _integrals_DL_finite(a,b,c,L)
        np.testing.assert_almost_equal(I30(a,b,c) , I30L  , 7)
        np.testing.assert_almost_equal(tI50(a,b,c), tI50L, 7)
        np.testing.assert_almost_equal(tI51(a,b,c), tI51L, 7)
        np.testing.assert_almost_equal(tI52(a,b,c), tI52L, 7)
        # --- On Axis, z>0 (a==b/c)
        y = 0; z = 1 ; m = 0.3; x = z*m
        a = np.sqrt(x**2+y**2+z**2); b = m*x+z ; c = np.sqrt(m**2+1)
        a = np.sqrt(x**2+y**2+z**2); b = m*x+z ; c = np.sqrt(m**2+1)
        I30L, tI50L, tI51L, tI52L = _integrals_DL_finite(a,b,c,L)
        #np.testing.assert_almost_equal(I30(a,b,c) , I30L  , 7)
        #np.testing.assert_almost_equal(tI50(a,b,c), tI50L, 7)
        #np.testing.assert_almost_equal(tI51(a,b,c), tI51L, 7)
        #np.testing.assert_almost_equal(tI52(a,b,c), tI52L, 7)


    def test_Doublet_example(self):
        # ---- Singularity check
        u=doublet_u(0,0,0,[0,0,100])
        np.testing.assert_almost_equal(u,(0,0,0))
        # --- Random point check
        u     = doublet_u(1,2,3,[700,-800,900],x0=[-4,-5,6])
        u_ref = (-0.16496, -0.043617, -0.039940)
        np.testing.assert_almost_equal(u,u_ref,5)

    def test_Doublet_polar(self):
        # ---- Singularity check
        u=doublet_u_polar(0,0,0,100)
        np.testing.assert_almost_equal(u,(0,0))
        # ---- Random points check
        u     = doublet_u_polar(1,  3,     100,  z0=6)
        u_ref = doublet_u      (1,0,3,[0,0,100], x0=[0,0,6])
        np.testing.assert_almost_equal(u,(u_ref[0],u_ref[2]))
        u     = doublet_u_polar( 1,  -3,     100,  z0=6)
        u_ref = doublet_u      ( 1,0,-3,[0,0,100], x0=[0,0,6])
        np.testing.assert_almost_equal(u,(u_ref[0],u_ref[2]))

        u     = doublet_u_polar( 0.011,  -0.326,  100    ,  z0=0)
        u_ref = doublet_u      ( 0.011,0,-0.326,[0,0,100], x0=[0,0,0])
        np.testing.assert_almost_equal(u,(u_ref[0],u_ref[2]))


    def test_Doublet_axis(self):
        gamma_t = -10; R = 5
        m       = 2
        dmzdz   = np.pi * R**2 * gamma_t

        #--- Induction zone on axis
        z=np.linspace(-10*R,-1*R,51)
        y=0*z
        x=0*z
        urp,   uzp = doublet_line_polar_u(x, z, dmzdz)
        uri,_ ,uzi = doublet_line_u(x, y, z, dmzdz, m=0, L=np.inf)
        urL,_ ,uzL = doublet_line_u(x, y, z, dmzdz, m=0, L=40*R)
        np.testing.assert_almost_equal(uzp, uzi, 6)
        np.testing.assert_almost_equal(urp, uri, 6)
        np.testing.assert_almost_equal(uzi, uzL, 3)
        np.testing.assert_almost_equal(uri, urL, 3)

        # --- Far field on axis
        z=np.linspace(-100*R,-30*R,51)
        y=0*z
        x=0*z
        uzc = gamma_t/2*(1+z/np.sqrt(z**2 + R**2))
        _,_ ,uzi = doublet_line_u(x, y, z, dmzdz, m=0, L=np.inf)
        _,_ ,uzL = doublet_line_u(x, y, z, dmzdz, m=0, L=100*R)
        np.testing.assert_almost_equal(uzi, uzc, 5)
        np.testing.assert_almost_equal(uzL, uzc, 4)

        #--- Skewed line - on axis - comparison finite infinite 
        # Compare axisformula of a doublet with the one of a cylinder (should match at infinite)
        # TODO axis formula : v_z=-Gamma/(2) *( 1 + z / sqrt(R^2+z^2))
        # gamma_t/2 * (1 + z / np.sqrt(z** 2 + R**2))
        z=np.linspace(-10*R,10*R,31)
        x=m*z
        y=0*z
        uri,_ ,uzi = doublet_line_u(x, y, z, dmzdz, m=m, L=np.inf)
        urL,_ ,uzL = doublet_line_u(x, y, z, dmzdz, m=m, L=100*R)
        np.testing.assert_almost_equal(uzi, uzL, 5)
        np.testing.assert_almost_equal(uri, urL, 4)

        #--- Skewed, on line, enforced to be zero
        z=np.linspace(0*R,10*R,31)
        x=m*z
        y=0*z
        uri,_ ,uzi = doublet_line_u(x, y, z, dmzdz, m=m, L=np.inf)
        urL,_ ,uzL = doublet_line_u(x, y, z, dmzdz, m=m, L=10*R)
        np.testing.assert_almost_equal(uzi, z*0.0, 5)
        np.testing.assert_almost_equal(uri, z*0.0, 5)
        np.testing.assert_almost_equal(uzL, z*0.0, 5)
        np.testing.assert_almost_equal(urL, z*0.0, 5)
#         import matplotlib.pyplot as plt
#         fig,ax = plt.subplots(1,1)
#         ax.plot(z/R, uzL, 'o',label='uzL')
#         ax.plot(z/R, uzi, 'k+',label='uz')
#         ax.plot(z/R, uri,'d', label='uri')
#         ax.plot(z/R, urL,'^', label='urL')
#         ax.set_xlabel('')
#         ax.set_ylabel('')
#         ax.legend()
#         plt.show()

    def test_Doublet_line(self):
        gamma_t = -10
        R       = 10
        dmzdz   = np.pi * R**2 * gamma_t

        # --- Singularity 0,0 (should be infinity)
        u=doublet_line_polar_u  (0,0,dmzdz)
        np.testing.assert_almost_equal(u,(0,0))

        # --- Formula on the z axis from analytical and numerical integration
        zcp = np.linspace(-10,-1,50)*R
        rcp = zcp*0 
        urn,uzn=doublet_line_polar_u_num(rcp,zcp,dmzdz,0,1000,10000)
        urt,uzt=doublet_line_polar_u    (rcp,zcp,dmzdz)

        # --- Enforced 0 on positive z (NOTE: a solution in fact exist
#         u=doublet_line_polar_u  (0,0,dmzdz)
#         np.testing.assert_almost_equal(u,(0,0))


        # ---- Plot
        #uzc = gamma_t/2*(1+zcp/np.sqrt(zcp**2 + R**2))
        #import matplotlib.pyplot as plt
        #fig,ax = plt.subplots(1,1)
        #ax.plot(zcp/R, urt    , label='urt')
        #ax.plot(zcp/R, urn ,'--'   , label='urn')
        #ax.plot(zcp/R, uzt    , label='uzt')
        #ax.plot(zcp/R, uzn ,'--'   , label='uzn')
        #print(zcp)
        #ax.plot(zcp/R, uzc ,'k.'   , label='uzn')
        #ax.set_xlabel('')
        #ax.set_ylabel('')
        #ax.legend()
        #plt.show()

        np.testing.assert_almost_equal(urn,urt,1)
        np.testing.assert_almost_equal(uzn,uzt,1)


if __name__ == "__main__":
#     TestDoublet().test_Doublet_polar()
#     TestDoublet().test_Integrals()
    unittest.main()
