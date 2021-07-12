"""
References:
    [1] E. Branlard, M. Gaunaa - Cylindrical vortex wake model: skewed cylinder, application to yawed or tilted rotors - Wind Energy, 2015
    [2] E. Branlard - Wind Turbine Aerodynamics and Vorticity Based Method, Springer, 2017
"""
#--- Legacy python 2.7
from __future__ import division
from __future__ import print_function
# --- General
import unittest
import numpy as np
import numpy.matlib
# --- Local
try:
    from .elliptic import ellipticPiCarlson, ellipe, ellipk
    from .VortexLine import vl_semiinf_u
except:
    from elliptic import ellipticPiCarlson, ellipe, ellipk
    from VortexLine import vl_semiinf_u

# --------------------------------------------------------------------------------}
# --- Helper function 
# --------------------------------------------------------------------------------{
def skew_components(u_x,u_z,m):
    coschi = 1/np.sqrt(1+m**2)
    sinchi = m/np.sqrt(1+m**2)
    u_zeta = u_z * coschi + u_x * sinchi
    u_xi = - u_z * sinchi + u_x * coschi
    return u_zeta,u_xi
def polar_components(u_x,u_y,vpsi):
    u_r  =  np.multiply(u_x,np.cos(vpsi)) + np.multiply(u_y,np.sin(vpsi))
    u_psi= -np.multiply(u_x,np.sin(vpsi)) + np.multiply(u_y,np.cos(vpsi))
    return u_r,u_psi

# --------------------------------------------------------------------------------}
# --- Core functions, polar coordinates inputs
# --------------------------------------------------------------------------------{
def svc_tang_u_polar(vr,vpsi,vz,gamma_t=-1,R=1,m=0,ntheta=180,polar_out=False):
    """ Induced velocity from a skewed semi infinite cylinder of tangential vorticity.
    Takes polar coordinates as inputs, returns velocity either in Cartesian (default) or polar.
    The cylinder axis is defined by x=m.z, m=tan(chi). The rotor is in the plane z=0.
    The algorithm loops over the control points and performs the integration over theta
    INPUTS:
       vr,vpsi,vz : flat list of control points in polar coordinates
       gamma_t    : tangential vorticity of the vortex sheet (circulation per unit of length oriented along psi). (for WT rotating positively along psi , gamma psi is negative)
       R          : radius of cylinder
       m =tan(chi): tangent of wake skew angle
       ntheta    : number of points used for integration
    Reference: [1,2]"""
    EPSILON_AXIS=1e-7; # relative threshold for using axis formula
    vtheta = np.pi/2 + np.linspace(0, 2*np.pi, ntheta)
    # Flattening
    shape_in=vr.shape
    vr   = np.asarray(vr).ravel()
    vpsi = np.asarray(vpsi).ravel()
    vz   = np.asarray(vz).ravel()
    # Constants of theta
    c   = 1 + m**2
    bz  = R * m * np.cos(vtheta)
    u_z = np.zeros(vr.shape)
    if polar_out:
        u_r   = np.zeros(vr.shape)
        u_psi = np.zeros(vr.shape)
        # ---- Loop on all control points to find velocity
        for i,(r,psi,z) in enumerate(zip(vr,vpsi,vz)):
            # Functions of theta in the integrand
            a = R**2 + r** 2 + z**2 - 2*R*r*np.cos(vtheta - psi)
            b = 2 * m * R * np.cos(vtheta) - 2 * m * r * np.cos(psi) - 2 * z
            ap, bp = R * z * np.sin(vtheta - psi), -R * np.sin(vtheta - psi)
            ar, br = R * z * np.cos(vtheta - psi), -R * np.cos(vtheta - psi)
            az = R * (R - r * np.cos(vtheta - psi))
            D = 2*gamma_t/(4*np.pi)/(np.multiply(np.sqrt(a),(2 * np.sqrt(a * c)+ b)))
            # Integrations
            u_r[i]    = np.trapz((ar * np.sqrt(c)+ np.multiply(br,np.sqrt(a)))*D, vtheta)
            u_psi[i]  = np.trapz((ap * np.sqrt(c)+ np.multiply(bp,np.sqrt(a)))*D, vtheta)
            u_z[i]    = np.trapz((az * np.sqrt(c)+ np.multiply(bz,np.sqrt(a)))*D, vtheta)
        # Reshaping to desired shape
        u_r   =  u_r.reshape(shape_in)   
        u_psi =  u_psi.reshape(shape_in)   
        u_z   =  u_z.reshape(shape_in)   
        return u_r,u_psi,u_z
    else:
#         print('>>>>>> HACK')
#         bx, by = -R * np.cos(vtheta), -R * np.sin(vtheta)
#         u_x   = np.zeros(vr.shape)
#         u_y   = np.zeros(vr.shape)
#         # ---- Loop on all control points to find velocity
#         for i,(r,psi,z) in enumerate(zip(vr,vpsi,vz)):
#             # Functions of theta in the integrand
#             a = R**2 + r** 2 + z**2 - 2*R*r*np.cos(vtheta - psi)
#             b = 2 * m * R * np.cos(vtheta) - 2 * m * r * np.cos(psi) - 2 * z
#             ax, ay = R * z * np.cos(vtheta),  R * z * np.sin(vtheta)
#             az = R * (R - r * np.cos(vtheta - psi))
#             #D = 2*gamma_t/(4*np.pi)/(np.multiply(np.sqrt(a),(2 * np.sqrt(a * c)+ b)))
#             D = -4*gamma_t/(np.sqrt(c)*4*np.pi)/(4*a*c-b**2)
#             # Integrations
#             u_x[i]    = np.trapz((b*bx-2*ax*c)*D, vtheta)
#             u_y[i]    = np.trapz((b*by-2*ay*c)*D, vtheta)
#             u_z[i]    = np.trapz((b*bz-2*az*c)*D, vtheta)


        bx, by = -R * np.cos(vtheta), -R * np.sin(vtheta)
        u_x   = np.zeros(vr.shape)
        u_y   = np.zeros(vr.shape)
        # ---- Loop on all control points to find velocity
        for i,(r,psi,z) in enumerate(zip(vr,vpsi,vz)):
            # Functions of theta in the integrand
            a = R**2 + r** 2 + z**2 - 2*R*r*np.cos(vtheta - psi)
            b = 2 * m * R * np.cos(vtheta) - 2 * m * r * np.cos(psi) - 2 * z
            ax, ay = R * z * np.cos(vtheta),  R * z * np.sin(vtheta)
            az = R * (R - r * np.cos(vtheta - psi))
            D = 2*gamma_t/(4*np.pi)/(np.multiply(np.sqrt(a),(2 * np.sqrt(a * c)+ b)))
            # Integrations
            u_x[i]    = np.trapz((ax * np.sqrt(c)+ np.multiply(bx,np.sqrt(a)))*D, vtheta)
            u_y[i]    = np.trapz((ay * np.sqrt(c)+ np.multiply(by,np.sqrt(a)))*D, vtheta)
            u_z[i]    = np.trapz((az * np.sqrt(c)+ np.multiply(bz,np.sqrt(a)))*D, vtheta)

        # Reshaping to desired shape
        u_x   =  u_x.reshape(shape_in)   
        u_y   =  u_y.reshape(shape_in)   
        u_z   =  u_z.reshape(shape_in)   
        return u_x,u_y,u_z


def svc_longi_u_polar(vr,vpsi,vz,gamma_l=-1,R=1,m=0,ntheta=180,polar_out=False):
    """ Raw function, not intended to be exported. 
    Induced velocity from a skewed semi infinite cylinder of longitudinal vorticity.
    Takes polar coordinates as inputs, returns velocity either in Cartesian (default) or polar.
    The cylinder axis is defined by x=m.z, m=tan(chi). The rotor is in the plane z=0.
    INPUTS:
       vr,vpsi,vz : control points in polar coordinates, may be of any shape
       gamma_t    : tangential vorticity of the vortex sheet (circulation per unit of length oriented along psi). (for WT rotating positively along psi , gamma psi is negative)
       R          : radius of cylinder
       m =tan(chi): tangent of wake skew angle
       ntheta    : number of points used for integration
    Reference: [1,2]"""
    EPSILON_AXIS=1e-7; # relative threshold for using axis formula
    vtheta = np.linspace(0,2 * np.pi,ntheta) + np.pi / ntheta
    # Flattening, and dimensionless!
    shape_in=vr.shape
    vr   = np.asarray(vr/R).ravel()
    vpsi = np.asarray(vpsi).ravel()
    vz   = np.asarray(vz/R).ravel()
    u_z = np.zeros(vr.shape)
    if polar_out:
        u_r = np.zeros(vr.shape)
        u_psi = np.zeros(vr.shape)
        for i,(r,psi,z) in enumerate(zip(vr,vpsi,vz)):
            Den1 = np.sqrt(1 + r**2 + z**2 - 2*r* np.cos(vtheta - psi))
            Den2 = - z + m * np.cos(vtheta) + np.sqrt(1 + m ** 2) * np.sqrt(1 + r ** 2 + z ** 2 - 2 * r * np.cos(vtheta - psi)) - m * r * np.cos(psi)
            DenInv = gamma_l/(4*np.pi)/np.multiply(Den1,Den2)
            u_r[i]   = np.trapz((  - m*z*np.sin(psi) + np.sin(vtheta-psi))*DenInv,vtheta)
            u_psi[i] = np.trapz((r - m*z*np.cos(psi) - np.cos(vtheta-psi))*DenInv,vtheta)
            u_z[i]   = np.trapz(m * (-np.sin(vtheta) + r*np.sin(psi))     *DenInv,vtheta)
        # Reshaping to input shape
        u_r   =  u_psi.reshape(shape_in) 
        u_psi =  u_psi.reshape(shape_in) 
        u_z   =  u_z.reshape(shape_in)   
        return (u_r,u_psi,u_z)
    else:
        u_x = np.zeros(vr.shape)
        u_y = np.zeros(vr.shape)
        for i,(r,psi,z) in enumerate(zip(vr,vpsi,vz)):
            Den1 = np.sqrt(1 + r**2 + z**2 - 2*r* np.cos(vtheta - psi))
            Den2 = - z + m * np.cos(vtheta) + np.sqrt(1 + m ** 2) * np.sqrt(1 + r ** 2 + z ** 2 - 2 * r * np.cos(vtheta - psi)) - m * r * np.cos(psi)
            DenInv = gamma_l/(4*np.pi)/np.multiply(Den1,Den2)
            u_x[i]   = np.trapz(        (np.sin(vtheta) - r*np.sin(psi))  *DenInv,vtheta)
            u_y[i]   = np.trapz((- m*z - np.cos(vtheta) + r*np.cos(psi))  *DenInv,vtheta)
            u_z[i]   = np.trapz(m * (-np.sin(vtheta) + r*np.sin(psi))     *DenInv,vtheta)
        # Reshaping to input shape
        u_x   =  u_x.reshape(shape_in)   
        u_y   =  u_y.reshape(shape_in)   
        u_z   =  u_z.reshape(shape_in)   
        return (u_x,u_y,u_z)

def svc_root_u_polar(vr,vpsi,vz,Gamma_r=-1,m=0,polar_out=False):
    """
    Induced velocity from a skewed root vortex
    Takes polar coordinates as inputs, returns velocity either in Cartesian (default) or polar.
    The cylinder axis is defined by x=m.z, m=tan(chi). The rotor is in the plane z=0.
    INPUTS:
       vr,vpsi,vz : control points in polar coordinates, may be of any shape
       Gamma_r    : Root vortex circulation, negative for a wind turbine
       m =tan(chi): tangent of wake skew angle
    Reference: [1,2]"""
    EPSILON_AXIS=1e-7; # relative threshold for using axis formula
    chi = np.arctan(m)
    if Gamma_r==0:
        return vr*0,vr*0,vr*0
    # Flattening
    shape_in=vr.shape
    vr   = np.asarray(vr).ravel()
    vpsi = np.asarray(vpsi).ravel()
    vz   = np.asarray(vz).ravel()
    # --- Computes ux, uy, uz, upsi
    u_z = np.zeros(vr.shape)
    if (m == 0):
        u_psi = np.multiply(Gamma_r/(4*np.pi*vr), (1+vz/np.sqrt(vr** 2 + vz**2)))
        u_x   = -np.sin(vpsi)*u_psi
        u_y   =  np.cos(vpsi)*u_psi
    else:
        if (np.max(np.abs(vz)) > 0):
            # need to use general formula
            u_x = np.zeros(vr.shape)
            u_y = np.zeros(vr.shape)
            e = np.array([np.sin(chi),0,np.cos(chi)])
            for i,(r,psi,z) in enumerate(zip(vr,vpsi,vz)):
                u_x[i],u_y[i],u_z[i]= vl_semiinf_u(r*np.cos(psi),r*np.sin(psi),z,e[0],e[1],e[2],Gamma_r,visc_model=0,t=0)
            u_psi = - u_x*np.sin(vpsi) + u_y* np.cos(vpsi)
        else:
            # rotor plane analytical (see Yaw article)
            u_psi = np.zeros(vr.shape)
            coschi = 1 / np.sqrt(1 + m ** 2)
            sinchi = m / np.sqrt(1 + m ** 2)
            Iz = vr > (EPSILON_AXIS)
            bnIz = np.logical_not(Iz)
            u_z  [Iz] = np.multiply(Gamma_r/(4*np.pi*vr[Iz]), 1.0/(1-np.cos(vpsi[Iz])*sinchi)*sinchi*np.sin(vpsi[Iz]))
            u_psi[Iz] = np.multiply(Gamma_r/(4*np.pi*vr[Iz]), 1.0/(1-np.cos(vpsi[Iz])*sinchi)*coschi)
            u_z  [bnIz] =0
            u_psi[bnIz] =0
            u_x = -np.sin(vpsi)*u_psi
            u_y =  np.cos(vpsi)*u_psi
    # Reshaping to input shape
    u_z   =  u_z.reshape(shape_in) 
    if polar_out:
        u_r   = u_x * np.cos(vpsi) + u_y * np.sin(vpsi)
        u_psi = u_psi.reshape(shape_in) 
        return (u_r,u_psi,u_z)
    else:
        u_x   =  u_x.reshape(shape_in)   
        u_y   =  u_y.reshape(shape_in)   
    return (u_x,u_y,u_z)

def svc_u_polar(vr,vpsi,vz,gamma_t,gamma_l,Gamma_r,R=1,m=0,ntheta=180,polar_out=False):
    """ Induced velocities from a skewed semi infinite cylinder with:
       - tangential vorticity gamma_t
       - longitudinal vorticity gamma_l
       - a root vortex, Gamma_r
    """
    u1 ,u2 ,u3  = svc_longi_u_polar(vr,vpsi,vz,gamma_l,R=R,m=m,ntheta=ntheta,polar_out=False)
    u1t,u2t,u3t = svc_tang_u_polar (vr,vpsi,vz,gamma_t,R=R,m=m,ntheta=ntheta,polar_out=False)
    u1 += u1t
    u2 += u2t
    u3 += u3t
    u1t,u2t,u3t = svc_root_u_polar (vr,vpsi,vz,Gamma_r    ,m=m,           polar_out=False)
    u1 += u1t
    u2 += u2t
    u3 += u3t
    return u1,u2,u3

# --------------------------------------------------------------------------------}
# --- Main functions with Cartesian inputs
# --------------------------------------------------------------------------------{
def svc_longi_u(Xcp,Ycp,Zcp,gamma_l=-1,R=1,m=0,ntheta=180,polar_out=False):
    """ Induced velocity from a skewed semi infinite cylinder of longitudinal vorticity.
    The cylinder axis is defined by x=m.z, m=tan(chi). The rotor is in the plane z=0.
    INPUTS:
       Xcp,Ycp,Zcp: vector or matrix of control points Cartesian Coordinates
       gamma_l    : longitudinal vorticity of the vortex sheet (circulation per unit of length oriented along zeta), negative for a WT
       R          : radius of cylinder
       m =tan(chi): tangent of wake skew angle
       ntheta     : number of points used for integration
    Reference: [1,2]"""
    vr, vpsi = np.sqrt(Xcp**2+Ycp**2), np.arctan2(Ycp,Xcp) # polar coords
    u1,u2,u3=svc_longi_u_polar(vr,vpsi,Zcp,gamma_l,R,m,ntheta,polar_out=polar_out)
    return u1,u2,u3 # ux,uy,uz OR ur,upsi,uz

def svc_tang_u(Xcp,Ycp,Zcp,gamma_t=-1,R=1,m=0,ntheta=180,polar_out=False):
    """ Induced velocity from a skewed semi infinite cylinder of tangential vorticity.
    The cylinder axis is defined by x=m.z, m=tan(chi). The rotor is in the plane z=0.
    INPUTS:
       Xcp,Ycp,Zcp: vector or matrix of control points Cartesian Coordinates
       gamma_t    : tangential vorticity of the vortex sheet (circulation per unit of length oriented along psi), negative for a WT
       R          : radius of cylinder
       m =tan(chi): tangent of wake skew angle
       ntheta     : number of points used for integration
    Reference: [1,2]"""
    vr, vpsi = np.sqrt(Xcp**2+Ycp**2), np.arctan2(Ycp,Xcp) # polar coords
    u1,u2,u3 = svc_tang_u_polar(vr,vpsi,Zcp,gamma_t,R,m,ntheta,polar_out=polar_out)
    return u1,u2,u3 # ux,uy,uz OR ur,upsi,uz

def svc_root_u(Xcp,Ycp,Zcp,Gamma_r=-1,m=0,polar_out=False):
    """ Induced velocity from a skewed root vortex.
    The root vortex axis is defined by x=m.z, m=tan(chi). The rotor is in the plane z=0.
    INPUTS:
       Xcp,Ycp,Zcp: vector or matrix of control points Cartesian Coordinates
       Gamma_r    : Root vortex circulation, negative for a wind turbine
       m =tan(chi): tangent of wake skew angle
       ntheta     : number of points used for integration
    Reference: [1,2]"""
    vr, vpsi = np.sqrt(Xcp**2+Ycp**2), np.arctan2(Ycp,Xcp) # polar coords
    u1,u2,u3 = svc_root_u_polar(vr,vpsi,Zcp,Gamma_r,m,polar_out=polar_out)
    return u1,u2,u3 # ux,uy,uz OR ur,upsi,uz

def svcs_tang_u(Xcp,Ycp,Zcp,gamma_t,R,m,Xcyl,Ycyl,Zcyl,ntheta=180, Ground=False):
    """ 
    Computes the velocity field for nCyl*nr cylinders, extending along z:
        nCyl: number of main cylinders
        nr  : number of concentric cylinders within a main cylinder 
    INPUTS: 
        Xcp,Ycp,Zcp: cartesian coordinates of control points where the velocity field is not be computed
        gamma_t: array of size (nCyl,nr), distribution of gamma for each cylinder as function of radius
        R      : array of size (nCyl,nr), 
        m      : array of size (nCyl,nr), 
        Xcyl,Ycyl,Zcyl: array of size nCyl) giving the center of the rotor
        Ground: boolean, True if ground effect is to be accounted for
    All inputs (except Ground) should be numpy arrays
    """ 
    Xcp=np.asarray(Xcp)
    Ycp=np.asarray(Ycp)
    Zcp=np.asarray(Zcp)
    ux = np.zeros(Xcp.shape)
    uy = np.zeros(Xcp.shape)
    uz = np.zeros(Xcp.shape)
    nCyl,nr = R.shape
    print('Tang.  (skewed)   ',end='')
    for i in np.arange(nCyl):
        Xcp0,Ycp0,Zcp0=Xcp-Xcyl[i],Ycp-Ycyl[i],Zcp-Zcyl[i]
        if Ground:
            YcpMirror = Ycp0+2*Ycyl[i]
            Ylist = [Ycp0,YcpMirror]
        else:
            Ylist = [Ycp0]
        for iy,Y in enumerate(Ylist):
            for j in np.arange(nr):
                if iy==0:
                    print('.',end='')
                else:
                    print('m',end='')
                if np.abs(gamma_t[i,j]) > 0:
                    ux1,uy1,uz1 = svc_tang_u(Xcp0,Y,Zcp0,gamma_t[i,j],R[i,j],m[i,j],ntheta=ntheta,polar_out=False)
                    ux = ux + ux1
                    uy = uy + uy1
                    uz = uz + uz1
    print('')
    return ux,uy,uz

def svcs_longi_u(Xcp,Ycp,Zcp,gamma_l,R,m,Xcyl,Ycyl,Zcyl,ntheta=180,Ground=False):
    """ See svcs_tang_u """ 
    Xcp=np.asarray(Xcp)
    Ycp=np.asarray(Ycp)
    Zcp=np.asarray(Zcp)
    ux = np.zeros(Xcp.shape)
    uy = np.zeros(Xcp.shape)
    uz = np.zeros(Xcp.shape)
    nCyl,nr = R.shape
    print('Longi. (skewed)   ',end='')
    for i in np.arange(nCyl):
        Xcp0,Ycp0,Zcp0=Xcp-Xcyl[i],Ycp-Ycyl[i],Zcp-Zcyl[i]
        if Ground:
            YcpMirror = Ycp0+2*Ycyl[i]
            Ylist = [Ycp0,YcpMirror]
        else:
            Ylist = [Ycp0]
        for iy,Y in enumerate(Ylist):
            for j in np.arange(nr):
                if iy==0:
                    print('.',end='')
                else:
                    print('m',end='')
                if np.abs(gamma_l[i,j]) > 0:
                    ux1,uy1,uz1 = svc_longi_u(Xcp0,Ycp0,Zcp0,gamma_l[i,j],R[i,j],m[i,j],ntheta=ntheta,polar_out=False)
                    ux = ux + ux1
                    uy = uy + uy1
                    uz = uz + uz1
    print('')
    return ux,uy,uz



# --------------------------------------------------------------------------------}
# --- Rewrite of Matlab functions, legacy
# --------------------------------------------------------------------------------{
def fV_Trailed(vr,vpsi,vz,m,gamma_longi,ntheta,nout=7):
    """ See Yaw article for notations and coordinate system
    Return induced velocity by an infinite number of trailed vortices (semi-infinite lines whose starting points lay on the rotor circle)
    """
    u_x,u_y,u_z = svc_longi_u_polar(vr,vpsi,vz,gamma_longi,R=1,m=m,ntheta=ntheta,polar_out=False)
    if nout==1:
        return u_z
    u_zeta,u_xi = skew_components(u_x,u_z,m)
    u_r,u_psi   = polar_components(u_x,u_y,vpsi)
    outputs=(u_z,u_psi,u_x,u_y,u_zeta,u_xi,u_r)
    return outputs[:nout]

def fV_Tangential(vr,vpsi,vz,m,gamma_t,ntheta,nout=7): 
    """ This function is purely for backward compatibility with Matlab scripts"""
    u_x,u_y,u_z = svc_tang_u_polar(vr,vpsi,vz,gamma_t,R=1,m=m,ntheta=ntheta,polar_out=False)
    if nout==1:
        return u_z
    u_zeta,u_xi=skew_components(u_x,u_z,m)
    u_r,u_psi =polar_components(u_x,u_y,vpsi)
    outputs=(u_z,u_psi,u_x,u_y,u_zeta,u_xi,u_r)
    return outputs[:nout]

def fV_Root(vr,vpsi,vz, m =0, Gamma_r=-1,nout=1):
    """ Return induced velocity by the root vortex
    Coordinate system is true polar coordinates, with convention of Yaw article
    """
    u_x,u_y,u_z= svc_root_u_polar(vr,vpsi,vz,Gamma_r=Gamma_r,m=m,polar_out=False)
    if nout==1:
        return u_z
    u_zeta,u_xi = skew_components(u_x,u_z,m)
    u_r,u_psi   = polar_components(u_x,u_y,vpsi)
    outputs=(u_z,u_psi,u_x,u_y,u_zeta,u_xi,u_r)
    return outputs[:nout]

# --------------------------------------------------------------------------------}
# --- Rotor plane flow expansions 
# --------------------------------------------------------------------------------{
def fKxit(vr,m):
    """ Returns Kxit according to yaw article . vr is in [0;1], m=tan(chi)"""
    EPSILON_AXIS=1e-7; # relative threshold for using axis formula
    fOye = 0.5 * (vr + 0.4 * vr ** 3 + 0.4 * vr ** 5)
    Kxit_num = np.zeros((1,len(vr)))
    k2 = ((1 - vr) ** 2) / ((1 + vr) ** 2)
    m1 = (np.sqrt(1 + m ** 2) + np.sqrt(vr ** 2 + m ** 2)) / (1 + vr)
    m2 = (np.sqrt(1 + m ** 2) - np.sqrt(vr ** 2 + m ** 2)) / (1 + vr)
    b1 = m1 ** 2 - 1
    b2 = 1 - m2 ** 2
    j2 = 1 - k2
    kr2 = ellipk(vr ** 2)
    Pi1 = ellipticPiCarlson(- b1,j2)
    Pi2 = ellipticPiCarlson(b2,j2)
    Kxit=np.zeros(vr.shape)
    if (m == 0):
        k2 = 4 * vr / ((vr + 1) ** 2)
        k = np.sqrt(k2)
        K = ellipk(k2)
        E = ellipe(k2)
        b1 = (vr) > (EPSILON_AXIS) 
        b0 = np.logical_not(b1)
        Kxit[b1] = np.multiply(1/(np.pi)*np.sqrt(1.0/vr[b1]),(np.multiply((2 - k2[b1]) / k[b1],K[b1]) - np.multiply(2.0/k[b1],E[b1])))
        Kxit[b0] = 0
    else:
        b1 = (vr) > (EPSILON_AXIS) 
        b0 = np.logical_not(b1)
        Kxit[b1] = np.multiply(2*(1+m**2)*vr[b1]/(m**2*np.pi),kr2[b1]) - np.multiply(np.multiply(vr[b1],(vr[b1] + 1))*np.sqrt(m ** 2 + 1)/(2*m**2*np.pi*np.sqrt(m**2+vr[b1]**2)),(np.multiply((b1[b1]+j2[b1]),Pi1[b1]) + np.multiply((b2[b1]-j2[b1]),Pi2[b1])))
        Kxit[b0] = 0
    # See yaw article
    chi = np.arctan(m)
    vtheta = np.linspace(0,np.pi / 2,1000)
    Kxit_num=np.zeros(vr.shape)
    for ir in np.arange(len(vr)):
        r = vr[ir]
        Kxit_num[ir] = 2 * r / np.pi * np.trapz(np.sin(2 * vtheta) ** 2.0 / (np.multiply(np.sqrt((1 + r) ** 2 - 4 * r * np.cos(vtheta) ** 2),((r - np.cos(2 * vtheta)) ** 2 * np.cos(chi) ** 2 + np.sin(2 * vtheta) ** 2))), vtheta)
    
    return Kxit,Kxit_num,fOye

def fKzt(r,m,nout=2):
    """ Returns Kzt according to yaw article """
    fOye = 0.5 * (r + 0.4 * r ** 3 + 0.4 * r ** 5)
    vr = r
    Kzt = np.zeros(vr.shape)
    Kztnum = np.zeros(vr.shape)
    if m == 0:
        raise Exception('Not intended for m==0')
    k2 = ((1 - r) ** 2) / ((1 + r) ** 2)
    m1 = (np.sqrt(1 + m ** 2) + np.sqrt(r ** 2 + m ** 2)) / (1 + r)
    m2 = (np.sqrt(1 + m ** 2) - np.sqrt(r ** 2 + m ** 2)) / (1 + r)
    b1 = m1 ** 2 - 1
    b2 = 1 - m2 ** 2
    j2 = 1 - k2
    kr2 = ellipk(r ** 2)
    Pi1 = ellipticPiCarlson(- b1,j2)
    Pi2 = ellipticPiCarlson(b2,j2)
    Kzt = np.multiply(2 * np.sqrt(1 + m ** 2) * r / (m * np.pi),kr2) - np.multiply(np.multiply(r,(r + 1)) / (2 * m * np.pi * np.sqrt(m ** 2 + r ** 2)),(np.multiply((b1 + j2),Pi1) + np.multiply((b2 - j2),Pi2)))
    # Coleman formula B.5 term 3 and 4  !!!!! Note the minus sign added
    vtheta = np.linspace(0,np.pi,1000)
    for ir,r in enumerate(vr):
        Kztnum[ir] = - 1 / (np.pi) * r * np.sqrt(1 + m ** 2) / m * np.trapz(- 1.0 / (np.sqrt(1 + r ** 2 - 2 * r * np.cos(vtheta))) + np.sqrt(1 - 2 * r * np.cos(vtheta) + r ** 2) / (1 + r ** 2 - 2 * r * np.cos(vtheta) + m ** 2 * np.sin(vtheta) ** 2),vtheta)

    if nout==1:
        return Kzt
    elif nout <=3:
        outputs=(Kzt,Kztnum,fOye)
    elif nout > 3:
        Kztnum2 = np.zeros(vr.shape)
        Kztnum3 = np.zeros(vr.shape)
        Kztnum4 = np.zeros(vr.shape)
        # My formula Alternative form1
        vtheta = np.linspace(0,np.pi,1000)
        for ir,r in enumerate(vr):
            Kztnum2[ir] = r * m * np.sqrt(1 + m ** 2) / np.pi * np.trapz(np.sin(vtheta) ** 2.0 / (np.multiply(np.sqrt(1 + r ** 2 - 2 * r * np.cos(vtheta)),(1 + r ** 2 - 2 * r * np.cos(vtheta) + m ** 2 * np.sin(vtheta) ** 2))),vtheta)
        # My formula Alternative form3 (WEIRD RESULTS !!!!!!!!!!!!)
        vtheta = np.linspace(0,np.pi / 2,1000)
        for ir,r in enumerate(vr):
            Kztnum3[ir] = 2 * r * np.sqrt(1 + m ** 2) * m / np.pi * np.trapz(np.sin(2 * vtheta) ** 2.0 / (np.multiply(np.sqrt((1 + r) ** 2 - 4 * r * np.cos(vtheta) ** 2),((1 + r) ** 2 - 4 * r * np.cos(vtheta) ** 2 + m ** 2 * np.sin(2 * vtheta) ** 2))),vtheta)
        # My formula Alternative form4
        vtheta = np.linspace(0,np.pi / 2,1000)
        for ir,r in enumerate(vr):
            Kztnum4[ir] = 2 * r / np.pi * 1 * np.sin(chi) * np.trapz(np.sin(2 * vtheta) ** 2.0 / (np.multiply(np.sqrt((1 + r) ** 2 - 4 * r * np.cos(vtheta) ** 2),((r - np.cos(2 * vtheta)) ** 2 * np.cos(chi) ** 2 + np.sin(2 * vtheta) ** 2))),vtheta)
        outputs=(Kzt,Kztnum,fOye,Kztnum2,Kztnum3,Kztnum4)
    return outputs[:nout]
    
    




# --------------------------------------------------------------------------------}
# --- TEST 
# --------------------------------------------------------------------------------{
class TestSkewedCylinder(unittest.TestCase):
    def test_SVC_rotor(self):
        """ """
        """ See paragraph "Properties on the rotor disk" of [1] """
        # data
        gamma_t,R,chi = -5, 10, 30*np.pi/180
        m=np.tan(chi) # tan(chi) 
        eps=10**-1 *R
        # --- At rotor center (see also next test, stronger)
        u_x,u_y,u_z = svc_tang_u(0,0,0,gamma_t,R,m)
        u_zeta,u_xi=skew_components(u_x,u_z,m)
        uz0=gamma_t/2
        np.testing.assert_almost_equal(u_x        ,np.tan(chi/2)*uz0      ,decimal=7)
        np.testing.assert_almost_equal(u_z        ,uz0 ,decimal=7)
        np.testing.assert_almost_equal(u_zeta     ,uz0 ,decimal=7)

        # --- At psi=pi/2 (i.e. x=0), z=0 (Eq 9 from [1]), ux,uz,uzeta,uxi constant!
        y=np.linspace(0,R-eps,4)
        x=y*0
        z=y*0
        u_x,u_y,u_z=svc_tang_u(x,y,z,gamma_t,R,m)
        u_zeta,u_xi=skew_components(u_x,u_z,m)
        uz0=np.asarray([gamma_t/2]*len(x))
        np.testing.assert_almost_equal(u_zeta     ,uz0                    ,decimal=7)
        np.testing.assert_almost_equal(u_z        ,uz0                    ,decimal=7)
        np.testing.assert_almost_equal(u_xi/u_zeta,[-np.tan(chi/2)]*len(x),decimal=7)
        np.testing.assert_almost_equal(u_x /u_z   ,[ np.tan(chi/2)]*len(x),decimal=7)
        np.testing.assert_almost_equal(u_x        ,uz0*np.tan(chi/2)      ,decimal=7)

        # --- Component zeta over the entire plane is g_t/2 (Eq 10 from [1])
        vR,vPsi = np.meshgrid(np.linspace(0,R-eps,5), np.linspace(0,2*np.pi,12))
        x=vR*np.cos(vPsi)
        y=vR*np.sin(vPsi)
        z=x*0
        u_x,u_y,u_z=svc_tang_u(x,y,z,gamma_t,R,m)
        u_zeta,u_xi=skew_components(u_x,u_z,m)
        uz0=gamma_t/2
        np.testing.assert_almost_equal(u_zeta , uz0 ,decimal=5)

        # --- Plane y=0 (anti-)symmetry - x,z,zeta,xi: symmetric - y: anti-symmetric
        x,y = np.meshgrid(np.linspace(-R/3,R/3,3), [-R/2,R/2] )
        z=x*0
        u_x,u_y,u_z=svc_tang_u(x,y,z,gamma_t,R,m)
        u_zeta,u_xi=skew_components(u_x,u_z,m)
        np.testing.assert_almost_equal(u_x   [0,:], u_x   [1,:])
        np.testing.assert_almost_equal(u_z   [0,:], u_z   [1,:])
        np.testing.assert_almost_equal(u_zeta[0,:], u_zeta[1,:])
        np.testing.assert_almost_equal(u_xi  [0,:], u_xi  [1,:])
        np.testing.assert_almost_equal(u_y   [0,:],-u_y   [1,:]) # anti-symmetric

        # --- Radial anti-symmetry of components x,y,z about their origin value
        r0   = R/2 # cannot do negative r here since definition of r is positive
        psi0 = np.linspace(0,np.pi,10)
        vPsi = np.array([psi0 , psi0+np.pi] ).T
        x=r0*np.cos(vPsi)
        y=r0*np.sin(vPsi)
        z=x*0
        u_x,u_y,u_z     =svc_tang_u(x,y,z,gamma_t,R,m)
        u_zeta,u_xi=skew_components(u_x,u_z,m)
        u_x0,u_y0,u_z0=svc_tang_u(0,0,0,gamma_t,R,m)
        u_zeta0,u_xi0=skew_components(u_x0,u_z0,m)
        np.testing.assert_almost_equal(u_x   [:,0]+u_x   [:,1], 2*u_x0   )
        np.testing.assert_almost_equal(u_y   [:,0]+u_y   [:,1], 2*u_y0   )
        np.testing.assert_almost_equal(u_z   [:,0]+u_z   [:,1], 2*u_z0   )
        np.testing.assert_almost_equal(u_zeta[:,0]+u_zeta[:,1], 2*u_zeta0)
        np.testing.assert_almost_equal(u_xi  [:,0]+u_xi  [:,1], 2*u_xi0  )


    def test_SVC_farwake(self):
        """ """
        """ See paragraph "Properties on the rotor disk" of [1] """
        # data
        gamma_t,R,chi = -5, 10, 30*np.pi/180
        m=np.tan(chi) # tan(chi) 
        eps=10**-1 *R
        z0=1000*R # Far wake
        # --- At rotor center (see also next test, stronger)
        #u_x,u_y,u_z=svc_tang_u(0,0,0,gamma_t,R,m)
        #uz0=gamma_t/2
        #np.testing.assert_almost_equal(u_x        ,np.tan(chi/2)*uz0      ,decimal=7)
        #np.testing.assert_almost_equal(u_z        ,uz0 ,decimal=7)
        #np.testing.assert_almost_equal(u_zeta     ,uz0 ,decimal=7)
        # --- Component zeta over the entire plane is g_t/2 (Eq 10 from [1])
        vR,vTheta = np.meshgrid(np.linspace(0,R-eps,5), np.linspace(0,2*np.pi,12))
        x=vR*np.cos(vTheta)+z0*m
        y=vR*np.sin(vTheta)
        z=x*0 + z0
        u_x,u_y,u_z=svc_tang_u(x,y,z,gamma_t,R,m)
        u_zeta,u_xi=skew_components(u_x,u_z,m)
        np.testing.assert_almost_equal(u_zeta , gamma_t               , decimal=5)
        np.testing.assert_almost_equal(u_xi  , -gamma_t*np.tan(chi/2) , decimal=5)
        np.testing.assert_almost_equal(u_z   ,  gamma_t               , decimal=5)
        np.testing.assert_almost_equal(u_x   ,  gamma_t*np.tan(chi/2) , decimal=5)
        
        #print('ux',u_x)
        #print('uy',u_y)
        #print('uz',u_z)
        #print('uzeta',u_zeta)
        #print('uxi',u_xi)

#     def test_singularities(self):
#         # TODO!
# 
#     def test_regularization(self):
#         #TODO
# 
#     def test_multirotor(self):
#         #TODO

    def test_SVC_rings(self):
        # Test that induction is close to the one obtained from a series of rings
        try:
            from .VortexRing import rings_u
        except:
            try:
                from wiz.VortexRing import rings_u
            except:
                from VortexRing import rings_u

        # Parameters
        chi = 30*np.pi/180
        m=np.tan(chi)
        gamma_t, R= -1, 10
        eps=10**-6 *R
        # Parameters for rings
        nRings      = 1000
        z_max       = 20*2*R
        Zr          = np.linspace(0,z_max,nRings)
        dzeta       = (Zr[1]-Zr[0])/np.cos(chi)
        vGamma_r    = Zr*0 + gamma_t*dzeta
        vR_r        = Zr*0 + R
        Xr          = m*Zr
        Yr          = 0*Zr

        def compare(x,y,z,dec):
            ux,uy,uz       = svc_tang_u(x,y,z,gamma_t,R,m, polar_out=False)
            ux_r,uy_r,uz_r = rings_u(x,y,z,vGamma_r,vR_r,Xr,Yr,Zr,polar_out = False)
            np.testing.assert_almost_equal(ux,ux_r,decimal=dec)
            np.testing.assert_almost_equal(uy,uy_r,decimal=dec)
            np.testing.assert_almost_equal(uz,uz_r,decimal=dec)
            return ux,uy,uz,ux_r,uy_r,uz_r
        # --- test on rotor
        x0=np.linspace(-2*R,2*R,40)
        x,y,z=x0,x0*0,x0*0
        b=np.abs(np.sqrt((x-z*m)**2)-R)>0.1*R
        x,y,z=x[b],y[b],z[b]
        ux,uy,uz,ux_r,uy_r,uz_r=compare(x,y,z,1)
        # --- test at -R downstream
        x,y,z=x0,x0*0,x0*0-R
        b=np.abs(np.sqrt((x-z*m)**2)-R)>0.1*R
        x,y,z=x[b],y[b],z[b]
        ux,uy,uz,ux_r,uy_r,uz_r=compare(x,y,z,2)
        # --- test at +R upstream
        x,y,z=x0,x0*0,x0*0+R
        b=np.abs(np.sqrt((x-z*m)**2)-R)>0.1*R
        x,y,z=x[b],y[b],z[b]
        ux,uy,uz,ux_r,uy_r,uz_r=compare(x,y,z,2)

        #import matplotlib.pyplot as plt
        #plt.figure()
        #plt.plot(x,ux)
        #plt.plot(x,ux_r)
        #plt.figure()
        #plt.plot(x,uy)
        #plt.plot(x,uy_r)
        #plt.figure()
        #plt.plot(x,uz)
        #plt.plot(x,uz_r)
        #plt.show()

if __name__ == "__main__":
#     TestCylinder().test_singularities()
    unittest.main()
