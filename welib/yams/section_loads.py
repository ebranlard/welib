import numpy as np
import numpy as np
try:
    from scipy.integrate import cumulative_trapezoid 
except:
    from scipy.integrate import cumtrapz as cumulative_trapezoid
try:
    from numpy import trapezoid
except:
    from numpy import trapz as trapezoid


def beamSectionLoads1D(z, p, Ftop=0, Mtop=0, s=1, F_lumped=None, method='plin'):
    r"""
    Integrate section loads along a beam based on inline loads, lumped loads and top load.
    All the loads are assumed to be in the same direction.
    NOTE: loads in other directions will actually affect the section loads in the current direction.
          To account for the contributions from other components, see beamSectionLoads

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
        r""" 
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
            Fsec[i]=trapezoid(pabove, zabove) + Ftop + np.sum(Fabove)
            Msec[i]=trapezoid(pabove*(zabove-zcur), zabove)+Ftop*(z[-1]-zcur) + np.sum(Fabove*(zabove-zcur))
        Msec+=Mtop

    elif method=='cumtrapz':
        # NOTE: might not work as well when lumped forces are present
        zn = z[-1::-1]  # flip z so it goes from top to bottom for cumtrapz
        Fsec[:]  = Ftop
        Msec[:]  = Mtop
        Fsec      += np.cumsum(F_lumped[-1::-1])[-1::-1]
        Fsec[:-1] +=  - cumulative_trapezoid(p[-1::-1], zn)[-1::-1]    
        Msec[:-1] +=  - cumulative_trapezoid(Fsec[-1::-1], zn)[-1::-1] 

    else:
        raise NotImplementedError()

    return Fsec, Msec


def beamSectionLoads3D(p_ext, F_top, M_top, s_span, m, 
                       U=None, V=None, K=None, a_struct=None, 
        M_lumped=None, m_hydro=None, a_ext=None, F_ext_lumped=None, 
        corrections=1,
        bSelfWeight=False,
        main_axis='z',
        nWARNS=[0],
        debug=False
        ):
    """ 
    Determine section loads along a beam by integration of the external forces (distributed or lumped),
    inertia loads, and tip loads

    INPUTS:
     - p_ext: loads along the beam (without inertia), shape: (3 x n) [N/m]
     - F_top: beam top force , shape:(3) [N]
     - M_top: beam top moment, shape:(3) [N/m]
     - s_span: coordinates along the beam, shape:(n) [m]
     - m     : mass per length of the beam, shape:(n) [kg/m]
     - m_hydro: hydrodynamic added mass, shape: (n) [kg/m]
     - U     : Displacements along the span, shape:(3,n) [m]
     - V     : Slopes along the span       , shape:(3,n) [m/m]
     - K     : Curvature along the span    , shape:(3,n) [m/m/m]

     - M_lumped: lumped masses along the span: shape(n) [kg] (zero where no mass)
     - F_lumped: lumped forces along the span: shape(n) [N]  (zero where no force)

     - a_struct : structural acceleration of the beam, shape: (3,n)
     - a_ext : external linear-acceleration vector, (0,0,-g) for gravity

    OUTPUTS:
     - F_sec: section forces ,  shape:(3,len(s_span))
     - M_sec: section moments,  shape:(3,len(s_span)) 
     - out  : dictionary with additional info
    
    """
    # Main dimensions
    nSpan     = len(s_span)
    # --- Default outputs
    F_sec = np.zeros((3,nSpan))
    M_sec = np.zeros((3,nSpan))
    outD = dict()

    # --- Default values
    if m_hydro is None:
        m_hydro  = np.zeros(nSpan)           # added mass, only on wet surface of structure
    if M_lumped is None:
        M_lumped = np.zeros(nSpan)
    if F_ext_lumped is None:
        F_ext_lumped=np.zeros(nSpan)
    if a_struct is None:
        a_struct = np.zeros((3,nSpan))
    if a_ext is None:
        a_ext = [0, 0, 0]

    # Linear Translation, Velocity, Acceleration
    if U is None:
        U = np.zeros((3,nSpan))
    if V is None:
        print('[WARN] yams: flexibility: beamSectionLoads3D: computing V as gradient U')
        V = np.zeros((3,nSpan)) 
        V[0,:] = np.gradient(U[0,:],  s_span, edge_order=2)
        V[1,:] = np.gradient(U[1,:],  s_span, edge_order=2)
        V[2,:] = np.gradient(U[2,:],  s_span, edge_order=2)
        #from welib.mesh.gradient import gradient_regular
        #V[0,:] = gradient_regular(U[0,:], dx=s_span[1]-s_span[0], order=4)
        #V[1,:] = gradient_regular(U[1,:], dx=s_span[1]-s_span[0], order=4)
        #V[2,:] = gradient_regular(U[2,:], dx=s_span[1]-s_span[0], order=4)
    if K is None:
        if nWARNS[0]<5:
            print('[WARN] yams: flexibility: beamSectionLoads3D: computing K as gradient V')
            nWARNS[0]+=1
        K = np.zeros((3,nSpan)) 
        K[0,:] = np.gradient(V[0,:],  s_span, edge_order=2)
        K[1,:] = np.gradient(V[1,:],  s_span, edge_order=2)
        K[2,:] = np.gradient(V[2,:],  s_span, edge_order=2)

    # --- Acceleration
    a_tot = a_struct.copy()
    # Typically, a_ext is gravity
    # TODO Body root acceleration!
    a_tot[0,:] -= a_ext[0]
    a_tot[1,:] -= a_ext[1]
    a_tot[2,:] -= a_ext[2]

    # --- Inertial loads
    m_struct = m
    m_tot = m_struct + m_hydro
    p_inertia     =   m_struct * a_tot 
    p_inertia_st  =   m_struct * a_struct
    p_inertia_ext = - m_struct * a_ext[:, np.newaxis]
    F_inertia_lumped = M_lumped * a_tot

    p_inertia_hydro =   m_hydro * a_struct 

    # --- Total loads from external forces and inertia (inline and lumped)
    p_all        = p_ext        - p_inertia_st - p_inertia_ext -p_inertia_hydro
    F_lumped_all = F_ext_lumped - F_inertia_lumped

    # --- Axial force 
    p_corr = np.zeros((3,nSpan))
    p_x = np.zeros(nSpan) 
    pax = np.zeros(nSpan) 
    Pax = np.zeros(nSpan)  # Cumulative axial force \int_z^L pax dz
    if main_axis=='z':
        if bSelfWeight:
            pax_SW = - p_inertia[2,:] #   - m * (g + zddot)
            Pax_SW  = fcumtrapzlr(s_span, pax_SW)
            pax = pax + pax_SW # TODO lumped forces
            Pax = Pax + Pax_SW # TODO lumped forces
            #print(m)
            #print(np.sum(m))
            #print(pax)
            #print(np.sum(pax))
            #print(Pax)
            #Pax[:] += F_top[2] # Not so sure here
            #pax[-1] += F_top[2]
            p_x +=  K[0,:] * Pax
            p_x += - V[0,:] * pax
            p_corr[0,:] +=p_x
    # print('TODO TODO NEW Beam Section TODO Include')
    #    quadratic velocity terms
    #    Coriolis-type coupling
    #    axial–bending coupling due to large deflection kinematics
    #    Reference tools keep at least first-order convective inertia, which produces small oscillations at wave frequency.
    #    This is usually small but exactly the kind of “missing wiggles” you describe.


    p_all += p_corr

    # TODO self-weight correction
#     # FT=fcumtrapzlr(s_span,m);
#     FT = - sciint.cumtrapz( m[-1::-1], s_span[-1::-1],)[-1::-1] 
#     FT = np.concatenate((FT,[0]))
#     if V_tot is None: 
#         raise Exception('Please provide Vtot for axial correction'); end
#     if main_axis=='x':
#         Mxt[0,1]=+trapzs(V_tot[2,:]*FT) # m15
#         Mxt[0,2]=-trapzs(V_tot[1,:]*FT) # m16
#     else:
#         # TODO TODO TODO VERIFY ME
#         Mxt[2,0]=+trapzs(V_tot[1,:]*FT) # m15
#         Mxt[2,1]=-trapzs(V_tot[0,:]*FT) # m16
#     # --- Axial force 
#     Pacc    = np.zeros(nSpan) 
#     # TopMass contribution to Pacc
#     if bMtop:
#         Pacc_MT = -Mtop * gravity*np.ones(nSpan)
#         Pacc=Pacc+Pacc_MT
#     if bSelfWeight:
#         Pacc_SW  = fcumtrapzlr(s_span, -m * gravity)
#         Pacc=Pacc+Pacc_SW
#     if bRot:
#         Pacc_Rot = fcumtrapzlr(s_span,  m * Omega**2 * s_span)
#         Pacc=Pacc+Pacc_Rot
#     # Method 2
#     KKCorr = np.zeros((nf,nf))
#     for i in range(0,nf):
#         for j in range(0,nf):
#             #xx=trapz(s_span, Pacc .* PhiV{i}(1,:).* o.PhiV{j}(1,:));
#             if main_axis=='x':
#                 yy=trapzs(Pacc * dU[i][1,:] * dU[j][1,:])
#                 zz=trapzs(Pacc * dU[i][2,:] * dU[j][2,:])
#                 KKCorr[i,j]=yy+zz
#             elif main_axis=='z':
#                 xx=trapzs(Pacc * dU[i][0,:] * dU[j][0,:])
#                 yy=trapzs(Pacc * dU[i][1,:] * dU[j][1,:])
#                 KKCorr[i,j]=yy+xx

    # --- Intermediate storage
    outD['m']               = m
    outD['m_tot']           = m_tot
    outD['m_hydro']         = m_hydro
    outD['a_struct']        = a_struct
    outD['a_ext']           = a_ext
    outD['p_inertia']       = p_inertia
    outD['p_inertia_st']    = p_inertia_st
    outD['p_inertia_ext']   = p_inertia_ext
    outD['p_inertia_hydro'] = p_inertia_hydro
    outD['p_ext']           = p_ext
    outD['p_corr']          = p_corr
    outD['p_all']           = p_all
    outD['F_lumped_all']    = F_lumped_all


    # --- Section Loads
    z  = s_span-s_span[0]
    zn = z[-1::-1]  # flip z so it goes from top to bottom for cumtrapz
    # Bending momemts 
    F_sec[0,:], M_sec[1,:] = beamSectionLoads1D(z, p_all[0,:], F_top[0], M_top[1], s=1,  F_lumped = F_lumped_all[0,:])
    F_sec[1,:], M_sec[0,:] = beamSectionLoads1D(z, p_all[1,:], F_top[1], M_top[0], s=-1, F_lumped = F_lumped_all[1,:])
    # Axial force
    F_sec[2,:-1] =- cumulative_trapezoid(p_all[2, -1::-1], zn)[-1::-1] # NOTE: mostly m*acc, can use FXG
    F_sec[2,:] += F_top[2] 
    # Torsion moment
    M_sec[2,:] += M_top[2]  # TODO integrate external torsions - torsional inertias and contributions from sectionn loads due to lever arm of deflection

    # Additional forces and moments from top loads due to deflections (ExtraLeverArm)
    if corrections>=1:
        F_sec[0,:] += -F_top[2] * V[0,:] # Fx = Fz v_y 
        F_sec[1,:] +=  F_top[2] * V[1,:] # Fy = Fz v_x  # TODO check sign
        dx = U[0,-1] - U[0,:]
        dy = U[1,-1] - U[1,:]
        M_sec[1,:] += -F_top[2] * dx # My =-Fz dx 
        M_sec[0,:] += +F_top[2] * dy # Mx = Fz dy 
        M_sec[2,:] +=  F_top[1]*dx - F_top[0]*dy # Mz = Fy dx - Fx dy 


    # Torsion correction
    if corrections>=2:
        M_sec[2,1:] +=- V[1,1:]*M_sec[0,1:]-V[0,1:]*M_sec[1,1:] # Mx = - Vy Mx - Vx My # TODO check sign
    
    if debug:
        print('Debug in Flexibility beamSectionLoads3D')
        import pdb; pdb.set_trace()

    # KEEP ME: M_y approximation
    #M_sec[1,0] = F_top[1]*z[-1] + M_top[1] # approximation

    return F_sec, M_sec, outD



def beamSectionLoadsFromShapeFunctions(x, xd, xdd, p_ext, F_top, M_top, s_span, PhiU, PhiV, m, 
        M_lumped=None, m_hydro=None, a_ext=None, F_ext_lumped=None, corrections=1, PhiK=None, debug=False):
    """ 
    Compute section loads along a beam represented by shape functions
    INPUTS:
     - x, xd, xdd : elastic motion associated with the shape functions
                 array-like of shape nf
     - p_ext: loads along the beam (without inertia), shape: (3 x n) [N/m]

     - PhiU, PhiV, PhiK: Deflections, slopes, curvature of shape functions (nf x 3 x n)  
    
    """
    # Main dimensions
    shapeDisp = PhiU[0].shape
    nf        = len(PhiU)
    # Linear Translation, Velocity, Acceleration
    U        = np.zeros(shapeDisp)
    V        = np.zeros(shapeDisp)
    #v_struct = np.zeros(shapeDisp)
    a_struct = np.zeros(shapeDisp)

    # Compute elastic deformation, slope curvature:
    for j in np.arange(nf):
        U         += x  [j] * PhiU[j] # Deflections
        V         += x  [j] * PhiV[j] # Slopes
        #v_struct  += xd [j] * PhiU[j]
        a_struct  += xdd[j] * PhiU[j] # TODO base motion
    K = None
    if PhiK is not None:
        K = np.zeros(shapeDisp)
        for j in np.arange(nf):
            K += x[j] * PhiK[j] # Deflections

    return beamSectionLoads3D(p_ext=p_ext, F_top=F_top, M_top=M_top, s_span=s_span, m=m, U=U, V=V, K=K, a_struct=a_struct, 
            M_lumped=M_lumped, m_hydro=m_hydro, a_ext=a_ext, F_ext_lumped=F_ext_lumped, corrections=corrections, debug=debug)

#     def moments(KF):
#         WT=KF.WT2
#         z_test = fastlib.ED_TwrGag(WT.ED) - WT.ED['TowerBsHt']
#         EI     = np.interp(z_test, WT.Twr.s_span, WT.Twr.EI[0,:])
#         kappa  = np.interp(z_test, WT.Twr.s_span, WT.Twr.PhiK[0][0,:])
#         qx    = KF.X_hat[KF.iX['ut1']]
#         KF.M_sim = [qx*EI[i]*kappa[i]/1000 for i in range(len(z_test))]                 # in [kNm]
#         KF.M_ref=[]
#         for i in range(len(z_test)):
#             try:
#                 val=KF.df['TwHt{:d}MLyt_[kN-m]'.format(i+1)].values
#             except:
#                 try:
#                     val=KF.df['TwHt{:d}MLyt'.format(i+1)].values
#                 except:
#                    val=KF.time*0
#             KF.M_ref.append(val)
#         return KF.M_sim, KF.M_ref
