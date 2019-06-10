# -*- coding: utf-8 -*-
""" 
Reference:
     E.Branlard: Wind turbine Aerodynamics and Vorticity Based Method - Chapter 10, Springer, 2017
"""
from __future__ import print_function,division
import numpy as np
from numpy import pi, cos, exp, sqrt, sin, arctan2, arccos
from scipy.interpolate import interp1d

import matplotlib.pyplot as plt
import pdb

np.seterr(all='raise')

def fAeroCoeffWrap(fPolars, alpha, phi, bAIDrag=True, bTIDrag=True):
    """Tabulated airfoil data interpolation
        Inputs
        ----------
        Polars: interpolant function for each alpha
        alpha: Angle Of Attack [rad]
        phi  : flow angle  [rad]

        Outputs
        ----------
        Cl,Cd  : lift and drag coefficients
        cn,ct  : normal and tangential coefficients (for load computation)
        cnForAI: normal and tangential coefficient  (for induction computation)
    """
    alpha[alpha<-pi] += 2*pi
    alpha[alpha> pi] -= 2*pi
    Cl = np.zeros(alpha.shape)
    Cd = np.zeros(alpha.shape)
    for i,(fPolar,alph) in enumerate(zip(fPolars,alpha)):
        ClCdCm = fPolar(alph)
        Cl[i], Cd[i] = ClCdCm[0], ClCdCm[1]
    # --- Normal and tangential
    cn = Cl * cos(phi) + Cd * sin(phi)
    ct = Cl * sin(phi) - Cd * cos(phi)
    if (bAIDrag):
        cnForAI = cn
    else:
        cnForAI = Cl*cos(phi) # cnNoDrag
    if (bTIDrag):
        ctForTI = ct
    else:
        ctForTI =  Cl * sin(phi) # ctNoDrag
    return Cl, Cd, cn, ct, cnForAI, ctForTI

def fInductionCoefficients(a_last, Vrel_norm, V0, F, cnForAI, ctForTI,
        lambda_r, sigma, phi, relaxation=0.4,bSwirl=True):
    """Compute the induction coefficients

        Inputs
        ----------
        a_last    : last iteration axial induction factor
        Vrel_norm : normed relative velocity
        V0        : free stream velocity
        F         : total loss
        cnForAI   : normal force coefficient
        ctForTI   : tangential force coefficient
        lambda_r  : speed ratio distribution
        sigma     : blade solidity
        phi       : flow angle [deg]
        relaxation: relaxation factor in axial induction factor
        bSwirl    : swirl flow model enabled / disabled

        Outputs
        ----------
        a: axial induction factor
        aprime: tangential induction factor
        Ct: local thrust coefficient
    """
    # --- Default a and CT
    a = 1. / ((4.*F*sin(phi)**2)/(sigma*(cnForAI+10**-8))+1) # NOTE simgularity avoided
    # CT=(1-a_last).^2.*sigma.*CnForAI./((sind(phi)).^2)
    Ct = Vrel_norm**2 * sigma * cnForAI/(V0**2)  # that's a CT loc
    # --- Hight thrust correction
    # Glauert correction
    ac = 0.3
    bHigh = a > ac
    fg = 0.25*(5.-3.*a[bHigh])
    a[bHigh] = Ct[bHigh]/(4.*F[bHigh]*(1.-fg*a[bHigh]))
    #a_high=0.5*(2+K*(1-2*ac)-sqrt((K*(1-2*ac)+2)^2+4*(K*ac^2-1)));
    # --- Relaxation
    a = a*relaxation + (1.-relaxation)*a_last

    # --- Swirl
    if bSwirl is True:
        #aprime = 1/((4*F*sin(phi)*cos(phi))/(sigma*ctForTI+10**-8)-1)
        # HAWC2 method:
        aprime = (Vrel_norm**2*ctForTI*sigma)/(4.*(1.-a)*V0**2*lambda_r)
    else:
        aprime = a * 0.
    # Bounding values for safety
    aprime = np.clip(aprime,-1,1.5) 
    a      = np.clip(a     ,-1,1.5)
    Ct     = np.clip(Ct    ,-1,3)
    return a, aprime, Ct


class MiniBEM_Outputs:
    pass

def MiniBEM(Omega,pitch,V0,xdot,u_turb,
        nB, cone, r, chord, twist, polars, # Rotor
        rho=1.225,KinVisc=15.68*10**-6,    # Environment
        nItMax=100, aTol=10**-6, bTipLoss=True, bHubLoss=False, bAIDrag=True, bTIDrag=True, bSwirl=True, relaxation=0.4, a_init=None, ap_init=None):
    """ Run the BEM main loop
        Inputs:
        -------
        Omega [rpm]:
        pitch [deg]:
        twist [deg]:
        cone  [deg]:
        r     [m]  : from rhub to R
        chord [m]  :
        polars     : nSpan matrices  

        Outputs
        ----------
        BEM : class with attributes, such as BEM.r, BEM.a, BEM.Power
    """
    VHubHeight = V0
    # --- Converting units
    fulltwist = (twist+pitch) *pi/180    # [rad]
    Omega    = Omega*2*pi/60 # [rad/s]
    # --- Derived params
    rhub, R  = r[0], r[-1]
    # Computing a dr, such that sum(dr)=R-rhub
    dr    = np.diff(r) 
    MidPointAfter = np.concatenate((  r[0:-1]+dr/2 , [R] ))
    MidPointBefore= np.concatenate(( [r[0]] ,  r[1:]-dr/2))
    dr    = MidPointAfter-MidPointBefore
    cCone    = cos(cone*pi/180.)
    sigma    = chord * nB / (2.0 * pi * r * cCone)
    lambda_r = Omega * r * cCone/ V0
    # Creating interpolation functions for each polar, now in rad!
    fPolars = [interp1d(p[:,0]*pi/180,p[:,1:],axis=0) for p in polars]
    # Initializing outputs
    if a_init is None:
        a_init = np.ones((len(r)))*0.2
    if ap_init is None:
        ap_init = np.ones(len(r))*0.01
    # --- Vectorized BEM algorithm
    # Radial inputs: a_init,ap_init,r,chord,fulltwist,sigma,lambda_r,fPolars
    a, aprime  = a_init, ap_init
    for i in np.arange(nItMax):
        # --------------------------------------------------------------------------------
        # --- Step 0: Relative wind
        # --------------------------------------------------------------------------------
        # --------------------------------------------------------------------------------
        # --- Step 1: Wind Components
        # --------------------------------------------------------------------------------
        Ut = Omega * r * (1. + aprime)
        Un = V0 * (1. - a) - xdot + u_turb
        Vrel_norm = np.sqrt(Un** 2 + Ut** 2)
        # --------------------------------------------------------------------------------
        # --- Step 2: Flow Angle
        # --------------------------------------------------------------------------------
        phi = arctan2(Un, Ut) # flow angle [rad]
        # --------------------------------------------------------------------------------
        # --- Tip loss
        # --------------------------------------------------------------------------------
        Ftip = np.ones((len(r)))
        Fhub = np.ones((len(r)))
        IOK=sin(phi)>0.01
        try:
            if bTipLoss:
                # Glauert tip correction
                Ftip[IOK] = 2/pi*arccos(exp(-nB/2*(R-r[IOK])/(r[IOK]*sin(phi[IOK]))))
            if bHubLoss:
                # Prandtl hub loss correction
                Fhub[IOK] = 2/pi*arccos(exp(-nB/2*(r[IOK]-rhub)/(rhub*sin(phi[IOK]))));
        except:
            raise
#             pdb.set_trace()
        F=Ftip*Fhub;
        F[F<=0]=0.5 # To avoid singularities
        # --------------------------------------------------------------------------------
        # --- Step 3: Angle of attack
        # --------------------------------------------------------------------------------
        alpha = phi - fulltwist # [rad], contains pitch
        # --------------------------------------------------------------------------------
        # --- Step 4: Profile Data
        # --------------------------------------------------------------------------------
        Cl, Cd, cn, ct, cnForAI, ctForTI = fAeroCoeffWrap(fPolars, alpha, phi, bAIDrag, bTIDrag)
        # --------------------------------------------------------------------------------
        # --- Step 5: Induction Coefficients
        # --------------------------------------------------------------------------------
        # Storing last values
        a_last      = a
        aprime_last = aprime
        a, aprime, CT_loc = fInductionCoefficients(a_last,Vrel_norm,V0, F, cnForAI, ctForTI,
                                               lambda_r, sigma, phi, relaxation, bSwirl)

        if (i > 3 and (np.mean(np.abs(a-a_last)) + np.mean(np.abs(aprime - aprime_last))) < aTol):  # used to be on alpha
            break
    nIt = i + 1
    if i == nItMax-1:
        print('Maximum iterations reached : Omega=%.2f V0=%.2f' % (Omega,V0))
    #print('Converged: V0=%2.f om=%5.2f pit=%3.1f nIt=%d' % (V0,Omega,pitch, nIt))
    # --------------------------------------------------------------------------------
    # --- Step 6: Outputs
    # --------------------------------------------------------------------------------
    BEM=MiniBEM_Outputs();
    BEM.a,BEM.aprime,BEM.phi,BEM.Cl,BEM.Cd,BEM.Un,BEM.Ut,BEM.Vrel,BEM.F,BEM.nIt = a,aprime,phi,Cl,Cd,Un,Ut,Vrel_norm,F,nIt
    # L = 0.5 * rho * Vrel_norm ** 2 * chord[e] * Cl
    # D = 0.5 * rho * Vrel_norm ** 2 * chord[e] * Cd
    # Radial quantities (recomputed since thought as derived outputs)
    BEM.cn = BEM.Cl * cos(BEM.phi) + BEM.Cd * sin(BEM.phi)
    BEM.ct = BEM.Cl * sin(BEM.phi) - BEM.Cd * cos(BEM.phi)
    BEM.Pn    = 0.5 * rho * BEM.Vrel**2 * chord * BEM.cn   # [N/m]
    BEM.Pt    = 0.5 * rho * BEM.Vrel**2 * chord * BEM.ct   # [N/m] 
    BEM.alpha = (BEM.phi - fulltwist)*180/pi               # [deg]
    BEM.phi   = BEM.phi*180/pi                             # [deg]
    BEM.Re    = BEM.Vrel * chord / KinVisc / 10**6  # Reynolds number in Millions
    BEM.Gamma = 0.5 * BEM.Vrel * chord * BEM.Cl   # Circulation [m^2/s]
    # Radial quantities, "dr" formulation
    BEM.ThrLoc   = dr * BEM.Pn * cCone
    BEM.ThrLocLn = BEM.Pn * cCone
    BEM.Ct       = nB * BEM.ThrLoc / (0.5 * rho * VHubHeight** 2 * (2*pi * r * cCone * dr))
    BEM.TqLoc    = dr * r * BEM.Pt * cCone
    BEM.TqLocLn  = r * BEM.Pt * cCone
    BEM.Cq      = nB * BEM.TqLoc / (0.5 * rho * VHubHeight** 2 * (2*pi * r * cCone)) * dr * r * cCone
    BEM.Cp      = BEM.Cq*lambda_r
    # --- Integral quantities
    BEM.Torque = nB * np.trapz(r * (BEM.Pt * cCone), r)  # Rotor shaft torque [N]
    BEM.Thrust = nB * np.trapz(     BEM.Pn * cCone, r)   # Rotor shaft thrust [N]
    BEM.Flap   = np.trapz( (BEM.Pn * cCone) * (r - rhub), r)      # Flap moment at blade root [Nm]
    BEM.Edge   = np.trapz(  BEM.Pt * (r * cCone) * (r - rhub), r) # Edge moment at blade root [Nm]
    BEM.Power = Omega * BEM.Torque
    BEM.CP = BEM.Power  / (0.5 * rho * V0**3 * pi * R**2)
    BEM.CT = BEM.Thrust / (0.5 * rho * V0**2 * pi * R**2)
    BEM.CQ = BEM.Torque / (0.5 * rho * V0**2 * pi * R**3)
    BEM.r=r
    BEM.R=R
    BEM.uia = V0 * BEM.a
    BEM.uit = Omega * r * BEM.aprime
    BEM.Omega = Omega
    BEM.Pitch = pitch
    return BEM


def WriteRadialFile(BEM,filename):
    header='r_[m] a_[-] a_prime_[-] Ct_[-] Cq_[-] Cp_[-] cn_[-] ct_[-] phi_[deg] alpha_[deg] Cl_[-] Cd_[-] Pn_[N/m] Pt_[N/m] Vrel_[m/s] Un_[m/s] Ut_[m/s] F_[-] Re_[-] Gamma_[m^2/s] uia_[m/s] uit_[m/s]'
    header=' '.join(['{:14s}'.format(s) for s in header.split()])
    M=np.column_stack((BEM.r,BEM.a,BEM.aprime,BEM.Ct,BEM.Cq,BEM.Cp,BEM.cn,BEM.ct,BEM.phi,BEM.alpha,BEM.Cl,BEM.Cd,BEM.Pn,BEM.Pt,BEM.Vrel,BEM.Un,BEM.Ut,BEM.F,BEM.Re,BEM.Gamma,BEM.uia,BEM.uit))
    np.savetxt(filename,M,header=header,fmt='%14.7e',comments='#')


def FASTFile2MiniBEM(FASTFileName):
    import weio
    F=weio.FASTInputDeck(FASTFileName,readlist=['Aero','ED'])

    rho     = F.Aero['AirDens']
    KinVisc = F.Aero['KinVisc']

    nB   =  F.ED['NumBl']
    cone = -F.ED['PreCone(1)']
    r     = F.Aero.Bld1['BldAeroNodes'][:,0] + F.ED['HubRad']
    chord = F.Aero.Bld1['BldAeroNodes'][:,-2] 
    twist = F.Aero.Bld1['BldAeroNodes'][:,-3]
    polars=[]
    ProfileID=F.Aero.Bld1['BldAeroNodes'][:,-1].astype(int)
    for ipolar in  ProfileID:
        polars.append(F.Aero.AF[ipolar-1]['AFCoeff'])
    return nB,cone,r,chord,twist,polars,rho,KinVisc

if __name__=="__main__":
    import weio
    import pandas as pd
    from pybra.tictoc import Timer

    nB,cone,r,chord,twist,polars,rho,KinVisc = FASTFile2MiniBEM('OpenFAST_model/SNLV27_OF2.fst')

    df=weio.read('Simulated_data/RotorPerformances_OF2.csv').toDataFrame();
    Pref      = df['GenPower_[kW]']
    Tref      = df['AeroThrust_[kN]']
    CPref     = df['AeroCp_[-]']
    CTref     = df['AeroCt_[-]']
    FlapRef   = df['BldRootFlapM_[kN.m]']
    Omega_ref = df['RotSpeed_[rpm]']
    V0_ref    = df['WS_[m/s]']

#     fig=plt.figure(); ax=fig.add_subplot(111)
#     fig2=plt.figure(); ax2=fig2.add_subplot(111)
    with Timer('Total'):
        a0  = None
        ap0 = None
        for i in np.arange(len(df)):
            V0        = V0_ref[i]
            Omega     = Omega_ref[i]
            pitch=2     #[deg]
            xdot=0      #[m/s]
            u_turb=0    #[m/s]
            with Timer():
                BEM=MiniBEM(Omega,pitch,V0,xdot,u_turb,
                            nB,cone,r,chord,twist,polars,
                            rho=rho,KinVisc=KinVisc,bTIDrag=False,bAIDrag=True,
                            a_init =a0,
                            ap_init=ap0
                            )
                a0  = BEM.a
                ap0 = BEM.aprime

            filename='BEM_ws{}_radial.csv'.format(V0)
            WriteRadialFile(BEM,filename)

#             ax.plot(V0,BEM.Power*0.9/1000 ,'ro')
#             ax.plot(V0,BEM.Thrust/1000    ,'go')
#             ax.plot(V0,BEM.Flap/1000      ,'bo')
#             ax2.plot(V0,BEM.CP   ,'ro')
#             ax2.plot(V0,BEM.CT   ,'go')
#     ax.plot (V0_ref,Pref               ,'-',label = 'Power ref')
#     ax.plot (V0_ref,Tref               ,'-',label = 'Thurst ref')
#     ax.plot (V0_ref,FlapRef               ,'-',label = 'Flap ref')
#     ax2.plot(V0_ref,CPref   ,'-',label = 'CP ref')
#     ax2.plot(V0_ref,CTref   ,'-',label = 'CT ref')
#     ax.legend()
#     ax2.legend()
#     fig=plt.figure(); ax=fig.add_subplot(111)
#     ax.plot(BEM.r ,BEM.a     ,label='a'   )
#     ax.plot(BEM.r ,BEM.aprime,label='ap'  )
#     ax.legend()
#     ax2.legend()
#     plt.show()
#     print(F.Aero.AD['BldAeroNodes'])
#     for afName in F.Aero['AFNames']:
#         fpath = os.path.join(modeldir, afName.strip('"'))
#         print(fpath)

