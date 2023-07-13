""" 
Set of high thrust corrections

NOTE: ac=0.4 -> Ctc = (1-(1-2*ac)**2) = 0.96


Conventions (see [3])
-----------------------
NOTE: with skew (tilt or yaw), different conventions exist for the dimensionless coefficients:
- a, axial induction factor (what velocity is it multipled with)
- CT, thrust coefficient (what velocity is used to nondimensional the thrust)
different conventions exist

We write:
 - U0 the norm of the wind vector
 - Un the component of the wind normal to the actuator disk
(without skew, U0=Un)
We define two different thrust coefficients depending on the reference velocity:
 - CT,U0 = T / (1/2 rho A U0^2)
 - CT,Un = T / (1/2 rho A U0^2)

Similarly, we define two axial induction factors: 
 - Wn = - a_U0  U0  
 - Wn = - a_Un  Un 
where Wn is the inducted velocity component normal to the rotor plane.
(NOTE: there are further considerations when the disk is moving)

In this script we assume some kind of consistency between both definitions and define the following conventions:
 - "CtUn_aUn":  means "CT" is CT,Un and  "a" is  a_Un
 - "CtU0_aU0":  means "CT" is CT,U0 and  "a" is  a_U0


The default convention of this script is  "CtUn_aUn".
If there is no skew, the convention does not matter as both collapse to the same definitions.




References:
    [1] Branlard (2017) Wind turbine aerodynamics and vorticity-based methods. Springer
    [2] Madsen, Larsen, Pirrung, Li, Zahle (2020) Implementation of the blade element momentum model
on a polar grid and its aeroelastic load impact. Wind Energ. Sci.
    [3] Branlard (2024) 3d BEM To be published


"""
import numpy as np


MaxTanChi0 = 100  # maximum absolute value allowed for tan(chi0), an arbitary large number


def a_Ct(Ct, a=None, F=None, chi=0, CT=None, method='AeroDyn', outputMore=False, convention='CtUn_aUn'):
    """ 
    High thrust corrections of the form: 
       a=a(Ct)
   or 
       a=a(Ct,a)
    INPUTS:
      - Ct: Local thrust coefficient
      - a: induction factor
      - F : tip-loss factor
      - method
      - CT: mean CT (needed for method 'HAWC2'
      - convention: only matters if skew!=0 
            see top of the script for some documenation. Important when skew present!
            either 'CtUn_aUn' or 'CTU0_aU0'
    """


    Ct, scalar = to_array(Ct)
    F = array_like(F, Ct, 1)
    out={}
    # -------------------------------------------------------
    # --- a = a(Ct)
    # -------------------------------------------------------
    if method=='MomentumTheory':
        assertZero(chi) # Not suitable for skew
        a = np.zeros(Ct.shape)
        b = Ct<=1
        a[b] = 1/2*(1-np.sqrt(1-Ct[b]))
        a[~b] = np.nan

    elif method=='AeroDyn15': # Very close to Glauert Empirical, Very close to Buhl
        assertZero(chi) # Not suitable for skew
        Ct[Ct>2]  = 2
        Ct[Ct<-2] = -2
        Ic        = Ct/F>0.96 # Correction
        In        = np.logical_not(Ic) # Normal
        a = np.zeros(Ct.shape)
        a[Ic]     = 0.1432+np.sqrt(-0.55106+0.6427*Ct[Ic]/F[Ic])
        a[In]     = 0.5*(1-np.sqrt(1-Ct[In]/F[In]))
    elif method=='AeroDyn14':
        assertZero(chi) # Not suitable for skew
        Ic        = Ct/F>0.96 # Correction
        In        = np.logical_not(Ic) # Normal
        a = np.zeros(Ct.shape)
        a[Ic]     = (18*F[Ic]-20-3*np.sqrt(Ct[Ic]*(50-36*F[Ic]) + 12*F[Ic]*(3*F[Ic]-4) ) ) / (36*F[Ic] - 50)
        a[In]     = 0.5*(1-np.sqrt(1-Ct[In]/F[In])) # Baseline

    elif method=='Buhl':
        fa_Ct = a_Ct_numInverse(lambda a: Ct_a(a, method='Buhl', F=F[0])) # TODO F
        a = fa_Ct(Ct)

    elif method=='Glauert': 
        fa_Ct = a_Ct_numInverse(lambda a: Ct_a(a, method='Glauert', F=F[0])) # TODO F
        a = fa_Ct(Ct)

    elif method=='Spera':
        fa_Ct = a_Ct_numInverse(lambda a: Ct_a(a, method='Spera', F=F[0])) # TODO F
        a = fa_Ct(Ct)

    elif method=='GlauertEmpirical':
        assertZero(chi) # Not suitable for skew
        Ic    = Ct/F> 0.96  # Correction
        In    = np.logical_not(Ic) # Normal
        a=np.zeros(Ct.shape)
        a[Ic] = 1/F[Ic]*(0.143+np.sqrt(0.0203-0.6427*(0.889-Ct[Ic])))
        a[In] = 0.5*(1-np.sqrt(1-Ct[In]/F[In]))
    elif method=='WEHandbook':
        assertZero(chi) # Not suitable for skew
        Ic = Ct>0.96            # Correction
        In = np.logical_not(Ic) # Normal
        a=np.zeros(Ct.shape)
        a[Ic]   = 1/F[Ic]*(0.143 + np.sqrt( 0.0203-0.6427 *(0.889-Ct[Ic] ) ))
        a[In] = 0.5*(1-np.sqrt(1-Ct[In]/F[In]))
    elif method=='HAWC2':
        # NOTE: this is "CTU0_aU0" convention
        #k = [-0.0017077,0.251163,0.0544955,0.0892074]
        # Numerical values from [2]
        k = [0.0       ,0.2460  ,0.0586,   0.0883]
        Ct = Ct/F
        if chi==0:
            # Convention does not matter here
            ka = 1
            a = ka*(k[3]*Ct**3+k[2]*Ct**2+k[1]*Ct+k[0])
        else:
            Ca=np.array([[-0.5136 , 0.4438 , -0.1640],
                         [ 2.1735 ,-2.6145 ,  0.8646],
                         [-2.0705 , 2.1667 , -0.6481 ]])
            if CT is None:
                raise Exception('CT needs to be provided for HAWC2 method')
            if not scalar:
                CT= np.asarray(CT)
                if len(CT)>1:
                    raise Exception('CT must be of length 1')

            if convention.startswith('CtUn'):
                CTU0=CT*np.cos(chi)**2
                CtU0=Ct*np.cos(chi)**2
            elif convention.startswith('CtU0'):
                CTU0=CT
                CtU0=Ct
            else:
                raise NotImplementedError()
            CTU0 = min(CTU0,0.9)
            th = np.array([chi, chi**2, chi**3])
            ka1,ka2,ka3 = Ca.dot(th)
            ka = ka3 * CTU0**3 + ka2 * CTU0**2  + ka1 * CTU0  +1
            out['ka1'] = ka1
            out['ka2'] = ka2
            out['ka3'] = ka3
            out['ka']  = ka
            aU0 = ka*(k[3]*CtU0**3+k[2]*CtU0**2+k[1]*CtU0+k[0])
            if convention.find('aU0')>1:
                a = aU0 
            elif convention.find('aUn')>1:
                a = aU0/np.cos(chi)
            else:
                raise NotImplementedError()

    # -------------------------------------------------------
    # --- a = a(Ct,a) TODO TODO TODO might need rethinking..
    # with a = 1/2*(1-sqrt(1-Ct))
    # -------------------------------------------------------
    elif method=='Glauert_CTa': # see [1]
        assertZero(chi) # Not suitable for skew
        ac = ac_val(method=method)
        #Ic =Ct/F> 1-(1-(2*ac))**2  # Correction
        Ic = a>ac                  # Correction
        fg=0.25*(5-3*a[Ic])
        a[Ic] = Ct[Ic]/(4*F[Ic]*(1-fg*a[Ic]))
        #In = np.logical_not(Ic)    # Normal
        #a[In] = 0.5*(1-np.sqrt(1-Ct[In]/F[In]))
    elif method=='Spera_CTa':
        assertZero(chi) # Not suitable for skew
        ac = ac_val(method=method)
        Ic    = a>ac
        fgs   = ac/a[Ic]*(2-ac/a[Ic])
        a[Ic] = Ct[Ic]/(4*F[Ic]*(1-fgs*a[Iac]))
        #In = np.logical_not(Ic) # Normal
    elif method=='Shen_CTa':
        assertZero(chi) # Not suitable for skew
        ac = ac_val(method=method)
        Ic    = a>ac
        a[Ic] = (Ct[Ic]/4*F[Ic]-F[Ic]*ac**2)/(1-2*ac*F[Ic]);
    else:
        raise NotImplementedError('High-thrust correction method for `a-Ct` function unknown:'+method)

    if scalar:
        a=a[0]

    if outputMore:
        return a, out
    else:
        return a


def Ct_a(a, F=None, method='Glauert', ac=None, chi=0, convention='CtUn_aUn'):
    """ 
    High thrust corrections of the form: Ct = Ct(a)
        see a_Ct
    - convention: only matters if skew!=0 
          see top of the script for some documenation. Important when skew present!
          either 'CtUn_aUn' or 'CTU0_aU0'
    """
    a, scalar = to_array(a)
    F = array_like(F, a, 1)
    Ct = 4*a*F*(1-a)


    if method=='MomentumTheory':
        assertZero(chi) # Not suitable for skew
        Ct = 4*a*F*(1-a)

    elif method=='MomentumGlauertSkew':
        if convention.find('aU0')>1:
            a = a/np.cos(chi)       # input is aU0
        elif convention.find('aUn')>1:
            a = a  # input is aUn
        else:
            raise NotImplementedError()

        Ct = 4*a*F * np.sqrt((1-a)**2 + np.tan(chi)**2)

        if convention.startswith('CtUn'):
            Ct=Ct
        elif convention.startswith('CtU0'):
            Ct=Ct*np.cos(chi)**2
        else:
            raise NotImplementedError()


    elif method=='MomentumGlauertSkewCos':
        # NOTE: this is the exactly the same as momentum Glauert skew
        #       But the epxression if for "CTU0" and "aU0"
        if convention.find('aU0')>1:
            a0 = a               # input is already aU0
        elif convention.find('aUn')>1:
            a0 = a*np.cos(chi)   # input is aUn, convert to a0
        else:
            raise NotImplementedError()

        CtU0 = 4*a0*F * np.sqrt((1+a0**2 -2*a0 * np.cos(chi)))

        if convention.startswith('CtUn'):
            Ct=CtU0/np.cos(chi)**2
        elif convention.startswith('CtU0'):
            Ct=CtU0
        else:
            raise NotImplementedError()


    elif method=='Glauert': # Glauert hight thrust  # see [1]
        assertZero(chi) # Not suitable for skew
        if ac is None:
            ac = ac_val(method=method)
        Ic     = a>ac
        Ct[Ic] =4*a[Ic]*F[Ic]*(1-1/4*(5-3*a[Ic])*a[Ic])

    elif method=='Spera':
        assertZero(chi) # Not suitable for skew
        if ac is None:
            ac = ac_val(method=method)
        Ic     = a>ac
        Ct[Ic] =4*F[Ic]*(ac**2 + (1-2*ac)*a[Ic])

    elif method=='Leishman':
        assertZero(chi) # Not suitable for skew
        # NOTE: doesn't seem to work
        c0 = -3.402; c1 = 6.489; c2 = -3.167; c3 = 0.695; c4 = -0.0562
        if ac is None:
            ac = ac_val(method=method)
        Ic     = a>ac
        Ct[Ic] =F[Ic]*(c0 + c1*a[Ic] + c2*a[Ic]**2 + c3*a[Ic]**3 + c4*a[Ic]**4)

    elif method=='Buhl':
        assertZero(chi) # Not suitable for skew
        if ac is None:
            ac=ac_val(method=method)
        Ic = a>ac
        c0 = 8/9;
        c1 = 4*F[Ic]-40/9;
        c2 = 50/9-4*F[Ic];
        Ct[Ic] =c0 + c1*a[Ic] + c2*a[Ic]**2


    elif method=='BladedNoCorr':
        # TODO TODO MIGHT BE CONVENTION U0 
        raise NotImplementedError('TODO FIGURE OUT CONVENTION HERE')
        if convention=='CtU0_aU0':
            ac = ac_val(chi, method=method, convention='aU0')
            for ia, (aa,FF) in enumerate(zip(a,F)):
                if aa<ac:
                    Ct[ia] = 4*aa*FF*(1-aa)
                else:
                    c = getEmpiricalCoefficients(chi=chi, F=FF, method=method)
                    Ct[ia] =c[0] + c[1]*aa + c[2]*aa**2
           #if (k <= k0 ) then
           #    if (VxCorrected2 > 0.0) then
           #        a = k/(k+1.0)
           #    else
           #        a = k/(k-1.0)
           #    end if
           #    H = 1.0_R8Ki
           #else
           #    call axialInductionFromEmpiricalThrust( effectiveYaw, phi, k, F, a, H, MomentumCorr )
           #endif
        else:
            raise NotImplementedError()

    elif method=='BladedCorr':
        if convention=='CtUn_aUn':
            ac = ac_val(chi, method=method, convention='aUn')
            for ia, (aa,FF) in enumerate(zip(a,F)):
                if aa<ac:
                    Ct[ia] = 4*aa*FF*np.sqrt((1-aa)**2 + np.tan(chi)**2)
                else:
                    c = getEmpiricalCoefficients(chi=chi, F=FF, method=method)
                    Ct[ia] =c[0] + c[1]*aa + c[2]*aa**2
        else:
            raise NotImplementedError()

    elif method=='Branlard':
        if convention=='CtUn_aUn':
            ac = ac_val(chi, method=method, convention='aUn')
            for ia, (aa,FF) in enumerate(zip(a,F)):
                if aa<ac:
                    Ct[ia] = 4*aa*FF*np.sqrt((1-aa)**2 + np.tan(chi)**2)
                else:
                    c = getEmpiricalCoefficients(chi=chi, F=FF, method=method)
                    Ct[ia] =c[0] + c[1]*aa + c[2]*aa**2
        else:
            raise NotImplementedError()

    # --- Numerical inverses NOTE: only works for F scalar
    elif method=='HAWC2':
        fCt_a = Ct_a_numInverse(lambda ct: a_Ct(ct, chi=chi, CT=ct, F=F[0], method=method, convention=convention))
        Ct = fCt_a(a)
    else:
        raise NotImplementedError('High-thrust correction method for `Ct_a` function unknown:'+method)
    if scalar:
        Ct = Ct[0]
    return Ct

def ac_val(chi=0, method='Bladed', convention='aUn'):
    """ Return the critical value of the axial induction factor above which the high-thrust
        correction should be used 
    """

    if method in ['Glauert','Spera', 'Shen_CTa','Leishman']: # see [1]
        # For chi=0
        ac_n = 1/3  / np.cos(chi)
    elif method=='Glauert_CTa': # see [1]
        # For chi=0
        ac_n = 0.3  / np.cos(chi)
    elif method=='Spera_CTa':
        ac_n  = 0.34  / np.cos(chi)
    elif method=='Buhl':
        ac_n = 0.4 / np.cos(chi)
    elif method in ['Bladed','BladedCorr']:
        # a0 = 0.5*cos(45) ~= 0.3535
        ac_n = 0.5*np.cos(45*np.pi/180)/np.cos(chi)
        # NOTE: for polynomial continuation we don't want ac to get close to 1
        ac_n = np.clip(ac_n, 0, 0.5)
    elif method=='Branlard':
        ac_n = 0.35/ np.cos(chi)
        ac_n = np.clip(ac_n, 0, 0.5)
    else:
        raise NotImplementedError('Method: {}'.format(method))

    if convention=='aUn':
        return ac_n 
    else:
        return ac_n * np.cos(chi)

def ak_lim(chi, method='Bladed', convention='aUn'):
    a0 = ac_val(chi, method=method, convention=convention)
    k0 = a0/(1-a0)
    tan_chi0 = np.clip(np.tan(chi),  -MaxTanChi0, MaxTanChi0)
    if convention=='aUn':
        k0_lim = k0*np.sqrt(1+(tan_chi0/(1-a0))**2)
    else:
        raise NotImplementedError()
    return a0, k0, k0_lim, tan_chi0




def getEmpiricalCoefficients(chi=0, F=1, method='BladedCorr'):
    """ 
    Compute the coefficients of a second order polynomial that extends the Momenutm relationship CT(a) 
    above a value a>ac. The continuation is done such that the slope and value at a=a_c match 
    the momentum relation. The last constraint is the value of CT at a=1. 
    Currently a hard-coded model is used for the value at at=1.
    The polynomial is:
       CT(a) = c0 + c1*a + c2*a2    a>ac
    obtained with the constraints:
       CT(a_c)     = CT_c
       CT(1)       = CT_1
       dCT/da(a_c) = s_c

        CT = 4*a*(1-a)*F = c0 + c1*a + c2*a^2 + ...

    """
    if hasattr(F, '__len__'):
        raise Exception('F should be scalar')
    if hasattr(chi, '__len__'):
        raise Exception('chi should be scalar')

    if method=='Branlard':
        # Empirical CT = c2*a^2 + c1*a + c0 for a > ac
        #ac = ac_val(chi)
        ac = ac_val(chi, method=method)
        denom = (1-ac)**2
        tchi2 = np.clip(np.tan(chi), -MaxTanChi0, MaxTanChi0)**2  # tan(chi)**2
        # Value and slope at a=ac
        CT_c = 4*F*ac * np.sqrt( (1-ac)**2 + tchi2 ) 
        s_c = 4*F*(1-3*ac+2*ac**2+tchi2)/np.sqrt( (1-ac)**2 + tchi2 ) 
        s_c = max(s_c, 1)  # We don't want slopes<1
        # Empirial value of CT at a=1
        #CT_1 = max(1, 1.6822 +  -2.5176*(1-1/np.cos(chi)**1.35) )
        #CT_1= max(1, 1.86 +  -3.5*(1-1/np.cos(chi)**1.25))
        #CT_1 = -6.5470946189396155*(1-1/np.cos(chi)**0.8074219429690564)+2.016788590000479
        #CT_1 = Ct_a(1, chi=chi, method='HAWC2')
        CT_1 =  2 + 2.0*np.tan(chi)**1.6
        CT_1 =  2 + 2.113*np.tan(chi)**1.527
        CT_1 =  max(CT_1, CT_c + s_c * (1-ac) + 0.001 ) #Make sure c2>0
        c0, c1, c2 = polynomialOrder2Continuation(ac, s_c, CT_c, CT_1)
        if (CT_1-CT_c -s_c * (1-ac))<0:
            print('>>>', c2, chi)
        if c2<0: 
            print('===', c2, chi)
        #c0, c1 = linearContinuation(ac, s_c, CT_c); c2=0

    elif method=='BladedCorr':
        # Continuation of Glauert Skew Momentum    CT= 4 a F sqrt( (1-a)^2 + tan(chi)^2 ) 
        # Empirical CT = c2*a^2 + c1*a + c0 for a > ac
        ac = ac_val(chi, method=method)
        denom = (1-ac)**2
        tchi2 = np.clip(np.tan(chi), -MaxTanChi0, MaxTanChi0)**2  # tan(chi)**2
        # Value and slope at a=ac
        CT_c = 4*F*ac * np.sqrt( (1-ac)**2 + tchi2 ) 
        s_c = 4*F*(1-3*ac+2*ac**2+tchi2)/np.sqrt( (1-ac)**2 + tchi2 ) 
        # Empirial value of CT at a=1
        CT_1 = max(1, np.sqrt(((-0.64755/(np.cos(chi)*np.cos(chi)) - 0.8509/np.cos(chi) + 3.4984)*F)**2 + 16*tchi2))
        c0, c1, c2 = polynomialOrder2Continuation(ac, s_c, CT_c, CT_1)

    elif method=='BladedNoCorr':
        # Continuation of Glauert Momentum    CT= 4 a F (1-a)
        # Empirical CT = 4*a*(1-a)*F = c2*a^2 + c1*a + c0 for a > a0
        ac = ac_val(chi, method=method)
        CT_1 = max(1, (-0.64755/(np.cos(chi)*np.cos(chi)) - 0.8509/np.cos(chi) + 3.4984)*F )
        CT_c = 4*F*ac * (1-ac)  # CT(ac)
        s_c  = 4*F*(1-2*ac)     # dCT/da(ac) (slope)
        c0, c1, c2 = polynomialOrder2Continuation(ac, s_c, CT_c, CT_1)
    else:
        raise NotImplementedError('High-thrust polynomial coefficients method unknown: '+method)
    return [c0, c1, c2]


def a_Ct_numInverse(fCta):
    """ take a Ct(a) function and make it a a(Ct) function"""
    from  scipy.optimize import minimize_scalar
    def a_Ct_local(Ct): 
        # TODO add F here to interface
        Ct, scalar = to_array(Ct)
        a = []
        for ct in Ct:
            minimize = lambda a : abs(fCta(a) - ct) # NOTE: abs value
            res = minimize_scalar(minimize, bounds=(0,1), method='bounded') # a between 0 and 1 for now..  
            a.append(res.x)
        if scalar:
            return a[0]
        else:
            return a

    return a_Ct_local

def Ct_a_numInverse(faCt):
    """ take a a(Ct) function and make it a Ct(a) function"""
    from  scipy.optimize import minimize_scalar
    def Ct_a_local(a): 
        # TODO add F here to interface
        a, scalar = to_array(a)
        Ct = []
        for aa in a:
            minimize = lambda ct : abs(faCt(ct) - aa) # NOTE: abs value
            res = minimize_scalar(minimize, bounds=(0,25), method='bounded') # Ct between 0 and 8 for now..  
            Ct.append(res.x)
        if scalar:
            return Ct[0]
        else:
            return Ct
    return Ct_a_local





def a_k(k, chi=0, phi=0.1, F=None, method='Buhl', outputMore=False):
    """
    Return a = a(k) 
    INPUTS:
      - k:
      - phi: flow angle mostly used for sign (positive or negative)
      - F: tiploss ac_valtor
    """
    k, scalar = to_array(k)
    F = array_like(F, k, 1)

    A = 4*k*F
    a = np.zeros(k.shape)

    if method=='MomentumTheory':
        if chi!=0:
            raise Exception('Momentum Theory `a_k` only for chi=0')
        Ip = k>-1
        # TODO might need a switch based on "Vx" or phi
        a[Ip]  = k[ Ip]/(1+k[ Ip]) 
        a[~Ip] = k[~Ip]/(k[~Ip]-1) 

    elif method=='Buhl' or method=='AeroDyn15':
        if chi!=0:
            raise Exception('Buhl or AeroDyn15 `a_k` only for chi=0')

        InductionLimit = 1000000
        isValid=np.asarray([True]*len(k))
        for ik, kk in enumerate(k):
            #momentumRegion = (phi > 0.0_ReKi .and. Vx >= 0.0_ReKi) .or. (phi < 0.0_ReKi .and. Vx < 0.0_ReKi)
            momentumRegion=True # TODO TODO TODO 
            if momentumRegion:
                # Used in AeroDyn15
                if kk<=2/3: # momentum state for a < 0.4
                     if abs(kk - -1)<1e-6: #if ( EqualRealNos(k,-1.0_R8Ki) ) then
                         a[ik] = InductionLimit * np.sign(1+kk)
                     else:
                         a[ik] =  kk / (1+kk)
                     if (kk<-1):# k < -1 cannot be a solution in momentum region (this is equivalent to a>1.0)
                        isValid[ik] = False
                else: # Glauert(Buhl) correction for 0.4<a<1
                    T = A[ik]/2
                    FF = F[ik]
                    g1 = T - (10/9-FF)
                    g2 = T - ( 4/3-FF)*FF # Should always be >0
                    g3 = T - (25/9-2*FF)
                       
                    if abs(g3)<1e-6:
                        a[ik] =  1 - 1/(2*np.sqrt(g2))
                    else:
                        a[ik] = (g1 - np.sqrt(np.abs(g2)))/g3

                    #b = abs(g3) <1e-6
                    #a[b]  = 1 - 1/(2*np.sqrt(g2[b]))
                    #a[~b] = (g1[~b] - np.sqrt(np.abs(g2[~b])))/g3[~b]
            else:
                 if abs(kk -1)<1e-6: #if ( EqualRealNos(k,-1.0_R8Ki) ) then
                     a[ik] = InductionLimit
                     isValid[ik] = False
                 else:
                     a[ik] = kk/(kk-1)
                 # TODO I think should this be kk>=1
                 if kk<=1: #  k <= 1 cannot be a solution in propeller brake region (this is equivalent to a<1.0)
                     isValid[ik] = False
        a[~isValid]=np.nan


    elif method=='Thrust': # TODO rename

        for ik,kk in enumerate(k):
            a[ik], _ = axialInductionFromEmpiricalThrust(kk, chi, phi=phi, F=F[ik])

    elif method=='MomentumGlauertSkewRoots' or method=='MomentumGlauertSkewRootsWithHT':
        highThrustCorr =  method=='MomentumGlauertSkewRootsWithHT'
        cs=np.zeros((3,len(k)))
        roots = np.zeros((4,len(k)))
        for ik,kk in enumerate(k):
            a[ik],cs[:,ik], roots[:,ik] = axialInductionFromGlauertMomentum(kk, chi=chi, phi=phi, F=F[ik], highThrustCorr=highThrustCorr)
        if outputMore:
            return a, cs, roots
    else:
        raise NotImplementedError('High-thrust correction method for `a_k` function unknown:'+method)

    return a



def axialInductionFromGlauertMomentum(k, chi, phi=0.1, F=1, highThrustCorr=True):
    """ 
    Solve for `a` by equating blade element theory (BET) and momentum theory (MT) thrust
   
    At low loading, |k|<kc, Glauert's skew momentum theory is used:
   
           BET     =     MT
         (1-a)^2 k = a sqrt((1-a)^2 + tan(chi)^2)
   
    Which, when squared, leads to the fourth order polynomial:
   
        (1-k^2)a^4 + (4k^2-2)a^3 + (-6k^2 + tan^2\chi +1)a^2 + 4k^2 a - k^2 = 0 
   
    At high loading, |k|>kc, a hight thrust correction (2nd order polynomial) is used for "MT"

    Introduced to match OpenFAST interface

    INPUTS:
     - k: load factor [-]

    OUTPUTS:
     - a: axial induction ac_valtor, scalar [-]
     - c: c0, c1, c2

    """
    c=[0,0,0]
    ac, k0, kc, tan_chi0 = ak_lim(chi, method='Bladed', convention='aUn')
    if not highThrustCorr:
        kc = np.inf
    if abs(k)<=kc:
        c11 = tan_chi0**2
        c12 = k**2
        coeffs = (1-c12, 4*c12-2,  1+c11 -6*c12, 4*c12, -c12)
        roots = np.roots(coeffs) # Note: the roots are always real as far as I can tell
        #bReal = np.abs(np.imag(roots))<1e-10 
        #roots = np.sort(np.real(roots[bReal]))
        roots = np.sort(np.real(roots))
        if phi>=0:
            if roots[0]<0:
                # Will happen when k \in [0,1], we chose the solution of a in [0,1]
                a =  roots[1] # 
            else:
                a =  roots[0]
        else:
            a =  min( roots[0], roots[1])
    else:
        roots=None
        a, c = axialInductionFromEmpiricalThrust(k, chi=chi, phi=phi, F=F, momentumCorr=True)
    return a, c, roots

# 
def axialInductionFromEmpiricalThrust(k, chi, phi, F, momentumCorr=True, H=1):
    """ 

    Method where the induced velocity are factored out of the BET thrust
    and the HT is a second order polynomial

    CTs are:
          CT_BT = 4kF (1-a^2)   CT defined using Vxp
          CT_HT = c2 a^2 + c1 a  + c0 

    Equate them:
        (A-c2)a^2 - (2A +c1) a  + (1-c0) =0  A = 4kf

    Optional, square it:
      (A^2-c_2^2)a^4 + (-4A^2 - 2c_1 c_2) a^3 + (6A^2 - 2c_0 c_2 - c_1^2)a^2 + (-4A^2 -2c_0 c_1) a + (A^2-c_0^2) = 0

    Below are two methods, one where the equation above is solved for directly, one where it's squared first

    """
    # Coefficients for empirical CT = c2*a^2 + c1*a + c0 for a > a0
    if momentumCorr:
        [c0,c1,c2] = getEmpiricalCoefficients(chi=chi, F=F, method='BladedCorr')
    else:
        [c0,c1,c2] = getEmpiricalCoefficients(chi=chi, F=F, method='BladedNoCorr')
    A = 4*F*k

    # --- Solve for axial induction
    if not momentumCorr: 
        # Using un-squared version:
        #    (A-c2)a^2 - (2A +c1) a  + (1-c0) =0
        y1 = 2*A + c1
        y2 = 4*A*(c2+c1+c0) + c1*c1 - 4*c0*c2 
        y3 = 2*(A-c2)
        if abs(y3)<1e-16:
            a = 1 - 1/(2*np.sqrt(y2))
        else:
            if phi>=0:
                a = ( y1 - np.sqrt(y2) ) / y3
            else:
                a = ( y1 + np.sqrt(y2) ) / y3
        H=1.0 # TODO
#        if ((axInd>a0(chi0)).AND.(axInd<=1.0)) then
#           H = (4.0*axInd*(1.0-axInd)*F)/(c0+c1*axInd+c2*axInd*axInd)
#        elseif (axInd>1.0) then
#           H = (-4.0*axInd*(1.0-axInd)*F)/(c0+c1*axInd+c2*axInd*axInd)
#        else
#           H = 1.0
    else:
        # Using squared version:
        #    (A^2-c_2^2)a^4 + (-4A^2 - 2c_1 c_2) a^3 + (6A^2 - 2c_0 c_2 - c_1^2)a^2 + (-4A^2 -2c_0 c_1) a + (A^2-c_0^2) = 0
        A2 = A**2 
        coeffs = (A2 - c2*c2, -4*A2-2*c1*c2, 6*A2 -2*c0*c2 -c1*c1, -4*A2 - 2*c0*c1, A2 -c0*c0)
        roots = np.roots(coeffs)     # NOTE: those roots are often complex
        bReal = np.abs(np.imag(roots))<1e-10 # So we neglec the imaginary solutions below
        roots = np.sort(np.real(roots[bReal]))
        if phi>=0:
            if roots[0]<0:
                # Will happen when k \in [0,1], we chose the solution of a in [0,1]
                a =  roots[1] # 
            elif roots[1]<1:
                a =  roots[1]
            else:
                a =  roots[0]
        else:
            #a =  min( roots[1], roots[2])
            a =  roots[1]

#        if (phi >= 0.0) then
#            if (real(roots(0))<0.0_R8Ki) then
#                axInd = real(roots(1))
#            elseif (real(roots(1))<1.0_R8Ki) then
#                axInd = real(roots(1))
#            else
#                axInd = real(roots(0))
#            endif
#        else
#            axInd = real(roots(1))



#        tan_chi0 = min(MaxTanChi0, max(-MaxTanChi0, tan(chi0)))
#        if (equalrealnos(axInd,1.0_R8Ki)) then
#            H = 0       
#        elseif ((axInd>a0(chi)).AND.(axInd<=1.0)) then
#            H = 4.0_R8Ki*axInd*(1.0_R8Ki-axInd)*F*sqrt(1 + (tan_chi0/(1.0_R8Ki-axInd)*F)**2)/sqrt((c0+c1*axInd+c2*axInd*axInd)**2 + (4.0_R8Ki*axInd*tan_chi)**2)
#            ! Alternatively following implemention can be used but it keeps H from approaching zero as a -> 1
#            !H = (4.0_R8Ki*axInd*sqrt(((1.0_R8Ki-axInd)*F)**2 + tan(chi)**2))/sqrt((c0+c1*axInd+c2*axInd*axInd)**2 + (4.0_R8Ki*axInd*tan(chi))**2)           
#        elseif (axInd>1.0) then
#            H = -4.0_R8Ki*axInd*(1.0_R8Ki-axInd)*F*sqrt(1 + (tan_chi/(1.0_R8Ki-axInd)*F)**2)/sqrt((c0+c1*axInd+c2*axInd*axInd)**2 + (4.0_R8Ki*axInd*tan_chi0)**2)
#            ! Alternatively following implemention can be used but it keeps H from approaching zero as a -> 1
#            !H = -(4.0_R8Ki*axInd*sqrt(((1.0_R8Ki-axInd)*F)**2 + tan(chi0)**2))/sqrt((c0+c1*axInd+c2*axInd*axInd)**2 + (4.0_R8Ki*axInd*tan(chi))**2)
#        else
#            H = 1.0
#        if (k<0.0) then
#            H = 1.0
#
    return a, (c0,c1,c2)

def k_a(a, chi=0, F=None, method='Glauert'):
    """ Returns k = k(a) """
    a, scalar = to_array(a)
    F = array_like(F, a, 1)


    if method=='MomentumTheory':

        k= a/(1-a)

    elif method=='MomentumGlauertSkew':
        tan_chi0 = np.clip(np.tan(chi),  -MaxTanChi0, MaxTanChi0)
        k = a/(1-a)**2 * np.sqrt((1-a)**2 + tan_chi0**2)

    else:
        # Try a CT-a relation
        Ct = Ct_a(a, F=F, method=method, chi=chi)
        k = Ct/(4*a*F)
    return k


# --------------------------------------------------------------------------------}
# --- Helper functions 
# --------------------------------------------------------------------------------{
def to_array(a) :
    if not hasattr(a, '__len__'):
        scalar=True
        a = np.atleast_1d(a)
    else:
        scalar=False
        a = np.asarray(a)
    return a, scalar

def array_like(F, a, Fdefault):
    if F is None:
        F = Fdefault * np.ones(a.shape)
    else :
        if not hasattr(F, '__len__'):
            F = F * np.ones(a.shape)
        else:
            if len(F)!=len(a):
                raise Exception('F and a must have the same dimension')
    return F

def assertZero(x):
    if x!=0:
        raise Exception('Should be zero')


def polynomialOrder2Continuation(a_c, s_c, CT_c, CT_1):
    """ 
    return polynomial coeffs for a C1-continuous second order polynomial such that:
       CT = c0 + c1*a + c2*a2 
       CT(a_c)     = CT_c
       CT(1)       = CT_1
       dCT/da(a_c) = s_c
    """
    denom = (1-a_c)**2
    c0=(CT_1*a_c**2 - 2*CT_c*a_c + CT_c + a_c**2*s_c - a_c*s_c)/denom
    c1=(-2*CT_1*a_c + 2*CT_c*a_c - a_c**2*s_c + s_c)/denom
    c2=(CT_1 - CT_c + a_c*s_c - s_c)/denom
    return c0, c1, c2

def linearContinuation(a_c, s_c, CT_c):
    """ 
    return linear coeffs such that:
       CT = c1*a + c2*a2 
       CT(a_c)     = CT_c
       dCT/da(a_c) = s_c
    """
    c0 = CT_c - s_c*a_c 
    c1 = s_c
    return c0, c1




# --------------------------------------------------------------------------------}
# --- Test / debug 
# --------------------------------------------------------------------------------{
def main_plot_a_k():
    import matplotlib.pyplot as plt
    k = np.linspace(-2,2,50)

    fig,ax = plt.subplots(1,1)
    # Functions that depend on a only
    ax.plot(k, a_k(k, method='MomentumTheory'),'k-' ,label = 'Momentum theory')
    ax.plot(k, a_k(k, method='Buhl'),'--' ,label = 'Buhl')
    ax.set_xlim([-2,2])
    ax.set_ylim([-2,2])
    ax.legend()
    plt.show()

def main_plot():
    import matplotlib.pyplot as plt
    Ct=np.linspace(0,2,50)
    a =np.linspace(0,1,50)
    Ct_MT = 4*a*(1-a)

    fig,ax = plt.subplots(1,1)

    # Functions that depend on a only
    ax.plot(a ,Ct_MT,'k-' ,label = 'MomentumTheory'          )
#     ax.plot(a ,Ct_a(a,method='Glauert'),'-' ,label = 'Glauert (ac=1/3)')
#     ax.plot(a ,Ct_a(a,method='Spera')  ,'.' ,label = 'Spera (ac=0.3)')
    chi0=0
    ax.plot(a ,Ct_a(a,method='BladedCorr'  , chi=chi0),'o' ,label = 'BladedCorr')
    ax.plot(a ,Ct_a(a,method='BladedNoCorr', chi=chi0),'+' ,label = 'BladedNoCorr')
    chi0=30*np.pi/180
    ax.plot(a ,Ct_a(a,method='BladedCorr'  , chi=chi0),'o' ,label = 'BladedCorr')
    ax.plot(a ,Ct_a(a,method='BladedNoCorr', chi=chi0),'+' ,label = 'BladedNoCorr')
    #ax.plot(a ,Ct_a(a,method='Leishman')  ,'--' ,label = 'Leishman') # Buggy
    #ax.plot(a ,Ct_a(a,method='Buhl')  ,'+' ,label = 'Buhl') # Same as AeroDyn15 & 14

    # Inverting functions that depend on a to make them depend on Ct
    #a_Ct_Gl   = a_Ct_numInverse(lambda a: Ct_a(a, method='Glauert'))
    #a_Ct_Buhl = a_Ct_numInverse(lambda a: Ct_a(a, method='Buhl'))

#     # Functions that depend on Ct only
#     ax.plot(a_Ct(Ct,method = 'AeroDyn15'       ),Ct,'-' ,label = 'AeroDyn15'        )  # Same as Buhl
#     #ax.plot(a_Ct(Ct,method = 'AeroDyn14'       ),Ct,'d' ,label = 'AeroDyn14'        ) # Same as Buhl
#     #ax.plot(a_Ct_Buhl(Ct                       ),Ct,'d' ,label = 'Buhl'        )
#     ax.plot(a_Ct(Ct,method = 'HAWC2'           ),Ct,'--',label = 'HAWC2'            )
#     ax.plot(a_Ct(Ct,method = 'WEHandbook'      ),Ct,':' ,label = 'Handbook'         )
#     ax.plot(a_Ct(Ct,method = 'GlauertEmpirical'),Ct,'-.',label = 'Glauert Empirical')
    ax.set_xlabel('a [-]')
    ax.set_ylabel('Ct [-]')
    ax.set_xlim([0,1])
    ax.set_ylim([0,2])
    ax.legend()

    plt.show()  


def main_plot_coeffs():
    from welib.tools.colors import fColrs
    chi=np.linspace(0,50,100)*np.pi/180 # NOTE: might need abs

    C_corr   = np.zeros((3,len(chi)))
    C_noco = np.zeros((3,len(chi)))
    AC       = np.zeros(len(chi))

    for ic, chi in enumerate(chi):
        AC[ic] = ac_val(chi)
        C_corr[:,ic] = getEmpiricalCoefficients(chi=chi, F=1, method='BladedCorr')
        C_noco[:,ic] = getEmpiricalCoefficients(chi=chi, F=1, method='BladedNoCorr')

    chi*=180/np.pi 

    fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
    fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
    ax.plot(chi, C_corr[0,:],'-'  ,c=fColrs(1)   , label='c0 corr')
    ax.plot(chi, C_noco[0,:],'--' ,c=fColrs(1)   , label='c0 no corr')
    ax.plot(chi, C_corr[1,:],'-'  ,c=fColrs(2)   , label='c1 corr')
    ax.plot(chi, C_noco[1,:],'--' ,c=fColrs(2)   , label='c1 no corr')
    ax.plot(chi, C_corr[2,:],'-'  ,c=fColrs(3)   , label='c2 corr')
    ax.plot(chi, C_noco[2,:],'--' ,c=fColrs(3)   , label='c2 no corr')
    ax.plot(chi, AC    , label='a_c')
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.legend()


    # ---
    chi=np.linspace(0,50,5)*np.pi/180 # NOTE: might need abs
    a = np.linspace(0,1, 100)
    Ct = a*0

    fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
    fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)

    for ic, chi in enumerate(chi):
        Ct = a*0
        AC[ic] = ac_val(chi)
        c_corr = getEmpiricalCoefficients(chi=chi, F=1, method='BladedCorr')
        c_noco = getEmpiricalCoefficients(chi=chi, F=1, method='BladedNoCorr')

        Ct = Ct_a(a, F=None, method='BladedNoCorr', ac=None, chi=chi)
#         b=a>a0
#         aa= a[b]
#         Ct[b] = c[0]

        ax.plot(a, Ct, label='chi = {}'.format(chi*180/np.pi))
    ax.set_xlabel('a')
    ax.set_ylabel('Ct')
    ax.legend()
    plt.show()



def compareAD() :
    import pandas as pd
    Ct = np.linspace(0, 2, 500)
    a  = np.linspace(0, 1, 500)
    CtBuhl = Ct_a(a, method='Buhl')
    a_Ct_Buhl = a_Ct_numInverse(lambda a: Ct_a(a, method='Buhl'))
    a_AD   = a_Ct(Ct,method = 'AeroDyn15' )
    a_AD14 = a_Ct(Ct,method = 'AeroDyn14' )
    a_Buhl = a_Ct_Buhl(Ct )

    df = pd.DataFrame(data=np.column_stack((Ct,a_AD,a_AD14,a_Buhl)), columns=['Ct', 'a_AD','a_AD14','a_Buhl'])
    df.to_csv('AD_HighThrust.csv',index=False)

    fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
    fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
    ax.plot(Ct, a_AD , label='AD')
    ax.plot(Ct, a_AD14 , label='AD14')
    ax.plot(Ct, a_Buhl , label='Buhl')
    ax.set_xlabel('Ct')
    ax.set_ylabel('a')
    ax.legend()
    plt.show()


# # TODO TODO
# function [ a Ct] = fCorrectionHighThrust(CTCorrection,a, Cn ,phi, a_last,sigma, F,F, Ct )
#     switch CTCorrection
#             %%% Glauert correction
#         case 'Glauert'
#             ac=0.3;
#             Iac=a>ac;
#             A=sigma.*Cn./sind(phi)^2;
#             for i=1:length(a)
#                 a(i)=fzero(@(aa) -A(i)+aa*(4*F(i)+2*A(i))+aa.^2*(-5*F(i)-A(i))+3*F(i)*aa.^3    ,[0 1]);
#             end
#             %%% Glauert Exact correction
#         case 'GlauertExact'
#             ac=0.3;
#             error();
#             if a>ac
#                 A=sigma(e)*Cn/sind(phi)^2;
#                 asolutions=fGlauertSolutions(F,A);
#                 a=asolutions(whichmin(abs(asolutions-a_last)));
#             end
#             %%% Glauert correction REQUIRES RELAXATION
#         case 'Spera'
#             ac=0.34;
#             Iac=a>ac;
#             K=4*F(Iac).*(sind(phi(Iac))).^2./(sigma(Iac).*Cn(Iac));
#             a(Iac)=0.5*(2+K*(1-2*ac)-sqrt((K*(1-2*ac)+2 ).^2 + 4*(K*ac^2-1)    )  );
#             %Spera correction REQUIRES RELAXATION
# function [ as ] = fGlauertSolutions(F,A )
#     a1 = (-sqrt(3)*sqrt(-1)/2-1/2)*(sqrt(368*F^3+12*A*F^2+87*A^2*F-8*A^3) ...
#         /(2*3^(7/2)*F^(3/2)) ...
#         +(-290*F^3+15*A*F^2-24*A^2*F+2*A^3)/(1458*F^3)) ...
#         ^(1/3) ...
#         +(sqrt(3)*sqrt(-1)/2-1/2)*(-11*F^2-8*A*F+A^2) ...
#         /(81*F^2 ...
#         *(sqrt(368*F^3+12*A*F^2+87*A^2*F-8*A^3)/(2*3^(7/2)*F^(3/2)) ...
#         +(-290*F^3+15*A*F^2-24*A^2*F+2*A^3)/(1458*F^3)) ...
#         ^(1/3))+(5*F+A)/(9*F);
# 
#     a2 = (sqrt(3)*sqrt(-1)/2-1/2)*(sqrt(368*F^3+12*A*F^2+87*A^2*F-8*A^3) ...
#         /(2*3^(7/2)*F^(3/2)) ...
#         +(-290*F^3+15*A*F^2-24*A^2*F+2*A^3)/(1458*F^3)) ...
#         ^(1/3) ...
#         +(-sqrt(3)*sqrt(-1)/2-1/2)*(-11*F^2-8*A*F+A^2) ...
#         /(81*F^2 ...
#         *(sqrt(368*F^3+12*A*F^2+87*A^2*F-8*A^3)/(2*3^(7/2)*F^(3/2)) ...
#         +(-290*F^3+15*A*F^2-24*A^2*F+2*A^3)/(1458*F^3)) ...
#         ^(1/3))+(5*F+A)/(9*F);
#     a3 = (sqrt((368*F^3+12*A*F^2+87*A^2*F-8*A^3)/F)/(2*3^(7/2)*F) ...
#         +(-290*F^3+15*A*F^2-24*A^2*F+2*A^3)/(1458*F^3))^(1/3) ...
#         +(-11*F^2-8*A*F+A^2)/ ...
#         (81*F^2 *(sqrt((368*F^3+12*A*F^2+87*A^2*F-8*A^3)/F) ...
#         /(2*3^(7/2)*F) ...
#         +(-290*F^3+15*A*F^2-24*A^2*F+2*A^3) ...
#         /(1458*F^3)) ...
#         ^(1/3))+(5*F+A)/(9*F);
#     as=[real(a1) real(a2) real(a3)];
# 



if __name__ == '__main__':
    from  scipy.optimize import minimize_scalar
    import matplotlib.pyplot as plt
    main_plot_coeffs()
#     compareAD()
    main_plot() 
#     main_plot_a_k() 
