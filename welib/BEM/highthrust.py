""" 
Set of high thrust corrections

NOTE: ac=0.4 -> Ctc = (1-(1-2*ac)**2) = 0.96

"""
import numpy as np


def a_Ct(Ct, a=None, F=None, method='AeroDyn'):
    """ 
    High thrust corrections of the form: 
       a=a(Ct)
   or 
       a=a(Ct,a)
    INPUTS:
        Ct: Local thrust coefficient
        a: induction factor
        F : tip-loss factor
    """
    if F is None:
        F=np.ones(Ct.shape)
    # -------------------------------------------------------
    # --- a = a(Ct)
    # -------------------------------------------------------
    if method=='AeroDyn15': # Very close to Glauert Empirical, Very close to Buhl
        Ct[Ct>2]  = 2
        Ct[Ct<-2] = -2
        Ic        = Ct/F>0.96 # Correction
        In        = np.logical_not(Ic) # Normal
        a = np.zeros(Ct.shape)
        a[Ic]     = 0.1432+np.sqrt(-0.55106+0.6427*Ct[Ic]/F[Ic])
        a[In]     = 0.5*(1-np.sqrt(1-Ct[In]/F[In]))
    elif method=='AeroDyn14':
        Ic        = Ct/F>0.96 # Correction
        In        = np.logical_not(Ic) # Normal
        a = np.zeros(Ct.shape)
        a[Ic]     = (18*F[Ic]-20-3*np.sqrt(Ct[Ic]*(50-36*F[Ic]) + 12*F[Ic]*(3*F[Ic]-4) ) ) / (36*F[Ic] - 50)
        a[In]     = 0.5*(1-np.sqrt(1-Ct[In]/F[In])) # Baseline

    elif method=='GlauertEmpirical':
        Ic    = Ct/F> 0.96  # Correction
        In    = np.logical_not(Ic) # Normal
        a=np.zeros(Ct.shape)
        a[Ic] = 1/F[Ic]*(0.143+np.sqrt(0.0203-0.6427*(0.889-Ct[Ic])))
        a[In] = 0.5*(1-np.sqrt(1-Ct[In]/F[In]))
    elif method=='WEHandbook':
        Ic = Ct>0.96            # Correction
        In = np.logical_not(Ic) # Normal
        a=np.zeros(Ct.shape)
        a[Ic]   = 1/F[Ic]*(0.143 + np.sqrt( 0.0203-0.6427 *(0.889-Ct[Ic] ) ))
        a[In] = 0.5*(1-np.sqrt(1-Ct[In]/F[In]))
    elif method=='HAWC2':
        #k = [-0.0017077,0.251163,0.0544955,0.0892074]
        k = [0.0       ,0.2460  ,0.0586,   0.0883]
        Ct = Ct/F
        a = k[3]*Ct**3+k[2]*Ct**2+k[1]*Ct+k[0]
    # -------------------------------------------------------
    # --- a = a(Ct,a)
    # -------------------------------------------------------
    elif method=='Glauert_CTa':
       ac=0.3;
       #Ic =Ct/F> 1-(1-(2*ac))**2  # Correction
       Ic = a>ac                  # Correction
       fg=0.25*(5-3*a[Ic])
       a[Ic] = Ct[Ic]/(4*F[Ic]*(1-fg*a[Ic]))
       #In = np.logical_not(Ic)    # Normal
       #a[In] = 0.5*(1-np.sqrt(1-Ct[In]/F[In]))
    elif method=='Spera_CTa':
        ac    = 0.34
        Ic    = a>ac
        fgs   = ac/a[Ic]*(2-ac/a[Ic])
        a[Ic] = Ct[Ic]/(4*F[Ic]*(1-fgs*a[Iac]))
        #In = np.logical_not(Ic) # Normal
    elif method=='Shen_CTa':
        ac    = 1/3
        Ic    = a>ac
        a[Ic] = (Ct[Ic]/4*F[Ic]-F[Ic]*ac**2)/(1-2*ac*F[Ic]);
    else:
        raise NotImplementedError('High-thrust correction method for `a-Ct` function unknown:'+method)
    return a


def Ct_a(a, F=None, method='Glauert', ac=None, chi0=0):
    """ 
    High thrust corrections of the form: Ct = Ct(a)
        see a_Ct
    """
    if not hasattr(a, '__len__'):
        scalar=True
        a = np.atleast_1d(a)
    else:
        scalar=False
        a = np.asarray(a)
    if F is None:
        F=np.ones(a.shape)
    Ct     = 4*a*F*(1-a)
    if method=='Glauert':
        if ac is None:
            ac = 1/3
        Ic     = a>ac
        Ct[Ic] =4*a[Ic]*F[Ic]*(1-1/4*(5-3*a[Ic])*a[Ic])
    elif method=='Spera':
        if ac is None:
            ac = 1/3
        Ic     = a>ac
        Ct[Ic] =4*F[Ic]*(ac**2 + (1-2*ac)*a[Ic])

    elif method=='Leishman':
        # NOTE: doesn't seem to work
        c0 = -3.402; c1 = 6.489; c2 = -3.167; c3 = 0.695; c4 = -0.0562
        ac = 1/3 # <<< TODO
        Ic     = a>ac
        Ct[Ic] =F[Ic]*(c0 + c1*a[Ic] + c2*a[Ic]**2 + c3*a[Ic]**3 + c4*a[Ic]**4)

    elif method=='Buhl':
        ac = 0.4
        Ic = a>ac
        c0 = 8/9;
        c1 = 4*F[Ic]-40/9;
        c2 = 50/9-4*F[Ic];
        Ct[Ic] =c0 + c1*a[Ic] + c2*a[Ic]**2


    elif method=='BladedNoCorr':
        a0_local = a0(chi0)
        for ia, (aa,FF) in enumerate(zip(a,F)):
            if aa<a0_local:
                Ct[ia] = 4*aa*FF*(1-aa)
            else:
                c = getEmpiricalCoefficients(chi0=chi0, F=FF, method=method)
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

    elif method=='BladedCorr':
        #c = getEmpiricalCoefficients(chi0=chi0, F=F, method=method)
        #Ct =c[0] + c[1]*a + c[2]*a**2

        a0_local = a0(chi0)
        for ia, (aa,FF) in enumerate(zip(a,F)):
            if aa<a0_local:
                Ct[ia] = 4*aa*FF*(1-aa)
            else:
                c = getEmpiricalCoefficients(chi0=chi0, F=FF, method=method)
                Ct[ia] =c[0] + c[1]*aa + c[2]*aa**2


    else:
        raise NotImplementedError('High-thrust correction method for `Ct_a` function unknown:'+method)
    if scalar:
        Ct = Ct[0]
    return Ct

def a0(chi):
   return np.clip(0.5*np.cos(45*np.pi/180)/np.cos(chi), 0, 0.5)


def getEmpiricalCoefficients(chi0=0, F=1, method='BladedCorr'):
    """ Return polynomial coefficients 
        CT = 4*a*(1-a)*F = c0 + c1*a + c2*a^2 + ...
    """
    if hasattr(F, '__len__'):
        raise Exception('F should be scalar')
    if hasattr(chi0, '__len__'):
        raise Exception('chi0 should be scalar')

    MaxTanChi0 = 100  # maximum absolute value allowed for tan(chi0), an arbitary large number
    if method=='BladedCorr':
        # Empirical CT = 4*a*(1-a)*F = c2*a^2 + c1*a + c0 for a > a0
        a0_local = a0(chi0)
        denom = (a0_local**2 - 2*a0_local + 1)
        temp2 = np.clip(np.tan(chi0), -MaxTanChi0, MaxTanChi0)**2
        temp1 = np.sqrt((a0_local-1)**2 +temp2)
        
        CTata1 = np.sqrt(((-0.64755/(np.cos(chi0)*np.cos(chi0)) - 0.8509/np.cos(chi0) + 3.4984)*F)**2 + 16*temp2)
        CTata1 = max(1, CTata1)
        
        c2 = (CTata1 - 4*F/temp1 + 16*F*a0_local/temp1 - 4*F*a0_local*temp1 - 4*temp2*F/temp1 - 20*F*(a0_local**2)/temp1 + 8*F*(a0_local**3)/temp1 + 4*temp2*F*a0_local/temp1 ) /denom
        c1 = 2*( 2*F/temp1 - a0_local*CTata1 - 6*F*a0_local/temp1 + 2*temp2*F/temp1 + 2*F*(a0_local**2)/temp1 + 4*F*(a0_local**2)*temp1 + 6*F*(a0_local**3)/temp1 - 4*F*(a0_local**4)/temp1 - 2*temp2*F*(a0_local**2)/temp1 )/denom
        c0 = a0_local*( a0_local*CTata1 - 4*F/temp1 + 4*F*temp1 + 16*F*a0_local/temp1 - 8*F*a0_local*temp1 - 4*temp2*F/temp1 - 20*F*(a0_local**2)/temp1 + 8*F*(a0_local**3)/temp1 + 4*temp2*F*a0_local/temp1 )/denom
        return [c0, c1, c2]
     

    elif method=='BladedNoCorr':
        # Empirical CT = 4*a*(1-a)*F = c2*a^2 + c1*a + c0 for a > a0
        a0_local = a0(chi0)
        denom = (a0_local**2 - 2*a0_local + 1)

        CTata1 = (-0.64755/(np.cos(chi0)*np.cos(chi0)) - 0.8509/np.cos(chi0) + 3.4984)*F      
        CTata1 = max( 1, CTata1 )       

        c2 =  (-4*F*a0_local**2 + 8*F*a0_local - 4*F + CTata1)/denom    
        c1 = 2*(2*F*a0_local**2 - CTata1*a0_local - 4*F*a0_local  + 2*F)/denom    
        c0 = CTata1*(a0_local**2)/denom

        return [c0, c1, c2]
    else:
        raise NotImplementedError('High-thrust polynomial coefficients method unknown: '+method)

   



def a_Ct_numInverse(fCta):
    """ take a Ct(a) function and make it a a(Ct) function"""
    from  scipy.optimize import minimize_scalar
    def a_Ct_local(Ct): 
        if not hasattr(Ct, '__len__'):
            scalar=True
            Ct = np.atleast_1d(Ct)
        else:
            scalar=False
        Ct = np.asarray(Ct)
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


def a_k(k, F=None, method='Buhl'):
    """
    Return a = a(k) 
    """

    if not hasattr(k, '__len__'):
        scalar=True
        k = np.atleast_1d(k)
    else:
        scalar=False
    if F is None:
        F=np.ones(k.shape)

    A = 4*k*F
    a = np.zeros(k.shape)

    if method=='Momentum theory':
        Ip = k>-1
        # TODO might need a switch based on "Vx" or phi
        a[Ip]  = k[ Ip]/(1+k[ Ip]) 
        a[~Ip] = k[~Ip]/(k[~Ip]-1) 

    elif method=='Buhl' or method=='AeroDyn15':

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



    elif method=='Thrust':
        c = getEmpiricalCoefficients(chi0=chi0, F=F, method=method)

#        y1 = 2.0*A + c1
#        y2 = 4.0*A*(c2+c1+c0) + c1*c1 - 4.0*c0*c2 
#        y3 = 2.0*(A-c2)
#        if ( EqualRealNos( y3, 0.0_R8Ki ) ) then
#           axInd = 1.0 - 1.0/(2.0*SQRT(y2))
#        else
#           if (phi>=0.0) then
#              axInd = ( y1 - SQRT(y2) ) / y3
#           else
#              axInd = ( y1 + SQRT(y2) ) / y3
#           end if
#        end if

    else:
        raise NotImplementedError('High-thrust correction method for `a_k` function unknown:'+method)

    return a

def main_plot_a_k():
    import matplotlib.pyplot as plt
    k = np.linspace(-2,2,50)

    fig,ax = plt.subplots(1,1)
    # Functions that depend on a only
    ax.plot(k, a_k(k, method='Momentum theory'),'k-' ,label = 'Momentum theory')
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
    ax.plot(a ,Ct_MT,'k-' ,label = 'Momentum theory'          )
#     ax.plot(a ,Ct_a(a,method='Glauert'),'-' ,label = 'Glauert (ac=1/3)')
#     ax.plot(a ,Ct_a(a,method='Spera')  ,'.' ,label = 'Spera (ac=0.3)')
    chi0=0
    ax.plot(a ,Ct_a(a,method='BladedCorr'  , chi0=chi0),'o' ,label = 'BladedCorr')
    ax.plot(a ,Ct_a(a,method='BladedNoCorr', chi0=chi0),'+' ,label = 'BladedNoCorr')
    chi0=30*np.pi/180
    ax.plot(a ,Ct_a(a,method='BladedCorr'  , chi0=chi0),'o' ,label = 'BladedCorr')
    ax.plot(a ,Ct_a(a,method='BladedNoCorr', chi0=chi0),'+' ,label = 'BladedNoCorr')
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
    A0       = np.zeros(len(chi))

    for ic, chi0 in enumerate(chi):
        A0[ic] = a0(chi0)
        C_corr[:,ic] = getEmpiricalCoefficients(chi0=chi0, F=1, method='BladedCorr')
        C_noco[:,ic] = getEmpiricalCoefficients(chi0=chi0, F=1, method='BladedNoCorr')

    chi*=180/np.pi 

    fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
    fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
    ax.plot(chi, C_corr[0,:],'-'  ,c=fColrs(1)   , label='c0 corr')
    ax.plot(chi, C_noco[0,:],'--' ,c=fColrs(1)   , label='c0 no corr')
    ax.plot(chi, C_corr[1,:],'-'  ,c=fColrs(2)   , label='c1 corr')
    ax.plot(chi, C_noco[1,:],'--' ,c=fColrs(2)   , label='c1 no corr')
    ax.plot(chi, C_corr[2,:],'-'  ,c=fColrs(3)   , label='c2 corr')
    ax.plot(chi, C_noco[2,:],'--' ,c=fColrs(3)   , label='c2 no corr')
    ax.plot(chi, A0    , label='A0')
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.legend()


    # ---
    chi=np.linspace(0,50,5)*np.pi/180 # NOTE: might need abs
    a = np.linspace(0,1, 100)
    Ct = a*0

    fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
    fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)

    for ic, chi0 in enumerate(chi):
        Ct = a*0
        A0[ic] = a0(chi0)
        c_corr = getEmpiricalCoefficients(chi0=chi0, F=1, method='BladedCorr')
        c_noco = getEmpiricalCoefficients(chi0=chi0, F=1, method='BladedNoCorr')

        Ct = Ct_a(a, F=None, method='BladedNoCorr', ac=None, chi0=chi0)
#         b=a>a0
#         aa= a[b]
#         Ct[b] = c[0]

        ax.plot(a, Ct, label='chi = {}'.format(chi0*180/np.pi))
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
