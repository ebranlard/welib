""" 
Set of high thrust corrections

NOTE: ac=0.4 -> Ctc = (1-(1-2*ac)**2) = 0.96

"""
import numpy as np


def a_Ct(Ct, F=None, method='AeroDyn'):
    """ 
    High thrust corrections of the form: a=a(Ct)
    INPUTS:
        Ct: Local thrust coefficient
        F : tip-loss factor
    """
    a=np.zeros(Ct.shape)
    if F is None:
        F=np.ones(Ct.shape)

    if method=='AeroDyn': # Very close to Glauert Empirical
        Ct[Ct>2]  = 2
        Ct[Ct<-2] = -2
        Ic        = Ct/F>0.96 # Correction
        In        = np.logical_not(Ic) # Normal
        a[Ic]     = 0.1432+np.sqrt(-0.55106+0.6427*Ct[Ic]/F[Ic])
        a[In]     = 0.5*(1-np.sqrt(1-Ct[In]/F[In]))
    elif method=='GlauertEmpirical':
        Ic    = Ct/F> 0.96  # Correction
        In    = np.logical_not(Ic) # Normal
        a[Ic] = 1/F[Ic]*(0.143+np.sqrt(0.0203-0.6427*(0.889-Ct[Ic])))
        a[In] = 0.5*(1-np.sqrt(1-Ct[In]/F[In]))
    elif method=='Handbook':
        Ic = Ct>0.96            # Correction
        In = np.logical_not(Ic) # Normal
        a[Ic]   = 1/F[Ic]*(0.143 + np.sqrt( 0.0203-0.6427 *(0.889-Ct[Ic] ) ))
        a[In] = 0.5*(1-np.sqrt(1-Ct[In]/F[In]))
    elif method=='HAWC2':
        #k = [-0.0017077,0.251163,0.0544955,0.0892074]
        k = [0.0       ,0.2460  ,0.0586,   0.0883]
        Ct = Ct/F
        a = k[3]*Ct**3+k[2]*Ct**2+k[1]*Ct+k[0]
    else:
        raise NotImplementedError('High correcton method '+method)
    return a

def a_Ct_a(Ct, a, F=None, method='Glauert_CTa'):
    """ 
    High thrust corrections of the form: a = a(Ct,a)
        see a_Ct
    """
    if F is None:
        F=np.ones(Ct.shape)
    # These need "a", TODO, invert them
    if method=='Glauert_CTa':
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
        raise NotImplementedError('High correcton method '+method)
    return a

def Ct_a(a, F=None, method='Glauert', ac=None):
    """ 
    High thrust corrections of the form: Ct = Ct(a)
        see a_Ct
    """
    if F is None:
        F=np.ones(a.shape)
    if method=='Glauert':
        Ct     = 4*a*F*(1-a)
        if ac is None:
            ac = 1/3
        Ic     = a>ac
        Ct[Ic] =4*a[Ic]*F[Ic]*(1-1/4*(5-3*a[Ic])*a[Ic])
    elif method=='Spera':
        Ct     = 4*a*F*(1-a)
        if ac is None:
            ac = 1/3
        Ic     = a>ac
        Ct[Ic] =4*F[Ic]*(ac**2 + (1-2*ac)*a[Ic])

    else:
        raise NotImplementedError('High correcton method '+method)
    return Ct

def main_plot():
    import matplotlib.pyplot as plt
    Ct=np.linspace(0,2,50)
    a =np.linspace(0,1,50)
    Ct_MT = 4*a*(1-a)

    fig,ax = plt.subplots(1,1)

    # Functions that depend on a only
    ax.plot(a ,Ct_MT,'k-' ,label = 'Momentum theory'          )
    ax.plot(a ,Ct_a(a,method='Glauert'),'-' ,label = 'Glauert (ac=1/3)')
    ax.plot(a ,Ct_a(a,method='Spera')  ,'.' ,label = 'Spera (ac=0.3)')
    # Functions that depend on Ct only
    ax.plot(a_Ct(Ct,method = 'AeroDyn'         ),Ct,'-' ,label = 'AeroDyn'          )
    ax.plot(a_Ct(Ct,method = 'HAWC2'           ),Ct,'--',label = 'HAWC2'            )
    ax.plot(a_Ct(Ct,method = 'Handbook'        ),Ct,':' ,label = 'Handbook'         )
    ax.plot(a_Ct(Ct,method = 'GlauertEmpirical'),Ct,'-.',label = 'Glauert Empirical')
    ax.set_xlabel('a [-]')
    ax.set_ylabel('Ct [-]')
    ax.set_xlim([0,1])
    ax.set_ylim([0,2])
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
    main_plot() 
