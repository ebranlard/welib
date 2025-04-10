""" 
Ideal actuator disk theory results.

TODO merge with:

- Planforms for ideal rotors
       See welib.BEM.idealrotors
       See welib.BEM.examples/Example_IDealRotor
- 


References:
   [1]  Branlard (2017) Wind turbine aerodynamics and vorticity-based methods, Springer
   [2]  Manwell, et al. (2023) Wind Energy Explained, 3rd Edition, Wiley

"""
import numpy as np
from scipy.optimize import fsolve
from welib.BEM.idealrotors import *
try:
    from numpy import trapezoid
except:
    from numpy import trapz as trapezoid

# --------------------------------------------------------------------------------}
# --- Maximum power extraction from stream tube theory (STT)
# --------------------------------------------------------------------------------{
def ADMTO_inductions(lambda_, method='fzero'):
    """
    Actuator disc momentum theory optimal (ADMTO) induction factors
    See e.g. [1] Section 9.5.4
    """
    def eq_a_lambda_optim(a, lambda_r):
         return 16*a**3 - 24*a**2 + a*(9-3*lambda_r**2) - 1 + lambda_r**2 # [1] Eq. 9.85

    lambda_ = np.asarray(lambda_)
    a       = np.zeros_like(lambda_)
    ap      = np.zeros_like(lambda_)

    bValid = lambda_>0
    bZero  = lambda_==0

    if method=='fzero':
        for i,lambda_r in enumerate(lambda_):
            if lambda_r==0:
                a[i] = 1/4
            else:
                a[i] = fsolve(eq_a_lambda_optim, 0.33, args=(lambda_r))[0]

    elif method =='analytical':
        a[bValid] = 1/2 * (1 - np.sqrt(1+lambda_[bValid]**2) * np.sin ( 1/3*np.arctan(1/lambda_[bValid]) )  ) # [1] Eq. 9.86
    elif method =='analytical2':
        phi   = 2./3. * np.arctan(1./(lambda_[bValid]))
        a[bValid] = 1/(1+ np.sin(phi)**2/(  (1-np.cos(phi))*np.cos(phi)   ) )

    else:
        raise NotImplementedError()

    # Limits:
    a[bZero] = 1/4
    ap[bValid] = (1-3*a[bValid])/(4*a[bValid]-1) # [1] Eq. 9.84, only true for optimal
    #ap= 1/2*(-1+np.sqrt(1+4/lambda_r**2 * a * (1-a)))    # Always true 
    ap[bZero]  = np.nan


    return a, ap

def ADMTO_CP(lambda_, method='analytical'):
    """
    Actuator disc momentum theory optimal (ADMTO) induction factors
    See e.g. [2] Chap 3 
    """
    lambda_ = np.asarray(lambda_)
    CP      = np.zeros_like(lambda_)
    a, ap = ADMTO_inductions(lambda_, method='analytical')

    if method =='analytical':
        def dCPmax(x): # [2] Chapter 3, eq 43
            return  (64/5*x**5 + 72*x**4 + 124*x**3 + 38*x**2 - 63*x-12*np.log(x) - 4/x)
        for i, (l,a2) in enumerate(zip(lambda_,a)):
            CP[i] = 8/(729*l**2) * ( dCPmax(x=1/4) - dCPmax(x=1-3*a2)  )
    elif method == 'num_int':
        n=100000
        for i, l in enumerate(lambda_):
            lamb_ = np.linspace(0.00001, l , n)
            a_, ap_ = ADMTO_inductions(lamb_, method='analytical')
            if l>0:
                #CT[i]=8/(l**2) * trapezoid(a_  *(1-a_) * lamb_   , lamb_) # [1] Eq. 9.80
                CP[i]=8/(l**2) * trapezoid(ap_ *(1-a_) * lamb_**3, lamb_)  # [1] Eq. 9.82
    elif method == 'Wilson':
        # Results given by WilsonLissaman 1974 p57
        lambda_=np.array([0.5,      1, 1.5 ,   2 ,  2.5,  5  ,7.5  ,10])  #
        a, ap = ADMTO_inductions(lambda_, method='analytical')
        CP  =np.array([0.288,0.416,0.480,0.512,0.532,0.570,0.582,0.593])
    elif method == 'Hansen':
        # Results found in Hansen p39
        lambda_=np.array([0.5,      1, 1.5 ,   2 ,  2.5,  5  ,7.5  ,10])  #
        a, ap = ADMTO_inductions(lambda_, method='analytical')
        CP  =np.array([0.486,0.703,0.811,0.865,0.899,0.963,0.983,0.987])*16/27
    return CP, a



def CompareTheories():
    lambda_=np.array([0.5,      1, 1.5 ,   2 ,  2.5,  5  ,7.5  ,10])  #
    CP, a    = ADMTO_CP(vlambdaRef, method='analytical')
    CPnum, _ = ADMTO_CP(vlambdaRef, method='num_int')


    print(np.abs(CP-CPnum))
    print('CP Methods agree: ',np.all(np.abs(CP-CPnum)<1e-10))

    for l,cp,a2,cpWL,cpM,cpnum in zip(vlambdaRef, CP, a, CPWilson, CPMartin, CPnum):
        print('{:4.1f}\t{:6.4f}\t{:7.5f}\t{:7.5f}\t{:6.3f}\t{:6.3f}'.format(l,a2,cp,cpnum,cpWL, cpM))
    print(a)
    print(CPnum)



if __name__ == '__main__':
    

    from welib.tools.figure import *
    setFigureFont(14)



    vlambdaRef=np.array([0.5,      1, 1.5 ,   2 ,  2.5,  5  ,7.5  ,10])  #
    CPWilson  =np.array([0.288,0.416,0.480,0.512,0.532,0.570,0.582,0.593])#   WilsonLissaman 1974 p57
    CPMartin  =np.array([0.486,0.703,0.811,0.865,0.899,0.963,0.983,0.987])*16/27# Hansen p39


    CP, a    = ADMTO_CP(vlambdaRef)
    CPnum, _ = ADMTO_CP(vlambdaRef, method='num_int')


    print(np.abs(CP-CPnum))
    print('CP Methods agree: ',np.all(np.abs(CP-CPnum)<1e-10))

    for l,cp,a2,cpWL,cpM,cpnum in zip(vlambdaRef, CP, a, CPWilson, CPMartin, CPnum):
        print('{:4.1f}\t{:6.4f}\t{:7.5f}\t{:7.5f}\t{:6.3f}\t{:6.3f}'.format(l,a2,cp,cpnum,cpWL, cpM))
    print(a)
    print(CPnum)


    lambda_ = np.linspace(0.01,10,1000)
    CP, a    = ADMTO_CP(lambda_)

    fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,3.6)) # (6.4,4.8)
    fig.subplots_adjust(left=0.12, right=0.97, top=0.96, bottom=0.16, hspace=0.10, wspace=0.13)
    ax.plot(lambda_, 16/27+lambda_*0, 'k-'   , label='Betz limit (no wake rotation)')
    ax.plot(lambda_, CP             , 'k--'  , label='Including wake rotation')
    ax.set_xlabel(r'Tip speed ratio, $\lambda$ [-]')
    ax.set_ylabel(r'Power coefficient, $C_P$ [-]')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_ylim([0,0.6])
    ax.set_xlim([0,10.0])
    ax.set_yticks([0,0.1,0.2,0.3,0.4,0.5,0.6])

    bbox2 = dict(boxstyle='square,pad=-0.1', fc='w', ec='none')
    ax.annotate( r'$C_P=\frac{16}{27}$' , xy=(0.78, 16/27+0.005), xycoords='data', ha="center", va="center", xytext=(0.78, 0.52), textcoords='data', arrowprops=dict(arrowstyle="->"),bbox=bbox2, fontsize=14)

    ax.legend(frameon=False, loc='center right')
    lw=0.9
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(lw)
    ax.tick_params(length=6, width=lw)
    plt.show()



# 
# %% and yet another way of doing it, the deterministic way
# tic()
# disp('Third method...')
# a3=linspace(0.25+10^-9,1/3,n3); % why do i start at 0.25, because of the expression for aprime
# aprime3=(1-3*a3)./(4*a3-1);
# vlambda3=sqrt(a3.*(1-a3)./(aprime3.*(1+aprime3)));
# [vlambda3 Isort]=sort(vlambda3);% in fact not needed...
# a3=a3(Isort);
# aprime3=aprime3(Isort);
# toc()
# 
