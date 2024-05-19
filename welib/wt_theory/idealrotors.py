""" 
Ideal actuator disk theory results

- Planforms for ideal rotors, See welib/BEM/examples/Example_IDealRotor
- 


References:
   [1]  Branlard (2017) Wind turbine aerodynamics and vorticity-based methods, Springer
   [2]  Manwell, et al. (2023) Wind Energy Explained, 3rd Edition, Wiley

"""
import numpy as np
from scipy.optimize import fsolve
from scipy.integrate import cumtrapz

from welib.BEM.idealrotors import *

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
                a[i] = fsolve(eq_a_lambda_optim, 0.33, args=(lambda_r))

    elif method =='analytical':
        a[bValid] = 1/2 * (1 - np.sqrt(1+lambda_[bValid]**2) * np.sin ( 1/3*np.arctan(1/lambda_[bValid]) )  ) # [1] Eq. 9.86
    else:
        raise NotImplementedError()

    # Limits:
    a[bZero] = 1/4

    ap[bValid] = (1-3*a[bValid])/(4*a[bValid]-1) # [1] Eq. 9.84
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
                #CT[i]=8/(l**2) * np.trapz(a_  *(1-a_) * lamb_   , lamb_) # [1] Eq. 9.80
                CP[i]=8/(l**2) * np.trapz(ap_ *(1-a_) * lamb_**3, lamb_)  # [1] Eq. 9.82
    return CP, a


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
# 
# 
#     TSR = 7.5
#     lam = np.linspace(0.001,TSR,100)
#     a, ap = ADMTO_inductions(lam, method='fzero')
#     a2,ap2 = ADMTO_inductions(lam, method='analytical')
#     print('Methods agree: ',np.all(np.abs(a-a2)<1e-13))
# 
#     fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,3.9)) # (6.4,4.8)
#     fig.subplots_adjust(left=0.11, right=0.97, top=0.85, bottom=0.15, hspace=0.10, wspace=0.13)
#     ax.plot(lam/TSR, a  , 'k-'   , label='a')
#     ax.plot(lam/TSR, ap , 'k--'  , label='a\'')
#     ax.set_xlabel('Non-dimensional blade radius, r/R')
#     ax.set_ylabel('Induction factor [-]')
#     ax.spines['right'].set_visible(False)
#     ax.spines['top'].set_visible(False)
#     ax.set_ylim([0,0.4])
#     ax.set_xlim([0,1.0])
#     ax.set_yticks([0,0.1,0.2,0.3,0.4])
#     ax.legend(frameon=False, loc='center right')
#     lw=0.9
#     for axis in ['top','bottom','left','right']:
#         ax.spines[axis].set_linewidth(lw)
#     ax.tick_params(length=6, width=lw)
#     plt.show()
# 

# % Used to generate figures 9.8 of [1]
# % Reference:
# %  [1]  Branlard 2017 Wind turbine aerodynamics and vorticity-based methods
# % 
# %% Initialization
# clear all; close all; clc; % addpath()
# restoredefaultpath;
# addpath(genpath('C:/Config/path/MatlabPath/libs/')) % http://github.com/ebranlard/matlab-path.git
# % I explain this a bit in my book section Maximum Power extraction 1D momentum+rotation
# R=100;
# nr=200;
# n2=10000; % for second method
# n3=100000; % for third method
# 
# 
# 
# options = optimset('TolX',1e-11);
# %% Overview of the function whose zeroes are sought
# % figure, hold all
# % for x=vlambdaRef
# %     va=linspace(0.25,1/3,100);
# %     plot(va,16*va.^3-24*va.^2+va.*(9-3*x^2)-1+x^2)
# % end
# % plot(va,va*0,'k')
# % ylim([-0.01 0.01])
# % 
# 
# %% That's one way of doing it
# tic()
# disp('First method...')
# vlambda=unique(sort([vlambdaRef 10.^-(9:-1:1) 0.1:0.5:1 1:15]));
# vr=linspace(0,R,nr);
# a=zeros(1,nr);
# CP=zeros(1,length(vlambda));
# CT=zeros(1,length(vlambda));
# for il=1:length(vlambda)    
#     lambda=vlambda(il);
#     lambda_r=vr/R*lambda;  
#     for e = 1:nr
#         x=lambda_r(e);
#         a(e)=fzero(@(a) 16*a^3-24*a^2+a*(9-3*x^2)-1+x^2,0.3,options); %martin 8.1 combined from 4.32 and 4.38
#     end
#     aprime=(1-3*a)./(4*a-1);
#     CP(il)=8/lambda^2*trapz(lambda_r,aprime.*(1-a).*lambda_r.^3);    % martin 4.30
#     CT(il)=8/lambda^2*trapz(lambda_r,a.*(1-a).*lambda_r);    % See my book simplified 2D with rotation
# end
# toc()
# 
# %% Another way of doing it
# tic()
# disp('Second method...')
# vlambda2=unique(sort([vlambdaRef linspace(0,max(vlambda),n2)]));
# nn=length(vlambda2);
# a2=zeros(1,nn);
# aprime2=zeros(1,nn);
# for il=1:nn
#     x=vlambda2(il);
#     a2(il)=fzero(@(a) 16*a^3-24*a^2+a*(9-3*x^2)-1+x^2,0.3,options);  %martin 8.1 combined from 4.32 and 4.38
# end
# aprime2=(1-3*a2)./(4*a2-1);
# toc()
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
# 
# figure,hold all, 
# plot(vlambda2,aprime2.*(1-a2).*vlambda2.^3)
# plot(vlambda3,aprime3(Isort).*(1-a3(Isort)).*vlambda3.^3)
# 
# % CP
# figure,hold on,grid,box
# plot(vlambda,vlambda*0+16/27 ,'k--')
# plot(vlambda2,CP2 ,'k-')
# % plot(vlambda,CP,'b+')
# % plot(vlambda3,CP3,'r--')
# % plot(vlambdaRef,CPRef,'ko')
# legend('Betz limit','Ideal rotor with wake rotation')
# % ,'Same but different way of thinking it','Same but deterministic','Wilson Lissaman',0)
# xlim([0 10])
# ylim([0 0.62])
# xlabel('Tip speed ratio \lambda [-]')
# ylabel('C_P [-]')
# title('MomentumTheoryActuatorDiskOptimalCP')
# 
# 
# 
# % CT
# figure,hold on,grid,box
# plot(vlambda,vlambda*0+8/9 ,'k--')
# plot(vlambda2,CT2 ,'k-')
# % plot(vlambda,CT,'b+')
# % plot(vlambda3,CT3,'r--')
# legend('Betz limit','Ideal rotor with wake rotation')
# % ,'Same but different way of thinking it','Same but deterministic',0)
# xlim([0 10])
# xlabel('Tip speed ratio \lambda [-]')
# ylabel('C_T [.]')
# title('MomentumTheoryActuatorDiskOptimalCT')
# 
# % a
# figure,hold on,grid,box
# plot(vlambda,vlambda*0+1/3 ,'k--')
# plot(vlambda2,a2 ,'k-')
# % plot(vlambda3,a3,'r--')
# legend('Betz limit','Ideal rotor with wake rotation')
# % ,'Same but deterministic',0)
# xlim([0 10])
# ylim([0 0.35])
# xlabel('Local tip speed ratio \lambda_r [-]')
# ylabel('Axial induction factor a [-]')
# title('MomentumTheoryActuatorDiskOptimala')
# 
# % aprime
# figure,hold on,grid,box
# plot(vlambda2,aprime2 ,'k-')
# % plot(vlambda3,aprime3,'r--')
# ylim([0 0.1])
# legend('Ideal rotor with wake rotation')
# % ,'Same but deterministic',0)
# xlim([0 10])
# title('MomentumTheoryActuatorDiskOptimalaprime')
# xlabel('Local tip speed ratio \lambda_r [-]')
# ylabel('Tangential induction factor a'' [-]')
# 
# 
# %% 
# I =whichvalue(vlambda2,vlambdaRef);
# I2=whichvalue(vlambda,vlambdaRef);
# I3=whichvalue(vlambda3,vlambdaRef);
# 
# [vlambda3(I3)' CP2(I)' CP(I2)' CP3(I3)' CPRef' CPMartin']
# 
# fprintf('\n\nlambda\t CP2\t CPref CPMartin\n');
# fprintf('%.3f\t %.4f\t %.4f %.3f\n', [vlambda2(I)' CP2(I)' CPRef' CPMartin']');
# %% Cp
# vlambdaAD=linspace(0,30,151);
# vaAD=interp1(vlambda3(1:end-2),a3(1:end-2),vlambdaAD);
# vaprimeAD=interp1(vlambda3(1:end-2),aprime3(1:end-2),vlambdaAD);
# vCPAD=interp1(vlambda3(1:end-2),CP3(1:end-2),vlambdaAD);
# vCTAD=interp1(vlambda3(1:end-2),CT3(1:end-2),vlambdaAD);
# vCPAD(1)=0;
# vCTAD(1)=0;
# vaAD(1)=0.25;
# 
# 
# 
# % save([./MomentumTheoryActuatorDisk.mat'],'vlambdaAD','vCPAD','vCTAD','vaAD','vaprimeAD');
# 
# %%
# function [a aprime phi c beta]=getOptimizedParameters(Profile,alpha_d,lambda,B,R,r,x)
#     Cd_opt=Profile.Cd(whichvalue(Profile.alpha,alpha_d));
#     Cl_opt=Profile.Cl(whichvalue(Profile.alpha,alpha_d));
#     a=zeros(1,length(x));
#     for i=1:length(x)
#         a(i)=fzero(@(aa) getRelationAX(aa,x(i)),0.3);
#     end    
#     aprime=(1-3*a)./(4*a-1);
#     phi=atan( (1-a)./((1+aprime).*x) );  %rad
#     Cn=Cl_opt*cos(phi)+Cd_opt*sin(phi);
#     f=B/2*(R-r)./(r.*sin(phi));
#     F=2/pi*acos(exp(-f));
#     c=(8*pi*R.*F.*a.*x.*(sin(phi)).^2)./((1-a)*B*lambda.*Cn);
#     beta= ((phi*180/pi)-alpha_d);
# end

