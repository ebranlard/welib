import numpy as np
import matplotlib.pyplot as plt
# Local 

from vortilib.elements.VortexHelix import *


def Params(a=1/3,U0=10,Lambda=6,R=100,nB=3):
    h         = 2*np.pi*R*(1-a)/Lambda
    CT        = 4*a*(1-a)
    Gamma     = CT*np.pi*R*U0/(nB*Lambda) # NOTE assume aprime = 0, large tip-sped ratio 
    return  a,U0,R,Lambda,nB,h,CT,Gamma

def main():
    import matplotlib.pyplot as plt
#     import warnings
#     warnings.filterwarnings('error')

    # --------------------------------------------------------------------------------}
    # --- Radial survey on the lifting line 
    # --------------------------------------------------------------------------------{
    a,U0,R,Lambda,nB,h,CT,Gamma_B=Params()
    Gamma_tot = Gamma_B*nB
    psi_blade = 0*np.pi/2
    method    = 'wrench'
    bWT       = True # Convention WT or Propeller
    bSemi     = True # infinite helix or at the rotor
    vr        = np.linspace(0.1,2,100)*R
    Xcp       = vr
    Ycp       = Xcp*0
    Zcp       = Xcp*0
    Omega = U0*Lambda/R

    # --- Axial induction on the lifting line 
    fig,ax = plt.subplots(1,1)
    ax.plot([0,1],[-a,-a],'k--',label='Momentum T. (B=infty)',lw=3)
    ax.plot([1,2],[ 0, 0],'k--',lw=3)

    for nB in [1, 3, 10, 500]:
        ur,ut,uz = vh_u(Xcp,Ycp,Zcp,Gamma_tot/nB,R,h,psih=psi_blade,nB=nB,bWT=bWT,method=method,bSemi=bSemi)
        ax.plot(Xcp/R,uz/U0, label='B = {:d}'.format(nB))

    ax.set_xlabel('r/R [-]')
    ax.set_ylabel('uz/U0 [-]')
    ax.set_ylim([-2*a, 1/2*a])
    ax.set_title('Axial induced velocity on the lifting line')
    ax.legend()

    # --- Tangential induction on the lifting line 

    fig,ax = plt.subplots(1,1)
    #ax.plot([0,1],[-a,-a],'k--',label='Momentum T. (B=infty)',lw=3)
    ut_th_inf=np.zeros(vr.shape)
    ut_th_inf[vr<R]=0
    ut_th_inf[vr>R]=Gamma_tot/(2*np.pi*vr[vr>R])/2
    ax.plot(vr/R,ut_th_inf,'k--',label='Vortex cylinder (B=infty)',lw=3)

    for nB in [1, 3, 10, 500]:
        ur,ut,uz = vh_u(Xcp,Ycp,Zcp,Gamma_tot/nB,R,h,psih=psi_blade,nB=nB,bWT=bWT,method=method,bSemi=bSemi)
        ax.plot(vr/R,ut, label='B = {:d}'.format(nB))

    #a,U0,R,Lambda,nB,h,CT,Gamma_B=Params(Lambda=0.1,nB=10)
    #h=1000
    #nB=10
    #ur,ut,uz = vh_u(Xcp,Ycp,Zcp,Gamma_tot/nB,R,h,psih=psi_blade,nB=nB,bWT=bWT,method=method,bSemi=bSemi)
    #ax.plot(vr/R,ut, label='B= {:d} - Lambda = {:.1f}'.format(nB,Lambda))


    ax.set_xlabel('r/R [-]')
    ax.set_ylabel('ut [m/s]')
    ax.set_ylim([-1, 1])
    ax.set_title('Tangential induced velocity on the lifting line')
    ax.legend()




    plt.show()



if __name__ == '__main__':
    main()
