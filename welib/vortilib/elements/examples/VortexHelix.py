import numpy as np
import matplotlib.pyplot as plt
# Local 

from welib.vortilib.elements.VortexHelix import *


def Params(a=1/3,U0=10,Lambda=6,R=100,nB=3):
    h         = 2*np.pi*R*(1-a)/Lambda
    CT        = 4*a*(1-a)
    Gamma     = CT*np.pi*R*U0/(nB*Lambda) # NOTE assume aprime = 0, large tip-sped ratio 
    return  a,U0,R,Lambda,nB,h,CT,Gamma

def main():
    import matplotlib.pyplot as plt

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

    fig,axes = plt.subplots(1, 2, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
    fig.subplots_adjust(left=0.12, right=0.95, top=0.91, bottom=0.11, hspace=0.20, wspace=0.39)

    # --- Axial induction on the lifting line 
    ax = axes[0]
    ax.plot([0,1],[-a,-a],'k--',label='Momentum T. (B=infty)',lw=3.5)
    ax.plot([1,2],[ 0, 0],'k--',lw=3)

    for nB in [1, 3, 10, 500]:
        ur,ut,uz = vh_u(Xcp,Ycp,Zcp,Gamma_tot/nB,R,h,psih=psi_blade,nB=nB,bWT=bWT,method=method,bSemi=bSemi)
        ax.plot(Xcp/R,uz/U0, label='B = {:d}'.format(nB))

    ax.set_xlabel(r'$r/R$ [-]')
    ax.set_ylabel(r'$u_z/U_0$ [-]')
    ax.set_ylim([-2*a, 1/2*a])
    ax.tick_params(direction='in')
    #ax.legend()

    # --- Tangential induction on the lifting line 
    ax = axes[1]
    #ax.plot([0,1],[-a,-a],'k--',label='Momentum T. (B=infty)',lw=3)
    ut_th_inf=np.zeros(vr.shape)
    ut_th_inf[vr<R]=0
    ut_th_inf[vr>R]=Gamma_tot/(2*np.pi*vr[vr>R])/2
    ax.plot(vr/R,ut_th_inf,'k--',label='Vortex cylinder (B=infty)',lw=3.5)

    for nB in [1, 3, 10, 500]:
        ur,ut,uz = vh_u(Xcp,Ycp,Zcp,Gamma_tot/nB,R,h,psih=psi_blade,nB=nB,bWT=bWT,method=method,bSemi=bSemi)
        ax.plot(vr/R,ut, label='B = {:d}'.format(nB))

    #a,U0,R,Lambda,nB,h,CT,Gamma_B=Params(Lambda=0.1,nB=10)
    #h=1000
    #nB=10
    #ur,ut,uz = vh_u(Xcp,Ycp,Zcp,Gamma_tot/nB,R,h,psih=psi_blade,nB=nB,bWT=bWT,method=method,bSemi=bSemi)
    #ax.plot(vr/R,ut, label='B= {:d} - Lambda = {:.1f}'.format(nB,Lambda))


    ax.set_xlabel('r/R [-]')
    ax.set_ylabel(r'$u_t$ [m/s]')
    ax.set_ylim([-1, 1])
    ax.tick_params(direction='in')
    ax.legend()

    fig.suptitle('Vortilib - Vortex helix lifting line velocity')



if __name__ == '__main__':
    main()
    plt.show()
if __name__=="__export__":
    main()
    from welib.tools.repo import export_figs_callback
    export_figs_callback(__file__)
