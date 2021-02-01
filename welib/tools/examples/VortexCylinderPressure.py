import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# Local 
from vortilib.elements.VortexCylinder import vc_tang_u
from pybra.fields    import *
from pybra.functions import *


def setPressure(p,X,R,x0,r0,p0):
    vxx=np.unique(X.ravel())
    vrr=np.unique(R.ravel())
    ir=np.argmin(np.abs(vrr-r0))
    ix=np.argmin(np.abs(vxx-x0))
    print('Setting pressure at x={:.2f} r={:.2f} ({:.2f} {:.2f}): {:.5f}'.format(vxx[ix],vrr[ir],x0,r0,p0))
    p = p - p[ir,ix] + p0
    return p

def getPressure(p,X,R,x0,r0):
    vxx=np.unique(X.ravel())
    vrr=np.unique(R.ravel())
    ir=np.argmin(np.abs(vrr-r0))
    ix=np.argmin(np.abs(vxx-x0))
    p0=p[ir,ix]
    print('Getting pressure at x={:.2f} r={:.2f} ({:.2f} {:.2f}): {:.5f}'.format(vxx[ix],vrr[ir],x0,r0,p0))
    return p0

def pressure_fourquadrants(vr,vx,ur,ux,xoff=0.1,roff=0.1):
    """ 
    1  | 2
    3  | 4
    """
    # acceleration
    ar, ax = vdel_cylindrical_axi( (ur,ux), (ur,ux), vr, vx)
    # Pressure gradient
    dpdr= - rho* ar # NOTE: no force since we avoid the AD
    dpdx= - rho* ax # NOTE: no force since we avoid the AD
    pBern = -rho*( Ef + (ur**2+ux**2)/2 )

    Xcp,Rcp = np.meshgrid(vx,vr)  # shape nr x nx
    b1 = np.logical_and( Xcp<0-xoff, Rcp>1+roff )
    b2 = np.logical_and( Xcp>0+xoff, Rcp>1+roff )
    b3 = np.logical_and( Xcp<0-xoff, Rcp<1-roff )
    b4 = np.logical_and( Xcp>0+xoff, Rcp<1-roff )
    Ir=np.array(np.arange(len(vr)))
    Ix=np.array(np.arange(len(vx)))
    br1 = Ir[vr>1+roff]
    br2 = Ir[vr>1+roff]
    br3 = Ir[vr<1-roff]
    br4 = Ir[vr<1-roff]
    bx1 = Ix[vx<0-xoff]
    bx2 = Ix[vx>0+xoff]
    bx3 = Ix[vx<0-xoff]
    bx4 = Ix[vx>0+xoff]
    X1, R1 = Xcp[np.ix_(br1,bx1)],Rcp[np.ix_(br1,bx1)]
    X2, R2 = Xcp[np.ix_(br2,bx2)],Rcp[np.ix_(br2,bx2)]
    X3, R3 = Xcp[np.ix_(br3,bx3)],Rcp[np.ix_(br3,bx3)]
    X4, R4 = Xcp[np.ix_(br4,bx4)],Rcp[np.ix_(br4,bx4)]

    p1 = intgrad((dpdx[np.ix_(br1,bx1)],dpdr[np.ix_(br1,bx1)]), (vx[bx1],vr[br1]), constant=0, tol=1e-16)
    p2 = intgrad((dpdx[np.ix_(br2,bx2)],dpdr[np.ix_(br2,bx2)]), (vx[bx2],vr[br2]), constant=0, tol=1e-16)
    p3 = intgrad((dpdx[np.ix_(br3,bx3)],dpdr[np.ix_(br3,bx3)]), (vx[bx3],vr[br3]), constant=0, tol=1e-16)
    p4 = intgrad((dpdx[np.ix_(br4,bx4)],dpdr[np.ix_(br4,bx4)]), (vx[bx4],vr[br4]), constant=0, tol=1e-16)


    # Normalizing pressure..
    rmax=np.max(vr)
    xmax=np.max(vx)
    xmin=np.min(vx)
    p0=getPressure(p1, X1, R1, -xoff, rmax)

    p0=0
    p1    = setPressure(p1, X1, R1     , -xoff, rmax, p0)
    p2    = setPressure(p2, X2, R2     ,  xoff, rmax, p0)
    pBern = setPressure(pBern, Xcp, Rcp, -xoff, rmax, p0)


    p31   = getPressure(p1, X1, R1     , xmin, 1+roff)
    p3    = setPressure(p3, X3, R3     , xmin, 1-roff, p31)
    p42   = getPressure(p2, X2, R2     , xmax, 1+roff)
    p42   = getPressure(pBern, Xcp, Rcp, xmax, 1+roff)
    p4    = setPressure(p4, X4, R4     , xmax, 1-roff, p42) # <<< This is the critical part


    X=np.block([[X3,X4],[X1,X2]])
    R=np.block([[R3,R4],[R1,R2]])
    p=np.block([[p3,p4],[p1,p2]])



    plevels=np.linspace(-0.6*DeltaP,0.6*DeltaP,20)
    fig,axes = plt.subplots(2,2, figsize=(9,5))
#     im=axes[0,0].contourf(X1, R1, p1, plevels)
#     im=axes[0,0].contourf(X2, R2, p2, plevels)
#     im=axes[0,0].contourf(X3, R3, p3, plevels)
#     im=axes[0,0].contourf(X4, R4, p4, plevels)
    im=axes[0,0].contourf(X, R, p, plevels)
    cb=fig.colorbar(im, ax=axes[0,0])
    im=axes[0,1].contourf(Xcp,Rcp, pBern, plevels)
    cb=fig.colorbar(im, ax=axes[0,1])

    axes[1,0].plot(vx,pBern[1,:]  , '-', label='Bern')
#     axes[1,0].plot(vx,p[1,:], '--', label='p')
    axes[1,0].plot(vx[bx3],p3[1,:], '--', label='p3')
    axes[1,0].plot(vx[bx4],p4[1,:], '--', label='p4')
    axes[1,0].legend()

    
def main():
    nx=250
    nr=200
    CT=0.4

    rho     = 1
    U0      = 1
    R       = 1
    A       = np.pi*R**2
    gamma_t = -U0*(1-np.sqrt(1-CT))
    DeltaP  = CT * 1/2*rho*U0**2
    print('DeltaP:', DeltaP)

    # All
    vx      = np.linspace(-3, 6, nx); dx=vx[1]-vx[0]
    vr      = np.linspace(0 ,2, nr); dr=vr[1]-vr[0]

    # Upstream - OK
    # vx      = np.linspace(-3,-0.1, nx); dx=vx[1]-vx[0]
    # vr      = np.linspace(0 ,2, nr); dr=vr[1]-vr[0]
    # # 
    # # # Donwstream (issue with forces along r)
    # vx      = np.linspace(0.1,3, nx); dx=vx[1]-vx[0]
    # vr      = np.linspace(0 ,2, nr); dr=vr[1]-vr[0]
    # 
    # Donwstream-lower - OK
    # vx      = np.linspace(0.1, 3  , nx); dx=vx[1]-vx[0]
    # vr      = np.linspace(0  ,0.99, nr); dr=vr[1]-vr[0]

    # # Donwstream-upper - OK
    # vx      = np.linspace(0.1, 3  , nx); dx=vx[1]-vx[0]
    # vr      = np.linspace(1.1 ,1.99, nr); dr=vr[1]-vr[0]
    #  upper - OK
    # vx      = np.linspace(-3, 3  , nx); dx=vx[1]-vx[0]
    # vr      = np.linspace(1.1 ,1.99, nr); dr=vr[1]-vr[0]

    Xcp,Rcp = np.meshgrid(vx,vr)  # shape nr x nx
    ur, ux  = vc_tang_u(Rcp,Rcp*0,Xcp,gamma_t = gamma_t)
    ur=ur*np.sign(Rcp)
    ux=ux+1

    Ef = np.zeros(ur.shape)
    fx = np.zeros(ur.shape)
    om = np.zeros(ur.shape)
    bCylinder= np.abs(Rcp)<=R 
    nCell=2
    bDisk = np.logical_and( Xcp>=0 , Xcp<nCell*dx ) # 2cell extent
    bDisk = np.logical_and( bDisk, bCylinder)
    bHeavyM= np.logical_and(Xcp<0, bCylinder) 
    bHeavyP= np.logical_and(Xcp>0, bCylinder) 



    epsilonx=20*dx
    epsilonr=2*dx
    print('Epsilon x,r',epsilonx, epsilonr)

    # --- Vorticity analytical and from numerical gradient
    om_theta = gamma_t # Note /dr included in delta function
    omega    = om_theta * delta(Rcp-R, epsilon=epsilonr) * Pi(Xcp, epsilon=epsilonx)
    omega_num = curl_cylindrical_axi((ur,ux), vr, vx)

    # --- Force, analytical
    fx =  -DeltaP * delta(Xcp, epsilon=epsilonx) * Pi(R-np.abs(Rcp), epsilon=epsilonr)
    fr =  fx*0

    # --- Force Potential, analytical and numerical
    Ef[bHeavyM]=0
    Ef[bHeavyP]=DeltaP # Note: force is negative, pointing along -x
    dEdx= -fx
    dEdr= -fr
    E_num = intgrad((dEdx,dEdr), (vx,vr), constant=0, tol=1e-16)
    E_num = setPressure(E_num  , Xcp, Rcp,-3, 0, 0)

    # --- Force, gradient of potential
    frb, fxb = np.gradient(-E_num, vr, vx)
    frc, fxc = np.gradient(-Ef   , vr, vx)

    # --- Verity that integration over disk gives DeltaP
    i=np.argmin(np.abs(vr))
    fx_int  = np.trapz(fx[i,:], vx)
    fxb_int = np.trapz(fxb[i,:], vx)
    print('Fx integration:', fx_int, fxb_int, 'DeltaP:',DeltaP)

    # --- Verity that vorticity integration across cylinder gives gamma_t
    i=np.argmin(np.abs(vx-1.0))
    bRp=vr>0
    om_int      = np.trapz(omega[bRp,i], vr[bRp])
    om_int_num  = np.trapz(omega_num[bRp,i], vr[bRp])
    print('omega integration:', om_int, om_int_num, 'gamma_t:',gamma_t)



    print(Xcp.shape)
    print(ur.shape)


    # --- Compare different input variables
    levOm = np.linspace(1.1*gamma_t/dr,0, 30)
    levE  = np.linspace(-DeltaP,DeltaP, 30)
    fig,axes = plt.subplots(4,2, figsize=(9,7))
    fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.50, wspace=0.20)
    ax=axes[0,0];im=ax.contourf(Xcp, Rcp, omega    ,levOm); cb=fig.colorbar(im, ax=ax);ax.set_title('omega')
    ax=axes[0,1];im=ax.contourf(Xcp, Rcp, omega_num,levOm); cb=fig.colorbar(im, ax=ax);ax.set_title('omega num')
    ax=axes[1,0];im=ax.contourf(Xcp, Rcp, fx             ); cb=fig.colorbar(im, ax=ax);ax.set_title('fx')
    ax=axes[1,1];im=ax.contourf(Xcp, Rcp, fxb            ); cb=fig.colorbar(im, ax=ax);ax.set_title('fx=d(Enum)/dx')
    ax=axes[2,0];im=ax.contourf(Xcp, Rcp, fr             ); cb=fig.colorbar(im, ax=ax);ax.set_title('fr')
    ax=axes[2,1];im=ax.contourf(Xcp, Rcp, frb            ); cb=fig.colorbar(im, ax=ax);ax.set_title('fr=d(Enum)/dr')
    ax=axes[3,0];im=ax.contourf(Xcp, Rcp, Ef    , levE   ); cb=fig.colorbar(im, ax=ax);ax.set_title('E (Heaviside)')
    ax=axes[3,1];im=ax.contourf(Xcp, Rcp, E_num , levE   ); cb=fig.colorbar(im, ax=ax);ax.set_title('E_num=intgrad(f)')

    # --- Integrate pressure on four quadrants
    pressure_fourquadrants(vr,vx,ur,ux,xoff=0.000,roff=0.000)

    # --- Verifying if u GradU = grad(u**2/2) - u x omega
    # ar, ax = vdel_cylindrical_axi( (ur,ux), (ur,ux), vr, vx)
    # 
    # u_x_om = (-ux*omega_num, ur*omega_num)
    # u2 = nr**2+ux**2
    # du2dr, du2dx = np.gradient(u2/2, vr, vx)
    # ar_bis = du2dr- u_x_om[0]
    # ax_bis = du2dx- u_x_om[1]
    # fig,axes = plt.subplots(2,2, figsize=(9,7))
    # # im=axes[0,0].contourf(Xcp, Rcp, omega    ); cb=fig.colorbar(im, ax=axes[0,0])
    # # im=axes[0,1].contourf(Xcp, Rcp, omega_num); cb=fig.colorbar(im, ax=axes[0,1])
    # im=axes[0,0].contourf(Xcp, Rcp, ar    ); cb=fig.colorbar(im, ax=axes[0,0])
    # im=axes[0,1].contourf(Xcp, Rcp, ar_bis); cb=fig.colorbar(im, ax=axes[0,1])
    # im=axes[1,0].contourf(Xcp, Rcp, ax    ); cb=fig.colorbar(im, ax=axes[1,0])
    # im=axes[1,1].contourf(Xcp, Rcp, ax_bis); cb=fig.colorbar(im, ax=axes[1,1])
    # plt.show()


    # --- Integrate on specified domain
    Acc_r, Acc_x = vdel_cylindrical_axi( (ur,ux), (ur,ux), vr, vx)
    dpdr   = fr   - rho* Acc_r
    dpdx   = fx   - rho* Acc_x
    p      = intgrad((dpdx,dpdr), (vx,vr), constant= 0, tol = 1e-16)

    pBern = -rho*(Ef    + (ur**2+ux**2)/2 )
    # pBern = -rho*(E_num + (ur**2+ux**2)/2 )

    # Setting pressure gage
    p0=0
    rmax=np.max(vr)
    xmax=np.max(vx)
    if rmax>1:
        p     = setPressure(p    , Xcp, Rcp, 0, rmax, p0)
        pBern = setPressure(pBern, Xcp, Rcp, 0, rmax, p0)
    elif xmax>0:
        p     = setPressure(p    , Xcp, Rcp, xmax, 0, p0)
        pBern = setPressure(pBern, Xcp, Rcp, xmax, 0, p0)

    # dpdr3, dpdx3 = np.gradient(pBern, vr, vx)
    # # p= intgrad((dpdx3,dpdr3), (vx,vr), constant=pBern[0,0], tol=1e-16)
    # # x,r    = np.meshgrid(xp,rp)       
    # # fr,fx  = np.gradient(f, rp, xp)      
    # plevels=np.linspace(-0.5,0.5,9)
    plevels=np.linspace(-0.6*DeltaP,0.6*DeltaP,20)
    # # plevels=10
    fig,axes = plt.subplots(2,2, figsize=(9,5))
    im=axes[0,0].contourf(Xcp,Rcp, p, plevels); cb=fig.colorbar(im, ax=axes[0])
    im=axes[0,1].contourf(Xcp,Rcp, pBern, plevels); cb=fig.colorbar(im, ax=axes[1])
    axes[1,0].plot(vx,pBern[1,:], '-', label='Bern')
    axes[1,0].plot(vx,p    [1,:], '--', label='p')
    axes[1,0].legend()


    # im=ax.contourf(Xcp,Rcp, ur**2+(1+ux)**2)
    # im=axes[0].contourf(Xcp,Rcp, dpdx)
    # im=axes[0].contourf(Xcp,Rcp, dpdr)
    # im=axes[0].contourf(Xcp,Rcp, p3-np.mean(p3))
    # im=axes[1].contourf(Xcp,Rcp, dpdx2)
    # im=axes[1].contourf(Xcp,Rcp, dpdr2)
    # im=axes[1].contourf(Xcp,Rcp, dpdx3)
    # im=axes[1].contourf(Xcp,Rcp, dpdr3)
    # fig,axe = plt.subplots()
    # im=axe.contourf(Xcp,Rcp, dpdx)
    # # im=axe.contourf(Xcp,Rcp, ur**2+(1+ux)**2)
    # # im=axe.contourf(Xcp,Rcp, fxb)
    # # im=axe.contourf(Xcp,Rcp, Ef)
    # cb=fig.colorbar(im)
    # # 

if __name__=="__main__":
    plt.show()
if __name__=="__test__":
    pass
if __name__=="__export__":
    pass
    #from welib.tools.repo import export_figs_callback
    #export_figs_callback(__file__)
