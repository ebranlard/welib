"""
Trailed (longitudinal) Vorticity from a skewed vortex cylinder.
This test reproduces the plots of the section 3.3 of the following reference:
    [1] E. Branlard, M. Gaunaa - Cylindrical vortex wake model: skewed cylinder, application to yawed or tilted rotors - Wind Energy, 2015
""" 

## Parameters
# --- General
import unittest
import matplotlib.pyplot as plt
import numpy as np
from numpy import pi
# --- Local
from welib.vortilib.elements.VortexCylinderSkewed import fV_Trailed

def pol2cart(rho,phi):
    x= rho * np.cos(phi)
    y= rho * np.sin(phi)
    return (x,y)


def main(test=False):
    R               = 1
    gamma_longi     = 1
    bBW             = 1
    nr              = 100
    ntheta          = 601
    nz              = 100
    epsilon         = 10 ** - 3

    # --- Y-Component Azimuthal survey and Model !!!
    chi = 30 * pi / 180
    m = np.tan(chi)
    npsi = 360
    vpsi = np.linspace(- pi,pi,npsi)
    vr = np.array([0.2,0.5,0.8])
    nr = len(vr)
    vz0 = 0
    gamma_longi = 1
    VR,VPSI,VZ0 = np.meshgrid(vr,vpsi,vz0)
    uiz_trailed,v_psi,v_x,v_y,v_chi,v_xi,v_r = fV_Trailed(VR,VPSI,VZ0,m,gamma_longi,ntheta,nout=7)
    uy0 = gamma_longi * (np.sqrt(1 + m ** 2) - 1) / (2 * m)
    uy0 - gamma_longi / 2 * np.tan(chi / 2)
    # Kx,Ky0,Kz,Ky2 = fKl(vr,m)
    # Gl,Fl,Gl_fit,Fl_fit = fKlApprox(vr)
    plt.figure()
    for ir,r in enumerate(vr):
        plt.plot(vpsi*180/pi,np.squeeze(v_y[:,ir,:])/uy0,label='$r/R=%2.1f$'%r)
    plt.legend()
    plt.xlabel(r'Azimuthal position $\psi$ [deg]')
    plt.ylabel(r'Vertical induced velocity $u_{y,l} / u_{y,0}$')
    plt.xlim(np.array([- 180,180]))
    # set(gca,'XTick',np.array([- 180,- 90,0,90,180]))
    plt.title('TrailedOnlyUyRadial_chi{:03d}'.format(int(np.round(np.arctan(m) * 180 / pi))))


    # ## Polyfit of Myl and Ayl (G_l and F_l)
    # plt.figure()
    # hold('all')
    # grid('on')
    # box('on')
    # vr = np.linspace(0.0,1.0,100)
    # Mylapprox,Aylapprox = fKlApprox(vr)
    # options = optimset('TolX',0.01)
    # cM = fminsearch(lambda c = None: fPolyOrder135(c,vr,Mylapprox),np.array([0.1,0.1,0.1]),options)
    # __,__,pM = fPolyOrder135(cM,vr,Mylapprox)
    # cA = fminsearch(lambda c = None: fPolyOrder135(c,vr,Aylapprox),np.array([0.1,0.1,0.1]),options)
    # __,__,pA = fPolyOrder135(cA,vr,Aylapprox)
    # # # SIMPLEST
    # cM = np.array([0.1,- 0.0,0.7])
    # cA = np.array([0.1,- 0.2,0.6])
    # __,__,pM = fPolyOrder135(cM,vr,Mylapprox)
    # __,__,pA = fPolyOrder135(cA,vr,Aylapprox)
    # iSlope = 40
    # slopeG = Mylapprox(iSlope) / vr(iSlope)
    # slopeF = Aylapprox(iSlope) / vr(iSlope)
    # plt.plot(vr,Mylapprox,'k-')
    # plt.plot(vr,Aylapprox,'-','Color',fColrs(2,2,bBW))
    # plt.plot(vr,pM,'k--')
    # plt.plot(vr,pA,'--','Color',fColrs(2,2,bBW))
    # # plot(vr,slopeG*vr,'k-.')
    # # plot(vr,slopeF*vr,'-.','Color',fColrs(2,2,bBW))
    # plt.xlabel('r/R')
    # plt.ylabel('Expansion functions')
    # plt.legend('G_l','F_l','Polynomial fit',0)
    # plt.title('ExpansionFunctionTrailed')



    ## --- Z-Component Azimuthal survey
    plt.figure()
    gamma_longi = 1
    npsi        = 360
    vpsi        = np.linspace(- pi,pi,npsi)
    vr          = np.array([0.2,0.5,0.8])
    vz0         = 0
    chi         = 30 * pi / 180
    m           = np.tan(chi)
    VR,VPSI,VZ0 = np.meshgrid(vr,vpsi,vz0)
    uiz_trailed,__,uix_trailed = fV_Trailed(VR,VPSI,VZ0,m,gamma_longi,ntheta,nout=3)
    uy0 = gamma_longi * (np.sqrt(1 + m ** 2) - 1) / (2 * m)
    # Gl,Fl = fKlApprox(vr)
    for ir,r in enumerate(vr):
        plt.plot(vpsi*180/pi,np.squeeze(uiz_trailed[:,ir,:]), label='$r/R=%2.1f$'%r)
    plt.legend()
    plt.xlabel(r'Azimuthal position $\psi$ [deg]')
    plt.ylabel(r'Axial induced velocity $u_{z,l} / u_{y,0}$')
    plt.xlim(np.array([- 180,180]))
    # set(gca,'XTick',np.array([- 180,- 90,0,90,180]))
    plt.title('TrailedOnlyUzRadial_chi%03d'%np.round(np.arctan(m) * 180 / pi))



    ## Quiver plot in rotor plane
    # vchi=[0 30 60]*pi/180;
    gamma_longi = 1
    chi         = 30*pi/180
    m           = np.tan(chi)
    ntheta      = 600
    npsi        = 12
    nr          = 8
    vpsi = np.linspace(0,2 * pi,npsi + 1)
    vpsi = vpsi[:-1]
    vr = np.linspace(0,R,nr)/R
    vr = vr[np.abs((np.abs(vr) - 1)) > epsilon]
    vz0 = 0
    VR,VPSI,VZ0 = np.meshgrid(vr,vpsi,vz0)
    uz,uipsi,ux,uy = fV_Trailed(VR,VPSI,VZ0,m,gamma_longi,ntheta,nout=4)
    X = np.multiply(VR,np.cos(VPSI))
    Y = np.multiply(VR,np.sin(VPSI))
    scale = 9
    uy0 = gamma_longi * (np.sqrt(1 + m ** 2) - 1) / (2 * m)
    plt.figure()
    plt.quiver(np.squeeze(X),np.squeeze(Y),np.squeeze(ux) / uy0 / scale,np.squeeze(uy) / uy0 / scale,scale=1,angles='xy',scale_units='xy',headwidth=2,headlength=3,headaxislength=3)
    # ArrowLength = 0.05
    # ArrowAngle = 20
    # vectarrowb(np.array([1.2,0.5]),np.array([1.2,0.5 + 1 / scale]),1,'k',1.5,ArrowLength,ArrowAngle)
    # text(1.16,0.5 + 1 / (2 * scale),'= 1','HorizontalAlignment','Left','VerticalAlignment','Middle')
    vtheta = np.linspace(0,2 * pi,100)
    plt.plot(R * np.cos(vtheta),R * np.sin(vtheta),'k')
    plt.axis('square')
    plt.xlim(np.array([  1.3,- 1.3])) # NOTE: flipped
    plt.ylim(np.array([- 1.3,1.3]))
    plt.title('TrailedOnlyQuiver_chi%03d'%np.round(np.arctan(m) * 180 / pi))
    plt.xlabel('x/R')
    plt.ylabel('y/R')
    # 
    ## Surface plot of axial velocity
    # chi         = 30*pi/180
    # m           = np.tan(chi)
    #     ntheta = 600
    #     if (bBW == 1):
    #         # Control Points
    #         npsi = 30
    #         nr = 20
    #         #   npsi=20; nr=10;
    #         vpsi = np.linspace(0,2 * pi,npsi + 1)
    #         vr = np.linspace(0,R,nr) / R
    #         vr = vr(np.abs((np.abs(vr) - 1)) > epsilon)
    #         #   vr=vr([2 8:length(vr)])
    #         vr = np.array([0.1052632,0.3157895,0.5789474,0.7894737,0.8947368,0.9573684])
    #     else:
    #         # Control Points
    #         npsi = 50
    #         nr = 30
    #         #   npsi=20; nr=10;
    #         vpsi = np.linspace(0,2 * pi,npsi + 1)
    #         vr = np.linspace(0,R,nr) / R
    #         vr = vr(np.abs((np.abs(vr) - 1)) > epsilon)
    #     vz0 = 0
    #     VR,VPSI,VZ0 = np.meshgrid(vr,vpsi,vz0)
    #     X = np.multiply(VR,np.cos(VPSI))
    #     Y = np.multiply(VR,np.sin(VPSI))
    #     # Ui
    #     uz,uipsi,ux,uy = fV_Trailed(VR,VPSI,VZ0,m,gamma_longi,ntheta,nout=4)
    #     plt.figure
    #     plt.axis('square')
    #     vtheta = np.linspace(0,2 * pi,100)
    #     plt.plot(R * np.cos(vtheta),R * np.sin(vtheta),'k','LineWidth',1.5)
    #     zz = 0.1
    #     #$plot3(np.array([0,zz * np.tan(chi)]),np.array([0,0]),np.array([0,zz]),'k--','LineWidth',1.5)
    #     #   plot3([1 -1],[0 0],[-tan(chi/2) tan(chi/2)],'k-.','LineWidth',1.5)
    #     #plt.surf(np.squeeze(X),np.squeeze(Y),np.squeeze(uz),'EdgeColor','k')
    #     plt.axis('square')
    #     plt.title('TrailedOnlyContourVz_chi%03d'%np.round(np.arctan(m) * 180 / pi))
    #     plt.xlabel('x/R')
    #     plt.ylabel('y/R')
    #     plt.zlabel('z/R and -K_{z,l}')
    #     plt.legend('Rotor Plane','Wake axis')
    # 
    # 
    # # set(gcf,'color','w');
    # # setFigureTight(0);
    # # setFigureRePosition(0);
    # # setFigureDoNothing(1);
    # 
    ## Contour plot of axial velocity - Computation
    n_azimuth   = 180
    n_radial    = 35
    gamma_longi = 1
    uy0         = gamma_longi * (np.sqrt(1 + m ** 2) - 1) / (2 * m)
    vpsi        = np.linspace(0,360,n_azimuth) * pi / 180
    vr          = np.linspace(0.0,0.99,n_radial) * R
    vz0         = 0
    chi         = 30*pi/180
    m           = np.tan(chi)
    VR,VPSI,VZ0 = np.meshgrid(vr,vpsi,vz0)
    VR   = VR.squeeze()
    VPSI = VPSI.squeeze()
    VZ0  = VZ0.squeeze()
    Xplane,Yplane = pol2cart(VR,VPSI)
    # Ui
    uz = fV_Trailed(VR,VPSI,VZ0,m,gamma_longi,ntheta,nout=1)
    uz[Yplane==0] = 0
    #
    plt.figure()
    ContourVal = np.array([0.5,2,5,10])
    ContourVal =np.unique(np.sort(np.array(np.concatenate([- ContourVal,ContourVal]))))
    CS = plt.contour(Xplane/R, Yplane/R, 100*uz/uy0, levels=ContourVal)
    manual_locations=[(-0.608,0.603),(-0.48,0.48),(-0.338,0.34),(-0.196,0.196),(0.488,0.49),(0.482,-0.481),(-0.465,-0.468)]
    lbs=plt.clabel(CS,inline=1,fontsize=10,manual=manual_locations)
    rots=[-45.0,-45.0,-45.0,-45.0,-315.0,-45.0,-315.0]
    for l,r in zip(lbs,rots):
        l.set_rotation(r)
    plt.axis('equal')
    vtheta = np.linspace(0,2 * pi,100)
    plt.plot(R * np.cos(vtheta),R * np.sin(vtheta),'k')
    plt.plot(np.array([- R,R]),np.array([0,0]))
    plt.plot(np.array([0,0]),np.array([- R,R]))
    plt.axis('square')
    plt.xlim(np.array([  1.3,-1.3]))# NOTE: flipped
    plt.ylim(np.array([- 1.3,1.3]))
    plt.title('LongiOnlyContourAxial_chi%03d'%np.round(np.arctan(m) * 180 / pi))
    plt.xlabel('x/R')
    plt.ylabel('y/R')


    if not test:
        plt.show()
    else:
        plt.close('all')

class Test(unittest.TestCase):
    def test_Article_Skew_TrailedVorticity(self):
        import sys
        if sys.version_info >= (3, 0):
            main(test=True)
        else:
            print('Test skipped due to travis display error')

if __name__ == "__main__":
    main()
