"""
Tangential Vorticity from a skewed vortex cylinder.
This test reproduces the plots of the section 3.1 of the following reference:
    [1] E. Branlard, M. Gaunaa - Cylindrical vortex wake model: skewed cylinder, application to yawed or tilted rotors - Wind Energy, 2015
""" 


#--- Legacy python 2.7
# --- General
import unittest
import matplotlib.pyplot as plt
import numpy as np
from numpy import pi
# --- Local
from welib.vortilib.elements.VortexCylinderSkewed import fV_Tangential, fKxit, fKzt
# 
# --------------------------------------------------------------------------------}
# --- Helper functions 
# --------------------------------------------------------------------------------{
def pol2cart(rho,phi):
    x= rho * np.cos(phi)
    y= rho * np.sin(phi)
    return (x,y)


def rsquare(y,f, c = True): 
    """ Compute coefficient of determination of data fit model and RMSE
    [r2 rmse] = rsquare(y,f)
    [r2 rmse] = rsquare(y,f,c)
    RSQUARE computes the coefficient of determination (R-square) value from
    actual data Y and model data F. The code uses a general version of
    R-square, based on comparing the variability of the estimation errors
    with the variability of the original values. RSQUARE also outputs the
    root mean squared error (RMSE) for the user's convenience.
    Note: RSQUARE ignores comparisons involving NaN values.
    INPUTS
      Y       : Actual data
      F       : Model fit
    
    # OPTION
      C       : Constant term in model
                R-square may be a questionable measure of fit when no
              constant term is included in the model.
      [DEFAULT] TRUE : Use traditional R-square computation
               FALSE : Uses alternate R-square computation for model
                     without constant term [R2 = 1 - NORM(Y-F)/NORM(Y)]
    # OUTPUT
      R2      : Coefficient of determination
      RMSE    : Root mean squared error """
    # Compare inputs
    if not np.all(y.shape == f.shape) :
        raise Exception('Y and F must be the same size')
    # Check for NaN
    tmp = np.logical_not(np.logical_or(np.isnan(y),np.isnan(f))) 
    y = y[tmp]
    f = f[tmp]
    if c:
        r2 = max(0,1-np.sum((y-f)**2)/np.sum((y-np.mean(y))** 2))
    else:
        r2 = 1 - np.sum((y - f) ** 2) / np.sum((y) ** 2)
        if r2 < 0:
            import warnings
            warnings.warn('Consider adding a constant term to your model')
            r2 = 0
    rmse = np.sqrt(np.mean((y - f) ** 2))
    return r2,rmse

# --------------------------------------------------------------------------------}
# --- Main 
# --------------------------------------------------------------------------------{
def main(test=False):
    # Parameters (!!!moved to a cell-wise pseudo standalone, so paramters are overriden)
    R               = 1
    gamma_rings     = - 1
    bBW             = 1
    npsi            = 3
    nr              = 100
    ntheta          = 600
    nz              = 100
    bNormalizeRotor = 0
    epsilon         = 10 ** - 3
    m               = np.tan(30 * pi / 180)
    # derived

    ## --- Quiver plot in rotor plane
    vchi = np.array([30]) * pi / 180
    for ichi in np.arange(len(vchi)):
        chi = vchi[ichi]
        m = np.tan(chi)
        ntheta = 600
        # Control Points
        npsi = 12
        nr = 8
        vpsi = np.linspace(0,2*pi,npsi+1)
        vpsi = vpsi[:-1]
        vr   = np.linspace(0,R,nr) / R
        vr   = vr[np.abs((np.abs(vr) - 1)) > epsilon]
        vz0 = 0
        VR,VPSI,VZ0 = np.meshgrid(vr,vpsi,vz0)
        X = np.multiply(VR,np.cos(VPSI))
        Y = np.multiply(VR,np.sin(VPSI))
        scale = 2
        # Ui
        uz,uipsi,ux,uy = fV_Tangential(VR,VPSI,VZ0,m,gamma_rings,ntheta,nout=4)
        plt.figure
        plt.quiver(np.squeeze(X),np.squeeze(Y),np.squeeze(ux)/scale,np.squeeze(uy)/scale,scale=1,angles='xy',scale_units='xy',headwidth=2,headlength=3,headaxislength=3)
        #ArrowLength = 0.12
        #ArrowAngle = 20
        #vectarrowb(np.array([1.2,0.5]),np.array([1.2,1]),1,'k',1.5,ArrowLength,ArrowAngle)
        #plt.text(1.16,0.75,'= 1','HorizontalAlignment','Left','VerticalAlignment','Middle')
        vtheta = np.linspace(0,2*pi,100)
        plt.plot(R * np.cos(vtheta),R * np.sin(vtheta),'k')
        plt.axis('square')
        plt.xlim(np.array([1.3,-1.3])) # NOTE FLIPPED!
        plt.ylim(np.array([- 1.3,1.3]))
        plt.title('CastlesQuiver_chi{:03d}'.format(int(np.round(np.arctan(m)*180/pi))))
        plt.xlabel('$x/R$')
        plt.ylabel('$y/R$')


    ## --- Side view only
    LIM         = 1.58
    vchi        = np.array([30]) * pi / 180
    uz0         = - 0.5
    gamma_rings = - 1
    chi         = vchi[0]
    m           = np.tan(chi)
    ntheta      = 600
    vr_low      = np.linspace(- 0.95,0.95,20)
    vr_high     = np.linspace(- 0.99,0.99,100)
    vpsi        = 0
    vz0         = 0
    # low res
    VR,VPSI,VZ0 = np.meshgrid(vr_low,vpsi,vz0)
    uz = fV_Tangential(VR,VPSI,VZ0,m,gamma_rings,ntheta,nout=1)
    # high res
    VR,VPSI,VZ0 = np.meshgrid(vr_high,vpsi,vz0)
    uz_high = fV_Tangential(VR,VPSI,VZ0,m,gamma_rings,ntheta,nout=1)
    plt.figure()
    # Rotor Plane
    plt.plot(np.array([- 1,1]),np.array([0,0]),'k')
    # Tan chi/2 slope
    plt.plot(np.array([1,- 1]) * LIM,np.array([- np.tan(chi / 2),np.tan(chi / 2)]) * LIM,'k-.','LineWidth',1.5)
    # Induction envelope
    plt.plot(vr_high,(uz_high.ravel()-uz0) / np.abs(uz0))
    # Inductions
    X = np.array([vr_low])
    Y = np.array([0 * vr_low])
    scale = 1
    LW = 1
#     try:
    plt.quiver(X.ravel(),Y.ravel() - uz.ravel()/ scale,0*X.ravel(),uz.ravel()/scale, scale=1,angles='xy',scale_units='xy',headwidth=2,headlength=3,headaxislength=3)
#     except:
#         print('Quiver failed')
    # Restoring Moment
    # ArrowLength = 0.22
    # ArrowAngle = 30
    # MomentStart = 0.3
    # LW = 1
    # fDrawRotArrow(MomentStart,220,140,LW,'k',2,ArrowLength / 2,ArrowAngle)
    # Free stream
    # alpha = 90 - 18
    # X0 = 0
    # Y0 = - 0.0
    # U0 = 1.5
    # ArrowLength = 0.12
    # ArrowAngle = 20
    # LW = 1
    # vectarrowb(np.array([X0 - U0 * np.cos(np.pi/180*alpha),Y0 - U0 * np.sin(np.pi/180*alpha)]),np.array([X0,Y0]),scale,'k',1.5,ArrowLength,ArrowAngle)
    # Wake Axis
    plt.plot(np.array([0,np.tan(chi)]) * LIM,np.array([0,1]) * LIM,'k--',linewidth=1.5)
    # Wake Envelop
    plt.plot(np.array([0,np.tan(chi)]) * LIM - 1,np.array([0,1]) * LIM,'k-',linewidth=1.5)
    plt.plot(np.array([0,np.tan(chi)]) * LIM + 1,np.array([0,1]) * LIM,'k-',linewidth=1.5)
    # text(0.637,1.218,'Wake axis','FontSize',11)
    # text(0.32,- 0.746,'Wind','FontSize',11)
    # text(1.538,- 0.122,'tan(\chi/2) slope','FontSize',11)
    # text(- 0.628,0.855,'Envelope','Color',fColrs(1),'FontSize',11)
    # text(- 0.408,- 0.128,'Rotor Plane','FontWeight','bold','FontSize',11)
    # text(- 0.59,- 0.493,'(Moment)','FontSize',11)
    # text(0.5,0.676,'Inductions','Color',fColrs(1),'FontSize',11)
    # Annotation for moment
    # plt.plot(np.array([- 0.544,- MomentStart]),np.array([- 0.433,- 0.2]),'k','LineWidth',0.5)
    # Legend for arrows
    # alpha = 90 - 18
    # X0 = 0
    # Y0 = - 0.0
    # U0 = 1.5
    # vectarrowb(np.array([1.35,- 1.5]),np.array([1.35,- 1]),scale,'k',1.5,ArrowLength,ArrowAngle)
    # text(1.32,- 1.26,'= 0.5','HorizontalAlignment','Left','VerticalAlignment','Middle','FontSize',11)
    plt.xlabel('$x/R$')
    plt.ylabel('$z/R$ and dimensionless velocity')
    plt.axis('square')
    plt.xlim(np.array([LIM,-LIM])) # NOTE: flipped
    plt.ylim(np.array([- LIM,LIM]))
    plt.title('TangentialOnlySideViewVz_chi{:03d}'.format(int(np.round(np.arctan(m)*180/pi))))
    # view(np.array([- 180,- 90]))

    ## --- Contour plot of axial velocity - Computation
    n_azimuth     = 180
    n_radial      = 25
    vpsi          = np.linspace(0,360,n_azimuth) * pi / 180
    vr            = np.linspace(0.0,0.99,n_radial) * R
    vz0           = 0
    chi           = 30 * pi / 180
    m             = np.tan(chi)
    VR,VPSI,VZ0   = np.meshgrid(vr,vpsi,vz0)
    VR   = VR.squeeze()
    VPSI = VPSI.squeeze()
    VZ0  = VZ0.squeeze()
    Xplane,Yplane = pol2cart(VR,VPSI)
    # Ui
    uz = fV_Tangential(VR,VPSI,VZ0,m,gamma_rings,ntheta,nout=1)
    uz=np.squeeze(uz)
    plt.figure()
    ContourVal = np.array([0,0.05,0.1,0.2])
    ContourVal = np.unique(np.sort(np.concatenate((-ContourVal,ContourVal)) - 0.5))
    manual_locations=[(-0.01,0.01),(0.34,0.01),(0.58,0.01),(0.86,0.01),(-0.36,0.01),(-0.62,0.01),(-0.86,0.01)]
    CS=plt.contour(Xplane/R, Yplane/R, uz, levels=ContourVal)
    plt.clabel(CS,inline=1,fontsize=10,manual=manual_locations)
    vtheta = np.linspace(0,2 * pi,100)
    plt.plot(R * np.cos(vtheta),R * np.sin(vtheta),'k')
    plt.axis('square')
    plt.xlim(np.array([  1.3,-1.3])) # NOTE: flipped
    plt.ylim(np.array([- 1.3,1.3]))
    plt.title('TangentialOnlyContourAxial_chi{:03d}'.format(int(np.round(np.arctan(m)*180/pi))))
    plt.xlabel('$x/R$')
    plt.ylabel('$y/R$')
    # --------------------------------------------------------------------------------
    # --- Z COMPONENT
    # --------------------------------------------------------------------------------
    # --- Z component Azimuthal survey
    m    = np.tan(30 * pi / 180)
    npsi = 3
    nr   = 200
    vpsi = np.linspace(0,pi / 2,npsi)
    vr   = np.linspace(- 2 * R,2 * R,nr) / R
    vr   = vr[np.abs((np.abs(vr) - 1)) > epsilon]
    vz0  = 0
    VR,VPSI,VZ0 = np.meshgrid(vr,vpsi,vz0)
    uiz = fV_Tangential(VR,VPSI,VZ0,m,gamma_rings,ntheta,nout=1)
    plt.figure()
    legds = []
    for ip in np.arange(len(vpsi)):
        plt.plot(vr,np.squeeze(uiz[ip,:,:]))
        legds.append('$\\psi=%3d$ deg.'%np.round(vpsi[ip] * 180 / pi))

    plt.legend(legds)
    plt.xlabel(r'Radial position $r/R$')
    plt.ylabel(r'Axial induced velocity $-u_{z,t}/\gamma_t$ ')
    plt.title('CastlesAzimuth_chi{:03d}'.format(int(np.round(np.arctan(m)*180/pi))))
    # --- Z component Radial survey
    npsi = 360
    vpsi = np.linspace(- pi,pi,npsi)
    vr = np.array([0.1,0.5,0.9])
    vz0 = 0
    VR,VPSI,VZ0 = np.meshgrid(vr,vpsi,vz0)
    uiz = fV_Tangential(VR,VPSI,VZ0,m,gamma_rings,ntheta,nout=1)
    plt.figure()
    legds = []
    for ip in np.arange(len(vr)):
        plt.plot(vpsi*180/pi,np.squeeze(uiz[:,ip,:]))
        legds.append('$r/R=%2.1f$'%vr[ip])
    plt.legend(legds)
    plt.xlabel(r'Azimuthal position $\psi$ [deg]')
    plt.ylabel(r'Axial induced velocity  $-u_{z,t}/\gamma_t$')
    plt.xlim(np.array([- 180,180]))
    plt.title('CastlesRadialchi{:03d}'.format(int(np.round(np.arctan(m)*180/pi))))

    # --- Flow expansion - Yaw article plot
    gamma_rings     = - 1
    nchi            = 3
    vchi            = np.linspace(0,90,nchi) * pi / 180
    nr              = 100
    vr              = np.linspace(0,1,nr)
    vr              = vr[:-1]
    nchi            = len(vchi)
    ntheta          = 600
    bNormalizeRotor = 0
    nz              = 0
    plt.figure()
    legds = []
    for ichi in np.arange(nchi):
        chi = vchi[ichi]
        m = np.tan(chi)
        Kxit,Kxit_num,fOye = fKxit(vr,m)
        plt.plot(vr,Kxit_num*(np.cos(chi/2)**2),'-')
        legds.append('$\\chi={:2d}$ deg.'.format(int(np.round(vchi[ichi] * 180 / pi))))

    legds.append('$F_t$ approx')#,np.round(vchi(ichi) * 180 / pi))
    plt.plot(vr,fOye,'k--')
    plt.legend(legds)
    plt.xlabel(r'Radial position $r/R$')
    plt.ylabel(r'Flow expansion $F_{t} (r,\chi)$')
    plt.title('FlowExpansion')



    # ## Just the influence function and the slope
    nchi = 3
    vchi = np.linspace(0,90,nchi) * pi / 180
    nr = 100
    vr = np.linspace(0,1,nr)
    vr = vr[:-1]
    plt.figure()
    legds = []
    slope=np.zeros(nchi)
    slope_th=np.zeros(nchi)
    for ichi in np.arange(nchi):
        chi = vchi[ichi]
        m = np.tan(chi)
        vpsi = 0
        vz0 = 0
        VR,VPSI,VZ0 = np.meshgrid(vr,vpsi,vz0)
        uiz = fV_Tangential(VR,VPSI,VZ0,m,gamma_rings,ntheta,nout=1).squeeze()
        slope[ichi]    = (uiz[2] - uiz[0]) / (vr[2] - vr[0]) / uiz[0]
        slope_th[ichi] = np.tan(chi/2)
        plt.plot(vr,np.squeeze(uiz)/uiz[0]-1)
        legds.append('$\\chi={:2d}$ deg.'.format(int(np.round(vchi[ichi] * 180 / pi))))

    for ichi in np.arange(nchi):
        plt.plot(vr,vr * slope_th[ichi],'--')

    legds.append('Slope at r=0')
    plt.legend(legds)
    plt.xlabel(r'Radial position r/R')
    plt.ylabel(r'$K_{z,t} (r,\chi)$')
    plt.title('CastlesKzChi')

    ## Z component - Is it a really cosine behavior? NO IT IS NOT - NOTE: NOT FINISHED
    vchi = np.linspace(0,pi / 2,15)
    vchi = vchi[1:-1]
    ntheta = 1000

    npsi = 180
    vpsi = np.linspace(- pi,pi,npsi)
    vr = np.array([0.4,0.7,0.95])
    vz0 = 0
    VR,VPSI,VZ0 = np.meshgrid(vr,vpsi,vz0)
    bBW = 1
    err_per_chi_r = np.zeros((len(vchi),len(vr)))
    gamma_rings=-1
    gamma_t = gamma_rings
    uz0 = gamma_t / 2
    for ichi,chi in enumerate(vchi):
        #print('.')
        # Radial survey - z component
        m = np.tan(chi)
        uz = fV_Tangential(VR,VPSI,VZ0,m,gamma_t,ntheta,nout=1).squeeze()
        K_r = uz[0,:]/uz0 - 1
        Kzt,Kztnum = fKzt(vr,m,nout=2)
        # Computing error
        for ir in np.arange(len(vr)):
            uz_model = uz0 * (1 + Kzt[ir] * np.cos(vpsi))
            err_per_chi_r[ichi,ir],_ = rsquare(uz[:,ir].ravel(),uz_model.ravel())
    err_per_chi_r[0,:] = 1
    bBW = 1
    plt.figure()
    legds = []
    for ir in np.arange(len(vr)):
        plt.plot(vchi*180/pi,err_per_chi_r[:,ir])
        legds.append('$r/R=%2.1f$'%vr[ir])
    
    plt.xlim(np.array([0,90]))
    plt.legend(legds)
    plt.xlabel(r'Skew angle $\chi$ [deg]')
    plt.ylabel(r'$R^2$')
    plt.title('ErrorSinModelForKzt')

    if not test:
        plt.show()



class Test(unittest.TestCase):
    def test_Article_Skew_TangentialVorticity(self):
        import sys
        if sys.version_info >= (3, 0):
            main(test=True)
        else:
            print('Test skipped due to travis display error')

if __name__ == "__main__":
    main()
