#!/usr/bin/env python
from __future__ import print_function
from math import pi, sin, cos, radians, degrees, tan, ceil, floor
import numpy as np
import copy
import os
# from scipy.interpolate import RectBivariateSpline



class Polar(object):
    """
    Defines section lift, drag, and pitching moment coefficients as a
    function of angle of attack at a particular Reynolds number.

    """

    def __init__(self, Re, alpha, cl, cd, cm):
        """Constructor

        Parameters
        ----------
        Re : float
            Reynolds number
        alpha : ndarray (deg)
            angle of attack
        cl : ndarray
            lift coefficient
        cd : ndarray
            drag coefficient
        cm : ndarray
            moment coefficient
        """

        self.Re = Re
        self.alpha = np.array(alpha)
        self.cl = np.array(cl)
        self.cd = np.array(cd)
        self.cm = np.array(cm)

    @classmethod
    def fromfile(cls,filename,fformat='delimited'):
        """Constructor based on a filename"""
        if fformat=='delimited':
            if not os.path.exists(filename):
                raise Exception('File not found:',filename)
            M=np.loadtxt(filename,comments=['#','!'])
            if M.shape[1]!=4:
                raise Exception('Only reading delimited files with 4 columns : alpha cl cd cm')
            alpha = M[:,0]
            cl    = M[:,1]
            cd    = M[:,2]
            cm    = M[:,3]
            Re    = np.nan
            return cls(Re,alpha,cl,cd,cm,)
        else:
            raise NotImplementedError('Format not implemented: {}'.format(fformat))

    def blend(self, other, weight):
        """Blend this polar with another one with the specified weighting

        Parameters
        ----------
        other : Polar
            another Polar object to blend with
        weight : float
            blending parameter between 0 and 1.  0 returns self, whereas 1 returns other.

        Returns
        -------
        polar : Polar
            a blended Polar

        """

        # generate merged set of angles of attack - get unique values
        alpha = np.union1d(self.alpha, other.alpha)

        # truncate (TODO: could also have option to just use one of the polars for values out of range)
        min_alpha = max(self.alpha.min(), other.alpha.min())
        max_alpha = min(self.alpha.max(), other.alpha.max())
        alpha = alpha[np.logical_and(alpha >= min_alpha, alpha <= max_alpha)]
        # alpha = np.array([a for a in alpha if a >= min_alpha and a <= max_alpha])

        # interpolate to new alpha
        cl1 = np.interp(alpha, self.alpha, self.cl)
        cl2 = np.interp(alpha, other.alpha, other.cl)
        cd1 = np.interp(alpha, self.alpha, self.cd)
        cd2 = np.interp(alpha, other.alpha, other.cd)
        cm1 = np.interp(alpha, self.alpha, self.cm)
        cm2 = np.interp(alpha, other.alpha, other.cm)

        # linearly blend
        Re = self.Re + weight*(other.Re-self.Re)
        cl = cl1 + weight*(cl2-cl1)
        cd = cd1 + weight*(cd2-cd1)
        cm = cm1 + weight*(cm2-cm1)

        return type(self)(Re, alpha, cl, cd, cm)



    def correction3D(self, r_over_R, chord_over_r, tsr, alpha_max_corr=30,
                     alpha_linear_min=-5, alpha_linear_max=5):
        """Applies 3-D corrections for rotating sections from the 2-D data.

        Parameters
        ----------
        r_over_R : float
            local radial position / rotor radius
        chord_over_r : float
            local chord length / local radial location
        tsr : float
            tip-speed ratio
        alpha_max_corr : float, optional (deg)
            maximum angle of attack to apply full correction
        alpha_linear_min : float, optional (deg)
            angle of attack where linear portion of lift curve slope begins
        alpha_linear_max : float, optional (deg)
            angle of attack where linear portion of lift curve slope ends

        Returns
        -------
        polar : Polar
            A new Polar object corrected for 3-D effects

        Notes
        -----
        The Du-Selig method :cite:`Du1998A-3-D-stall-del` is used to correct lift, and
        the Eggers method :cite:`Eggers-Jr2003An-assessment-o` is used to correct drag.


        """

        # rename and convert units for convenience
        alpha = np.radians(self.alpha)
        cl_2d = self.cl
        cd_2d = self.cd
        alpha_max_corr = radians(alpha_max_corr)
        alpha_linear_min = radians(alpha_linear_min)
        alpha_linear_max = radians(alpha_linear_max)

        # parameters in Du-Selig model
        a = 1
        b = 1
        d = 1
        lam = tsr/(1+tsr**2)**0.5  # modified tip speed ratio
        expon = d/lam/r_over_R

        # find linear region
        idx = np.logical_and(alpha >= alpha_linear_min,
                             alpha <= alpha_linear_max)
        p = np.polyfit(alpha[idx], cl_2d[idx], 1)
        m = p[0]
        alpha0 = -p[1]/m

        # correction factor
        fcl = 1.0/m*(1.6*chord_over_r/0.1267*(a-chord_over_r**expon)/(b+chord_over_r**expon)-1)

        # not sure where this adjustment comes from (besides AirfoilPrep spreadsheet of course)
        adj = ((pi/2-alpha)/(pi/2-alpha_max_corr))**2
        adj[alpha <= alpha_max_corr] = 1.0

        # Du-Selig correction for lift
        cl_linear = m*(alpha-alpha0)
        cl_3d = cl_2d + fcl*(cl_linear-cl_2d)*adj

        # Eggers 2003 correction for drag
        delta_cl = cl_3d-cl_2d

        delta_cd = delta_cl*(np.sin(alpha) - 0.12*np.cos(alpha))/(np.cos(alpha) + 0.12*np.sin(alpha))
        cd_3d = cd_2d + delta_cd

        return type(self)(self.Re, np.degrees(alpha), cl_3d, cd_3d, self.cm)



    def extrapolate(self, cdmax, AR=None, cdmin=0.001, nalpha=15):
        """Extrapolates force coefficients up to +/- 180 degrees using Viterna's method
        :cite:`Viterna1982Theoretical-and`.

        Parameters
        ----------
        cdmax : float
            maximum drag coefficient
        AR : float, optional
            aspect ratio = (rotor radius / chord_75% radius)
            if provided, cdmax is computed from AR
        cdmin: float, optional
            minimum drag coefficient.  used to prevent negative values that can sometimes occur
            with this extrapolation method
        nalpha: int, optional
            number of points to add in each segment of Viterna method

        Returns
        -------
        polar : Polar
            a new Polar object

        Notes
        -----
        If the current polar already supplies data beyond 90 degrees then
        this method cannot be used in its current form and will just return itself.

        If AR is provided, then the maximum drag coefficient is estimated as

        >>> cdmax = 1.11 + 0.018*AR


        """

        if cdmin < 0:
            raise Exception('cdmin cannot be < 0')

        # lift coefficient adjustment to account for assymetry
        cl_adj = 0.7

        # estimate CD max
        if AR is not None:
            cdmax = 1.11 + 0.018*AR
        self.cdmax = max(max(self.cd), cdmax)

        # extract matching info from ends
        alpha_high = radians(self.alpha[-1])
        cl_high = self.cl[-1]
        cd_high = self.cd[-1]
        cm_high = self.cm[-1]

        alpha_low = radians(self.alpha[0])
        cl_low = self.cl[0]
        cd_low = self.cd[0]

        if alpha_high > pi/2:
            raise Exception('alpha[-1] > pi/2')
            return self
        if alpha_low < -pi/2:
            raise Exception('alpha[0] < -pi/2')
            return self

        # parameters used in model
        sa = sin(alpha_high)
        ca = cos(alpha_high)
        self.A = (cl_high - self.cdmax*sa*ca)*sa/ca**2
        self.B = (cd_high - self.cdmax*sa*sa)/ca

        # alpha_high <-> 90
        alpha1 = np.linspace(alpha_high, pi/2, nalpha)
        alpha1 = alpha1[1:]  # remove first element so as not to duplicate when concatenating
        cl1, cd1 = self.__Viterna(alpha1, 1.0)

        # 90 <-> 180-alpha_high
        alpha2 = np.linspace(pi/2, pi-alpha_high, nalpha)
        alpha2 = alpha2[1:]
        cl2, cd2 = self.__Viterna(pi-alpha2, -cl_adj)

        # 180-alpha_high <-> 180
        alpha3 = np.linspace(pi-alpha_high, pi, nalpha)
        alpha3 = alpha3[1:]
        cl3, cd3 = self.__Viterna(pi-alpha3, 1.0)
        cl3 = (alpha3-pi)/alpha_high*cl_high*cl_adj  # override with linear variation

        if alpha_low <= -alpha_high:
            alpha4 = []
            cl4 = []
            cd4 = []
            alpha5max = alpha_low
        else:
            # -alpha_high <-> alpha_low
            # Note: this is done slightly differently than AirfoilPrep for better continuity
            alpha4 = np.linspace(-alpha_high, alpha_low, nalpha)
            alpha4 = alpha4[1:-2]  # also remove last element for concatenation for this case
            cl4 = -cl_high*cl_adj + (alpha4+alpha_high)/(alpha_low+alpha_high)*(cl_low+cl_high*cl_adj)
            cd4 = cd_low + (alpha4-alpha_low)/(-alpha_high-alpha_low)*(cd_high-cd_low)
            alpha5max = -alpha_high

        # -90 <-> -alpha_high
        alpha5 = np.linspace(-pi/2, alpha5max, nalpha)
        alpha5 = alpha5[1:]
        cl5, cd5 = self.__Viterna(-alpha5, -cl_adj)

        # -180+alpha_high <-> -90
        alpha6 = np.linspace(-pi+alpha_high, -pi/2, nalpha)
        alpha6 = alpha6[1:]
        cl6, cd6 = self.__Viterna(alpha6+pi, cl_adj)

        # -180 <-> -180 + alpha_high
        alpha7 = np.linspace(-pi, -pi+alpha_high, nalpha)
        cl7, cd7 = self.__Viterna(alpha7+pi, 1.0)
        cl7 = (alpha7+pi)/alpha_high*cl_high*cl_adj  # linear variation

        alpha = np.concatenate((alpha7, alpha6, alpha5, alpha4, np.radians(self.alpha), alpha1, alpha2, alpha3))
        cl = np.concatenate((cl7, cl6, cl5, cl4, self.cl, cl1, cl2, cl3))
        cd = np.concatenate((cd7, cd6, cd5, cd4, self.cd, cd1, cd2, cd3))

        cd = np.maximum(cd, cdmin)  # don't allow negative drag coefficients


        # Setup alpha and cm to be used in extrapolation
        cm1_alpha = floor(self.alpha[0] / 10.0) * 10.0
        cm2_alpha = ceil(self.alpha[-1] / 10.0) * 10.0
        alpha_num = abs(int((-180.0-cm1_alpha)/10.0 - 1))
        alpha_cm1 = np.linspace(-180.0, cm1_alpha, alpha_num)
        alpha_cm2 = np.linspace(cm2_alpha, 180.0, int((180.0-cm2_alpha)/10.0 + 1))
        alpha_cm = np.concatenate((alpha_cm1, self.alpha, alpha_cm2))  # Specific alpha values are needed for cm function to work
        cm1 = np.zeros(len(alpha_cm1))
        cm2 = np.zeros(len(alpha_cm2))
        cm_ext = np.concatenate((cm1, self.cm, cm2))
        if np.count_nonzero(self.cm) > 0:
            cmCoef = self.__CMCoeff(cl_high, cd_high, cm_high)  # get cm coefficient
            cl_cm = np.interp(alpha_cm, np.degrees(alpha), cl)  # get cl for applicable alphas
            cd_cm = np.interp(alpha_cm, np.degrees(alpha), cd)  # get cd for applicable alphas
            alpha_low_deg = self.alpha[0]
            alpha_high_deg = self.alpha[-1]
            for i in range(len(alpha_cm)):
                cm_new = self.__getCM(i, cmCoef, alpha_cm, cl_cm, cd_cm, alpha_low_deg, alpha_high_deg)
                if cm_new is None:
                    pass  # For when it reaches the range of cm's that the user provides
                else:
                    cm_ext[i] = cm_new
        cm = np.interp(np.degrees(alpha), alpha_cm, cm_ext)
        return type(self)(self.Re, np.degrees(alpha), cl, cd, cm)




    def __Viterna(self, alpha, cl_adj):
        """private method to perform Viterna extrapolation"""

        alpha = np.maximum(alpha, 0.0001)  # prevent divide by zero

        cl = self.cdmax/2*np.sin(2*alpha) + self.A*np.cos(alpha)**2/np.sin(alpha)
        cl = cl*cl_adj

        cd = self.cdmax*np.sin(alpha)**2 + self.B*np.cos(alpha)

        return cl, cd

    def __CMCoeff(self, cl_high, cd_high, cm_high):
        """private method to obtain CM0 and CMCoeff"""

        found_zero_lift = False

        for i in range(len(self.cm)-1):
            if abs(self.alpha[i]) < 20.0 and self.cl[i] <= 0 and self.cl[i+1] >= 0:
                p = -self.cl[i] / (self.cl[i + 1] - self.cl[i])
                cm0 = self.cm[i] + p * (self.cm[i+1] - self.cm[i])
                found_zero_lift = True
                break

        if not found_zero_lift:
            p = -self.cl[0] / (self.cl[1] - self.cl[0])
            cm0 = self.cm[0] + p * (self.cm[1] - self.cm[0])
        self.cm0 = cm0
        alpha_high = radians(self.alpha[-1])
        XM = (-cm_high + cm0) / (cl_high * cos(alpha_high) + cd_high * sin(alpha_high))
        cmCoef = (XM - 0.25) / tan((alpha_high - pi/2))
        return cmCoef

    def __getCM(self, i, cmCoef, alpha, cl_ext, cd_ext, alpha_low_deg, alpha_high_deg):
        """private method to extrapolate Cm"""

        cm_new = 0
        if alpha[i] >= alpha_low_deg and alpha[i] <= alpha_high_deg:
            return
        if alpha[i] > -165 and alpha[i] < 165:
            if abs(alpha[i]) < 0.01:
                cm_new = self.cm0
            else:
                if alpha[i] > 0:
                    x = cmCoef * tan(radians(alpha[i]) - pi/2) + 0.25
                    cm_new = self.cm0 - x * (cl_ext[i] * cos(radians(alpha[i])) + cd_ext[i] * sin(radians(alpha[i])))
                else:
                    x = cmCoef * tan(-radians(alpha[i]) - pi/2) + 0.25
                    cm_new = -(self.cm0 - x * (-cl_ext[i] * cos(-radians(alpha[i])) + cd_ext[i] * sin(-radians(alpha[i]))))
        else:
            if alpha[i] == 165:
                cm_new = -0.4
            elif alpha[i] == 170:
                cm_new = -0.5
            elif alpha[i] == 175:
                cm_new = -0.25
            elif alpha[i] == 180:
                cm_new = 0
            elif alpha[i] == -165:
                cm_new = 0.35
            elif alpha[i] == -170:
                cm_new = 0.4
            elif alpha[i] == -175:
                cm_new = 0.2
            elif alpha[i] == -180:
                cm_new = 0
            else:
                print("Angle encountered for which there is no CM table value "
                      "(near +/-180 deg). Program will stop.")
        return cm_new

    def unsteadyparam(self, alpha_linear_min=-5, alpha_linear_max=5):
        """compute unsteady aero parameters used in AeroDyn input file

        Parameters
        ----------
        alpha_linear_min : float, optional (deg)
            angle of attack where linear portion of lift curve slope begins
        alpha_linear_max : float, optional (deg)
            angle of attack where linear portion of lift curve slope ends

        Returns
        -------
        aerodynParam : tuple of floats
            (control setting, stall angle, alpha for 0 cn, cn slope,
            cn at stall+, cn at stall-, alpha for min CD, min(CD))

        """

        alpha = np.radians(self.alpha)
        cl = self.cl
        cd = self.cd

        alpha_linear_min = radians(alpha_linear_min)
        alpha_linear_max = radians(alpha_linear_max)

        cn = cl*np.cos(alpha) + cd*np.sin(alpha)

        # find linear region
        idx = np.logical_and(alpha >= alpha_linear_min,
                             alpha <= alpha_linear_max)

        # checks for inppropriate data (like cylinders)
        if len(idx) < 10 or len(np.unique(cl)) < 10:
            return 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,

        # linear fit
        p = np.polyfit(alpha[idx], cn[idx], 1)
        m = p[0]
        alpha0 = -p[1]/m

        # find cn at stall locations
        alphaUpper = np.radians(np.arange(40.0))
        alphaLower = np.radians(np.arange(5.0, -40.0, -1))
        cnUpper = np.interp(alphaUpper, alpha, cn)
        cnLower = np.interp(alphaLower, alpha, cn)
        cnLinearUpper = m*(alphaUpper - alpha0)
        cnLinearLower = m*(alphaLower - alpha0)
        deviation = 0.05  # threshold for cl in detecting stall

        alphaU = np.interp(deviation, cnLinearUpper-cnUpper, alphaUpper)
        alphaL = np.interp(deviation, cnLower-cnLinearLower, alphaLower)

        # compute cn at stall according to linear fit
        cnStallUpper = m*(alphaU-alpha0)
        cnStallLower = m*(alphaL-alpha0)

        # find min cd
        minIdx = cd.argmin()

        # return: control setting, stall angle, alpha for 0 cn, cn slope,
        #         cn at stall+, cn at stall-, alpha for min CD, min(CD)
        return (0.0, degrees(alphaU), degrees(alpha0), m,
                cnStallUpper, cnStallLower, alpha[minIdx], cd[minIdx])

    def plot(self):
        """plot cl/cd/cm polar

        Returns
        -------
        figs : list of figure handles

        """
        import matplotlib.pyplot as plt

        p = self

        figs = []

        # plot cl
        fig = plt.figure()
        figs.append(fig)
        ax = fig.add_subplot(111)
        plt.plot(p.alpha, p.cl, label='Re = ' + str(p.Re/1e6) + ' million')
        ax.set_xlabel('angle of attack (deg)')
        ax.set_ylabel('lift coefficient')
        ax.legend(loc='best')

        # plot cd
        fig = plt.figure()
        figs.append(fig)
        ax = fig.add_subplot(111)
        ax.plot(p.alpha, p.cd, label='Re = ' + str(p.Re/1e6) + ' million')
        ax.set_xlabel('angle of attack (deg)')
        ax.set_ylabel('drag coefficient')
        ax.legend(loc='best')

        # plot cm
        fig = plt.figure()
        figs.append(fig)
        ax = fig.add_subplot(111)
        ax.plot(p.alpha, p.cm, label='Re = ' + str(p.Re/1e6) + ' million')
        ax.set_xlabel('angle of attack (deg)')
        ax.set_ylabel('moment coefficient')
        ax.legend(loc='best')

        return figs


    def _alpha_window_in_bounds(self,window):
        """ Ensures that the window of alpha values is within the bounds of alpha"""
        IBef=np.where(self.alpha<=window[0])[0]
        if len(IBef)>0:
            im=IBef[-1]
        else:
            im=0
        IAft=np.where(self.alpha>=window[1])[0]
        if len(IAft)>0:
            ip=IAft[0]
        else:
            ip=len(self.alpha)-1
        window=[self.alpha[im], self.alpha[ip]]
        return window


    def alpha0(self,window=[-20,20]):
        """ Finds alpha0, angle of zero lift """
        # Constant case or only one value
        if np.all(self.cl == self.cl[0]) or len(self.cl)==1:
            if self.cl[0]==0:
                return 0
            else:
                return np.nan

        # Ensuring window is within our alpha values
        window = self._alpha_window_in_bounds(window)
        
        # Finding zero up-crossing within window
        iwindow=np.where((self.alpha>=window[0]) & (self.alpha<=window[1]))
        alpha = self.alpha[iwindow]
        cl    = self.cl[iwindow]
        alpha_zc,i_zc = _zero_crossings(x=alpha,y=cl,direction='up')
        if len(alpha_zc)>1:
            raise Exception('Cannot find alpha0, {} zero crossings of Cl in the range of alpha values: [{} {}] '.format(len(alpha_zc),window[0],window[1]))
        elif len(alpha_zc)==0:
            raise Exception('Cannot find alpha0, no zero crossing of Cl in the range of alpha values: [{} {}] '.format(window[0],window[1]))

        alpha0=alpha_zc[0]
        return alpha0

    def cl_linear_slope(self,window=None):
        """ Find slope of linear region """

        # finding our alpha0 is possible
        alpha0 = self.alpha0()
        try:
            alpha0 = self.alpha0()
        except:
            alpha0=None

        # Constant case or only one value
        if np.all(self.cl == self.cl[0]) or len(self.cl)==1:
            return 0,alpha0,[self.alpha[0],self.alpha[-1]],[self.alpha[0],self.alpha[-1]]


        if window is None:
            UserWindow=False
            if np.nanmin(self.cl)>0 or np.nanmax(self.cl)<0 or alpha0 is None:
                window =[self.alpha[0],self.alpha[-1]]
            else:
                # define a window around alpha0
                window = [alpha0-5, alpha0+20];
        else:
            UserWindow=True

        # Ensuring window is within our alpha values
        window = self._alpha_window_in_bounds(window)

        # Selecting range of values within window
        idx = np.where((self.alpha >= window[0]) & (self.alpha <= window[1]) & ~np.isnan(self.cl))[0]
        cl, alpha = self.cl[idx], self.alpha[idx]
        if not UserWindow:
            # Selecting within the min and max of this window
            imin=np.where(cl==np.min(cl))[0][-1]
            idx = np.arange(imin,np.argmax(cl)+1)
            cl, alpha = cl[idx], alpha[idx]
        if len(cl)<=0:
            raise Exception('Cannot find Cl slope')
        elif len(cl)<4:
            if alpha0 is not None:
                sl=np.linalg.lstsq(alpha.reshape((-1,1))-alpha0,cl.reshape((-1,1)),rcond=None)[0][0]
                return sl[0],alpha0,[alpha[0],alpha[-1]],[alpha[0],alpha[-1]]
            else:
                p = np.polyfit(alpha, cl, 1)
                return p[0],p[1],[alpha[0],alpha[-1]],[alpha[0],alpha[-1]]

        def find_linear_region(x,y,nMin,x0=None):
            """ Find a linear region by computing all possible slopes
                The objective function tries to minimize the error with the linear slope
                and maximize the length of the linear region.
                If x0 is provided, the function a*(x-x0) is fitted
            """
            if x0 is not None:
                x = x.reshape((-1,1))-x0
                y = y.reshape((-1,1))
            n=len(x)-nMin+1
            err = np.zeros((n,n))*np.nan
            slp = np.zeros((n,n))*np.nan
            off = np.zeros((n,n))*np.nan
            spn = np.zeros((n,n))*np.nan
            for iStart in range(n):
                for j in range(iStart,n):
                    iEnd=j+nMin
                    if x0 is not None:
                        sl=np.linalg.lstsq(x[iStart:iEnd],y[iStart:iEnd],rcond=None)[0][0]
                        slp[iStart,j]= sl
                        off[iStart,j]= alpha0
                        y_lin = x[iStart:iEnd] * sl
                    else:
                        coefs = np.polyfit(x[iStart:iEnd],y[iStart:iEnd],1)
                        slp[iStart,j]= coefs[0]
                        off[iStart,j]= -coefs[1]/coefs[0]
                        y_lin = x[iStart:iEnd] * coefs[0] + coefs[1]
                    err[iStart,j]= np.mean((y[iStart:iEnd]-y_lin)**2)
                    spn[iStart,j]= iEnd-iStart
            spn=1/(spn-nMin+1)
            err=(err)/(np.nanmax(err)) 
            obj = np.multiply(spn,err) 
            obj=err
            (iStart,j) = np.unravel_index(np.nanargmin(obj), obj.shape)
            iEnd=j+nMin-1 # note -1 since we return the index here
            return slp[iStart,j],off[iStart,j],iStart,iEnd

        # Performing minimization of slope
        if UserWindow:
            nMin = len(alpha) # using full window
        else:
            nMin = max(3,round(0.6*len(alpha))) # optimizing between length and linear error  
        slope,off,iStart,iEnd = find_linear_region(alpha,cl,nMin,alpha0)

        # Slope around alpha 0
        if alpha0 is not None:
            ia0 = (np.abs(alpha-alpha0)).argmin()
            slope_0 = (cl[ia0+1]-cl[ia0-1])/(alpha[ia0+1]-alpha[ia0-1])
            if abs(slope-slope_0)/slope_0*100>20:
                print('Warning: More than 20% error between estimated slope ({:.4f}) and the slope around alpha0 ({:.4f}). The window for the slope search ([{} {}]) is likely wrong.'.format(slope,slope_0,window[0],window[-1]))

#             print('slope0',slope_0)
#         print('slope ',slope,' Alpha range: {:.3f} {:.3f} - nLin {}  nMin {}  nMax {}'.format(alpha[iStart],alpha[iEnd],len(alpha[iStart:iEnd+1]),nMin,len(alpha)))

        p = np.polyfit(alpha, cl, 1)
#         print('slpol ',p[0],p[1],alpha0)


        # Making sure the linear fit is bigger than ClMax
        DeltaAlpha=alpha[iEnd]-alpha[iStart]
        DeltaClMax = ((alpha[-1]-off)*slope)-cl[-1]
        if DeltaClMax<0:
            print('NEED TO ADJUST SLOPE FOR CL MAX',DeltaClMax) 

        LinWindow    = [alpha[iStart],alpha[iEnd]]
        SearchWindow = [alpha[0],alpha[-1]]

        return slope,off,LinWindow,SearchWindow


def _zero_crossings(y,x=None,direction=None):
    """
      Find zero-crossing points in a discrete vector, using linear interpolation.

      direction: 'up' or 'down', to select only up-crossings or down-crossings

      returns: 
          x values xzc such that y(yzc)==0
          indexes izc, such that the zero is between y[izc] (excluded) and y[izc+1] (included)

      if direction is not provided, also returns:
              sign, equal to 1 for up crossing
    """
    if x is None:
        x=np.arange(len(y))

    if np.any((x[1:] - x[0:-1]) <= 0.0):
        raise Exception('x values need to be in ascending order')

    # Indices before zero-crossing
    iBef = np.where(y[1:]*y[0:-1] < 0.0)[0]
    
    # Find the zero crossing by linear interpolation
    xzc = x[iBef] - y[iBef] * (x[iBef+1] - x[iBef]) / (y[iBef+1] - y[iBef])
    
    # Selecting points that are exactly 0 and where neighbor change sign
    iZero = np.where(y == 0.0)[0]
    iZero = iZero[np.where((iZero > 0) & (iZero < x.size-1))]
    iZero = iZero[np.where(y[iZero-1]*y[iZero+1] < 0.0)]

    # Concatenate 
    xzc  = np.concatenate((xzc, x[iZero]))
    iBef = np.concatenate((iBef, iZero))

    # Sort
    iSort = np.argsort(xzc)
    xzc, iBef = xzc[iSort], iBef[iSort]

    # Return up-crossing, down crossing or both
    sign = np.sign(y[iBef+1]-y[iBef])
    if direction == 'up':
        I= np.where(sign==1)[0]
        return xzc[I],iBef[I]
    elif direction == 'down':
        I= np.where(sign==-1)[0]
        return xzc[I],iBef[I]
    elif direction is not None:
        raise Exception('Direction should be either `up` or `down`')
    return xzc, iBef, sign


if __name__ == '__main__':
    pass
