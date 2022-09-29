"""
Tools to work with OLAF the vortex code implemented in openfast
"""
import numpy as np



def OLAFParams(omega_rpm, U0, R, a=0.3, aScale=1.5,
          deltaPsiDeg=6, nPerRot=None,
          targetWakeLengthD=6,
          nNWrot=8, nFWrot=2, nFWrotFree=1,
          verbose=True, dt_glue_code=None):
    """ 
    Computes recommended time step and wake length based on the rotational speed in RPM

    INPUTS:
     - omega_rpm: rotational speed [RPM]
     - U0: mean wind speed [m/s]
     - R: rotor radius [m]

    OPTIONS FOR TIME STEP:
      - either:
         - deltaPsiDeg : target azimuthal discretization [deg]
              or
         - nPerRot     : number of time step per rotations.
                deltaPsiDeg  -  nPerRot
                     5            72    
                     6            60    
                     7            51.5  
                     8            45    
     - dt_glue_code: glue code time step. If provided, the time step of OLAF will be approximated
                     such that it is a multiple of the glue-code time step.

    OPTIONS FOR WAKE LENGTH:
     - a: average axial induction factor at the rotor [-]
     - aScale: scaling factor to estimate induction, such that the wake convection velocity is:
               Uc=U0(1-aScale*a)
     - targetWakeLengthD: target wake length in diameter [D]
     - nNWrot     : minimum number of near wake rotations
     - nFWrot     : minimum number of far wake rotations
     - nFWrotFree : minimum number of far wake rotations (free panels)

    """
    def myprint(*args, **kwargs):
        if verbose:
            print(*args, **kwargs)

    # Rotational period
    omega_rpm = np.asarray(omega_rpm)
    omega = omega_rpm*2*np.pi/60
    T = 2*np.pi/omega

    # Desired time step
    if nPerRot is not None:
        dt_wanted    = np.around(T/nPerRot,5)
        deltaPsiDeg  = np.around(omega*dt_wanted*180/np.pi ,2)
    else:
        dt_wanted    = np.around(deltaPsiDeg/(6*omega_rpm),5)
        nPerRot = int(2*np.pi /(deltaPsiDeg*np.pi/180))

    # Adapting desired time step based on glue code time step
    if dt_glue_code is not None:
        dt_rounded = round(dt_wanted/dt_glue_code)*dt_glue_code
        deltaPsiDeg2 = np.around(omega*dt_rounded *180/np.pi ,2)
        myprint('>>> To satisfy glue-code dt:')
        myprint('    Rounding dt   from {} to {}'.format(dt_wanted, dt_rounded    ))
        myprint('    Changing dpsi from {} to {}'.format(deltaPsiDeg, deltaPsiDeg2))
        dt_fvw   = dt_rounded
        deltaPsiDeg = deltaPsiDeg2
        nPerRot = int(2*np.pi /(deltaPsiDeg*np.pi/180))
    else:
        dt_fvw = dt_wanted

    # Wake length from mean wind speed
    targetWakeLength = targetWakeLengthD * 2 * R
    Uc = U0 * (1-aScale*a)
    nPanels_FromU0 = int(targetWakeLength / (Uc*dt_fvw))

    # Wake length from rotational speed and number of rotations
    nNWPanel_FromRot     = int(nNWrot*nPerRot)
    nFWPanel_FromRot     = int(nFWrot*nPerRot)
    nFWPanelFree_FromRot = int(nFWrotFree*nPerRot)
    nPanels_FromRot      = nNWPanel_FromRot +  nFWPanel_FromRot

    # Below we chose between criteria on number of rotation or donwstream distance
    # This can be adapted/improved
    myprint('Number of panels from wind speed and distance:{:15d}'.format(nPanels_FromU0))
    myprint('Number of panels from number of rotations    :{:15d}'.format(nPanels_FromRot))
    if nPanels_FromRot>nPanels_FromU0:
        # Criteria based on rotation wins: 
        myprint('[INFO] Using number of rotations to setup number of panels')
        nNWPanel     = nNWPanel_FromRot
        nFWPanel     = nFWPanel_FromRot
        nFWPanelFree = nFWPanelFree_FromRot
    else:
        myprint('[INFO] Using wind speed and distance to setup number of panels')
        # Wake distance wins, we keep the nFW from rot but increase nNW
        nPanels      = nPanels_FromU0
        nFWPanel     = nFWPanel_FromRot
        nFWPanelFree = nFWPanelFree_FromRot
        nNWPanel     = nPanels - nFWPanel # increase nNW

    # Recompute
    nNWrot     = nNWPanel    *dt_fvw/T
    nFWrot     = nFWPanel    *dt_fvw/T
    nFWrotFree = nFWPanelFree*dt_fvw/T

    wakeLengthRot = nPanels * dt_fvw/T
    wakeLengthEst = (nPanels * dt_fvw * Uc)/(2*R)


    # Transient time (twice the time to develop the full wake extent)
    # This is the minimum recommended time before convergence of the wake is expected 
    # (might be quite long)
    tMax = 2 * dt_fvw*nPanels

    if verbose:
        myprint('Wake extent - panels: {:d} - est. distance {:.1f}D - {:5.1f} rotations'.format(nPanels, wakeLengthEst, wakeLengthRot))
        myprint('Transient time: {:.6f} ({:.1f} rotations)'.format(tMax, tMax/T))
        myprint('')
        myprint('OLAF INPUT FILE:')
        myprint('----------------------- GENERAL OPTIONS ---------------------')
        myprint('{:15.6f} DT_FVW       (delta psi = {:5.1f})'.format(dt_fvw, deltaPsiDeg))
        myprint('--------------- WAKE EXTENT AND DISCRETIZATION --------------')
        myprint('{:15d} nNWPanel     ({:5.1f} rotations)'.format(nNWPanel, nNWrot))
        myprint('{:15d} nFWPanel     ({:5.1f} rotations) (previously called WakeLength)'.format(nFWPanel, nFWrot))
        myprint('{:15d} nFWPanelFree ({:5.1f} rotations) (previously called FreeWakeLength)'.format(nFWPanelFree, nFWrotFree))

    return dt_fvw, tMax, nNWPanel, nFWPanel, nFWPanelFree



def OLAFParamsRPM(omega_rpm, deltaPsiDeg=6, nNWrot=2, nFWrot=10, nFWrotFree=3, nPerRot=None,  verbose=True, dt_glue_code=None):
    """ 
    Computes recommended time step and wake length based on the rotational speed in RPM

    INPUTS:
     - omega_rpm: rotational speed in RPM
     - deltaPsiDeg : azimuthal discretization in deg
     - nNWrot : number of near wake rotations
     - nFWrot : total number of far wake rotations
     - nFWrotFree : number of far wake rotations that are free

        deltaPsiDeg  -  nPerRot
             5            72    
             6            60    
             7            51.5  
             8            45    
    """
    omega_rpm = np.asarray(omega_rpm)
    omega = omega_rpm*2*np.pi/60
    T = 2*np.pi/omega
    if nPerRot is not None:
        dt_wanted    = np.around(T/nPerRot,5)
        deltaPsiDeg  = np.around(omega*dt_wanted*180/np.pi ,2)
    else:
        dt_wanted    = np.around(deltaPsiDeg/(6*omega_rpm),5)
        nPerRot = int(2*np.pi /(deltaPsiDeg*np.pi/180))
    if dt_glue_code is not None:
        dt_rounded = round(dt_wanted/dt_glue_code)*dt_glue_code
        deltaPsiDeg2 = np.around(omega*dt_rounded *180/np.pi ,2)
        print('>>> To satisfy glue-code dt:')
        print('    Rounding dt   from {} to {}'.format(dt_wanted, dt_rounded    ))
        print('    Changing dpsi from {} to {}'.format(deltaPsiDeg, deltaPsiDeg2))
        dt_wanted   = dt_rounded
        deltaPsiDeg = deltaPsiDeg2
        nPerRot = int(2*np.pi /(deltaPsiDeg*np.pi/180))

    nNWPanel     = int(nNWrot*nPerRot)
    nFWPanel     = int(nFWrot*nPerRot)
    nFWPanelFree = int(nFWrotFree*nPerRot)

    if totalRot is None:
        totalRot = (nNWrot + nFWrot)*3 # going three-times through the entire wake

    tMax = dt_wanted*nPerRot*totalRot

    if verbose:
        print(dt_wanted              , '  dt')
        print(int      (nNWPanel    ), '  nNWPanel          ({} rotations)'.format(nNWrot))
        print(int      (nFWPanel    ), '  FarWakeLength     ({} rotations)'.format(nFWrot))
        print(int      (nFWPanelFree), '  FreeFarWakeLength ({} rotations)'.format(nFWrotFree))
        print(tMax              , '  Tmax ({} rotations)'.format(totalRot))

    return dt_wanted, tMax, nNWPanel, nFWPanel, nFWPanelFree






if __name__ == '__main__':
    OLAFParams(omega_rpm = 4.87558, U0=8, R=63, deltaPsiDeg=6, verbose=True)
