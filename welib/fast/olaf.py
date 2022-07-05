"""
Tools to work with OLAF the vortex code implemented in openfast
"""
import numpy as np


def OLAFParams(omega_rpm, deltaPsiDeg=6, nNWrot=2, nFWrot=10, nFWrotFree=3, nPerRot=None, totalRot=None, show=True, dt_glue_code=None):
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

    if show:
        print(dt_wanted              , '  dt')
        print(int      (nNWPanel    ), '  nNWPanel          ({} rotations)'.format(nNWrot))
        print(int      (nFWPanel    ), '  FarWakeLength     ({} rotations)'.format(nFWrot))
        print(int      (nFWPanelFree), '  FreeFarWakeLength ({} rotations)'.format(nFWrotFree))
        print(tMax              , '  Tmax ({} rotations)'.format(totalRot))

    return dt_wanted, tMax, nNWPanel, nFWPanel, nFWPanelFree


if __name__ == '__main__':
    OLAFParams(omega_rpm = 4.87558, deltaPsiDeg=6, show=True)
