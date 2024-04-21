"""

NACA4:
- First digit describing maximum camber as percentage of the chord.
- Second digit describing the distance of maximum camber from the airfoil leading edge in tenths of the chord.
- Last two digits describing maximum thickness of the airfoil as percent of the chord.[3]



References:
  [1] E. Jacobs, K. Ward, R. Pinkerton (1933) "The characteristics of 78 related airfoil sections from tests in the variable-density wind tunnel", NACA Report 460, 1933.

"""
import numpy as np


def naca_shape(digits, chord=1, n=151, thickTEZero=False, pitch=0, xrot=0.25):
    """ 
    INPUTS:
     - digits: 4 digits string, e.g. '0012'
     - thickTEZero: if true, forces the values at x=1 to be y=0
                 The original NACA equations gives a non zero thickness at the trailing edge
     - pitch: pitch angle, positive nose up (angle of attack) [rad]
     - xrot: center of rotation for pitch (in chord coordinates) [-]

    """
    if len(digits)!=4:
        raise NotImplementedError()

    maxCamb      = int(digits[0])
    d_LE_maxCamb = int(digits[1])*10
    t            = int(digits[2:4])/100  # Maximum thickness

    if maxCamb!=0:
        raise NotImplementedError()
    if d_LE_maxCamb!=0:
        raise NotImplementedError()

    # --- Symmetric airfoils
    x = np.linspace(0, 1, n)
    if not thickTEZero:
        # Original NACA equation - the TE is not sharp
        y = 5 * t * (0.2969*np.sqrt(x) + ((((- 0.1015 )*x + 0.2843 )*x - 0.3516)*x - 0.1260)*x)
        # NOTE: the TE point will be repeated
        xa = np.concatenate((x, np.flip( x     ,0)))
        ya = np.concatenate((y, np.flip(-y     ,0)))
    else:
        # Small modifications to ensure the thickness is zero at the TE
        y = 5 * t * (0.2969*np.sqrt(x) + ((((- 0.1036 )*x + 0.2843 )*x - 0.3516)*x - 0.1260)*x)
        # NOTE: for a sharp trailing edge, we avoid repetition of the TE point
        xa = np.concatenate((x, np.flip( x[:-1],0)))
        ya = np.concatenate((y, np.flip(-y[:-1],0)))


    # --- Rotate
    x = (xa-xrot)*np.cos(pitch) +        ya*np.sin(pitch) + xrot
    y =        ya*np.cos(pitch) - (xa-xrot)*np.sin(pitch)

    # --- Scale
    x*=chord
    y*=chord

    return x, y


if __name__ == '__main__':
    import matplotlib.pyplot as plt

    digits='0022'
    n = 151
    x,y= naca_shape(digits,n=n)
    xu= x[:n]
    xl= x[-1:n-1:-1]
    yu= y[:n]
    yl= y[-1:n-1:-1]
#     print(xu)
#     print(xl)
#     print(x)
#     print(y)

    dxu = xu-xu[-1]
    dyu = yu-yu[-1]
    dxl = xl-xl[-1]
    dyl = yl-yl[-1]

    au = np.arctan(dyu[:-1]/dxu[:-1])*180/np.pi
    al = np.arctan(dyl[:-1]/dxl[:-1])*180/np.pi
    #print('au',au[:-1])
    #print('al',al[:-1])
    alpha=np.mean(-au[-6:]+al[-6:])
    print(np.mean(-au[-6:]))
    print(np.mean( al[-6:]))
    print(alpha)
    fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
    fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
    ax.plot(xu[:-1],-au[:], '-', label='')
    ax.plot(xl[:-1], al[:], ':', label='')
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.legend()



    fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
    fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
    ax.plot(x, y, label='')

    x_=np.linspace(-1,0,10)
    ax.plot(x_+1,  np.tan(alpha/2*np.pi/180)*x_ + yl[-1], 'k--')
    ax.plot(x_+1, -np.tan(alpha/2*np.pi/180)*x_ + yu[-1], 'k--')



    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.set_title('{}'.format(digits))
    plt.axis ( 'equal' )

    #writeToFile('NACA{}.csv'.format(digits),x,y)

    plt.show()


