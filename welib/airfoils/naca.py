"""

NACA4:
- First digit describing maximum camber as percentage of the chord.
- Second digit describing the distance of maximum camber from the airfoil leading edge in tenths of the chord.
- Last two digits describing maximum thickness of the airfoil as percent of the chord.[3]



References:
  [1] E. Jacobs, K. Ward, R. Pinkerton (1933) "The characteristics of 78 related airfoil sections from tests in the variable-density wind tunnel", NACA Report 460, 1933.

"""
import numpy as np


def naca_shape(digits, chord=1, n=151, sharp=False, pitch=0, xrot=0.25):
    """ 
    INPUTS:
     - digits: 4 digits string, e.g. '0012'
     - n     : number of points on the 
     - sharp: if true, forces the values at x=1 to be y=0
                 The original NACA equations gives a non zero thickness at the trailing edge
     - pitch: pitch angle, positive nose up (angle of attack) [rad]
     - xrot: center of rotation for pitch (in chord coordinates) [-]

    """
    if len(digits)!=4:
        raise NotImplementedError()

    maxCamb      = int(digits[0])
    d_LE_maxCamb = int(digits[1])*10
    t            = int(digits[2:4])/100  # Maximum thickness

    m = maxCamb/100.0
    p = d_LE_maxCamb/100.0

    x = np.linspace(0, 1, n)

    if not sharp:
        # Original NACA equation - the TE is not blunt
        y_t = 5 * t * (0.2969*np.sqrt(x) + ((((- 0.1015 )*x + 0.2843 )*x - 0.3516)*x - 0.1260)*x)
    else:
        # Small modification to ensure the trailing edge is sharp
        y_t = 5 * t * (0.2969*np.sqrt(x) + ((((- 0.1036 )*x + 0.2843 )*x - 0.3516)*x - 0.1260)*x)

    if m == 0 or p == 0:
        yc = np.zeros_like(x)
        dyc = np.zeros_like(x)
    else:
        yc = np.zeros_like(x)
        dyc = np.zeros_like(x)
        idx = x <= p
        yc[idx] = (m/p**2)*(2*p*x[idx] - x[idx]**2)
        dyc[idx] = (2*m/p**2)*(p - x[idx])
        idx = ~idx
        yc[idx] = (m/(1-p)**2)*((1-2*p) + 2*p*x[idx] - x[idx]**2)
        dyc[idx] = (2*m/(1-p)**2)*(p - x[idx])

    theta = np.arctan(dyc)

    xu = x - y_t*np.sin(theta)
    yu = yc + y_t*np.cos(theta)
    xl = x + y_t*np.sin(theta)
    yl = yc - y_t*np.cos(theta)

    if not sharp:
        xa = np.concatenate((xu, np.flip(xl, 0)))
        ya = np.concatenate((yu, np.flip(yl, 0)))
    else:
        xa = np.concatenate((xu, np.flip(xl[:-1], 0)))
        ya = np.concatenate((yu, np.flip(yl[:-1], 0)))

    # --- Rotate
    x = (xa-xrot)*np.cos(pitch) +        ya*np.sin(pitch) + xrot
    y =        ya*np.cos(pitch) - (xa-xrot)*np.sin(pitch)

    # --- Scale
    x *= chord
    y *= chord

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


