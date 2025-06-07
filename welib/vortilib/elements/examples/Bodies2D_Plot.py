""" 
Velocity field about the rankine nose (point source + freestream)
Note: methods are stored in SourcePoint package

Reference: 
    [1] Branlard (2017) Wind turbine aerodynamics and vorticity-based methods, Springer

"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from welib.tools.curves import streamQuiver
from welib.vortilib.elements.SourcePoint import *
from welib.vortilib.elements.RankineBody import *
from welib.vortilib.elements.VortexCylinder2D import *

def main(whichPlots=[0]):
    # --- Compare Rankine nose and cylinder
    U0=1
    R=2
    Gamma=0
    alpha=0*np.pi/180

    Sigma = 2*U0*R
    x0 = Sigma/(2*np.pi*U0)

    # --- Cylinder
    x_cyl  = np.linspace(-R, R, 200)
    y_cyl  = vc_coords_yx(x_cyl, R)
    Cp_cyl = vc_Cpx(x_cyl, U0 = U0, R = R, Gamma = Gamma, alpha = alpha)

    # --- Rankine nose
    x_rn = np.linspace(-x0,6*R+x0-R, 200)
    y_rn = rn_coords_yx(x_rn, U0=U0, Sigma=Sigma)
    x_rn2 = x_rn+x0-R
    Cp_rn = rn_Cpx(x_rn, U0=U0, Sigma=Sigma)

    # --- Rankine oval
    Sigma, b = ro_params(R=2*R, t=R, U0=U0) #, method='ajkj')
    x_ro = np.linspace(-2*R, 2*R, 200)
    y_ro = ro_coord_yx(x_ro, b=b, U0=U0, Sigma=Sigma)
    x_ro2 = x_ro+2*R/2
    Cp_ro = ro_Cpx(x_ro, b=b, U0=U0, Sigma=Sigma)


    fig,axes = plt.subplots(2, 1, sharex=True, figsize=(6.4,4.8)) # (6.4,4.8)
    fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
    ax=axes[0]
    ax.plot(x_cyl, Cp_cyl, label='Cylinder')
    ax.plot(x_rn2, Cp_rn, label='Rankine nose')
    ax.plot(x_ro2, Cp_ro, label='Rankine oval')
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.legend()

    ax=axes[1]
    ax.plot(x_cyl, y_cyl, label='Cylinder')
    ax.plot(x_rn2, y_rn, label='Rankine nose')
    ax.plot(x_ro2, y_ro, label='Rankine oval')
    ax.set_aspect('equal','box')

#     ax.legend()


if __name__ == '__main__':
    main(whichPlots=[0])
    plt.show()
if __name__=="__test__":
    main()
    pass
if __name__=="__export__":
    main(whichPlots=[0])
    from welib.tools.repo import export_figs_callback
    export_figs_callback(__file__)
