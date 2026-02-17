"""
Compute velocity field in a cavity using discontinuous constant source panel method

Cavity looks as follows:

    ----------
    |         |
    |          outlet
    |         |
inlet         |
    |         |
    -----------
"""

import matplotlib.pyplot as plt
from welib.essentials import *

from welib.vortilib.panelcodes.panel_tools import *

import welib.vortilib.panelcodes.pointSources_discontinuous as dps
import welib.vortilib.panelcodes.constantSourcePanels_discontinuous as dcsp

from welib.CFD.flows2D import *    
from welib.tools.curves import curve_interp

np.set_printoptions(linewidth=300, precision=3)


def cavity_geometry(x1=0, x2=1, y1=0, y2=1, yin=0.1, yout=0.7, hin=0.5, hout=0.3,
                    nin=1, nout=1, ds=0.1, plot=False, eps=0.0):
    """ 
    INPUTS:
     - x1, x2: Horizontal limits of cavity
     - y1, y2: Vertical limits of cavity
     - hin, hout: heights of inlet and outlets
    """
    upper_0 = np.array([
            [x1, yin+hin],
            [x1, y2],
            [x2, y2],
            [x2, yout+hout]])
    lower_0 = np.array([
            [x2, yout],
            [x2, y1],
            [x1, y1],
            [x1, yin]])
    inlet_0  = np.array([[x1+eps, yin],      [x1+eps,yin+hin]])
    outlet_0 = np.array([[x2-eps, yout+hout],[x2-eps,yout]])

    # Remesh
    upper  = curve_interp(line=upper_0, ds = ds, keepOri = True)
    lower  = curve_interp(line=lower_0, ds = ds, keepOri = True)
    inlet  = curve_interp(line=inlet_0, n = nin+1)
    outlet = curve_interp(line=outlet_0, n = nout+1)

    if plot:
        fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8))
        fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
        ax.plot(upper_0[:,0], upper_0[:,1],'o')
        ax.plot(upper  [:,0], upper  [:,1],'+-', label='Upper')
        ax.plot(lower_0[:,0], lower_0[:,1],'o')
        ax.plot(lower  [:,0], lower  [:,1],'+-', label='Lower')
        ax.plot(inlet  [:,0], inlet  [:,1],'d--', label='Inlet' , ms=3)
        ax.plot(outlet [:,0], outlet [:,1],'d--', label='Outlet', ms=3)
        ax.set_xlabel('')
        ax.set_ylabel('')
        ax.legend()
        ax.set_aspect('equal','box')
    p = {'hin':hin, 'hout':hout}

#     upper = np.asarray(opposite_contour(upper[:,0], upper[:,1])).T
#     lower = np.asarray(opposite_contour(lower[:,0], lower[:,1])).T
    # TODO loop from the start in the right direction
    upper = upper[-1::-1,:]
    lower = lower[-1::-1,:]
    return upper, lower, inlet, outlet, p


if __name__ == '__main__':

    #upper, lower, inlet, outlet = cavity_geometry(plot=False, ds=0.01, yin=0.0, nin=5, nout=5)
    #upper, lower, inlet, outlet, p = cavity_geometry(plot=False, ds=1.0, yin=0.0,  nin=5, nout=5, yout=0, hin=1, hout=1)
    #upper, lower, inlet, outlet, p = cavity_geometry(plot=False, ds=0.05, yin=0.0, nin=1, nout=10, yout=0, hin=1, hout=1)
    upper, lower, inlet, outlet, p = cavity_geometry(plot=False, ds=0.05, yin=0.0, nin=1, nout=10, yout=0.5, hin=0.5, hout=0.3)
    SP_k = np.vstack([inlet, outlet])



    Vin = 10
    Vout = Vin*p['hin']/p['hout'] # Conservation of mass


    SP1_ki = inlet[0:-1,:]
    SP2_ki = inlet[1:,:]
    sigmas_ki = np.asarray([2*Vin]*(len(inlet)-1))
    SP1_ko = outlet[0:-1,:]
    SP2_ko = outlet[1:,:]
    sigmas_ko = np.asarray([-2*Vout]*(len(outlet)-1))


    SP1_k = np.vstack([inlet[0:-1,:], outlet[0:-1,:]])
    SP2_k = np.vstack([inlet[1:,:], outlet[1:,:]])
    sigmas_k = 2*np.concatenate([ [Vin]*(len(inlet)-1), [-Vout]*(len(outlet)-1)])
    print('Vout', Vin, Vout)

    #ax = plot_line(upper[:,0], upper[:,1] )
    #plot_line(lower[:,0], lower[:,1] , ax=ax)
    #fU = lambda X,Y: dcsp_u(X, Y, SP1_k, SP2_k, sigmas_k)
    fUcst = lambda X,Y: (X*0+Vin, X*0)
    fU = fUcst

    SP1 = np.vstack([upper[0:-1,:], lower[0:-1,:]])
    SP2 = np.vstack([upper[1:,:]  , lower[1:,:]])

#     print('SP_k\n', np.column_stack((SP1_k,SP2_k, sigmas_k.reshape((-1,1)))))
#     print('SP\n',   np.column_stack((SP1  ,SP2  )))


    n_hat, t_hat, mids, ds, _, _ = line_params2(SP1,SP2)


    # --- Panel method
    method='SP'
    method='DCSP'
    if method=='SP':
        out = dps.dPS_solve(SP1, SP2, fU=fU, offset=1)
        CP = out['CP']
        Vpnl = np.asarray( dps.PS_velocity(CP[:,0], CP[:,1], out['SP'], out['sigmas']))
        vel = lambda X, Y : dps.PS_velocity(X, Y, out['SP'], out['sigmas'])
        # TODO
        V_k1 = np.asarray( dcsp.dcsp_velocity(CP[:,0], CP[:,1], SP1_k, SP2_k, sigmas_k) ) # TODO

    else:
        out = dcsp.dCSP_solve(SP1, SP2, fU=fU)
        CP = out['CP']
        Vpnl = np.asarray( dcsp.dcsp_u(CP[:,0], CP[:,1], SP1, SP2, out['sigmas'], debug=False))
        vel = lambda X, Y : dcsp.dcsp_u(X, Y, SP1, SP2, out['sigmas'])

        V_k1 = np.asarray( dcsp.dcsp_u(CP[:,0], CP[:,1], SP1_k, SP2_k, sigmas_k) )
        V_k2 = np.asarray(  fU( CP[:,0], CP[:,1] ))

    Vtot = Vpnl+V_k1
#     printMat(CP,'CP')
#     printMat(Vtot,'Vtot')
#     printMat(Vpnl,'Vpnl')
#     printMat(V_k1,'V_k1')
#     printMat(V_k2,'V_k2')
#     printMat(out['sigmas'], 'sigmas')

    # --- Flow field
    fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8))
    fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)

#     ys = np.linspace(0,1,20)
#     xs = [0.5]*len(ys)
    ys = np.linspace(inlet[0,1],inlet[-1,1],21)
    xs = [0.01]*len(ys)
#     xs=None
#     ys=None
    streamStart=True
    streamStart=False

    X, Y, U, V =  flowfield2D(vel, xmin=0, xmax=1, ymin=0, ymax=1, nx=350, fU=fU, rel=False)
    #vel = lambda X, Y : (X*0,X*0)
#     X, Y, U, V =  flowfield2D(vel, xmin=-0.1, xmax=1.1, ymin=-0.1, ymax=1.1, nx=350, fU=fU, rel=False)
    ax =  flowfield2D_plot(X, Y, U, V, ax=ax, minVal=0, maxVal=max(Vin,Vout), bounded=False, rel=False, xs=xs, ys=ys, streamStart=streamStart)

#     ax.plot(upper  [:,0], upper  [:,1],'+-', label='Upper')
#     ax.plot(lower  [:,0], lower  [:,1],'+-', label='Lower')
    ax.plot(inlet  [:,0], inlet  [:,1],'d--', label='Inlet' , ms=3)
    ax.plot(outlet [:,0], outlet [:,1],'d--', label='Outlet', ms=3)

    
    #plot_line2(SP1, SP2, Uwall=Vtot.T, ntScale=0.1, UScale=0.1, nt=True, ax=ax)
    plot_line2(SP1, SP2, Uwall=out['Vwall'], ntScale=0.3, UScale=0.05, nt=True, ax=ax)

    ax.set_xlim([-0.1, 1.1])
    ax.set_ylim([-0.1, 1.1])

    ax.set_title('Source panel - Cavity')


    plt.show()

