""" 
Example using "yams numeric resursive"

One (uniform) beam is splitted into several sub-beams.
An eigenvalue analysis is run to look at the frequencies and modes, comparing them to the analytical
modes for a uniform beam.

Top mass 

"""

import numpy as np
import matplotlib.pyplot as plt
from welib.tools.clean_exceptions import *
import welib.beams.theory as bt
from welib.tools.strings import printVec
from welib.yams.yams_rec import YAMSRecUniformBeamBody, YAMSRecGroundBody, YAMSRecRigidBody
from welib.yams.utils import buildRigidBodyMassMatrix

def main(test=False, useNac=False):
    nBeams_twr  = 2 # Number of times we split the tower.
    nShapes_twr = 6 # Number of shapes per sub beam
    nModesPlot  = 6
    norm = 'tip' # 'max', 'tip'  Method to normalize modes at the end
    useAnalyticalShapeFunctions = False # <<< TODO implement this
    shapeFunctions='admissible'
    #shapeFunctions='masslessbeam'
    # shapeFunctions='Guyan'

    nSpan_twr         = 101
    bInit             = False # Use some default initial conditions (will change M&K and pos).
    Mtop              = 10.0e4     # Not sure if this works well for now
    bStiffening       = True  # Not tried yet
    bAxialCorr        = False # Not tried yet
    gravity = 9.81

    L_twr  = 100
    EI_twr = 2.817E+11
    m_twr  = 1.724E+03


    # --- Analytical modes
    x_th = np.linspace(0, L_twr, nSpan_twr*nBeams_twr);
    freq_th,_,U_th,_,_ = bt.UniformBeamBendingModes('unloaded-topmass-clamped-free', EI_twr, rho=1*m_twr, A=1, L=L_twr, x=x_th, Mtop=Mtop, nModes=nModesPlot, norm=norm)


    # --- Model Setup
    # Creating reference frame
    grd = YAMSRecGroundBody()
    # Bodies
    twrs=[]
    L_sub = L_twr/nBeams_twr  # sub beam length
    M_sub = L_sub*m_twr       # sub beam mass
    if not test:
        print('Beam     length {:10.1f}m, mass {:10.1f}kg'.format(L_twr, L_twr*m_twr))
        print('Sub Beam length {:10.1f}m, mass {:10.1f}kg'.format(L_sub, M_sub))
    for itwr in np.arange(nBeams_twr):
        Mtop_sub = Mtop + (nBeams_twr-itwr-1)*M_sub 
        if not test:
            print('Beam {} has {:10.1f}kg on top'.format(itwr+1, Mtop_sub))
        # Internal sub-beams are attached to another beam at their tip.
        # Use an interface-admissible basis there, and clamped-free only for the last segment.
        if shapeFunctions=='admissible':
            topBC = 'interface' if itwr < nBeams_twr-1 else 'free'
        else:
            topBC = 'free'
        concentrated_inertias = None
        if not useNac:
            if itwr==nBeams_twr-1:
                MM = buildRigidBodyMassMatrix(Mtop, (0,0,0), COG=(0,0,0))
                concentrated_inertias=[{'iNode':-1, 'MM':MM}]

        twr= YAMSRecUniformBeamBody('T'+str(itwr+1), nShapes_twr, nSpan_twr, L_sub, EI_twr, m_twr, 
                                    Mtop=Mtop_sub,                               # For stiffness
                                    bottomBC='clamped',
                                    topBC=topBC,
                                    shapeFunctions=shapeFunctions,
                                    concentrated_inertias=concentrated_inertias, # For Mass matrix
                                    bAxialCorr=bAxialCorr, main_axis='x', bStiffening=bStiffening, gravity=gravity)
        twrs.append(twr)


    if useNac:
        nac = YAMSRecRigidBody(name='Nac', mass=Mtop, J=0, rho_G=None, s_OP=None, r_O=None, R_b2g=None, sympy=False)


    # Connection between bodies
    grd.connectTo(twrs[0], Point=None, Type='Rigid')
    for twr_prev, twr_next in zip(twrs[0:-1], twrs[1:]):
        twr_prev.connectTo(twr_next, BodyPoint='LastPoint', Type='Rigid', RelOrientation = np.eye(3), OrientAfter=True)
    if useNac:
        twrs[-1].connectTo(nac, BodyPoint='LastPoint', Type='Rigid')



    # Init of DOF
    nDOF = grd.setupDOFIndex()
    q=np.zeros(nDOF)
    if bInit:
        q[:] = np.arange(len(q))+1 # Define some kind of initial conditions
    #     q[1] = 0


    # --- Apply a given DOF vector
    # Set degree of freedom value, and update kinematics associated with dofs
    grd.setDOF(q)

    # print(twrs[0])
    # print(twrs[1])
    # print(twrs[2])
    # print(grd)
    # print('>>q',q)

    # Retrieve global position of bodies (including mean line of flexible bodies)
    pos = grd._all_positions_global
    # x= pos[0,:]
    # print(x)
    # fig,ax = plt.subplots(1,1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
    # fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
    # ax.plot(pos[2,1:], pos[0,1:])



    # --- Full system matrices
    # MM = grd.M
    # KK = grd.K
    # DD = grd.D
    # print(MM)
    # print(KK)

    # --- Eigenvalue analysis 
    grd.setDOF(q*0) # <<<<< NOTE: unless the q_init was a good initial value, for now we use all zero 
    freq_d, zeta, Q, freq_0 = grd.eva()
    nMaxModes = min(4, len(freq_0))

    # --- Modes
    modes = grd.modes(norm=norm)
    nMaxModes = min(nModesPlot, len(modes))
    freq_err = np.abs(freq_0[:nMaxModes]-freq_th[:nMaxModes])/freq_th[:nMaxModes]


    
    if not test:
        printVec('Freq Num (Hz)', freq_0[:nMaxModes])
        printVec('Freq Th  (Hz)', freq_th[:nMaxModes])
        printVec('Freq Err  (%)', freq_err*100)

        fig,axes = plt.subplots(1,nMaxModes, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
        fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
        print(len( modes[0][0,1:]))
        for im, (m,ax) in enumerate(zip(modes,axes)):
            ax.plot(U_th[im,:], x_th, 'k--',label='Theory')
            # NOTE: first point is "reference/ground" frame (0,0,0)
            #ax.plot(m[1,1:]/m[1,-1], m[0,1:], label='y')
            ax.plot(m[2,1:]/m[2,-1], m[0,1:], label='z')
            ax.set_title('Mode {}'.format(im+1))
        ax.set_xlabel('')
        ax.set_ylabel('')
        ax.legend()

    return freq_err

if __name__ == '__main__':
    main(useNac=True)
    main(useNac=False)
    plt.show()

if __name__ == '__test__':
    freq_err = main(test=True, useNac=True)
    np.testing.assert_array_less(freq_err, 0.005) # 0.05%

    freq_err = main(test=True, useNac=False)
    np.testing.assert_array_less(freq_err, 0.005) # 0.05%


