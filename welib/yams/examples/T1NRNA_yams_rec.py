# -- IPython specific
# from IPython.display import display
# from IPython.core.interactiveshell import InteractiveShell
# InteractiveShell.ast_node_interactivity = "all"
# K_ShearTors=(M1_ST*M2_ST*M3_ST).expand()
import numpy as np
import sympy
from sympy import Matrix, Symbol, symbols
from sympy.printing import lambdarepr
from sympy import init_printing
from sympy import lambdify
from sympy import cos,sin, expand_trig
from sympy import trigsimp
from sympy import simplify

from welib.yams.yams_sympy import colvec, R_x, R_y, R_z, cross
from welib.yams.yams_sympy import GroundBody
from welib.yams.yams_sympy import BeamBody
from welib.yams.yams_sympy import RigidBody


def main():
    init_printing(use_unicode=False, wrap_line=False, no_global=True)
    #init_printing(wrap_line=False)
    # init_printing(use_latex='mathjax')
    display=lambda x: sympy.pprint(x, use_unicode=False,wrap_line=False)
    disp=lambda x: print(lambdarepr.lambdarepr(x))
    sep=lambda : print('--------')


    # ---
    L  = symbols('L')
    theta_yaw, theta_tilt   = symbols('theta_yaw,theta_tilt')
    alpha_x,alpha_y,alpha_z = symbols('alpha_x,alpha_y,alpha_z')
    ux1c,ux2c,ux3c,ux4c     = symbols('ux1c,ux2c,ux3c,ux4c')
    uy1c,uy2c,uy3c,uy4c     = symbols('uy1c,uy2c,uy3c,uy4c')
    uz1c,uz2c,uz3c,uz4c     = symbols('uz1c,uz2c,uz3c,uz4c')
    vx1c,vx2c,vx3c,vx4c     = symbols('vx1c,vx2c,vx3c,vx4c')
    vy1c,vy2c,vy3c,vy4c     = symbols('vy1c,vy2c,vy3c,vy4c')
    vz1c,vz2c,vz3c,vz4c     = symbols('vz1c,vz2c,vz3c,vz4c')

    rhoN_x,rhoN_y,rhoN_z = symbols('rhoN_x,rhoN_y,rhoN_z')  # Position of Nac center of gravity in N
    rNR_x,rNR_y,rNR_z = symbols('rNR_x,rNR_y,rNR_z')              # Position of rotor center in N
    # T_x,T_y,T_z = symbols('T_x,T_y,T_z')  #Thurst components in nacelle system
    T     = symbols('T')  #Thrust along the main shaft
    M_RNA = symbols('M_RNA')
    g     = symbols('g')

    subs={}
    subs.update({alpha_x:0})
    # subs.update({alpha_y:0})
    subs.update({alpha_z:0})
    # alpha_z:0})
    subs.update({ux1c:1,  ux2c:1, ux3c:1,  ux4c:1})
    subs.update({uz1c:1,  uz2c:1, uz3c:1,  uz4c:1})
    subs.update({uy1c:1,  uy2c:1, uy3c:1,  uy4c:1})
    subs.update({theta_yaw:0  })
    # subs.update({theta_tilt:0 })


    # --- Main parameters
    bTiltBeforeNac=False
    main_axis  = 'z'
    nD         = 1
    theta_yaw  = 0

    nShapes_twr=1
    nShapes_bld=0
    nDOF = nShapes_twr + nShapes_bld * 3
    q = np.zeros((nDOF,1))
    q[[0]]=0
    # q[[1]]=0.0
    # q[[2]]=0*np.pi/4.



    r_ET_inE  = colvec([0,0,0])
    r_TN_inT  = colvec([0,0,L])
    g_inE     = colvec([0,0,-g])
    r_NR_inN  = colvec([rNR_x,0,rNR_z])
    rho_N_inN = colvec([rhoN_x,0,rhoN_z])


    # --- Independent bodies
    Grd = GroundBody()
    Twr = BeamBody ('Twr',nShapes_twr, main_axis=main_axis, nD=nD)
    Nac = RigidBody('Nac'      ,0,0,0)

    # --- Connect bodies together
    if bTiltBeforeNac:
        # R_cn0 = R_z (theta_yaw) * R_y(theta_tilt)
        R_cn0 = R_z (theta_yaw) * R_y(theta_tilt)
        T_inN = colvec([T,0,0])
    else:
        R_cn0 = R_z (theta_yaw)
        T_inN = colvec([T*cos(theta_tilt),0,-T*sin(theta_tilt)])

    Grd.connectTo(Twr, Point=r_ET_inE, Type='Rigid')
    Twr.connectTo(Nac, Point=r_TN_inT, Type='Rigid', RelOrientation = R_cn0 )
    nq=Grd.setupDOFIndex(0);

    print('Number of DOFs: ')
    if nq!=len(q):
       print('>>> ',nq,len(q))
       raise Exception('Wrong number of dof')



    print('------------------ p=GROUND   i=TOWER --------------------------------------')
    Grd.updateChildrenKinematicsNonRecursive(q)
    print('------------------ p=TOWER   i=NACELLE --------------------------------------')
    Twr.updateChildrenKinematicsNonRecursive(q)


    print('------------------ TOWER --------------------------------------')
    print('B_T')
    display(Twr.B)
    print(np.array(Twr.B.subs(subs)))
    print('B_T_in_T')
    display(Twr.B_inB)
    print(np.array(Twr.B_inB.subs(subs)))
    print('BB_T_in_T')
    display(Twr.BB_inB)
    print(np.array(Twr.BB_inB.subs(subs)))

    print('------------------ NACELLE --------------------------------------')
    print('B_N')
    display(Nac.B)
    print(np.array(Nac.B.subs(subs)))
    print('B_N_in_N')
    display(Nac.B_inB)
    print(np.array(Nac.B_inB.subs(subs)))
    print('BB_N_in_N')
    display(Nac.BB_inB)
    print(np.array(Nac.BB_inB.subs(subs)))


    print('------------------ TOWER TOP FORCES IN EARTH--------------------------------------')
    print('Thrust in E')
    T_inE = Nac.R_0b*T_inN
    display(T_inE)
    print(np.array(T_inE.subs(subs)))

    print('Moment from thrust in E')
    r_NR_inE = Nac.R_0b*r_NR_inN
    MT_inE = Matrix(cross(r_NR_inE, T_inE))
    display(MT_inE)
    print(np.array(MT_inE.subs(subs)))


    W_inE= M_RNA*g_inE

    rho_N_inE = Nac.R_0b*rho_N_inN
    MW_inE = M_RNA*Matrix(cross(rho_N_inE , g_inE))

    print('Moment from weight in E')
    display(MW_inE)
    print(np.array(MW_inE.subs(subs)))


    print('FullForce in E')
    F_inE = W_inE + T_inE
    display(F_inE)
    print(np.array(F_inE.subs(subs)))
    print('Fullmoment in E at N')
    M_inE = MW_inE + MT_inE
    print(np.array(M_inE.subs(subs)))


    print('FullLoad in E at N')
    f_inE = Matrix(np.vstack((F_inE,M_inE)))
    print(np.array(f_inE.subs(subs)))


    print('')
    print('Fx in E<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< ')
    print(np.array(simplify(F_inE[0].subs(subs))))
    print('')
    print('My in E<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< ')
    print(np.array(simplify(M_inE[1].subs(subs))))
    print('')
    print('Fz in E<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< ')
    print(np.array(simplify(F_inE[2].subs(subs))))
    print('')


    print('Generalized force in E at N')
    GF_N_fromE = Nac.B.T * f_inE
    print(np.array(GF_N_fromE.subs(subs)))



    print('------------------ TOWER TOP FORCES IN NAC--------------------------------------')
    print('Thrust in N')
    display(T_inN)
    print(np.array(T_inN.subs(subs)))

    print('Moment from thrust in N')
    MT_inN = Matrix(cross(r_NR_inN, T_inN))
    display(MT_inN)
    print(np.array(MT_inN.subs(subs)))

    g_inN = Nac.R_0b.T * g_inE
    W_inN= M_RNA* g_inN

    MW_inN = M_RNA*Matrix(cross(rho_N_inN , g_inN))
    # 
    print('Moment from weight in N')
    display(MW_inN)
    print(np.array(MW_inN.subs(subs)))


    print('FullForce in N')
    F_inN = W_inN + T_inN
    display(F_inN)
    print(np.array(F_inN.subs(subs)))




    print('Fullmoment in N at N')
    M_inN = MW_inN + MT_inN
    print(np.array(M_inN.subs(subs)))


    print('FullLoad in E at N')
    f_inN = Matrix(np.vstack((F_inN,M_inN)))
    print(np.array(f_inE.subs(subs)))

    print('Generalized force in N at N')
    GF_N = Nac.B_inB.T * f_inN
    print(np.array(GF_N.subs(subs)))

    print('Generalized force in E at N')
    print(np.array(GF_N_fromE.subs(subs)))

    print('')
    print('Generalized force simplified')
    display(simplify(GF_N))
    print('')
    print('---------------------------')
    display(simplify(GF_N-GF_N_fromE))

    print('----------Fx in E----------')
    display(simplify(F_inE[0]).subs(subs))
    print('----------Fz in E----------')
    display(simplify(F_inE[2]).subs(subs))
    print('----------My in E----------')
    display(simplify(M_inE[1]).subs(subs))










    # display(Grd.R_bc)
    # sep()
    # display(Grd.Bhat_x_bc)
    # sep()
    # display(Grd.Bhat_t_bc)
    # sep()
    # 
    # display(Twr.R_bc)
    # sep()
    # display(Twr.Bhat_x_bc)
    # sep()
    # display(Twr.Bhat_t_bc)



    # M=Matrix([[],[],[]])
    # Matrix()
    # print(len(M))
    # print(M.shape)
    # display(M*r_T_inE)


    # r_T     = r_E+ r_TN_inE
    # 
    # # Grd.s_C_inB = r_T_inE              # NOTE: only one connection supported
    # # r_pi          = Grd.R_0p * Grd.s_C_inB
    # # 
    # # 
    # # Grd.Baug_B_inB=eye(4)
    # 
    # # --- Tower
    # 
    # # --- Tower nacelle connection
    # r_N
    # 
    # R_tc  = R_x(alpha_x) * R_y(alpha_y)  * R_z(alpha_z)  # R_pc
    # R_cn =  R_z (theta_yaw) * R_y(theta_tilt)            # R_ci0
    # R_TN = R_tc*R_cn
    # R_EN = R_ET * R_TN
    # 
    # display(R_TN)
    # display(R_EN)
    # 
    # Grd.connectTo(Twr,'Point',[0;0;0],'Type','Rigid');
    # Grd.connectTo(Twr,'Point',[0;0;0],'Type','Rigid');
    # Twr.connectTo(Nac,'Point','LastPoint','Type','Rigid');
    # Nac.connectTo(Sft,'Point',r_NS_inN,'Type','SphericalJoint','JointRotations',{'z'},'Orientation',fRotz(pi));


    # q        : [u_xf   , u_zf     , phi_y   , u_yf   , phi_z   , phi_x      , u_zt   , phi_yt   , u_yt    , phi_zt   , theta_yaw   , theta_tilt  ] $
    # JH_T : addcol(C0, C0, ey, C0,  R_y(phi_y) . ez , R_y(phi_y) . R_z(phi_z) . ex )$  /** 
    # JH_K : addcol(Ckob*ey, ey, -Ckob * R_y(alpha_y). ez ,  R_y(alpha_y) . ez , R_y(alpha_y) . R_z(alpha_z) . ex,  R_y(alpha_y) . R_z(alpha_z). R_x(theta_yaw). ey )

    # /* Jacobian Ksys - Translational, Hat and Rotational */ 
    # JT_K :  jacobian(r_K, q )$
    # JT_K_inK: (A_FK . JT_K) $
    # JR_K : addcol(JR_T, R_ET . JH_K)$
    # JR_K_inK: (A_FK . JR_K)$
    # J_K :     append(JT_K, JR_K)$
    # J_K_inK : append(JT_K_inK, JR_K_inK)$
    # A_GK: J_K_inK $


if __name__=="__main__":
    main()
    plt.show()
if __name__=="__test__":
    pass
if __name__=="__export__":
    pass
    #from welib.tools.repo import export_figs_callback
    #export_figs_callback(__file__)
