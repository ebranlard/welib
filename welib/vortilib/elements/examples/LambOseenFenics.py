"""
Solving Poisson equation for LambOseen vortex using the finite element method (fenics)

The boundary condition of the vector potential needs to be set on a circle to properly retrieve the solution 

Comparison is done with analytical formulae of velocity, vector potential and vorticity.

"""

import matplotlib.pyplot as plt
import numpy as np

from scipy.special import expi, gammainc
# try:
#     from mshr import *
#     from fenics import *
# except:
#     pass

def main():

    bDiskMesh=True

    dt=0.2; tmax=2.81; x_lim=4.0; t0=2; nu=5.e-2; Gamma=4*np.pi*nu*t0; 
    nx=200
    ny=nx;
    y_lim=x_lim*1
    if bDiskMesh:
        vx=np.linspace(0,x_lim,nx);
        vy=vx*0
    else:
        vx=np.linspace(-x_lim,x_lim,nx);
        vy=np.linspace(-y_lim,y_lim,ny);
    t=t0

    # --- LambOseen formulae
    def fU_LambOseen(X,Y,Gamma,t,nu):
        X=np.asarray(X)
        Y=np.asarray(Y)
        r      = np.sqrt(X**2+Y**2)
        theta  = np.arctan2(Y,X)
        utheta = Gamma/(2*np.pi*r)*(1-np.exp(- r**2 / (4*nu*t)))
        utheta=np.asarray(utheta)
        utheta[r==0]=0
        ux=-utheta*np.sin(theta)
        uy=utheta*np.cos(theta)
        return ux,uy,utheta

    def fPsi_LambOseen(X,Y,Gamma,t,nu,cst=0):
        gammaEulerMascheroni = 0.57721566490153286

        X=np.asarray(X)
        Y=np.asarray(Y)

        r      = np.sqrt(X**2+Y**2)
        k=1/(4*nu*t)
        psi= -Gamma/(4*np.pi)*( np.log(k*r**2) - expi(-k*r**2))  + cst
        psi=np.asarray(psi)
        psi[r==0] = Gamma/(4*np.pi)*gammaEulerMascheroni + cst
        return psi

    def fUmax_LambOseen(Gamma,t,nu):
        _,u,_= fU_LambOseen(2.24*np.sqrt(nu*t),0,Gamma,t,nu) # note approx
        return u

    def fOmega_LambOseen(X,Y,Gamma,t,nu):
        r=np.sqrt(X**2+Y**2);
        theta=np.arctan2(Y,X);
        return(Gamma/(4*np.pi*nu*t)*(np.exp(-r**2/(4*nu*t))));

    # --- Theory
    # X,Y=np.meshgrid(vx,vy)
    # u_th,v_th,ut_th=fU_LambOseen(X,Y,Gamma,t,nu)
    # Speed_th = np.sqrt(u_th**2 + v_th**2)
    # fig=plt.figure()
    # ax = fig.add_subplot(111)
    # # cf=ax.contourf(vx,vy,Speed_th)
    # qv=ax.quiver(X,Y,u_th,v_th,color='k')
    # ax.set_aspect('equal', 'box')
    # # fig.colorbar(cf)
    # plt.show()

    # --- Poisson solve
    if bDiskMesh:
        domain = Circle(Point(0, 0), x_lim)
        mesh = generate_mesh(domain, nx)
    #         plot(mesh)
    else:
        mesh = RectangleMesh(Point(-x_lim,-y_lim),Point(x_lim,y_lim),nx-1, ny-1)
    V    =       FunctionSpace(mesh, 'P', 1)

    # Define boundary condition
    u_D = Constant(0) #('{} + x[0]*x[0] + 2*x[1]*x[1]'.format(Gamma), degree=2)

    def boundary(x, on_boundary):
        return on_boundary

    bc = DirichletBC(V, u_D, boundary)

    # Define variational problem
    u = TrialFunction(V)
    v = TestFunction(V)
    #f = Constant(-6.0)
    f = Expression('{}*exp(-( pow(x[0],2) + pow(x[1],2)  )/{})'.format(Gamma/(4*pi*nu*t),4*nu*t), degree=2)
    a = dot(grad(u), grad(v))*dx
    L = f*v*dx

    # Compute solution
    u = Function(V)
    solve(a == L, u, bc)
    #solve(a == L, u, bc)


    psi=u

    curlpsi = project(curl(psi))
    # Plot solution and mesh
    # plot(u)
    #plt.colorbar()
    # plot(mesh)

    # --- Mesh values
    PSI = psi.compute_vertex_values(mesh)
    if not bDiskMesh:
        PSI = PSI.reshape(ny,nx)
    U = curlpsi.compute_vertex_values(mesh)
    XY  = mesh.coordinates()

    if bDiskMesh:
        UU = U.reshape((2,-1))
        print(U.shape)
        u= UU[0,:]
        v= UU[1,:]
        X = XY[:,0]
        Y = XY[:,1]
    else:
        U = U.reshape(2,ny,nx)
        u= U[0,:,:]
        v= U[1,:,:]
        XY = XY.reshape(ny,nx,2)
        X = XY[:,:,0]
        Y = XY[:,:,1]
    Speed = np.sqrt(u**2+v**2)
    R     = np.sqrt(X**2+Y**2)

    # --- Value on an axis (NASTY FOR LOOP)
    if bDiskMesh:
        psi_cst = fPsi_LambOseen(x_lim,0,Gamma,t,nu)
    else:
        psi_cst = fPsi_LambOseen(x_lim*np.sqrt(2),0,Gamma,t,nu)
    _,_,ut_axis_th= fU_LambOseen(vx,vx*0,Gamma,t,nu)
    psi_axis_th   = fPsi_LambOseen(vx,vx*0, Gamma, t,nu, -psi_cst )
    psi_axis_num = np.array([psi(    (x,0)) for x in vx])
    u_axis       = np.array([curlpsi((x,0)) for x in vx])
    omega_axis   = np.array([f(      (x,0)) for x in vx])
    ut_axis_num_grad = np.gradient(-psi_axis_th , vx)
    ut_axis_th_grad  = np.gradient(-psi_axis_num, vx)

    # fig=plt.figure()
    # ax = fig.add_subplot(111)
    # ax.plot(vx,omega_axis ,label='omega')
    # ax.plot(vx,u_axis[:,0],label='u'   )
    # ax.plot(vx,u_axis[:,1],label='v'   )
    # ax.plot(vx,psi_axis,label='psi'   )
    # ax.legend()


    # --- Theory
    u_th,v_th,ut_th = fU_LambOseen(X,Y,Gamma,t,nu)
    ut_max          = fUmax_LambOseen(Gamma,t,nu)
    psi_th          = fPsi_LambOseen(X,Y,Gamma,t,nu,-psi_cst)

    # --- Plot simulation/theory
    #fig=plt.figure()
    #ax = fig.add_subplot(111)
    ## cf=ax.contourf(vx,vy,PSI)
    #if not bDiskMesh:
    #    cf=ax.contourf(vx,vy,Speed)
    #qv=ax.quiver(X,Y,u,v,color='b')
    #qv=ax.quiver(X,Y,u_th,v_th,color='k')
    #ax.set_aspect('equal', 'box')
    # fig.colorbar(cf)

    # --- Comparison simulation/theory
    theta = np.arctan2(Y,X)                                ;
    ut = np.sign(Gamma)*np.sqrt(u**2 + v**2)
    # AbsErrU = np.abs( u-u_th)
    # AbsErrV = np.abs( v-v_th)
    # AbsErr  = np.abs((v-v_th)**2 + (u-u_th)**2)
    AbsErr  = np.abs((ut-ut_th))
    RelErr  = np.abs((ut-ut_th)/ut_th)


    fig=plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(vx        , ut_axis_th       , 'k-'  , label  = 'ut theory (on axis)')
    ax.plot(vx        , ut_axis_th_grad  , '.'   , label  = 'ut th (on axis gradient psi th))')
    ax.plot(vx        , ut_axis_num_grad , '--'  , label  = 'ut num (on gradient psi FE))')
    ax.plot(R.ravel() , ut_th.ravel()    , '.'   , label  = 'ut theory (on mesh)')
    ax.plot(R.ravel() , ut.ravel()       , '.'   , label = 'ut num (on mesh, FE)'    )
    ax.plot([0        , x_lim]           , [ut_max , ut_max] , 'k--')
    ax.legend()


    fig=plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(vx,psi_axis_th,'k-',          label = 'psi theory (on axis)')
    ax.plot(vx,psi_axis_num,'--',         label = 'psi num (on axis)')
    ax.plot(R.ravel(),psi_th.ravel(),'.', label = 'psi theory (on mesh)',markersize=5)
    ax.plot(R.ravel(),PSI.ravel   (),'.', label = 'psi num (on mesh FE)',markersize=4)
    ax.legend()


    fig=plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(R.ravel(),AbsErr.ravel(),'.',label = 'AbsErrUt'  )
    ax.plot(R.ravel(),RelErr.ravel(),'.',label = 'RelErrUt'  )
    # ax.plot(R,AbsErrV,'.',label='AbsErrV')
    ax.legend()


    plt.show()



    # Save solution to file in VTK format
    # vtkfile = File('poisson/solution.pvd')
    # vtkfile << u

    # # Compute error in L2 norm
    # error_L2 = errornorm(u_D, u, 'L2')
    # 
    # # Compute maximum error at vertices
    # vertex_values_u   =   u.compute_vertex_values(mesh)
    # vertex_values_u_D = u_D.compute_vertex_values(mesh)
    # error_max = np.max(np.abs(vertex_values_u_D - vertex_values_u))
    # 
    # # Print errors
    # print('error_L2  =', error_L2)
    # print('error_max =', error_max)

    # Hold plot
    plt.show()
