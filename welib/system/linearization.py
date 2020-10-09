import numpy as np
import copy


def numerical_jacobian(f, op, arg_number, deltas, *f_args):
    """
    Compute the jacobian of the function `f` at the operating point `op`
    with respect to its argument `arg_number` using the symmetric difference quotient method.

    example:
        f(x,u,p) where x,u are arrays

        partial f / partial x  = numerical_jacobian(f, dx, (x0,u0), 0, p) 
        partial f / partial u  = numerical_jacobian(f, du, (x0,u0), 1, p) 

    INPUTS:
        f: function of n arguments
        op: tuple of operating point values, e.g. (x0, y0) 
        arg_number: index of the argument of f (starting at 0) about which the jacobian needs to be computed
        deltas: array of numerical delta to be used to perform perturbations of the argument `arg_number` of the function. This array should have the same length as the input argument being perturbed
        *f_args: list of additional arguments required for the function f

    OUTPUTS:
       jac: jacobian, partial f/partial arg at op
    
    """
    if not isinstance(op,tuple):
        raise Exception('Operating point needs to be specified as a tuple')
    op     = list(op)
    f_args = list(f_args)

    # Number of states assumed from call at operating point
    f_op= f(*(op+f_args))
    nx=len(f_op)
    
    # Number of variables obtained from argument number 
    nj = len(op[arg_number])
    if nj!=len(deltas):
        raise Exception('Number of deltas ({}) different from dimension of operating point number {} ({}) '.format(len(deltas), arg_number, len(op[arg_number])))

    jac = np.zeros((nx,nj))


    for j in range(nj):
        # Positive and negative perturbation of operating point for arg
        op_p = copy.deepcopy(op)
        op_m = copy.deepcopy(op)
        op_p[arg_number][j] += deltas[j]
        op_m[arg_number][j] -= deltas[j]
        # Evaluation of the function at these points
        f_j_p= f(*(op_p+f_args)).flatten()
        f_j_m= f(*(op_m+f_args)).flatten()
        # Partial derivative using symmetric difference quotient
        jac[:,j] = (f_j_p-f_j_m) / (2*deltas[j])

    return jac




def linearize_function(F, xop, Iargs, delta_args,  *p):
    """ 
    Compute the jacobians of a vectorial function

    INPUTS:

    F: function with the following interface
          F(x0, x1, ..., xn, p0, p1)
       where each xi is a numpy array
    op: tuple of values for the operating point:
         (x0op, x1op, ..., xnop)

    Iargs: list of argument indices for which jacobians are to be computed
          for instance [2, 4,  n] implies jacobains for x2 x4 and xn

    delta_args: list of deltas used for finite differences, the list follows Iargs
          [ deltax2, deltax4] where each deltaxi is an array of the same size of xi

    p: optional list of arguments for the function F
         p0, p1 ...

    returns: list of jacobians
        [  partial F/partial xi  for i in Iargs ]
    
    """
    Jacs = []
    for iarg, deltas in zip(Iargs, delta_args):
        jac = numerical_jacobian(F, xop, iarg, deltas, *p)
        Jacs.append(jac)
    return Jacs

def linearize_explicit_system():
    pass



if __name__ == '__main__':
    from pybra.clean_exceptions import *

    def f_test(x, u):
        xdot=np.zeros(len(x))

        xdot[0] = x[0]**2 + x[2] + 3*u[0] 
        xdot[1] = x[1]**3            
        xdot[2] =             u[2]**3

        return xdot

    def f_testp(x, u, p):
        xdot=np.zeros(len(x))

        xdot[0] = x[0]**2 + x[2] + 3*u[0] 
        xdot[1] = x[1]**3                 + p[0]
        xdot[2] =             u[2]**3

        return xdot

    def jacx_test(x, u):
        jacx=np.zeros((x.shape[0],x.shape[0]))
        jacx[0,0] = 2*x[0];  jacx[0,2] = 1 
        jacx[1,1] = 3*x[1]**2; 
        return jacx

    def jacu_test(x, u):
        jacu=np.zeros((x.shape[0],u.shape[0]))
        jacu[0,0] = 3
        jacu[2,2] = 3*u[2]**2
        return jacu

    # Define an operating point
    p=np.zeros(2)
    x0=np.zeros((3,1))
    u0=np.zeros((3,1))
    x0[:,0]=[1,2,3]
    u0[:,0]=[10,20,30]
    dx=[0.01]*3
    du=[0.01]*3
    p[0]=100
    p[1]=2

    print(f_test (x0,u0))
    print(f_testp(x0,u0,p))
    print('---- Jac x')
    print(jacx_test(x0,u0))

    jacx_num = numerical_jacobian(f_test , dx, (x0,u0), 0)
    print(jacx_num)
    jacx_nump= numerical_jacobian(f_testp, dx, (x0,u0), 0, p)
    print(jacx_nump)

    print('---- Jac u')
    print(jacu_test(x0,u0))
    jacu_num = numerical_jacobian(f_test , du, (x0,u0), 1)
    print(jacu_num)
    jacu_nump= numerical_jacobian(f_testp, du, (x0,u0), 1, p)
    print(jacu_nump)




    
