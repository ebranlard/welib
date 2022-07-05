import numpy as np
import copy


# --------------------------------------------------------------------------------}
# --- Generic functions 
# --------------------------------------------------------------------------------{
def numerical_jacobian(f, op, arg_number, deltas, *f_args):
    """
    Compute the jacobian of the function `f` at the operating point `op`
    with respect to its argument `arg_number` using the symmetric difference quotient method.

    example:
        f(x,u,p) where x,u are arrays, p parameters array/dict, optional

        partial f / partial x  = numerical_jacobian(f, (x0,u0), dx, 0, p) 
        partial f / partial u  = numerical_jacobian(f, (x0,u0), du, 1, p) 

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

    # Convert op[arg_number] to array of floats
    deltas = np.asarray(deltas)
    op = copy.deepcopy(op)
    op[arg_number] = np.asarray(op[arg_number]).astype(float)
    #dtype_op   = op[arg_number].dtype
    #dtype_delta =deltas.dtype
    #if dtype_op!=dtype_op:
    #    raise Exception('Type of op ({}) and deltas ({}) are not the same'.format(dtype_op, dtype_delta))

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
          F(x, u, p)                 (typically)
       where each xi is a numpy array. 

    op: tuple of values for the operating point:
         (x0op, x1op, ..., xnop)
         (xop, uop)                (typically)

    Iargs: list of argument indices for which jacobians are to be computed
          for instance [2, 4,  n] implies jacobains for x2 x4 and xn
          typically: [0,1] gives "A" and "B" matrix

    delta_args: list of deltas used for finite differences, the list follows Iargs
          [ deltax2, deltax4] where each deltaxi is an array of the same size of xi

    p: optional list of arguments for the function F
         p0, p1 ...

    returns: list of jacobians
        [  partial F/partial xi  for i in Iargs ]

    example:
        def F(x,u,p):
            return A(p).dot(x) + B(p).dot(u)

        dx=[0.01]*3
        du=[0.01]*3
        A,B  = linearize_function(F, (xop, uop), [0,1], (dx, du), p)
    
    """
    Jacs = []
    for iarg, deltas in zip(Iargs, delta_args):
        jac = numerical_jacobian(F, xop, iarg, deltas, *p)
        Jacs.append(jac)
    return Jacs


def linearize_explicit_system():
    pass


# --------------------------------------------------------------------------------}
# --- Dedicated functions for usual cases 
# --------------------------------------------------------------------------------{
def linearize_Fx(F, x0, dx, *p):
    """ 
    Computes jacobian wrt to states and inputs for function F(x, p)

    INPUTS:
        F: function with the following interface: F(x, p)
            where x is an array, p is an optional array or dict
        x0: operating point values for x (array)
        dx: array of perturbations for finite differences (same size as x)
        p : optional arry or dict
    OUTPUTS:
       df_dx: Jacobian df/dx
    """
    return numerical_jacobian(F, (x0,), 0, dx, *p)
    #return linearize_function(F, (x0,), [0], (dx,), *p)[0]

def linearize_Fxu(F, x0, u0, dx, du, *p):
    """ 
    Computes jacobian wrt to states and inputs for function F(x, u, p)

    INPUTS:
        F: function with the following interface: F(x, u, p)
            where x is an array, u is an array, p is an optional array or dict
        x0: operating point values for x (array)
        u0: operating point values for u (array)
        dx: array of perturbations for finite differences (same size as x)
        du: array of perturbations for finite differences (same size as y)
        p : optional arry or dict
    OUTPUTS:
       df_dx, df_du: Jacobian df/dx df/du
    """
    return linearize_function(F, (x0, u0), [0,1], (dx, du), *p)





if __name__ == '__main__':
    # see tests/test_linearization.py
    pass



    
