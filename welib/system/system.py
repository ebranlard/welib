""" 
Handles a general system made of a state equation and an output equation.
The equations can be implicit or explicit non linear differential equations 

"""
import numpy as np
from numpy.linalg import inv
#from scipy.integrate import  solve_ivp #odeint

from .linearization import * 

# --------------------------------------------------------------------------------}
# --- Linear system
# --------------------------------------------------------------------------------{
class System():
    """ 
    Handles a general system made of a state equation and an optional output equation.
    The equations can be implicit or explicit non linear differential equations.

    General implicit form:
          0 = Fx(t,xdot,x,u,p)

    General explicit form:
          xdot= Fx(t,x,u,p)

    General output equation
          y = Y(t,x,u,p)  
    """

    def __init__(self, Fx, Y=None, interface='xup', param=None):
        """
        Fx: state equation function
            explicit form:
               xdot = Fx(t,x,..)
            implicit forms:
               Fx(t,xdot,x,...) = 0

        Y : output equation function
           y = Y(t,x,..)

        interface: string representing the interface of the equations
          'x'      :  Fx(t,x)            Y(t,x)            (explicit)
          'xp '    :  Fx(t,x,p)          Y(t,x,p)          (explicit)
          'xu '    :  Fx(t,x,u)          Y(t,x,u)          (explicit)
          'xup'    :  Fx(t,x,u,p)        Y(t,x,u,p)        (explicit)
          'xdotx'  :  Fx(t,xdot,x)       Y(t,x)            (implicit)
          'xdotxp' :  Fx(t,xdot,x,p)     Y(t,x,p)          (implicit)
          'xdotxu' :  Fx(t,xdot,x,u)     Y(t,x,u)          (implicit)
          'xdotxup':  Fx(t,xdot,x,u,p)   Y(t,x,u,p)        (implicit)

        """
        # TODO Need to handle different interface for Fx and Y!


        interface=interface.replace(',','').replace('q','x')

        self.Fx        = Fx
        self.Y         = Y
        self.interface = interface
        self.param     = param


        self.implicit  = interface.find('xdot')==0
        self.has_param = interface.find('p')>0 
        self.has_input = interface.find('u')>0
        # Being relaxed about parameters being provided
        #if self.has_param and param is None:
        #    raise Exception('Parameters needs to be provided since necessary for the function interface')

    @property
    def has_output(self):
        return self.Y is not None


    def linearize(self, op, dx, dxd=None, du=None, use_implicit=False):
        """ 
        Linearize the system

        if implicit:
            op = (t_op, xdot_op, x_op, u_op) (with u_op optional)
        if explicit:
            op = (t_op, x_op, u_op)          (with u_op optional)
        
        """
        if self.has_param and self.param is None:
            raise Exception('Parameters needs to be provided since necessary for the function interface')

        # Checks
        if not isinstance(op,tuple):
            raise Exception('Operating point needs to be specified as a tuple')

        if len(op)!=self.nArgs:
            raise Exception('Number of values of operating point ({}) does not match number of main argument of function ({}), with interface {} '.format(len(op),self.nArgs,self.interface))
        if self.has_input and du is None:
            raise Exception('du needs to be specicified for linearization')



        if self.implicit or use_implicit:
            # --- Linearization of implicit state equation

            if use_implicit and not self.implicit:
                # We need to add xdot0 to operating point
                xdot0 = self.Fx(*(op+(self.param,)))
                op_imp = (op[0], xdot0)+tuple(op[1:]) # (t, xdot0, x0, u0)
            else:
                op_imp=op

            if len(op_imp)!=self.nArgsImplicit:
                raise Exception('Number of values of implicit operating point ({}) does not match number of main argument of implicit function ({}), with interface {} '.format(len(op_imp),self.nArgsImplicit,self.interface))

            if dxd is None:
                raise Exception('delta_xdot needs to be specicified for linearization of implicit equations')

            F=self.implicit_function

            deltas = np.array([None, dxd, dx, du], dtype=object)     # (dt, dxd, dx, du)
            Iargs  = list(range(1,self.nArgsImplicit)) # [1,2] or [1,2,3]
            deltas = deltas[Iargs]

            # --- Compute necessary jacobians
            jacs= linearize_function(F, op_imp, Iargs, deltas, self.param)

            E  =   jacs[0]
            Ap = - jacs[1]
            Bp = - jacs[2]
            A = inv(E).dot(Ap)
            if self.has_input:
                Bp = - jacs[2]
                B = inv(E).dot(Bp)
            else:
                Bp= None
                B = None
        else:
            # --- Linearization of explicit equation
            deltas = np.array([None, dx, du], dtype=object)  # (dt, dx, du)
            Iargs  = list(range(1,self.nArgs)) # [1] or [1,2]
            deltas = deltas[Iargs]

            # --- Compute necessary jacobians
            if self.has_param:
                A, B = linearize_function(self.Fx, op, Iargs, deltas, self.param)
            else:
                A, B = linearize_function(self.Fx, op, Iargs, deltas)

        # --- Linearization of output equation
        if self.has_output:
            deltas = np.array([None, dx, du], dtype=object)  # (dt, dx, du)
            Iargs  = list(range(1,self.nArgs)) # [1] or [1,2]
            deltas = deltas[Iargs]
            if self.has_param:
                C, D = linearize_function(self.Y, op, Iargs, deltas, self.param)
            else:
                C, D = linearize_function(self.Y, op, Iargs, deltas)

            return A, B, C, D
        else:
            return A, B

    @property
    def implicit_function(self):
        """ return implicit state equation """
        if self.implicit:
            return self.Fx
        else:
            if self.interface =='x':
                return lambda t, xdot, x, u    : self.Fx(t,x)-xdot
            elif self.interface =='xu':
                return lambda t, xdot, x, u    : self.Fx(t,x,u)-xdot
            elif self.interface =='xp':
                return lambda t, xdot, x, p    : self.Fx(t,x,p)-xdot
            elif self.interface =='xup':
                return lambda t, xdot, x, u, p : self.Fx(t,x,u,p)-xdot
            else:
                raise Exception('Unsupported interface')


    @property
    def nArgs(self):
        """ Number of main arguments of Fx (t,x,u,xdot), not counting parameters """
        n=2 # x and t
        if self.has_input:
            n+=1 #u
        if self.implicit:
            n+=1 # xdot
        return n 

    @property
    def nArgsImplicit(self):
        """ Number of main arguments of the implicit form of Fx (t,x,u,xdot), not counting parameters """
        if self.implicit:
            return self.nArgs
        else:
            return self.nArgs+1

    @property
    def nArgsOutput(self):
        """ Number of main arguments of the output equation Fx (t,x,u), not counting parameters """
        return self.nArgsImplicit-1

    @property
    def nOuputs(self):
        pass

    def __repr__(self):
        s ='<{} object>:\n'.format(type(self).__name__)
        s+=' - interface: {}\n'.format(self.interface)
        s+=' - implicit:  {}\n'.format(self.implicit)
        s+=' - has_param: {}\n'.format(self.has_param)
        s+=' - has_input: {}\n'.format(self.has_input)
        return s

if __name__ == '__main__':

    def model(t, x, u, p):
        return p['A'].dot(x) + p['B'].dot(u)

    def output(t, x, u, p):
        return p['C'].dot(x) + p['D'].dot(u)

    nu = 3
    nx = 2
    p=dict()
    p['A']      = np.eye(nx)
    p['A'][0,1] = 3
    p['B']      = np.zeros((nx,nu))
    p['B'][0,0] = 1
    p['B'][0,1] = 2
    p['B'][0,2] = 3
    p['B'][1,2] = 12
    p['C']      = np.zeros((nx,nx))
    p['C'][0,0] = 1
    p['C'][1,0] = 2
    p['C'][1,1] = 5
    p['D']      = np.zeros((nx,nu))
    p['D'][1,0] = 1
    p['D'][1,1] = 2
    p['D'][1,2] = 3
    p['D'][0,1] = 12
    
    x0 = np.zeros((nx,1))
    u0 = np.zeros((nu,1))
    x0[0]=0.1
    x0[1]=0.2
    u0[0]=-1
    u0[1]=3
    xdot0=model(0,x0,u0,p) 
    print(xdot0)

    sys = System(Fx=model, Y=output, interface='xup', param = p )
    print(sys)
    F=sys.implicit_function
    print(F(0,xdot0,x0,u0,p))


    delta_x  = np.zeros(x0.shape)+0.001
    delta_xd = np.zeros(x0.shape)+0.001
    delta_u  = np.zeros(u0.shape)+0.001

    op=(0,x0,u0)
    A,B,C,D=sys.linearize(op, dx=delta_x, dxd=delta_xd, du=delta_u, use_implicit=True)
    print('A')
    print(A)
    print('B')
    print(B)
    print('C')
    print(C)
    print('D')
    print(D)

    A,B,C,D=sys.linearize(op, dx=delta_x, du=delta_u)
    print('A')
    print(A)
    print('B')
    print(B)
    print('C')
    print(C)
    print('D')
    print(D)




