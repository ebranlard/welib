""" 
Tools for yams_sympy
"""

import numpy as np
import os
import sympy
from sympy import latex, python, symarray, diff
from sympy import Symbol, Matrix, collect
from sympy import Function, DeferredVector, Derivative
import sympy.physics.mechanics as me
from sympy.physics.vector import dynamicsymbols
from sympy.physics.mechanics.functions import find_dynamicsymbols
from sympy import sin, cos, exp, sqrt


# --------------------------------------------------------------------------------}
# --- Basic math 
# --------------------------------------------------------------------------------{
def colvec(v): 
    return Matrix([[v[0]],[v[1]],[v[2]]])
def cross(V1,V2):
    return [V1[1]*V2[2]-V1[2]*V2[1], V1[2]*V2[0]-V1[0]*V2[2], (V1[0]*V2[1]-V1[1]*V2[0]) ]
def eye(n): 
    return Matrix( np.eye(n).astype(int) )


def skew(x):
    """ Returns the skew symmetric matrix M, such that: cross(x,v) = M v """
    #S = Matrix(np.zeros((3,3)).astype(int))
    if hasattr(x,'shape') and len(x.shape)==2:
        if x.shape[0]==3:
            return Matrix(np.array([[0, -x[2,0], x[1,0]],[x[2,0],0,-x[0,0]],[-x[1,0],x[0,0],0]]))
        else:
            raise Exception('fSkew expect a vector of size 3 or matrix of size 3x1, got {}'.format(x.shape))
    else:
        return Matrix(np.array([[0, -x[2], x[1]],[x[2],0,-x[0]],[-x[1],x[0],0]]))

# --------------------------------------------------------------------------------}
# --- PYDY VIZ related...
# --------------------------------------------------------------------------------{
try:
    import pydy.viz.scene 
    from pydy.viz.scene import Scene
    from pydy.viz.shapes import Shape
    MyScene=Scene
except:
    Scene=object
class MyScene(Scene):
    def __init__(self, reference_frame, origin, *visualization_frames,
                 **kwargs):
        Scene.__init__(self, reference_frame, origin, *visualization_frames,
                 **kwargs)
        self.pydy_directory = "yams-viz-resources"
    def create_static_html(self, overwrite=False, silent=False, prefix=None):
        import distutils
        import os
        cur_dir=os.getcwd()
        base = os.path.basename(cur_dir)
        while base== self.pydy_directory:
            parent = os.path.dirname(cur_dir)
            os.chdir(parent)
            print('>>>> CHANGING DIRECTORY')
            cur_dir=os.getcwd()
            base = os.path.basename(cur_dir)
        Scene.create_static_html(self, overwrite=overwrite, silent=silent, prefix=prefix)


def frame_viz(Origin, frame, l=1, r=0.08, name=''):
    """
    Generate visualization objects for a frame, i.e. three axes, represented by cylinders 
      x: red
      y: green
      z: blue

    """
    from pydy.viz.shapes import Cylinder
    from pydy.viz.visualization_frame import VisualizationFrame
    from sympy.physics.mechanics import Point

    #
    X_frame  = frame.orientnew('ffx', 'Axis', (-np.pi/2, frame.z) ) # Make y be x
    Z_frame  = frame.orientnew('ffz', 'Axis', (+np.pi/2, frame.x) ) # Make y be z
    
    X_shape   = Cylinder(radius=r, length=l, color='red') # Cylinder are along y
    Y_shape   = Cylinder(radius=r, length=l, color='green')
    Z_shape   = Cylinder(radius=r, length=l, color='blue')
    
    X_center=Point('X'); X_center.set_pos(Origin, l/2 * X_frame.y)
    Y_center=Point('Y'); Y_center.set_pos(Origin, l/2 *   frame.y)
    Z_center=Point('Z'); Z_center.set_pos(Origin, l/2 * Z_frame.y)
    X_viz_frame = VisualizationFrame(name+'_xaxis',X_frame, X_center, X_shape)
    Y_viz_frame = VisualizationFrame(name+'_yaxis',  frame, Y_center, Y_shape)
    Z_viz_frame = VisualizationFrame(name+'_zaxis',Z_frame, Z_center, Z_shape)
    return X_viz_frame, Y_viz_frame, Z_viz_frame

def body_viz(b, vo):
    """ 
    Return a "vizualization" of a body, based on options 
    For "pydy" viz for now
    
    - b: YAMSBody
    - vo: dictionary of options
    """
    from pydy.viz.shapes import Cylinder, Sphere, Cube, Box
    from pydy.viz.visualization_frame import VisualizationFrame
    from sympy.physics.mechanics import Point
    # NOTE: rememeber that cylinders are along "y" by default and about their center
    X_frame  = b.frame.orientnew('ffx', 'Axis', (-np.pi/2, b.frame.z) ) # Make y be x
    Z_frame  = b.frame.orientnew('ffz', 'Axis', (+np.pi/2, b.frame.x) ) # Make y be z
    vo = b.viz_opts
    if 'type' in vo:
        if vo['type']=='cylinder':
            l = vo['length']
            cyl = Cylinder(radius=vo['radius'], length=l, color=vo['color'])
            C_center=Point('C'); 
            if vo['normal'] == 'y':
                C_center.set_pos(b.origin, l/2 * b.frame.y)
                bod_viz = VisualizationFrame(b.name+'_body',b.frame, C_center, cyl)
            elif vo['normal'] == 'z':
                C_center.set_pos(b.origin, l/2 * Z_frame.y)
                bod_viz = VisualizationFrame(b.name+'_body',Z_frame, C_center, cyl)
            elif vo['normal'] == 'x':
                C_center.set_pos(b.origin, l/2 * X_frame.y)
                bod_viz = VisualizationFrame(b.name+'_body',X_frame, C_center, cyl)
            return [bod_viz]
        elif vo['type']=='three-blades':
            l = vo['radius']
            blade_shape = Cylinder(radius=0.03*vo['radius'] , length=l, color=vo['color'])
            # TODO might need to be thought again
            if vo['normal'] == 'z':
                B1_frame  = b.frame.orientnew('b1', 'Axis', (np.pi/2+0 ,        b.frame.z) ) # Y pointing along blade-Z
                B2_frame  = b.frame.orientnew('b2', 'Axis', (np.pi/2+2*np.pi/3, b.frame.z) )
                B3_frame  = b.frame.orientnew('b3', 'Axis', (np.pi/2+4*np.pi/3, b.frame.z) )
                B1_center=Point('B1'); B1_center.set_pos(b.origin, l/2 * B1_frame.y)
                B2_center=Point('B2'); B2_center.set_pos(b.origin, l/2 * B2_frame.y)
                B3_center=Point('B3'); B3_center.set_pos(b.origin, l/2 * B3_frame.y)
            B1_viz = VisualizationFrame(b.name+'_body_blade1',B1_frame, B1_center, blade_shape)
            B2_viz = VisualizationFrame(b.name+'_body_blade2',B2_frame, B2_center, blade_shape)
            B3_viz = VisualizationFrame(b.name+'_body_blade3',B3_frame, B3_center, blade_shape)
            return [B1_viz, B2_viz, B3_viz]
        elif vo['type']=='cube':
            cube = Cube(length=vo['length'], color=vo['color'], name=b.name+'_body')
            return [VisualizationFrame(b.name+'_body', b.frame, b.origin, cube)]
        elif vo['type']=='box':
            box = Box(width=vo['width'], height=vo['height'], depth=vo['depth'], color=vo['color'], name=b.name+'_body')
            return [VisualizationFrame(b.name+'_body', b.frame, b.origin, box)]
        else:
            raise NotImplementedError('Type, ',vo['type'])


# --------------------------------------------------------------------------------}
# --- Tools
# --------------------------------------------------------------------------------{
def exprHasFunction(expr):
    """ return True if a sympy expression contains a function"""
    if hasattr(expr, 'atoms'): 
        return len(expr.atoms(Function))>0
    else:
        return False


def subs_no_diff(expr, subslist):
    """ 
    Perform substitution in an expression, but not in the time derivatives
    Only works if Subslist is a simple lists of ( var, value) 


    see also: sympy.physics.mechanics.functions.msubs

    TODO extend to matrix

    """


    # Set mapping between time derivatives and dummy variables
    time=dynamicsymbols._t
    Dummys=symarray('DUM', len(subslist))
    DT2Dummy=[]
    for i, (var,b) in enumerate(subslist):
        if len(var.atoms())!=1:
            print(var)
            print(type(var))
            raise Exception('subs_no_diff only works for simple (atomic) substitution')
        if exprHasFunction(var):
            dvar = diff(var,time)
            DT2Dummy.append((dvar,Dummys[i]))
    Dummy2DT=[(b,a) for a,b in DT2Dummy]

    # Remove time derivative
    expr_clean = expr.subs(DT2Dummy)

    # Perform user substitution
    expr_new = expr_clean.subs(subslist)

    # Reinsert time derivatives
    return expr_new.subs(Dummy2DT)

def mycollect(expr, var_list, evaluate=True, **kwargs):
    """ Acts as collect but substitute the symbols with dummy symbols first so that it can work with partial derivatives. 
        Matrix expressions are also supported. 
    """
    if not hasattr(var_list, '__len__'):
        var_list=[var_list]
    # Mapping Var -> Dummy, and Dummy-> Var
    Dummies=symarray('DUM', len(var_list))
    Var2Dummy=[(var, Dummies[i]) for i,var in enumerate(var_list)]
    Dummy2Var=[(b,a) for a,b in Var2Dummy]
    # Replace var with dummies and apply collect
    expr = expr.expand().doit()
    expr = expr.subs(Var2Dummy)
    if hasattr(expr, '__len__'):
        expr = expr.applyfunc(lambda ij: collect(ij, Dummies, **kwargs))
    else:
        expr = collect(expr, Dummies, evaluate=evaluate, **kwargs)
    # Substitute back
    if evaluate:
        return expr.subs(Dummy2Var)
    d={}
    for k,v in expr.items():
        k=k.subs(Dummy2Var)
        v=v.subs(Dummy2Var)
        d[k]=v
    return d

def myjacobian(expr, var_list, value_list=None):
    """ Compute jacobian of expression, matrix or not. 
    Perform symbol substitution first to have support for derivatives

    J = [ \partial fi / \partial_xj ]_x0

    Inputs:
      - expr:  (f_i) expression to compute the jacobian of. Scalar, list or matrix,
      - var_list: (x_j) list of variables 
    Optional:
      - value_list: (x0): values at operating point (same length as var_list)
    
    Examples:
        x,y = symbols('x, y')
        a,b = symbols('a, b')
        f1  = a*x + b*y**2
        f2  = a*y + b*x**2
        f = Matrix([[f1],[f2]])
        J = myjacobian(f, [x,y])

    """
    if not hasattr(var_list, '__len__'):
        var_list=[var_list]
    # Mapping Var -> Dummy, and Dummy-> Var
    Dummies=symarray('DUM', len(var_list))
    Var2Dummy=[(var, Dummies[i]) for i,var in enumerate(var_list)]
    Dummy2Var=[(b,a) for a,b in Var2Dummy]
    if isinstance(expr, list):
        expr = Matrix(expr)
    expr = expr.expand().doit()
    expr = expr.subs(Var2Dummy)
    if hasattr(expr, '__len__'):
        jac = expr.jacobian(Dummies)
    else:
        jac = Matrix([expr]).jacobian(Dummies)
    jac = jac.subs(Dummy2Var)
    if value_list is not None:
        sub_list = [(var,val) for var,val in zip(var_list, value_list)]
        jac = jac.subs(sub_list)
    return jac

def linearize(expr, x0, order=1, sym=False, doSimplifyIfDeriv=True):
    """ 
    Return a Taylor expansion of the expression at the operating point x0
    INPUTS:
        expr: expression to be linearized
        x0: operating point [(x,0),(y,0),...]
        order: order of Taylor expansion
        sym: If sym is true and expr is an array of dimension 2, symmetry is assumed to avoid computing the linearization of symmetric elements
    OUTPUTS:
       linearized expression

    Basically computes:
      f(x,y) = f(x0,y0) + df/dx_0 (x-x0) + df/dy_0 (y-y0)

    if y=dx/dt, y is considered to be an independent variable (dxdt is substituted for a dummy variable y)


    Examples: 
        linearize(Function('f')(x,y), [(x, x0),(y, y0)], order=1)




    """
    time=dynamicsymbols._t
    x0_bkp = x0
    # --- First protect variables that are derivatives
    # TODO We really need to do something cleaner for that. Consider using:
    #    qd.is_Derivative
    #    qd.variables
    #    qd.variable_count
    #    qd.free_symbols
    #    qd._vwr_variables
    #    qd.expr
    DummysXD=symarray('DUMQD', len(x0))
    xd_2_dum=[]
    x0_new=[]
    for i,(q,q0) in enumerate(x0):
        if q.is_Derivative:
            xd_dum = DummysXD[i]
            x0_new.append((xd_dum, q0))
            xd_2_dum=[(q, xd_dum)]
        else:
            x0_new.append((q,q0))
    x0 = x0_new
    dum_2_xd = [(b,a) for a,b in xd_2_dum]

    # --- Then protect expression from derivatives
    # Mapping derivative -> Dummy
    Dummys=symarray('DUM', len(x0))
    DT2Dummy=[]
    for i,(v,x) in enumerate(x0):
        dv = diff(v,time)
        if dv.is_Derivative:
            DT2Dummy.append( (dv, Dummys[i])  )
    Dummy2DT=[(b,a) for a,b in DT2Dummy]

    # --- Remove Dummys
    RemoveDummy = dum_2_xd + Dummy2DT
    #print('RemoveDummy',RemoveDummy)
    #print('x0',x0)

    # Helper function, linearizes one expression (for matrices and lists..)
    def __linearize_one_expr(myexpr):
        # backup myexpr
        myexpr0=myexpr
        # Make operating point derivatives "dummy"
        myexpr = myexpr.subs(xd_2_dum)
        # Make other derivatives of states "dummy"
        myexpr = myexpr.subs(DT2Dummy)

        if myexpr.has(Derivative) and doSimplifyIfDeriv:
            print('[WARN] Yams sympy linearize: Expression still contains derivative. Behavior might not be the right one. To be safe, we are running simplify on the expression.')
            #NOTE: if expressions contains diff(-x,t) it won't work (for instance after a substitution). We need it simplified to -diff(x,t)
            # We simplify expression
            myexpr0=myexpr0.simplify()
            # Make operating point derivatives "dummy"
            myexpr = myexpr0.subs(xd_2_dum)
            # Make other derivatives of states "dummy"
            myexpr = myexpr.subs(DT2Dummy)
        # --- Order 0
        flin=myexpr.subs(x0)  # Order 0
        # --- Order 1
        if order>=1:
            df = [(myexpr.diff(v1).subs(x0))*(v1-v10) for v1,v10 in x0]
            df = sum(df)
            flin += df
        # --- Order 2
        if order>=2:
            df2 = [ (myexpr.diff(v1)).diff(v2).subs(x0)*(v1-v10)*(v2-v20) for v1,v10 in x0 for v2,v20 in x0]
            df2 = sum(df2)/2
            flin+=df2
        # --- Order 3
        if order>=3:
            df3 = [ (myexpr.diff(v1)).diff(v2).diff(v3).subs(x0)*(v1-v10)*(v2-v20)*(v3-v30) for v1,v10 in x0 for v2,v20 in x0 for v3,v30 in x0]
            df3 = sum(df3)/6
            flin+=df3

        if order>3:
            raise NotImplementedError('Higher order')
        return flin.subs(RemoveDummy)
    
    # --- Hanlde inputs of multiple dimension, scalar, vectors, matrices
    try:
        ndim=len(expr.shape)
    except:
        ndim=0

    expr_lin=expr*0

    if ndim==0:
        expr_lin = __linearize_one_expr(expr)

    elif ndim==1:
        for i in range(expr.shape[0]):
            expr_lin[i] = __linearize_one_expr(expr[i])

    elif ndim==2:
        if expr.shape[0] != expr.shape[1]:
            sym=False

        for i in range(expr.shape[0]):
            for j in range(expr.shape[1]):
                if (j<=i and sym) or (not sym):
                    expr_lin[i,j] = __linearize_one_expr(expr[i,j])
                if sym:
                    expr_lin[j,i] = expr_lin[i,j]
    else:
        raise NotImplementedError()

    return expr_lin



def smallAngleApprox(expr, angle_list, order=1, sym=True):
    """ 
    Perform small angle approximation of an expression by linearizing about the "0" operating point of each angle

    INPUTS:
        expr: expression to be linearized
        angle_list: list of angle to be considered small:  [phi,alpha,...]
        order: order of Taylor expansion
        sym: If sym is true and expr is an array of dimension 2, symmetry is assumed to avoid computing the linearization of symmetric elements
    OUTPUTS:
       linearized expression
    """
    # operating point
    x0=[(angle,0) for angle in angle_list]
    return linearize(expr, x0, order=order, sym=sym, doSimplifyIfDeriv=False)



# --------------------------------------------------------------------------------}
# --- Misc sympy utils 
# --------------------------------------------------------------------------------{
def cleantex(expr):
    """ clean a latex expression 
    
    Example:
        print( cleantex( vlatex(MM, symbol_names={Jxx:'J_{xx}'})))
    """
    s=latex(expr)
    D_rep={
        '\\operatorname{sin}':'\\sin',
        '\\operatorname{cos}':'\\cos',
        '\\left(\\theta\\right)':'\\theta',
        '\\left[\\begin{matrix}':'\\begin{bmatrix}\n', 
        '\\end{matrix}\\right]':'\n\\end{bmatrix}\n' ,
        '\\operatorname{q_{T1}}':'q_{T1}',
        '\\operatorname{q_{T2}}':'q_{T2}',
        '\operatorname{T_{a}}'  : 'T_a',
        '\\left(t \\right)'    : '(t)',
        '{(t)}': '(t)',
        '{d t}': '{dt}',
        'L_{T}':'L_T',
        'M_{N}':'M_N',
        '_{x}':'_x',
        '_{y}':'_y',
        '_{z}':'_z',
        'x(t)':'x',
        'y(t)':'y',
        'z(t)':'z',
        'T_a(t)': 'T_a',
        'g(t)'  : 'g',
        'phi_x(t)':'phi_x',
        'phi_y(t)':'phi_y',
        'phi_z(t)':'phi_z',
        'psi(t)'  :'psi',
        'q_{T1}(t)': 'q_{T1}',
        'q_{T2}(t)': 'q_{T2}',
        'frac{d}{dt} q_{T1}':'dot{q}_{T1}',
        'frac{d}{dt} \\phi_x':'dot{\\phi}_x',
        'frac{d}{dt} \\phi_y':'dot{\\phi}_y',
        'frac{d}{dt} \\phi_z':'dot{\\phi}_z',
        'frac{d}{dt} \\psi':'dot{\\psi}',
        'frac{d}{dt} q_{T2}':'dot{q}_{T2}',
#         '{RNA}':'{N}',
#         'RNA':'N',
        'x_{NG}^{2} + z_{NG}^{2}':'r_{NG}^2',
        'M_N x_{NG}^{2} + M_N z_{NG}^{2}':'M_N r_{NG}^{2}',
        'M_{R} x_{NR}^{2} + M_{R} z_{NR}^{2}':'M_{R} r_{NR}^{2}',
        '\\theta_{tilt}':'\\theta_t',
        '\\cos{\\left(\\theta_t \\right)}':'\\cos\\theta_t',
        '\\sin{\\left(\\theta_t \\right)}':'\\sin\\theta_t',
        '\\sin{\\left(\\phi_x \\right)}':'\\sin\\phi_x',
        '\\sin{\\left(\\phi_y \\right)}':'\\sin\\phi_y',
        '\\sin{\\left(\\phi_z \\right)}':'\\sin\\phi_z',
        '\\cos{\\left(\\phi_x \\right)}':'\\cos\\phi_x',
        '\\cos{\\left(\\phi_y \\right)}':'\\cos\\phi_y',
        '\\cos{\\left(\\phi_z \\right)}':'\\cos\\phi_z',
        '\\left(\\dot{\\phi}_x\\right)^{2}':'\dot{\\phi}_x^2',
        '\\left(\\dot{\\phi}_y\\right)^{2}':'\dot{\\phi}_y^2',
        '\\left(\\dot{\\phi}_z\\right)^{2}':'\dot{\\phi}_z^2',
        '\\\\':'\\\\ \n' 
        }
    for k in D_rep.keys():
        s=s.replace(k,D_rep[k])
    return s


def saveTex(expr, filename):
    print('Exporting to {}'.format(filename))
    with open(filename,'w') as f:
        f.write(cleantex(expr))
    print('Done')

def cleanPySmallMat(expr, varname='R', indent=0, replDict=None, noTimeDep=True): 
    """ Export a "small" matrix to python (row by row) """
    def cleanPyAtom(atom):
        s=repr(atom).replace(' ','')
        if replDict is not None:
            for k,v in replDict.items():
                if s.find(k)>=0:
                    s=s.replace(k,v)
        if noTimeDep:
            s=s.replace('(t)','',)
        return s
    s=''
    indent =''.join([' ']*indent)
    dims=expr.shape
    s+='{}{} = np.zeros(({},{}))\n'.format(indent,varname, dims[0],dims[1])
    for i in np.arange(dims[0]):
        s+='{}{}[{},:] = ['.format(indent,varname,i) + ','.join([cleanPyAtom(expr[i,j]) for j in np.arange(dims[1])]) +']\n' 
    return s

def cleanPySimple(expr, varname='R', indent=0, replDict=None, noTimeDep=False):
    """ 
    Clean a python sympy expression 
    """
    def cleanPyAtom(atom):
        s=repr(atom)
        # Replace
        if replDict is not None:
            for k,v in replDict.items():
                if s.find(k)>=0:
                    s=s.replace(k,v)
        s=s.replace(' ','')
        if noTimeDep:
            s=s.replace('(t)','',)
        return s
    try:
        dims=expr.shape
    except:
        dims=0
        return '{}{} = '.format(indent,varname) + cleanPyAtom(expr), list(set(parameters)), list(set(inputs)), sdofs

    indent =''.join([' ']*indent)
    s=''
    if len(dims)==1:
        s+='{}{} = np.zeros({})\n'.format(indent,varname,dims[0])
        for i in np.arange(dims[0]):
            s+='{}{}[{}] = '.format(indent,varname,i) + cleanPyAtom(expr[i]) +'\n'
    elif len(dims)==2:
        s+='{}{} = np.zeros(({},{}))\n'.format(indent,varname, dims[0],dims[1])
        for i in np.arange(dims[0]):
            for j in np.arange(dims[1]):
                s+='{}{}[{},{}] = '.format(indent,varname,i,j) + cleanPyAtom(expr[i,j])+'\n'
    return s



# TODO TODO TODO
def insertQQd(expr, dofs):
    #syms = [phi_x(t), phi_y(t), phi_z(t)]
    # See pydy.codegen.ode_function_generator _lambdaify
    subs = {}
    vec_name = 'q'

    q   = DeferredVector('q')
    qd  = DeferredVector('qd')
    for i, sym in enumerate(dofs):
        dsym = diff(sym, dynamicsymbols._t)
        subs[sym]  = q[i]
        subs[dsym] = qd[i]
    print(subs)
    if hasattr(expr, '__len__'):
        print(expr.shape)
        return Matrix([me.msubs(e, subs) for e in expr]).reshape(expr.shape[0],expr.shape[1])
    else:
        return me.msubs(expr, subs)


def cleanPy(expr, dofs=None, varname='R', indent=0, replDict=None, noTimeDep=False, method='subs'):
    """ 
    Clean a python sympy expression and perform replacements:
      - DOFs       -> q[i]
      - Velocities -> qd[i]
      - Constants ->  p['name']
      - Inputs    ->  u['name']
    INPUTS:
       - replDict: a dictionary of replacements, where the key define the expression to be replaced.
                the value is either a replacement string, or a tuple (for matrices), for instance::
                     replDict['M_B'] = 'Mass'
                     replDict['M_T'] = ('MM_T', [0,0])
                      
    """
    # list of parameters inputs and dofs
    parameters = []
    inputs     = []
    sdofs      = []
    sdofsDeriv = []
    parametersProblem = []
    if dofs is not None:
        for idof,dof in enumerate(dofs):
            s=repr(dof)
            sdofs.append(s)
            sdofsDeriv.append('Derivative({},t)'.format(s))
            #sdofsDeriv.append('Derivative({}(t),t)'.format(s))
        # sorting dofs by decreasing string length to avoid replacements (=> phi_y before y)
        IDOF = np.argsort([len(s) for s in sdofs])[-1::-1]
        dofs=np.array(dofs)[IDOF]


    d_dofs     = [d.diff(dynamicsymbols._t) for d in dofs]
    dd_dofs    = [d.diff(dynamicsymbols._t,2)   for d in dofs]
    all_dofs   = list(dofs) + list(d_dofs) + list(dd_dofs)
    all_dofs_s = [repr(s) for s in all_dofs] 



    def cleanPyAtom(atom, dofs=None):
        symbols     = list(atom.free_symbols)
        try:
            symbols.remove(Symbol('t'))
        except:
            pass
        symbols_s   = [repr(s) for s in symbols]
        dyn_symbols = find_dynamicsymbols(atom)
        dyn_symbols_u  = [symb for symb in dyn_symbols if repr(symb) not in all_dofs_s]
        functions   = list(atom.atoms(Function))
        functions   = [f for f in functions if f.args != (dynamicsymbols._t,) ] # remove functions of time only ("dyn_symbols")
        functions   = [f for f in functions if f.func not in [sin, cos, exp, sqrt] ] # remove usual math functions
        functions_symb = [f.func for f in functions]
        functions_args = [f.args for f in functions]
        #print('>>> symbols       ', symbols)
        #print('>>> dyn_symbols   ', dyn_symbols)
        #print('>>> dyn_symbols_u ', dyn_symbols_u)
        #print('>>> functions     ', functions)
        #print('>>> functions_symb', functions_symb)
        #print('>>> functions_args', functions_args)
        if method=='subs':
            # States: accelerations, derivatives, position
            for idof, dof, ddof, dddof in zip(IDOF, dofs, d_dofs, dd_dofs):
                sdof=repr(dof).replace('(t)','')
                atom = atom.subs([(dddof, 'qqdd__'+ str(idof) + '___')])
                atom = atom.subs([(ddof , 'qqd__' + str(idof) + '___')])
                atom = atom.subs([(dof  , 'qq__'  + str(idof) + '___')])
            # Functions (like inputs but have specific arguments)
            for symb, args in zip(functions_symb, functions_args):
                ssymb=repr(symb)
                atom = atom.subs([(symb  , 'uu__' + ssymb + '___')])
                inputs.append(ssymb)
            # Inputs
            for symb in dyn_symbols_u:
                ssymb=repr(symb).replace('(t)','')
                atom = atom.subs([(symb  , 'uu__' + ssymb + '__u')])
                inputs.append(ssymb)

            if replDict is not None:
                for k,v in replDict.items():
                    if k in symbols_s:
                        i = symbols_s.index(k)
                        if len(v)==2:
                            symbols.pop(i)
                            symbols_s.pop(i)
                        else:
                            print('yams: cleanPy, TODO replace dict: {} {}'.format(k,v) )
                            raise Exception()
            # Parameters
            for symb in symbols:
                ssymb=repr(symb)
                atom = atom.subs([(symb  , 'pp__' + ssymb + '__p')])
                parameters.append(ssymb)

            s= repr(atom).replace(' ','')
            s = s.replace('qqdd__','qdd[')
            s = s.replace('qqd__' ,'qd[')
            s = s.replace('qq__'  ,'q[')
            s = s.replace('uu__'  ,'u[\'')
            s = s.replace('pp__'  ,'p[\'')
            s = s.replace('___', ']')
            s = s.replace('__p', '\']')
            if noTimeDep:
                s = s.replace('__u', '\']')
            else:
                s = s.replace('__u', '\'](t,q,qd)')

            # User overrides
            if replDict is not None:
                for k,v in replDict.items():
                    if len(v)==2:
                        if s.find(k)>=0:
                            if v[1] is not None:
                                s=s.replace(k, 'p[\'{}\'][{}]'.format(v[0],','.join([str(ii) for ii in v[1]]) ) )
                                parameters.append(v[0])
                            else:
                                s=s.replace(k, 'p[\'{}\']'.format(v[0]))
                                parameters.append(v[0])
                    else:
                        if s.find(k)>=0:
                            s=s.replace(k,v)


        elif method=='string':
            s=repr(atom).replace(' ','')
            # Replace parameters provide in replDict, as priority
            if replDict is not None:
                for k,v in replDict.items():
                    if len(v)==2:
                        if s.find(k)>=0:
                            if v[1] is not None:
                                s=s.replace(k, 'p[\'{}\'][{}]'.format(v[0],','.join([str(ii) for ii in v[1]]) ) )
                                parameters.append(v[0])
                            else:
                                s=s.replace(k, 'p[\'{}\']'.format(v[0]))
                                parameters.append(v[0])
                    else:
                        if s.find(k)>=0:
                            s=s.replace(k,v)

            if dofs is not None:
                # time derivatives first!!! important
                for idof,dof in zip(IDOF,dofs):
                    sdof=repr(dof)
                    #s=s.replace('Derivative({}(t),t)'.format(sdof),'qd[{}]'.format(idof))
                    s=s.replace('Derivative({},t)'.format(sdof),'qd[{}]'.format(idof))
                # then Dof
                for idof,dof in zip(IDOF,dofs):
                    sdof=repr(dof)
                    s=s.replace(sdof+'(t)','q[{}]'.format(idof))
                    s=s.replace(sdof,'q[{}]'.format(idof))
            for symb in symbols:
                ssymb=repr(symb)
                if not any([p.find(ssymb)>0 for p in parameters]):
                    if s.find(ssymb)>=0:
                        s=s.replace(ssymb,'p[\'{}\']'.format(ssymb))
                        parameters.append(ssymb)
                else:
                    parametersProblem.append(ssymb)

            for symb in dyn_symbols:
                ssymb=repr(symb).replace(' ','')
                if ssymb not in sdofsDeriv and ssymb not in sdofs:
                    ssymb=ssymb.replace('(t)','')
                    if noTimeDep:
                        s=s.replace(ssymb+'(t)','u[\'{}\']'.format(ssymb)) # When linearizing, the "u" is a u0
                    else:
                        s=s.replace(ssymb,'u[\'{}\']'.format(ssymb))
                    inputs.append(ssymb)
        else:
            raise NotImplementedError()
        return s


    try:
        dims=expr.shape
    except:
        dims=0
        return '{}{} = '.format(indent,varname) + cleanPyAtom(expr,dofs=dofs), list(set(parameters)), list(set(inputs)), sdofs

    indent =''.join([' ']*indent)

    s=''
    if len(dims)==1:
        s+='{}{} = np.zeros({})\n'.format(indent,varname,dims[0])
        for i in np.arange(dims[0]):
            s+='{}{}[{}] = '.format(indent,varname,i) + cleanPyAtom(expr[i], dofs=dofs) +'\n'
    elif len(dims)==2:
        s+='{}{} = np.zeros(({},{}))\n'.format(indent,varname, dims[0],dims[1])
        for i in np.arange(dims[0]):
            for j in np.arange(dims[1]):
                s+='{}{}[{},{}] = '.format(indent,varname,i,j) + cleanPyAtom(expr[i,j], dofs=dofs)+'\n'
    # removing duplicates
    parameters = list(set(parameters))
    inputs     = list(set(inputs))
    #print('parameters',parameters)
    #print('inputs    ',inputs)
    #print('sdofs     ',sdofs)
    parameters.sort()
    inputs.sort()
    parametersProblem = list(set(parametersProblem))
    if len(parametersProblem)>0:
        print('>>> Parameters not replaced for python: ', parametersProblem)
    return s, parameters, inputs, sdofs
