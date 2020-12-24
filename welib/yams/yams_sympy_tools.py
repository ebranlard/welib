""" 
Tools for yams_sympy
"""

import numpy as np
from sympy import latex, python, symarray, diff
from sympy import Symbol
from sympy import Function
from sympy.physics.vector import dynamicsymbols
from sympy.physics.mechanics.functions import find_dynamicsymbols

def frame_viz(Origin, frame, l=1, r=0.08):
    """ """
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
    X_viz_frame = VisualizationFrame(X_frame, X_center, X_shape)
    Y_viz_frame = VisualizationFrame(  frame, Y_center, Y_shape)
    Z_viz_frame = VisualizationFrame(Z_frame, Z_center, Z_shape)
    return X_viz_frame, Y_viz_frame, Z_viz_frame


def exprHasFunction(expr):
    """ return True if a sympy expression contains a function"""
    if hasattr(expr, 'atoms'): 
        return len(expr.atoms(Function))>0
    else:
        return False


def subs_no_diff(expr, subslist):
    """ 
    Perform sustitution in an expression, but not in the time derivatives
    Only works if Subslist is a simple lists of ( var, value) 


    see also: sympy.physics.mechanics.functions.msubs


    """
    # Set mapping between time derivatives and dummy variables
    time=dynamicsymbols._t
    Dummys=symarray('DUM', len(subslist))
    DT2Dummy=[]
    for i, (var,b) in enumerate(subslist):
        if len(var.atoms())!=1:
            print(var)
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


def linearize(expr, x0, order=1, sym=False):
    """ 
    Return a Taylor expansion of the expression at the operating point x0
    x0=[(x,0),(y,0),...]
    """
    # --- First protect expression from derivatives
    time=dynamicsymbols._t
    # Mapping derivative -> Dummy
    Dummys=symarray('DUM', len(x0))
    DT2Dummy=[(diff(v,time),Dummys[i]) for i,(v,x) in enumerate(x0)]
    Dummy2DT=[(b,a) for a,b in DT2Dummy]

    # Helper function, linearizes one expression (for matrices and lists..)
    def __linearize_one_expr(myexpr):
        # Remove derivatives
        myexpr = myexpr.subs(DT2Dummy)
        flin=myexpr.subs(x0)  # Order 0
        if order==0:
            return flin.subs(Dummy2DT)
        
        df = [(myexpr.diff(v1).subs(x0))*(v1-v10) for v1,v10 in x0]
        flin+=sum(df)
        if order==1:
            return flin.subs(Dummy2DT)
        
        df2 = [ (myexpr.diff(v1)).diff(v2).subs(x0)*(v1-v10)*(v2-v20) for v1,v10 in x0 for v2,v20 in x0]
        flin+=sum(df2)/2
        if order==2:
            return flin.subs(Dummy2DT)
        
        df3 = [ (myexpr.diff(v1)).diff(v2).diff(v3).subs(x0)*(v1-v10)*(v2-v20)*(v3-v30) for v1,v10 in x0 for v2,v20 in x0 for v3,v30 in x0]
        flin+=sum(df3)/6
        if order==3:
            return flin.subs(Dummy2DT)
        raise NotImplementedError('Higher order')
    
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
    Perform small angle approximation of an expresion by linearizaing about the "0" operating point of each angle

    angle_list: list of angle to be considered small:  [phi,alpha,...]

    """
    # operating point
    x0=[(angle,0) for angle in angle_list]
    return linearize(expr, x0, order=order, sym=sym)



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



def cleanPy(expr, dofs=None, varname='R', indent=0, replDict=None, noTimeDep=False):
    """ 
    Clean a python sympy expression and perform replacements:
      - DOFs       -> q[i]
      - Velocities -> qd[i]
      - Constants ->  p['name']
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

    def cleanPyAtom(atom, dofs=None):
        symbols     = atom.free_symbols
        try:
            symbols.remove(Symbol('t'))
        except:
            pass
        dyn_symbols = find_dynamicsymbols(atom)
        #print('>>> symbols', symbols)
        #print('>>> Dynsymbols', dyn_symbols)
        s=repr(atom).replace(' ','')

        # Replace
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
    parameters.sort()
    inputs.sort()
    parametersProblem = list(set(parametersProblem))
    if len(parametersProblem)>0:
        print('>>> Parameters not replaced for python: ', parametersProblem)
    return s, parameters, inputs, sdofs
