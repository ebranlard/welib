""" 
Tools for yams_sympy
"""

import numpy as np
from sympy import latex, python

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



def linearize(expr, x0, order=1):
    """ 
    Return a Taylor expansion of the expression at the operating point x0
    x0=[(x,0),(y,0),...]
    """
    flin=expr.subs(x0)  # Order 0
    if order==0:
        return flin
    
    df = [(expr.diff(v1).subs(x0))*(v1-v10) for v1,v10 in x0]
    flin+=sum(df)
    if order==1:
        return flin
    
    df2 = [ (expr.diff(v1)).diff(v2).subs(x0)*(v1-v10)*(v2-v20) for v1,v10 in x0 for v2,v20 in x0]
    flin+=sum(df2)/2
    if order==2:
        return flin
    
    df3 = [ (expr.diff(v1)).diff(v2).diff(v3).subs(x0)*(v1-v10)*(v2-v20)*(v3-v30) for v1,v10 in x0 for v2,v20 in x0 for v3,v30 in x0]
    flin+=sum(df3)/6
    if order==3:
        return flin

    raise NotImplementedError('Higher order')


def smallAngleApprox(expr, angle_list, order=1, sym=True):
    """ 
    Perform small angle approximation of an expresion by linearizaing about the "0" operating point of each angle

    angle_list: list of angle to be considered small:  [phi,alpha,...]

    """

    # operating point
    x0=[(angle,0) for angle in angle_list]

    try:
        ndim=len(expr.shape)
    except:
        ndim=0

    expr_lin=expr*0

    if ndim==0:
        expr_lin = linearize(expr, x0, order=order)

    elif ndim==1:
        for i in range(expr.shape[0]):
            expr_lin[i] = linearize(expr[i], x0, order=1)

    elif ndim==2:
        if expr.shape[0] != expr.shape[1]:
            sym=False

        for i in range(expr.shape[0]):
            for j in range(expr.shape[1]):
                if (j<=i and sym) or (not sym):
                    expr_lin[i,j] = linearize(expr[i,j], x0, order=1)
                if sym:
                    expr_lin[j,i] = expr_lin[i,j]
    else:
        raise NotImplementedError()

    return expr_lin



# --------------------------------------------------------------------------------}
# --- Misc sympy utils 
# --------------------------------------------------------------------------------{
def cleantex(s):
    """ clean a latex expression 
    
    Example:
        print( cleantex( vlatex(MM, symbol_names={Jxx:'J_{xx}'})))
    """
    D_rep={
        '\\operatorname{sin}':'\\sin',
        '\\operatorname{cos}':'\\cos',
        '\\left(\\theta\\right)':'\\theta',
        '\\left[\\begin{matrix}':'\\begin{bmatrix}\n', 
        '\\end{matrix}\\right]':'\n\\end{bmatrix}\n' ,
        '\\operatorname{q_{T1}}':'q_{T1}',
        '\\left(t \\right)'    : '(t)',
        '{(t)}': '(t)',
        'L_{T}':'L_T',
        'M_{N}':'M_N',
        '_{x}':'_x',
        '_{y}':'_y',
        '_{z}':'_z',
        'phi_x(t)':'phi_x',
        'phi_y(t)':'phi_y',
        'phi_z(t)':'phi_z',
        'q_{T1}(t)': 'q_{T1}',
        '\\sin{\\left(\\phi_x \\right)}':'\\sin\\phi_x',
        '\\sin{\\left(\\phi_y \\right)}':'\\sin\\phi_y',
        '\\sin{\\left(\\phi_z \\right)}':'\\sin\\phi_z',
        '\\cos{\\left(\\phi_x \\right)}':'\\cos\\phi_x',
        '\\cos{\\left(\\phi_y \\right)}':'\\cos\\phi_y',
        '\\cos{\\left(\\phi_z \\right)}':'\\cos\\phi_z',
        '\\\\':'\\\\ \n' 
        }
    for k in D_rep.keys():
        s=s.replace(k,D_rep[k])
    return s


def saveTex(expr, filename):
    print('Exporting to {}'.format(filename))
    with open(filename,'w') as f:
        f.write(cleantex(latex(expr)))
    print('Done')
