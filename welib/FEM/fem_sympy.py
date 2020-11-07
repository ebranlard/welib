
import welib.FEM.frame2d as f2d
import welib.FEM.frame3d as f3d
import sympy
from sympy import symbols, Symbol, simplify
from sympy import diff, integrate 
from sympy import zeros

# --------------------------------------------------------------------------------}
# --- General  
# --------------------------------------------------------------------------------{
def stiffnessMatrixFromPot(Epot,u_dof):
    """ Return stiffness matrices from an experssion of potential energy and a list of DOF"""
    n=len(u_dof)
    Ke= zeros(n,n)
    for i in range(n):
        for j in range(n):
            Ke[i,j] = Epot.diff(u_dof[j]).diff(u_dof[i])
    return Ke

def stiffnessMatrixFromShapeFunctions(NN, x, L, EI):
    """
    For beam bending:
         ke = EI \int_0^L B^T B dx  where B = d^2N/dx^2 (curvature, or strain-displacement)
    TODO other
    """
    # --- "Slope"
    SS = NN.diff(x)
    # --- "Curvature" - Strain-displacement vector
    BB = NN.diff(x).diff(x)
    
    Ke = EI* integrate(BB.T * BB, (x, 0, L))
    
    return Ke

def massMatrixFromShapeFunctions(NN, x, L, m, coeffs=None):
    """
    Mass matrix for a Beam 

    If coefficients are provided, NN_rho is NN with the rows scaled
    Typical scalings are [m/L m/L m/L rx,m/L]
    
    """
    if coeffs is None:
        Me = m/L*integrate(NN.T * NN, (x, 0, L))
    else:
        if len(coeffs)!=NN.shape[0]:
            raise Exception('List of scaling coefficients need to be the same length as number of rows of NN')

        # Scaling the rows by the coeff provided
        NN_rho = NN.copy()
        for i,coeff in enumerate(coeffs):
            for j in range(NN.shape[1]):
                NN_rho[i,j] *= coeff

        Me = integrate(NN.T * NN_rho, (x, 0, L))
    return Me



# --------------------------------------------------------------------------------}
# --- Frame 2d 
# --------------------------------------------------------------------------------{
def frame2d_U_B(h,x,L, EI, u2=Symbol('u_2'), u3=Symbol('u_3'), u5=Symbol('u_5'), u6=Symbol('u_6')):
    """ Bending strain energy"""
    hpp=h(x,u2,u3, u5, u6, L).diff(x).diff(x)
    return EI/2*integrate( hpp * hpp, (x,0,L))

def frame2d_U_A(u, x, L, EA, u1=Symbol('u_1'), u4=Symbol('u_4')):
    """ Longitudinal strain energy"""
    up = u(x,u1,u4,L).diff(x)
    return EA/2*integrate( up * up, (x,0,L))

def frame2d_U_G(h,x,L,T, u2=Symbol('u_2'), u3=Symbol('u_3'), u5=Symbol('u_5'), u6=Symbol('u_6')):
    """ Potential energy due to axial loads """
    hp = h(x,u2,u3, u5, u6,L).diff(x)
    return T/2*integrate( hp*hp , (x,0,L))

def frame3d_U_B(fShape, x, L, EI):
    """ Bending-along x strain energy"""
    hpp=fShape.diff(x).diff(x)
    return EI/2*integrate( hpp * hpp, (x,0,L))

def frame3d_U_A(fShape, x, L, EA):
    """ Longitudinal strain energy"""
    up = fShape.diff(x)
    return EA/2*integrate( up * up, (x,0,L))


def frame3d_U_T(fShape, x, L, GK):
    """ Torsional strain energy"""
    up = fShape.diff(x)
    return GK/2*integrate( up * up, (x,0,L))

def frame3d_U_G(fShape,x,L,T):
    """ Potential energy due to axial loads """
    hp = fShape.diff(x)
    return T/2*integrate( hp*hp , (x,0,L))


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
        '\\\\':'\\\\ \n' 
        }
    for k in D_rep.keys():
        s=s.replace(k,D_rep[k])
    return s
