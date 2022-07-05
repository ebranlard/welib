"""

"Standard input data (SID)" used to define a flexible body.

The following references are used:

 [1]  Oskar Wallrapp, 1993, "Standard Input Data of Flexible Members in Multibody Systems"
      Advanced Multibody System Dynamics, 445-450

 [2]  Richard Schwertassek,  Oskar Wallrapp
      "Dynamik Flexibler Mehrkoerpersysteme : Methoden Der Mechanik 




"""

import os
import numpy as np
from welib.yams.utils import skew 

class Refmod(object):
    def __init__(self,nq):
        self.mass    = np.nan    # Mass of body
        self.nelastq = nq # Number of mode shapes
        self.ielastq = ['']*nq # Name of mode shapes
    def __repr__(self):
        s='refmod\n'
        s+='    mass           = {}\n'.format(self.mass)
        s+='    nelastq        = {}\n'.format(self.nelastq)
        s+='\n'.join(['    ielastq({}) = {}'.format(i+1,v) for i,v in enumerate(self.ielastq)])
        s+='\nend refmod\n'
        return s

class Taylor(object):
    """ Taylor modal expansion """
    def __init__(self, order, nrow, ncol, nq, nqn, structure=3, name=''):
        self.order     = order     # order of expansion
        self.nrow      = nrow      # number of rows
        self.ncol      = ncol      # number of columns
        self.nq        = nq        # dimension of q, number of shape functions
        self.nqn       = nqn       # dimension of higher order
        self.structure = structure # 3                                  = matrix, 2 = symmetric, 1 = diagonal
        self.M0=np.zeros((nrow,ncol))
        if self.order==1:
            self.M1=np.zeros((nrow,ncol,nq))
        self.name = name # name of taylor extension

    def __repr__(self):
        s='{}\n'.format(self.name)
        s+='    order        = {}\n'.format(self.order)
        s+='    nrow         = {}\n'.format(self.nrow )
        s+='    ncol         = {}\n'.format(self.ncol )
        s+='    nq           = {}\n'.format(self.nq   )
        s+='    nqn          = {}\n'.format(self.nqn  )
        s+='    structure    = {}\n'.format(self.structure)
        for i in np.arange(self.nrow):
            for j in np.arange(self.ncol):
                if self.structure==2 and j>i:
                    continue
                s+='    m0( {}, {})    = {:.6f}\n'.format(i+1,j+1,self.M0[i,j])
        if self.order==1:
            for i in np.arange(self.nrow):
                for j in np.arange(self.ncol):
                    for k in np.arange(self.nq):
                        s+='    m1( {}, {}, {})    = {:.6f}\n'.format(i+1,j+1,k+1,self.M1[i,j,k])
        s+='end {}\n'.format(self.name)
        return s

    def setorder(self,order):
        self.order=order
        if self.order==1:
            if not hasattr(self,'M1'):
                self.M1=np.zeros((self.nrow,self.ncol,self.nq))

class Node(object):
    """ """
    def __init__(self, nq, name=0, sigma=None):
        self.rframe = '' # name of reference frame
        self.origin = Taylor(1, 3, 1, nq, 0, name='origin') # position of node
        self.phi    = Taylor(0, 3, nq,nq, 0, name='phi'   ) # translational mode shapes
        self.psi    = Taylor(0, 3, nq,nq, 0, name='psi'   ) # rotational mode shapes
        self.APmod  = Taylor(1, 3, 3, nq, 0, name='AP'    ) # relative orientation matrix of node fixed frames
        if sigma is not None:
            self.sigma  = Taylor(1, 6, 1, nq, 0, name='sigma' )  # 
        else:
            self.sigma  = None
        self.name = name

    def __repr__(self):
        s='new node = {}\n'.format(self.name)
        s+='rframe = {}\n'.format(self.rframe)
        s+= str(self.origin)
        s+= str(self.phi   )
        s+= str(self.psi   )
        s+= str(self.APmod )
        if self.sigma is not None:
            s+= str(self.sigma )
        s+='end node\n'
        return s

class SID(object):
    """ """
    def __init__(self, nq, GrOrder=0):
        self.refmod = Refmod(nq=nq)
        self.frame  = []  # position, orientation, and mode shape matrices at each nodes
        self.mdCM   = Taylor(1,3 ,1   ,nq,0,   name = 'mdCM') # matrix of mass times position of center of mass
        self.J      = Taylor(1,3 ,3   ,nq,0,2, name = 'J' ) # mass moment of inertia
        self.Ct     = Taylor(1,nq,3   ,nq,0,   name = 'Ct') # coupling matrix with translational coordinates
        self.Cr     = Taylor(1,nq,3   ,nq,0,   name = 'Cr') # coupling matrix with rotational coordinates
        self.Me     = Taylor(0,nq,nq  ,0 ,0,2, name = 'Me') # Modal mass matrix
        self.Gr     = Taylor(GrOrder,3 ,3*nq,nq,0,   name = 'Gr') # Gyroscopic matrix for rotational coordinates
        self.Ge     = Taylor(0,nq,3*nq,0 ,0,   name = 'Ge') # Gyroscopic matrix for modal coordinates
        self.Oe     = Taylor(1,nq,6   ,nq,0,   name = 'Oe') # Centrifugal matrix for modal coordinates
        self.ksigma = Taylor(0,nq,1   ,nq,0,0, name = 'ksigma') # modal force vector of reference stresses
        self.Ke     = Taylor(0,nq,nq  ,0 ,0,2, name = 'Ke') # modal stiffness matrix due to deformations
        self.De     = Taylor(0,nq,nq  ,0 ,0,0, name = 'De') # modal dampoing matrix due to deformations

    def __repr__(self):
        s =''
        s+='part\n'
        s+='new modal\n'
        s+=str(self.refmod)+'\n'
        s+='new frame\n'
        for i,f in enumerate(self.frame):
            f.name=i+1
            s+=str(f)
        s+='end frame\n'
        s+=str(self.mdCM  )
        s+=str(self.J     )
        s+=str(self.Ct    )
        s+=str(self.Cr    )
        s+=str(self.Me    )
        s+=str(self.Gr    )
        s+=str(self.Ge    )
        s+=str(self.Oe    )
        s+=str(self.ksigma)
        s+=str(self.Ke    )
        s+=str(self.De    )
        s+='end modal\n'
        s+='end part\n'
        return s



# --------------------------------------------------------------------------------}
# ---  Converters
# --------------------------------------------------------------------------------{
def Beam2SID(xNodes, Imodes, m, Iy, Iz, G=None, Kv=None, A=None, E=None, phi=None):
    """ High level interface, converts a beam to SID"""
    from welib.FEM.fem_beam import cbeam_assembly_frame3dlin
    from welib.FEM.fem_beam import applyBC
    from welib.FEM.fem_beam import orthogonalizeModePair, normalize_to_last
    from welib.system.eva import eig

    # --- Assembling FEM beam model with frame3dlin
    MM, KK, xNodes, DCM, Elem2Nodes, Nodes2DOF, Elem2DOF = cbeam_assembly_frame3dlin(xNodes, m, Iy, Iz=Iz, A=A, Kv=Kv, E=E, G=G, phi=phi)

    # --- Constraints/ BC
    MMr, KKr, Tr,_,_ = applyBC(MM, KK, Elem2Nodes, Nodes2DOF, BC_root=[0,0,0,0,0,0], BC_tip=[1,1,1,1,1,1])
    iStart= 0; 

    # --- Selections, and orthogonlization of modes 
    # TODO orthogonalization options
    [Q, freq]= eig(KKr, MMr, freq_out=True)
    Imodes_all=[]
    for Im in Imodes:
        if type(Im) is tuple:
            if len(Im) == 2:
                Q[:,Im[0]],Q[:,Im[1]] = orthogonalizeModePair(Q[:,Im[0]],Q[:,Im[1]], iStart)
                Imodes_all.append(Im[0])
                Imodes_all.append(Im[1])
            else:
                raise Exception('Only tuples of length 2 accepted for orthogonalization')
        else:
            Imodes_all.append(Im)

    Q= normalize_to_last(Q, Imodes_all, iStart);
    # Selecting modes
    if len(Imodes)>0:
        Se= Tr.dot(Q[:, Imodes_all]) # nDOF_tot x nShapes
    else:
        Se= Tr # All

    # --- Export Modes
    #Mode1=np.column_stack([Se[0::6,0], Se[1::6,0], Se[2::6,0], Se[3::6,0], Se[4::6,0], Se[5::6,0]])
    #Mode2=np.column_stack([Se[0::6,1], Se[1::6,1], Se[2::6,1], Se[3::6,1], Se[4::6,1], Se[5::6,1]])
    #M=np.column_stack([xNodes[2,:], Mode1, Mode2])
    #np.savetxt('Modes.csv',M)

    sid = FEM2SID(xNodes, A, E, m, MM, KK, MMr, KKr, Tr, Se, DCM, Elem2Nodes, Nodes2DOF, Elem2DOF)

    # Additional info
    sid.Q    = Q
    sid.freq = freq
    return sid



def FEM2SID(xNodes, A, E, m, MM, KK, MMr, KKr, Tr, Se, DCM, Elem2Nodes, Nodes2DOF, Elem2DOF):
    """ Mid level interface - Convert FEM data to sid"""
    from welib.FEM.fem_beam import generalizedMassMatrix, shapeIntegrals
    from welib.FEM.fem_beam import geometricalStiffening
    from numpy.linalg import inv

    # --- Generalized mass matrix
    Mtt, J0, Mrt, Mgt, Mgr, Mgg, St, Sr= generalizedMassMatrix(xNodes, MM, Se)
    Ct0_ = (Tr.T).dot(MM).dot(St) # Mode mass matrix for all modes

    # --- Shape integrals
    C3, Kr, C4, KFom_ab, Kom, Kom0, Kom0_ = shapeIntegrals(xNodes, Nodes2DOF, Elem2Nodes, Elem2DOF, DCM, m, Se, Sr, Tr)

    # --- Stiffening terms
    Kinv= Tr.dot(inv(KKr)).dot(Tr.T);
    GKg= geometricalStiffening(xNodes, Kinv, Tr, Se, Nodes2DOF, Elem2Nodes, Elem2DOF, DCM, E, A, Kom0_, Ct0_)

    # --- Convert to SID 
    sid = FEMBeam2SID(Mtt, J0, Mrt, Mgt, Mgr, Mgg, KK, xNodes, DCM, Se, Kr, Kom0, Kom, C4=C4, GKg=GKg)

    # --- Additional info, shape functions and modes
    sid.C3    = C3
    sid.C4    = C4
    sid.Kr    = Kr
    sid.Kom   = Kom
    sid.Kom0  = Kom0
    sid.Kom0_ = Kom0_
    sid.Sr    = Sr
    sid.St    = St
    sid.Se    = Se
    sid.Mtt=Mtt
    sid.Mrt=Mrt
    sid.Mgt=Mgt
    sid.Mgr=Mgr
    sid.GKg=GKg

    return sid


# --- Populating SID
def FEMBeam2SID(Mtt, J0, Mrt, Mgt, Mgr, Mgg, KK, xNodes, DCM, Se, Kr, Kom0, Kom, C4=None, GKg=None):
    """ Low level interface, computes SID from a set of intermediate FEM variables"""

    # TODO TODO TODO call shapeIntegral2SID
    # 


    # TODO take MM, KK, Tr as inputs
    assert(xNodes.shape[0]==3)
    GKg = dict() if GKg is None else GKg
    nqk      = 6               # Number of DOF per nodes
    nNodes   = xNodes.shape[1]
    nDOF_tot = nNodes*nqk      # Total number of DOF without constraint (BC)

    iMaxDim = np.argmax(np.max(np.abs(xNodes),axis=1)-np.min(np.abs(xNodes),axis=1)) 

    # See [2] Table 6.9 S. 346
    nq  = Se.shape[1]
    sid=SID(nq=nq)
    sid.refmod.mass= Mtt[0,0]
    for i in np.arange(nq):
        sid.refmod.ielastq[i]= 'Eigen Mode {:4d}'.format(i)

    for i in np.arange(nNodes):
        f=Node(nq=nq)
        f.name= i
        f.rframe= 'body ref'
        Phi= np.zeros((3, nDOF_tot))
        Psi= np.zeros((3, nDOF_tot))
        Phi[:, i*nqk+0:i*nqk+3]= np.eye(3)
        Psi[:, i*nqk+3:i*nqk+6]= np.eye(3)

        # Origin
        f.origin.M0[:,0]= xNodes[:,i]
        for j in np.arange(nq):
            f.origin.M1[:, 0, j] = Phi.dot(Se[:, j]);
        # Node translations (6.46),(6.80), (6.337), (5.122), (6.415)
        f.phi.M0= Phi.dot(Se)
        # Stiffnening due to force
        if 'Fend' in GKg.keys():
            if i==nNodes-1:
                f.phi.setorder(1)
                f.phi.M1[iMaxDim, :, :]= GKg['Fend']

        # Node rotations (6.81), (6.337)
        f.psi.M0= Psi.dot(Se);
        if i==nNodes-1:
            f.APmod.M0= DCM[:,:,i-1] # orientation
        else:
            f.APmod.M0= DCM[:,:,i] # orientation
        for j in np.arange(nq):
            f.APmod.M1[:, :, j]= skew(f.psi.M0[:, j]);

        # --- Pretension
        #if i<nNodes-1:
        #  e= i;
        #else:
        #  e= i-1;
        #  Ke= T{e}' .dot(K{e}).dot(T{e})
        #  if i<nk:
        #      Ke= -Ke;
        #  for j= 1:nq:
        #      f.sigma.M1(:, 1, j)= Ke((i-1)*nqk+(1:6), :)*Se(:, j); % left sided nodal force except last node
        sid.frame.append(f)

    # --- mdCM
    # NOTE: Mrt (mc0) is anti symmetric and contains m*x_COG
    # The math below is to extract m*x_COG rom Mrt
    sid.mdCM.M0= np.array([
        np.sum(skew([0.5, 0 , 0  ])*Mrt), 
        np.sum(skew([0,  0.5, 0  ])*Mrt), 
        np.sum(skew([0,   0 , 0.5])*Mrt)
        ]).reshape(3,1)
    for j in np.arange(nq):
        sid.mdCM.M1[:, 0, j]= Mgt[j, :].T # Ct0 Mgt

    # --- J
    sid.J.M0= J0
    for j in np.arange(nq):
        sid.J.M1[:, :, j]= -C4[:, :, j] - C4[:, :, j].T;

    # --- Ct
    sid.Ct.M0= Mgt
    if 't_ax' in GKg.keys():
        sid.Ct.M1[:, iMaxDim, :]= GKg['t_ax'];
    # --- Cr
    sid.Cr.M0= Mgr
    sid.Cr.M1[:, 0, :]= Kr[0][:,:]; # nq x nq
    sid.Cr.M1[:, 1, :]= Kr[1][:,:]; 
    sid.Cr.M1[:, 2, :]= Kr[2][:,:]; 
    # --- Me
    sid.Me.M0= Mgg;

    # --- Gr 
    # see [2] (6.403) p. 339
    for j in np.arange(nq):
        sid.Gr.M0[0:3, 3*j:3*j+3]= -2*C4[:,:,j];  # columns concatenation 3 per shapes

    # --- Ge = 2 \int Phi^t phi_j~^t dm  = [2C5']
    # Ge Taylor(0,nq,3*nq,0 ,0) # Gyroscopic matrix for modal coordinates
    for j in np.arange(nq):
          # see [2] (6.405) p. 340 = 2*C5'
          M0j = 2*np.vstack([Kr[0,0:nq, j],Kr[1,0:nq, j],Kr[2,0:nq, j]])  # 3 x nq
          sid.Ge.M0[0:nq, 3*j:3*j+3]= M0j.T;  # columns concatenation 3 per shapes
    
    # --- Oe 
    # see [2] (6.407) 
    # M0: nq x 6
    # M1: nq x 6 x nq
    sid.Oe.M0= Kom0  # nq x 6
    sid.Oe.M1[:, 0, :]= Kom[0]   # nq x nq
    sid.Oe.M1[:, 1, :]= Kom[1]  
    sid.Oe.M1[:, 2, :]= Kom[2]  
    sid.Oe.M1[:, 3, :]= Kom[3]  
    sid.Oe.M1[:, 4, :]= Kom[4]  
    sid.Oe.M1[:, 5, :]= Kom[5]  
    if 'omxx' in GKg.keys():
        sid.Oe.M1[:, 0, :]+= GKg['omxx'];  # nq x nq
        sid.Oe.M1[:, 1, :]+= GKg['omyy'];
        sid.Oe.M1[:, 2, :]+= GKg['omzz'];
        sid.Oe.M1[:, 3, :]+= GKg['omxy'];
        sid.Oe.M1[:, 4, :]+= GKg['omxz'];
        sid.Oe.M1[:, 5, :]+= GKg['omyz'];

    # ---  sid.Ke 
    sid.Ke.M0= (Se.T).dot(KK).dot(Se);

    # --- De
    #sid.De.M0= [];
    return sid



# --- Populating SID
def shapeIntegrals2SID(p, consistency=''):
    """
    p: dictionnary with keys as returned by shapeIntegrals in flexibility

    """
    from welib.yams.utils import skew 

    nq     = p['C1'].shape[1]
    nNodes = p['s_G'].shape[1]
    sid=SID(nq=nq, GrOrder=1)
    sid.refmod.mass= p['mass']
    for i in np.arange(nq):
        sid.refmod.ielastq[i]= 'Eigen Mode {:4d}'.format(i+1)

    for i in np.arange(nNodes):
        f=Node(nq=nq)
        f.name= i
        f.rframe= 'body ref'
        # Origin
        f.origin.M0[:,0]= p['s_G'][:,i]
        for j in np.arange(nq):
            f.origin.M1[:, 0, j] = p['U'][j][:,i]
        # Phi - Node translations/shape functions (6.46),(6.80), (6.337), (5.122), (6.415)
        for j in np.arange(nq):
            f.phi.M0[:,j] = p['U'][j][:,i]
        ## Stiffnening due to force
        #if 'Fend' in GKg.keys():
        #    if i==nNodes-1:
        #        f.phi.setorder(1)
        #        f.phi.M1[iMaxDim, :, :]= GKg['Fend']
        # Psi - Node rotations/slope (6.81), (6.337)
        for j in np.arange(nq):
            f.psi.M0[:,j] = p['dU'][j][:,i]

        # Orientation
        f.APmod.M0 = np.eye(3) # keep original orientation for now E{i_elem};
        for j in np.arange(nq):
            f.APmod.M1[:, :, j]= skew(f.psi.M0[:, j]);

        # --- Pretension
        #f.sigma.M0 = np.zeros((6, 1)) # no pretension for now.
        sid.frame.append(f)

    # --- mdCM
    sid.mdCM.M0= p['mc0'].reshape(3,1)
    for j in np.arange(nq):
        sid.mdCM.M1[:, 0, j]= p['C1'][:, j]

    # --- J
    sid.J.M0= p['J0']
    #C4 is (nf,3,3) # NOTE: different convention that FEMBeam
    assert(p['C4'].shape==(nq,3,3))
    for j in np.arange(nq):
        sid.J.M1[:, :, j]= -p['C4'][j, :, :] - p['C4'][j, :, :].T;

    # --- Ct
    # NOTE: K0t is geometrical stiffening due to translational acceleration
    sid.Ct.M0= p['C1'].T
    sid.Ct.M1[:, 0, :] = p['K0t'][0,:,:] # Gkg['t_ax'], K0t is (3,nf,nf)
    sid.Ct.M1[:, 1, :] = p['K0t'][1,:,:]
    sid.Ct.M1[:, 2, :] = p['K0t'][2,:,:]

    # --- Cr
    # NOTE: K0r is geometrical stiffening due to rotational acceleration
    sid.Cr.M0= p['C2'].T
    #Kr is (3,nf,nf)
    sid.Cr.M1[:, 0, :]= p['Kr'][0][:,:]  + p['K0r'][0][:,:] # nq x nq
    sid.Cr.M1[:, 1, :]= p['Kr'][1][:,:]  + p['K0r'][1][:,:] 
    sid.Cr.M1[:, 2, :]= p['Kr'][2][:,:]  + p['K0r'][2][:,:] 

    # --- Me
    sid.Me.M0= p['Me']

    # --- Gr 
    # see [2] (6.403) p. 339
    for j in np.arange(nq):
        sid.Gr.M0[0:3, 3*j:3*j+3]= -2*p['C4'][j,:,:];  # columns concatenation 3 per shapes
    for j in np.arange(nq):
        for k in np.arange(nq):
            sid.Gr.M1[:, 3*j:3*j+3, k]= -2*p['C6'][k, j]


    # --- Ge = 2 \int Phi^t phi_j~^t dm  = [2C5']
    # Ge Taylor(0,nq,3*nq,0 ,0) # Gyroscopic matrix for modal coordinates
    # see [2] (6.405) p. 340 = 2*C5'
    for j in np.arange(nq):
          M0j = 2*np.vstack([p['Kr'][0,0:nq, j],p['Kr'][1,0:nq, j],p['Kr'][2,0:nq, j]])  # 3 x nq # TODO TODO CHECK THIS
          sid.Ge.M0[0:nq, 3*j:3*j+3]= M0j.T;  # columns concatenation 3 per shapes
    
    # --- Oe 
    # NOTE: K0omega is geometric stiffening due to centrifugal acceleration 
    # see [2] (6.407) 
    # M0: nq x 6
    # M1: nq x 6 x nq
    # Oe.M0 = Kom0
    for j in np.arange(nq):
        sid.Oe.M0[j, :]= [p['C4'][j,0,0], p['C4'][j,1,1], p['C4'][j,2,2], p['C4'][j,0,1]+p['C4'][j,1,0], p['C4'][j,1,2]+p['C4'][j,2,1], p['C4'][j,2,0]+p['C4'][j,0,2]]

    sid.Oe.M1_base = np.zeros(sid.Oe.M1.shape)
    sid.Oe.M1_geom = np.zeros(sid.Oe.M1.shape)

    sid.Oe.M1_base[:, 0, :]= p['Komega'][0, 0] 
    sid.Oe.M1_base[:, 1, :]= p['Komega'][1, 1] 
    sid.Oe.M1_base[:, 2, :]= p['Komega'][2, 2] 
    sid.Oe.M1_base[:, 3, :]= p['Komega'][0, 1] + p['Komega'][1, 0] 
    sid.Oe.M1_base[:, 4, :]= p['Komega'][1, 2] + p['Komega'][2, 1] 
    sid.Oe.M1_base[:, 5, :]= p['Komega'][0, 2] + p['Komega'][2, 0] 

    sid.Oe.M1_geom[:, 0, :]= p['K0omega'][0, 0] 
    sid.Oe.M1_geom[:, 1, :]= p['K0omega'][1, 1] 
    sid.Oe.M1_geom[:, 2, :]= p['K0omega'][2, 2] 
    sid.Oe.M1_geom[:, 3, :]= p['K0omega'][0, 1] + p['K0omega'][1, 0] 
    sid.Oe.M1_geom[:, 4, :]= p['K0omega'][1, 2] + p['K0omega'][2, 1] 
    sid.Oe.M1_geom[:, 5, :]= p['K0omega'][0, 2] + p['K0omega'][2, 0] 

    sid.Oe.M1 = sid.Oe.M1_base + sid.Oe.M1_geom 
    #sid.Oe.M1 = sid.Oe.M1_base 
    #sid.Oe.M1 = sid.Oe.M1_geom 

    # ---  sid.Ke 
    sid.Ke.M0= p['Ke']

    # --- De
    #sid.De.M0= [];

    # --- remove some off-diagonal terms to conform to FAST approach
    if consistency=='OpenFAST':
        sid.Me.M0 = np.diag(np.diag(sid.Me.M0))
        for a in [0,1,2]:
            sid.Ct.M1[:, a, :]= np.diag(np.diag(np.squeeze(sid.Ct.M1[:, a, :])))
        for m in np.arange(6):
            if m<3:
                sid.Oe.M1[:, m, :] += - p['K0omega'][m, m] + np.diag(np.diag(p['K0omega'][m, m]))
            else:
                c= m-3;
                d= np.mod(c+1, 3)
                sid.Oe.M1[:, m, :] += - p['K0omega'][c, d] - p['K0omega'][c, d].T + 2*np.diag(np.diag(p['K0omega'][c, d]))

    return sid


# --------------------------------------------------------------------------------}
# --- Beam Shape functions to SID 
# --------------------------------------------------------------------------------{
def BeamShapes2SID(s_G, s_span, m, EI, U, dU, ddU, int_method='trapz', damping=None, consistency='', shapeIntegrals=False):
    """ 
    Using generalized parameters as computed by GMBeam and GKBeam
    to populate a SID model

      int_method: 'trapz', 'OpenFAST'
    """
    # --- Method 1 using shape integrals
    if shapeIntegrals:
        from welib.yams.flexibility import shapeIntegrals

        p = shapeIntegrals(s_G, s_span, m, U, dU, ddU, method=int_method, EI=EI)
        sid = shapeIntegrals2SID(p, consistency=consistency)
        sid.p=p
        return sid

    # --- Method 2 relying on Generalized functions
    from welib.yams.flexibility import GMBeam, GKBeam, GKBeamStiffnening, GKBeamStiffneningSplit
    MM, IT = GMBeam(s_G, s_span, m, U, rot_terms=True, method=int_method, main_axis='z', M1=True) 
    Gr, Ge, Oe, Oe6 = IT['Gr'], IT['Ge'], IT['Oe'], IT['Oe6']

    KK  = GKBeam(s_span, EI, ddU, bOrth = False,                method = int_method)
    pKg = GKBeamStiffneningSplit(s_G, s_span, dU, m , main_axis='z', method = int_method)

    nq = U.shape[0]
    nNodes = U.shape[2]
    if int_method=='OpenFAST':
        I = np.arange(nNodes)[1:-1]
    else:
        I = np.arange(nNodes)

    sid=SID(nq=nq, GrOrder=1)

    sid.refmod.mass= MM[0,0]
    for i in np.arange(nq):
        sid.refmod.ielastq[i]= 'Eigen Mode {:4d}'.format(i+1)

    for ii, i in enumerate(I):
        f=Node(nq=nq)
        f.name= ii
        f.rframe= 'body ref'
        # Origin
        f.origin.M0[:,0]= s_G[:,i]
        for j in np.arange(nq):
            f.origin.M1[:, 0, j] = U[j][:,i]
        # Phi - Node translations/shape functions (6.46),(6.80), (6.337), (5.122), (6.415)
        for j in np.arange(nq):
            f.phi.M0[:,j] = U[j][:,i]
        ## Stiffnening due to force
        #if 'Fend' in GKg.keys():
        #    if i==nNodes-1:
        #        f.phi.setorder(1)
        #        f.phi.M1[iMaxDim, :, :]= GKg['Fend']
        # Psi - Node rotations/slope (6.81), (6.337)
        for j in np.arange(nq):
            f.psi.M0[:,j] = dU[j][:,i]

        # Orientation
        f.APmod.M0 = np.eye(3) # keep original orientation for now E{i_elem};
        for j in np.arange(nq):
            f.APmod.M1[:, :, j]= skew(f.psi.M0[:, j]);

        # --- Pretension
        #f.sigma.M0 = np.zeros((6, 1)) # no pretension for now.
        sid.frame.append(f)

    # --- mdCM
    sid.mdCM.M0= np.array([
        np.sum(skew([0.5, 0 , 0  ])*IT['Mxt'].T), 
        np.sum(skew([0,  0.5, 0  ])*IT['Mxt'].T), 
        np.sum(skew([0,   0 , 0.5])*IT['Mxt'].T)
        ]).reshape(3,1)
    for j in np.arange(nq):
        sid.mdCM.M1[:, 0, j]= IT['Mxg'][:, j].T # Ct0 Mgt
    #    sid.mdCM.M1[:, 0, j]= p['Ct'][j, :] # Ct0 Mgt # Verify

    # --- J
    sid.J.M0= IT['Mtt']
    for j in np.arange(nq):
        sid.J.M1[:, :, j]= IT['Mtt_M1'][:,:,j]
    #   sid.J.M1[:, :, j]= -C4[:, :, j] - C4[:, :, j].T;


    # --- Cr  -  Mtg  = \int [~s] Phi dm Or: Mrg, Cr^T
    sid.Cr.M0 = IT['Mtg'].T
    # TODO TODO TODO M1 term
    sid.Cr.M1[:,0,:] = IT['Mtg_M1'][0,:,:] # Kr[0][:,:]; # nq x nq # TODO TODO TODO WEIRD NO NEED FOR TRANSPOSE
    sid.Cr.M1[:,1,:] = IT['Mtg_M1'][1,:,:] # Kr[1][:,:]; 
    sid.Cr.M1[:,2,:] = IT['Mtg_M1'][2,:,:] # Kr[2][:,:]; 

    # --- Ct  -  Mxg = \int Phi dm     Or:  Psi , Ct^T
    # NOTE: Ct.M1 = K0t is geometrical stiffening due to translational acceleration/gravity
    sid.Ct.M0= IT['Mxg'].T
    #sid.Ct.M0= p['C1'].T
    sid.Ct.M1[:, 0, :] = pKg['K0t'][0,:,:] # Gkg['t_ax'], K0t is (3,nf,nf)
    sid.Ct.M1[:, 1, :] = pKg['K0t'][1,:,:]
    sid.Ct.M1[:, 2, :] = pKg['K0t'][2,:,:]

    # --- Gr 
    for j in np.arange(nq):
        sid.Gr.M0[0:3, 3*j:3*j+3]= IT['Gr'][j][:,:];  # columns concatenation 3 per shapes
    for j in np.arange(nq):
        for k in np.arange(nq):
            sid.Gr.M1[0:3, 3*j:3*j+3, k]= IT['Gr_M1'][j][:,:,k];  # columns concatenation 3 per shapes

    # --- Ge = 2 \int Phi^t phi_j~^t dm  = [2C5']
    # Ge Taylor(0,nq,3*nq,0 ,0) # Gyroscopic matrix for modal coordinates
    # IT['Ge'] = np.zeros((nf,nf,3))
    for j in np.arange(nq):
        sid.Ge.M0[0:nq, 3*j:3*j+3]= IT['Ge'][j,:,:] 
        # TODO revisit this
    
    # --- Oe 
    # see [2] (6.407) 
    # M0: nq x 6
    # M1: nq x 6 x nq
    sid.Oe.M0= IT['Oe6']  # nq x 6

    sid.Oe.M1_base = np.zeros(sid.Oe.M1.shape)
    sid.Oe.M1_geom = np.zeros(sid.Oe.M1.shape)

    sid.Oe.M1_base= IT['Oe6_M1']

    sid.Oe.M1_geom[:, 0, :]= pKg['K0omega'][0, 0] 
    sid.Oe.M1_geom[:, 1, :]= pKg['K0omega'][1, 1] 
    sid.Oe.M1_geom[:, 2, :]= pKg['K0omega'][2, 2] 
    sid.Oe.M1_geom[:, 3, :]= pKg['K0omega'][0, 1] + pKg['K0omega'][1, 0] 
    sid.Oe.M1_geom[:, 4, :]= pKg['K0omega'][1, 2] + pKg['K0omega'][2, 1] 
    sid.Oe.M1_geom[:, 5, :]= pKg['K0omega'][0, 2] + pKg['K0omega'][2, 0] 

    sid.Oe.M1 = sid.Oe.M1_base + sid.Oe.M1_geom 
    #sid.Oe.M1 = sid.Oe.M1_base 
    #sid.Oe.M1 = sid.Oe.M1_geom 

    ## --- Me, Ke, De
    sid.Me.M0= MM[6:,6:]
    sid.Ke.M0= KK[6:,6:]
    #sid.De.M0= p['De']

    # --- remove some off-diagonal terms to conform to FAST approach
    if consistency=='OpenFAST':
        sid.Me.M0 = np.diag(np.diag(sid.Me.M0))
        for a in [0,1,2]:
            sid.Ct.M1[:, a, :]= np.diag(np.diag(np.squeeze(sid.Ct.M1[:, a, :])))
        for m in np.arange(6):
            if m<3:
                sid.Oe.M1[:, m, :]+= - pKg['K0omega'][m, m] + np.diag(np.diag(pKg['K0omega'][m, m]))
            else:
                c= m-3;
                d= np.mod(c+1, 3)
                sid.Oe.M1[:, m, :]+= - pKg['K0omega'][c, d] - pKg['K0omega'][c, d].T + 2*np.diag(np.diag(pKg['K0omega'][c, d]))

    return sid


# --------------------------------------------------------------------------------}
# --- FAST Blade 2 SID
# --------------------------------------------------------------------------------{
def FASTBlade2SID(ed_file=None, Imodes_bld=[0,1,2], method='ShapeFunctions', startAtRoot=True, int_method='OpenFAST', consistency='OpenFAST', AdjBlMs=None):
    """ 
    Create a SID of the Blade from an OpenFAST ElastoDyn file
    """
    import welib.weio as weio
    from welib.fast.elastodyn import bladeParameters
    from welib.fast.elastodyn import bladeDerivedParameters


    if method=='FEM':

        parentDir=os.path.dirname(ed_file)
        ed = weio.read(ed_file)
        TipRad = ed['TipRad']
        HubRad = ed['HubRad']
        BldLen= TipRad-HubRad;
        BldFile = ed['BldFile(1)'].replace('"','')
        BldFile=os.path.join(parentDir,BldFile)
        bld = weio.FASTInputFile(BldFile).toDataFrame()
        z   = bld['BlFract_[-]']*BldLen + HubRad
        m   = bld['BMassDen_[kg/m]']               # mu
        EIy = bld['FlpStff_[Nm^2]']
        EIz = bld['EdgStff_[Nm^2]']               # TODO actually EIx, but FEM beams along x
        phi = bld['PitchAxis_[-]']/(180)*np.pi 
        # --- Derived parameters
        phi= np.concatenate(([0], np.diff(phi))) #% add phi_abs(1) to pitch angle
        A  = m*0 + 100       # Area
        Kv = m*0 + 100       # Saint Venant torsion
        E  = 214e9     # Young modulus  
        Iy = EIy/E
        Iz = EIz/E
        xNodes = np.zeros((3,len(m)))
        xNodes[2,:]=z
        sid = Beam2SID(xNodes, Imodes_bld, m, Iy, Iz, Kv=Kv, A=A, E=E, phi=phi)

    elif method=='ShapeIntegral':
        # Extract relevant blade data
        p = bladeParameters(ed_file, AdjBlMs=AdjBlMs)
        if not startAtRoot:
            p['s_G0'][2,:] += p['HubRad']
        # TODO downselect modes
        # Compute SID based on twisted shape functions
        sid = BeamShapes2SID(p['s_G0'], p['s_span'], p['m'], p['EI'], p['Ut'][Imodes_bld], p['Vt'][Imodes_bld], p['Kt'][Imodes_bld], int_method=int_method, damping=None, consistency=consistency, shapeIntegrals=True)


    elif method=='ShapeFunctions':
        p = bladeParameters(ed_file, AdjBlMs=AdjBlMs)
        #p = bladeDerivedParameters(p, inertiaAtBladeRoot=startAtRoot)
        if not startAtRoot:
            p['s_G0'][2,:] += p['HubRad']
        # TODO downselect modes
        # Compute SID based on twisted shape functions
        sid = BeamShapes2SID(p['s_G0'], p['s_span'], p['m'], p['EI'], p['Ut'][Imodes_bld], p['Vt'][Imodes_bld], p['Kt'][Imodes_bld], int_method=int_method, damping=None, consistency=consistency, shapeIntegrals=False)

    else:
        raise NotImplementedError(method)

    return sid

# --------------------------------------------------------------------------------}
# --- FAST Tower 2 SID 
# --------------------------------------------------------------------------------{
def FASTTower2SID(ed_file, method='ShapeFunctions', gravity=None, RotMass=None, Imodes_twr=[0,1,2,3], int_method='OpenFAST'):
    """ """
    import welib.weio as weio
    from welib.fast.elastodyn import rotorParameters
    from welib.fast.elastodyn import towerParameters

    if method=='FEM':
        # --- FEM method
        # Read data from ElastoDyn file
        parentDir=os.path.dirname(ed_file)
        ed = weio.read(ed_file)
        TowerHt = ed['TowerHt']
        TowerBs = ed['TowerBsHt']
        TwrFile = ed['TwrFile'].replace('"','')
        TwrFile=os.path.join(parentDir,TwrFile)
        twr = weio.FASTInputFile(TwrFile).toDataFrame()
        z   = twr['HtFract_[-]']*(TowerHt-TowerBs)
        m   = twr['TMassDen_[kg/m]']               # mu
        EIy = twr['TwFAStif_[Nm^2]']
        EIz = twr['TwSSStif_[Nm^2]']               # TODO actually EIx, but FEM beams along x
        # --- Create Beam FEM model
        # Derived parameters
        A  = m*0 + 100       # Area
        Kv = m*0 + 100       # Saint Venant torsion
        E  = 214e9     # Young modulus  
        Iy = EIy/E
        Iz = EIz/E
        xNodes = np.zeros((3,len(m)))
        xNodes[2,:]=z
        sid = Beam2SID(xNodes, Imodes_twr, m, Iy, Iz, Kv=Kv, A=A, E=E)

    elif method=='ShapeIntegral':
        from welib.yams.flexibility import shapeIntegrals
        # We need the blade parameters for now to get tower top mass...
        if RotMass is None:
            prot, pbld, phub = rotorParameters(ed_file, identicalBlades=False)
            RotMass = prot['RotMass']
        p = towerParameters(ed_file, RotMass=RotMass, gravity=gravity)
        # TODO downselect modes here
        sid = BeamShapes2SID(p['s_G0'], p['s_span'], p['m'], p['EI'], p['U'], p['V'], p['K'], int_method=int_method, damping=None, consistency='OpenFAST', shapeIntegrals=True)

    elif method=='ShapeFunctions':
        # We need the blade parameters for now to get tower top mass...
        if RotMass is None:
            prot, pbld, phub = rotorParameters(ed_file, identicalBlades=False)
            RotMass = prot['RotMass']
        p = towerParameters(ed_file, gravity, RotMass=RotMass)
        # TODO downselect modes here
        sid = BeamShapes2SID(p['s_G0'], p['s_span'], p['m'], p['EI'], p['U'], p['V'], p['K'], int_method=int_method, damping=None, consistency='OpenFAST', shapeIntegrals=False)

    else:
        raise NotImplementedError(method)
    return sid

def FAST2SID(ed_file, Imodes_twr=None, Imodes_bld=None, method='FEM', gravity=None):
    """ 
    method: 'FEM' or 'OpenFAST'
    """

    twr_sid=None
    if Imodes_twr is not None:
        twr_sid = FASTTower2SID(ed_file, Imodes_twr=Imodes_twr, method=method, gravity=gravity)

    bld_sid=None
    if Imodes_bld is not None:
        bld_sid = FASTBlade2SID(ed_file, Imodes_bld=Imodes_bld, method=method)

    return twr_sid, bld_sid

