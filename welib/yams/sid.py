"""

"Standard input data (SID)" used to define a flexible body.

The following references are used:

 [1]  Oskar Wallrapp, 1993, "Standard Input Data of Flexible Members in Multibody Systems"
      Advanced Multibody System Dynamics, 445-450

 [2]  Richard Schwertassek,  Oskar Wallrapp
      "Dynamik Flexibler Mehrkoerpersysteme : Methoden Der Mechanik 




"""

import numpy as np

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
    def __init__(self, nq, name=0):
        self.rframe = '' # name of reference frame
        self.origin = Taylor(1, 3, 1, nq, 0, name='origin') # position of node
        self.phi    = Taylor(0, 3, nq,nq, 0, name='phi'   ) # translational mode shapes
        self.psi    = Taylor(0, 3, nq,nq, 0, name='psi'   ) # rotational mode shapes
        self.APmod  = Taylor(1, 3, 3, nq, 0, name='AP'    ) # relative orientation matrix of node fixed frames
        self.sigma  = Taylor(1, 6, 1, nq, 0, name='sigma' )  # 
        self.name = name

    def __repr__(self):
        s='new node = {}\n'.format(self.name)
        s+='rframe = {}\n'.format(self.rframe)
        s+= str(self.origin)
        s+= str(self.phi   )
        s+= str(self.psi   )
        s+= str(self.APmod )
        s+= str(self.sigma )
        s+='end node\n'
        return s

class SID(object):
    """ """
    def __init__(self, nq):
        self.refmod = Refmod(nq=nq)
        self.frame  = []  # position, orientation, and mode shape matrices at each nodes
        self.mdCM   = Taylor(1,3 ,1   ,nq,0,   name = 'mdCM') # matrix of mass times position of center of mass
        self.J      = Taylor(1,3 ,3   ,nq,0,2, name = 'J' ) # mass moment of inertia
        self.Ct     = Taylor(1,nq,3   ,nq,0,   name = 'Ct') # coupling matrix with translational coordinates
        self.Cr     = Taylor(1,nq,3   ,nq,0,   name = 'Cr') # coupling matrix with rotational coordinates
        self.Me     = Taylor(0,nq,nq  ,0 ,0,2, name = 'Me') # Modal mass matrix
        self.Gr     = Taylor(0,3 ,3*nq,nq,0,   name = 'Gr') # Gyroscopic matrix for rotational coordinates
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
# --- Populating SID
def FEMBeam2SID(Mtt, J0, Mtr, Mgt, Mgr, Mgg, KK, Imodes, frequency, xNodes, DCM, Se, Kr, Kom0, Kom, C4=None, GKg=dict()):
    from welib.yams.utils import skew 
    # TODO take MM, KK, Tr as inputs
    assert(xNodes.shape[0]==3)
    nqk      = 6               # Number of DOF per nodes
    nNodes   = xNodes.shape[1]
    nDOF_tot = nNodes*nqk      # Total number of DOF without constraint (BC)

    iMaxDim = np.argmax(np.max(np.abs(xNodes),axis=1)-np.min(np.abs(xNodes),axis=1)) 

    # See [2] Table 6.9 S. 346
    nq = len(Imodes)
    sid=SID(nq=nq)
    sid.refmod.mass= Mtt[0,0]
    for i in np.arange(nq):
        imode=Imodes[i]
        sid.refmod.ielastq[i]= 'Eigen Mode {:4d}: {:.3f}Hz'.format(imode,frequency[imode])

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
    # NOTE: Mtr (mc0) is anti symmetric and contains m*x_COG
    # The math below is to extract m*x_COG rom Mtr
    sid.mdCM.M0= np.array([
        np.sum(skew([0.5, 0 , 0  ])*Mtr), 
        np.sum(skew([0,  0.5, 0  ])*Mtr), 
        np.sum(skew([0,   0 , 0.5])*Mtr)
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
    # Size nq x 6 x nq
    sid.Oe.M0= Kom0             
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
