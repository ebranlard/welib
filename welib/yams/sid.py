"""

"Standard input data (SID)" used to define a flexible body.

The format is defined in the following paper:
    Oskar Wallrapp, 1993, "Standard Input Data of Flexible Members in Multibody Systems"
    Advanced Multibody System Dynamics, 445-450



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
        s+=str(self.refmod)+'\n'
        s+='new frame\n'
        for i,f in enumerate(self.frame):
            f.name=i+1
            s+=str(f)+'\n'
        s+='end frame\n'
        s+=str(self.mdCM  )+'\n'
        s+=str(self.J     )+'\n'
        s+=str(self.Ct    )+'\n'
        s+=str(self.Cr    )+'\n'
        s+=str(self.Me    )+'\n'
        s+=str(self.Gr    )+'\n'
        s+=str(self.Ge    )+'\n'
        s+=str(self.Oe    )+'\n'
        s+=str(self.ksigma)+'\n'
        s+=str(self.Ke    )+'\n'
        s+=str(self.De    )+'\n'
        return s



# --------------------------------------------------------------------------------}
# ---  Converters
# --------------------------------------------------------------------------------{
# --- Populating SID
def FEMBeam2SID(Mtt, J0, Mgt, Mgr, Mgg, KK, Imodes, frequency, xNodes, DCM, Se, C4=None, Kom=None):
    from welib.yams.utils import skew 
    # TODO take MM, KK, Tr as inputs
    assert(xNodes.shape[0]==3)
    nqk      = 6               # Number of DOF per nodes
    nNodes   = xNodes.shape[1]
    nDOF_tot = nNodes*nqk      # Total number of DOF without constraint (BC)

    # Tabelle 6.9 S. 346
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
        # TODO, stiffnening due to force
        #         if i==nk
        #             f.Phi.order= 1;
        #             f.Phi.M1(1, :, :)= 0*K0Fend_ax; % stiffening due to force (6.46),(6.80), (6.337), (5.122), (6.415)
        #             f.Phi.M1(2, :, :)= 0*K0Fend_ax;
        #             f.Phi.M1(3, :, :)= 0*K0Fend_ax;
        #             f.Phi.M1(extentIdx, :, :)= K0Fend_ax;        
        #         end

        # Node rotations (6.81), (6.337)
        f.psi.M0= Psi.dot(Se);
        if i==nNodes-1:
            f.APmod.M0= DCM[:,:,i-1] # orientation
        else:
            f.APmod.M0= DCM[:,:,i] # orientation
        for j in np.arange(nq):
            f.APmod.M1[:, :, j]= skew(f.psi.M0[:, j]);

        # --- TODO Pretension
        #  f.sigma.M0= zeros(6, 1); % no pretension for now. Spannungen  % (6.413)
        #  if i<nk:
        #      e= i;
        #  else:
        #      e= i-1;
        #  Ke= T{e}' .dot(K{e}).dot(T{e})
        #  if i<nk:
        #      Ke= -Ke;
        #  for j= 1:nq:
        #      f.sigma.M1(:, 1, j)= Ke((i-1)*nqk+(1:6), :)*Se(:, j); % left sided nodal force except last node
        sid.frame.append(f)
    # sid.mdCM
    # TODO
#     sid.md.M0= [sum(sum(crossmat([0.5 0 0]).*mc0)) sum(sum(crossmat([0 0.5 0]).*mc0)) sum(sum(crossmat([0 0 0.5]).*mc0))]';
#     for j= 1:nq
#         sid.md.M1(:, 1, j)= Mgt[j, :]'; # Ct0 Mgt
#     end

    # --- J
    sid.J.M0= J0
    # TODO
    #for j in np.arange(nq):
    #    sid.J.M1[:, :, i]= -C4[:, :, i] - C4[:, :, i].T;

    # --- Ct
    sid.Ct.M0= Mgt
    # TODO
    #     sid.Ct.M1(:, 1, :)= 0*K0t_ax;
    #     sid.Ct.M1(:, 2, :)= 0*K0t_ax;
    #     sid.Ct.M1(:, 3, :)= 0*K0t_ax;
    #     sid.Ct.M1(:, extentIdx, :)= K0t_ax;
    # --- Cr
    sid.Cr.M0= Mgr
    # TODO
#     sid.Cr.M1(:, 1, :)= Kr{1}; % K0r= 0 because only axial stiffening.  Kr= Se*KFr?
#     sid.Cr.M1(:, 2, :)= Kr{2}; % K0r= 0 because only axial stiffening.  Kr= Se*KFr?
#     sid.Cr.M1(:, 3, :)= Kr{3}; % K0r= 0 because only axial stiffening.  Kr= Se*KFr?
    # --- Me
    sid.Me.M0= Mgg;

    # --- Gr TODO
    #sid.Gr.M0= -2*reshape(C4, 3, 3*nq); % (6.403) S. 339; or -2*KFr?
    # --- Ge TODO
    #     for i= 1:nq
    #         sid.Ge.M0(:, 3*(i-1)+(1:3))= 2*[Kr{1}(:, i) Kr{2}(:, i) Kr{3}(:, i)]; % (6.405) S. 340 = 2*C5'
    #     end
    # --- Oe (6.407) TODO
    #sid.Oe.M0= Kom0;
#     sid.Oe.M1(:, 1, :)= Kom{1}+K0omxx;
#     sid.Oe.M1(:, 2, :)= Kom{2}+K0omyy;
#     sid.Oe.M1(:, 3, :)= Kom{3}+K0omzz;
#     sid.Oe.M1(:, 4, :)= Kom{4}+K0omxy;
#     sid.Oe.M1(:, 5, :)= Kom{5}+K0omxz;
#     sid.Oe.M1(:, 6, :)= Kom{6}+K0omyz;
# 
    # ---  sid.Ke 
    sid.Ke.M0= (Se.T).dot(KK).dot(Se);
    #sid.ksigma.M0= [];
    # --- De
    #sid.De.M0= [];

    return sid
