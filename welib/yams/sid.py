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
def Beam2SID(xNodes, Imodes, m, Iy, Iz, G=None, Kv=None, A=None, E=None, phi=None):
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
    from welib.yams.utils import skew 
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

# --------------------------------------------------------------------------------}
# ---FAST 2 SID
# --------------------------------------------------------------------------------{
def FAST2SID(ed_file, Imodes_twr=None, Imodes_bld=None):
    import welib.weio as weio
    import os
    # --- Read data from ElastoDyn file
    parentDir=os.path.dirname(ed_file)
    ed = weio.read(ed_file)
#     edfile =
#     if twr:

    twr_sid=None
    if Imodes_twr is not None:
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

        twr_sid= Beam2SID(xNodes, Imodes_twr, m, Iy, Iz, Kv=Kv, A=A, E=E)

    bld_sid=None
    if Imodes_bld is not None:
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
        bld_sid= Beam2SID(xNodes, Imodes_bld, m, Iy, Iz, Kv=Kv, A=A, E=E, phi=phi)

    return twr_sid, bld_sid

