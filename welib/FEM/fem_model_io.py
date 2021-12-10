import numpy as np


def yaml_array(var, M, Fmt='{:15.6e}', comment=''):
    M = np.atleast_2d(M)
    if len(comment)>0:
        s='{}: # {} x {} {}\n'.format(var, M.shape[0], M.shape[1], comment)
    else:
        s='{}: # {} x {}\n'.format(var, M.shape[0], M.shape[1])
    if M.shape[0]==1:
        for l in M:
            s+= '  - [' + ','.join([Fmt.format(le) for le in l]) + ',]\n'
    else:
        for l in M:
            s+= '  - [' + ','.join([Fmt.format(le) for le in l]) + ']\n'
    s = s.replace('e+','E+').replace('e-','E-')
    return s

def toYAML(model, filename):
    with open(filename, 'w') as f:
        s=''
        s += '#____________________________________________________________________________________________________\n'
        s += '# RIGID BODY EQUIVALENT DATA\n'
        s += '#____________________________________________________________________________________________________\n'
        s0 = 'Mass: {:15.6e} # Total Mass\n'.format(model.M_O[0,0])
        s += s0.replace('e+','E+').replace('e-','E-')
        s0 = 'CM_point: [{:15.6e},{:15.6e},{:15.6e},] # Center of mass coordinates (Xcm,Ycm,Zcm)\n'.format(model.center_of_mass[0],model.center_of_mass[1],model.center_of_mass[2])
        s += s0.replace('e+','E+').replace('e-','E-')
        s0 = 'TP_point: [{:15.6e},{:15.6e},{:15.6e},] # Transition piece reference point\n'.format(model.TP[0],model.TP[1],model.TP[2])
        s += s0.replace('e+','E+').replace('e-','E-')
        s += yaml_array('MRB',  model.M_O,  comment = 'Rigid Body Equivalent Mass Matrix w.r.t. (0,0,0).')
        s += yaml_array('M_P' , model.M_TP, comment = 'Rigid Body Equivalent Mass Matrix w.r.t. TP Ref point')
        s += yaml_array('M_G' , model.M_G,  comment = 'Rigid Body Equivalent Mass Matrix w.r.t. CM (Xcm,Ycm,Zcm).')


        s += '#____________________________________________________________________________________________________\n'
        s += '# GUYAN MATRICES at the TP reference point\n'
        s += '#____________________________________________________________________________________________________\n'
        MBB = model.MM_CB[np.ix_(model.DOF_Leader_CB, model.DOF_Leader_CB)]
        KBB = model.KK_CB[np.ix_(model.DOF_Leader_CB, model.DOF_Leader_CB)]
        #CBB = model.KK_CB[np.ix_(model.DOF_Leader_CB, model.DOF_Leader_CB)]
        TI=model.TI
        MBBr = TI.T.dot(MBB).dot(TI)
        KBBr = TI.T.dot(KBB).dot(TI)
        MBBr[np.abs(MBBr)<1e-4] =0
        KBBr[np.abs(KBBr)<1e-4] =0
        s += yaml_array('KBBt' , KBBr,  comment = '')
        s += yaml_array('MBBt' , MBBr,  comment = '')
        s += 'CBBt: # 6 x 6 (user Guyan Damping + potential joint damping from CB-reduction)\n'

        s += '#____________________________________________________________________________________________________\n'
        s += '# SYSTEM FREQUENCIES\n'
        s += '#____________________________________________________________________________________________________\n'
        s += '#Eigenfrequencies [Hz] for full system, with reaction constraints (+ Soil K/M + SoilDyn K0) \n'
        s += yaml_array('Full_frequencies', model.freq_cr)
        s += '#Frequencies of Guyan modes [Hz]\n'
        s += yaml_array('GY_frequencies', model.f_G)
        s += '#Frequencies of Craig-Bampton modes [Hz]\n'
        s += yaml_array('CB_frequencies', model.f_CB)
        s += '#____________________________________________________________________________________________________\n'
        s += '# Internal FEM representation\n'
        s += '#____________________________________________________________________________________________________\n'
        s += 'nNodes_I: {:7d} # Number of Nodes: "interface" (I)\n'.format(len(model.interfaceNodes))
        s += 'nNodes_C: {:7d} # Number of Nodes: "reactions" (C)\n'.format(len(model.reactionNodes))
        s += 'nNodes_L: {:7d} # Number of Nodes: "internal"  (L)\n'.format(len(model.internalNodes))
        s += 'nNodes  : {:7d} # Number of Nodes: total   (I+C+L)\n'.format(len(model.Nodes))
        s += 'nDOF__B : {:7d} # Number of DOFs:             retained (__B)\n'.format(model.SD_IO_Vars['nDOF__B'])
        s += 'nDOF__L : {:7d} # Number of DOFs:             internal (__L)\n'.format(model.SD_IO_Vars['nDOF__L'])
        s += 'nDOF__F : {:7d} # Number of DOFs:             fixed    (__F)\n'.format(model.SD_IO_Vars['nDOF__F'])
        s += 'nDOF_red: {:7d} # Number of DOFs: total\n'                     .format(model.SD_IO_Vars['nDOF___'])
        s += yaml_array('Nodes_I', model.nodeID_py2SD([n.ID for n in model.interfaceNodes]), Fmt='{:7d}', comment='"interface" nodes"');
        s += yaml_array('Nodes_C', model.nodeID_py2SD([n.ID for n in model.reactionNodes ]), Fmt='{:7d}', comment='"reaction" nodes"');
        s += yaml_array('Nodes_L', model.nodeID_py2SD([n.ID for n in model.internalNodes ]), Fmt='{:7d}', comment='"internal" nodes"');
        s += yaml_array('DOF___B', np.array(model.DOF_Leader  )+1, Fmt='{:7d}',  comment='all         retained  DOFs');
        s += yaml_array('DOF___F', np.array(model.DOF_Fixed   )+1, Fmt='{:7d}',  comment='all         fixed     DOFs');
        s += yaml_array('DOF___L', np.array(model.DOF_Follower)+1, Fmt='{:7d}',  comment='all         internal  DOFs');
        s += '\n'
        s += '#Index map from DOF to nodes\n'
        s += '#     Node No.,  DOF/Node,   NodalDOF\n'
        s += 'DOF2Nodes: # {} x 3 (nDOFRed x 3, for each constrained DOF, col1: node index, col2: number of DOF, col3: DOF starting from 1)\n'.format(model.nDOF)
        DOF2Nodes = model.DOF2Nodes
        for l in model.DOF2Nodes:
            s +='  - [{:7d},{:7d},{:7d}] # {}\n'.format(l[1]+1, l[2], l[3], l[0]+1 )
        s += '#     Node_[#]          X_[m]           Y_[m]           Z_[m]       JType_[-]       JDirX_[-]       JDirY_[-]       JDirZ_[-]  JStff_[Nm/rad]\n'
        s += 'Nodes: # {} x 9\n'.format(len(model.Nodes))
        for n in model.Nodes:
            s += '  - [{:7d}.,{:15.3f},{:15.3f},{:15.3f},{:14d}.,   0.000000E+00,   0.000000E+00,   0.000000E+00,   0.000000E+00]\n'.format(model.nodeID_py2SD(n.ID), n.x, n.y, n.z, int(n.data['Type']) )
        s += '#    Elem_[#]    Node_1   Node_2   Prop_1   Prop_2     Type     Length_[m]      Area_[m^2]  Dens._[kg/m^3]        E_[N/m2]        G_[N/m2]       shear_[-]       Ixx_[m^4]       Iyy_[m^4]       Jzz_[m^4]          T0_[N]\n'
        s += 'Elements: # {} x 16\n'.format(len(model.Elements))
        for e in model.Elements:
            I = e.inertias
            s0='  - [{:7d}.,{:7d}.,{:7d}.,{:7d}.,{:7d}.,{:7d}.,{:15.3f},{:15.3f},{:15.3f},{:15.6e},{:15.6e},{:15.6e},{:15.6e},{:15.6e},{:15.6e},{:15.6e}]\n'.format(
                model.elemID_py2SD(e.ID), model.nodeID_py2SD(e.nodeIDs[0]), model.nodeID_py2SD(e.nodeIDs[1]), model.propID_py2SD(e.propIDs[0], e.propset), model.propID_py2SD(e.propIDs[1], e.propset), model.elemType_py2SD(e.data['Type']), 
                e.length, e.area, e.rho, e.E, e.G, e.kappa, I[0], I[1], I[2], e.T0)
            s += s0.replace('e+','E+').replace('e-','E-')
        s += '#____________________________________________________________________________________________________\n'
        s += '#User inputs\n'
        s += '\n'
        s += '#Number of properties (NProps):{:6d}\n'.format(len(model.NodePropertySets['Beam']))
        s += '#Prop No         YoungE         ShearG        MatDens          XsecD          XsecT\n'
        for ip,p in enumerate(model.NodePropertySets['Beam']):
            s0='#{:8d}{:15.6e}{:15.6e}{:15.6e}{:15.6e}{:15.6e}\n'.format(p.ID, p['E'],p['G'],p['rho'],p['D'],p['t'])
            s +=  s0.replace('e+','E+').replace('e-','E-')
        s +='\n'
        s += '#No. of Reaction DOFs:{:6d}\n'.format(len(model.SD_IO_Vars['IDC__']) )
        s += '#React. DOF_ID    BC\n'
        s += '\n'.join(['#{:10d}{:10s}'.format(idof+1,'     Fixed' ) for idof in model.SD_IO_Vars['IDC_F']])
        s += '\n'.join(['#{:10d}{:10s}'.format(idof+1,'     Free'  ) for idof in model.SD_IO_Vars['IDC_L']])
        s += '\n'.join(['#{:10d}{:10s}'.format(idof+1,'     Leader') for idof in model.SD_IO_Vars['IDC_B']])
        s += '\n\n'
        s += '#No. of Interface DOFs:{:6d}\n'.format(len(model.SD_IO_Vars['IDI__']))
        s += '#Interf. DOF_ID    BC\n'
        s += '\n'.join(['#{:10d}{:10s}'.format(idof+1,'     Fixed' ) for idof in model.SD_IO_Vars['IDI_F']])
        s += '\n'.join(['#{:10d}{:10s}'.format(idof+1,'     Leader') for idof in model.SD_IO_Vars['IDI_B']])
        s += '\n\n'
        s += '#Number of concentrated masses (NCMass):     0\n'
        s += '#JointCMas           Mass            JXX            JYY            JZZ            JXY            JXZ            JYZ           MCGX           MCGY           MCGZ\n'
        s += '\n'
        s += '#Number of members    18\n'
        s += '#Number of nodes per member:     2\n'
        s += '#Member I Joint1_ID Joint2_ID    Prop_I    Prop_J           Mass         Length     Node IDs...\n'
        s += '#       77        61        60        11        11   1.045888E+04   2.700000E+00       19    18\n'
        s += '#____________________________________________________________________________________________________\n'
        s += '#Direction Cosine Matrices for all Members: GLOBAL-2-LOCAL. No. of 3x3 matrices=    18\n'
        s += '#Member I        DC(1,1)        DC(1,2)        DC(1,3)        DC(2,1)        DC(2,2)        DC(2,3)        DC(3,1)        DC(3,2)        DC(3,3)\n'
        s += '#       77  1.000E+00  0.000E+00  0.000E+00  0.000E+00 -1.000E+00  0.000E+00  0.000E+00  0.000E+00 -1.000E+00\n'
        s += '#____________________________________________________________________________________________________\n'
        s += '#FEM Eigenvectors (114 x 108) [m or rad], full system with reaction constraints (+ Soil K/M + SoilDyn K0)\n'
        s += yaml_array('Full_Modes', model.Q_c)
        s += '#____________________________________________________________________________________________________\n'
        s += '#CB Matrices (PhiM,PhiR) (reaction constraints applied)\n'
        s += yaml_array('PhiM', model.Phi_CB[:,:model.nModesCB] ,comment='(CB modes)')
        s += yaml_array('PhiR', model.Phi_G,  comment='(Guyan modes)')
        s += '\n'
        s += '#____________________________________________________________________________________________________\n'
        s += '# ADDITIONAL DEBUGGING INFORMATION\n'
        s += '#____________________________________________________________________________________________________\n'
        s +=  ''
        e = model.Elements[0]
        rho=e.rho
        A = e.area
        L = e.length
        t= rho*A*L
        s0 = '{:15.6e}{:15.6e}{:15.6e}{:15.6e}{:15.6e}{:15.6e}{:15.6e}{:15.6e}{:15.6e}{:15.6e}{:15.6e}\n'.format(model.gravity,e.area, e.length, e.inertias[0], e.inertias[1], e.inertias[2], e.kappa, e.E, e.G, e.rho, t)
        s0 = s0.replace('e+','E+').replace('e-','E-')
        s += s0
        s += yaml_array('KeLocal' +str(), model.Elements[0].Ke(local=True))

        for ie,e in enumerate(model.Elements):
            s += yaml_array('DC' +str(ie+1), e.DCM.transpose())
            s += yaml_array('Ke' +str(ie+1), e.Ke())
            s += yaml_array('Me' +str(ie+1), e.Me())
            s += yaml_array('FGe'+str(ie+1), e.Fe_g(model.gravity))
            s += yaml_array('FCe'+str(ie+1), e.Fe_o())

            s += yaml_array('KeLocal' +str(ie+1), e.Ke(local=True))
            s += yaml_array('MeLocal' +str(ie+1), e.Me(local=True))
            s += yaml_array('FGeLocal'+str(ie+1), e.Fe_g(model.gravity, local=True))
            s += yaml_array('FCeLocal'+str(ie+1), e.Fe_o(local=True))

        s += '#____________________________________________________________________________________________________\n'
        e = model.Elements[0]
        s += yaml_array('Ke', e.Ke(local=True), comment='First element stiffness matrix'); # TODO not in local
        s += yaml_array('Me', e.Me(local=True), comment='First element mass matrix');
        s += yaml_array('FGe', e.Fe_g(model.gravity,local=True), comment='First element gravity vector');
        s += yaml_array('FCe', e.Fe_o(local=True), comment='First element cable pretension');
        s += '#____________________________________________________________________________________________________\n'
        s += '#FULL FEM K and M matrices. TOTAL FEM TDOFs:    {}\n'.format(model.MM_c.shape[0])
        s += yaml_array('K', model.KK_c, comment='Stiffness matrix');
        s += yaml_array('M', model.MM_c, comment='Mass matrix');
        s += '#____________________________________________________________________________________________________\n'
        s += '#Gravity and cable loads applied at each node of the system (before DOF elimination with T matrix)\n'
        s += yaml_array('FG', model.FF, comment='');
        print('>>> FG sum / gravity', np.sum(model.FF/model.gravity))
        s += '#____________________________________________________________________________________________________\n'
        s += '#Additional CB Matrices (MBB,MBM,KBB) (constraint applied)\n'
        s += yaml_array('MBB', model.MM_CB[np.ix_(model.DOF_Leader_CB, model.DOF_Leader_CB)], comment='');
        s += yaml_array('MBM', model.MM_CB[np.ix_(model.DOF_Leader_CB, model.DOF_Follower_CB[:model.nModesCB])], comment='');
        s += yaml_array('CMMdiag', model.CC_MM_CB, comment='(2 Zeta OmegaM)');
        s += yaml_array('KBB', model.KK_CB[np.ix_(model.DOF_Leader_CB, model.DOF_Leader_CB)], comment='');
        KMM = model.KK_CB[np.ix_(model.DOF_Follower_CB, model.DOF_Follower_CB)]
        s += yaml_array('KMM', np.diag(KMM), comment='(diagonal components, OmegaL^2)');
        s += yaml_array('KMMdiag', np.diag(KMM)[:model.nModesCB], comment='(diagonal components, OmegaL^2)');
        s += yaml_array('PhiL', model.Phi_CB, comment='');
        s += 'PhiLOm2-1: # 18 x 18 \n'
        s += 'KLL^-1: # 18 x 18 \n'
        s += '#____________________________________________________________________________________________________\n'
        s += yaml_array('T_red', model.T_c, Fmt = '{:9.2e}', comment='(Constraint elimination matrix)');
        s += 'AA: # 16 x 16 (State matrix dXdx)\n'
        s += 'BB: # 16 x 48 (State matrix dXdu)\n'
        s += 'CC: # 6 x 16 (State matrix dYdx)\n'
        s += 'DD: # 6 x 48 (State matrix dYdu)\n'
        s += '#____________________________________________________________________________________________________\n'
        s += yaml_array('TI', model.TI,  Fmt = '{:9.2e}',comment='(TP refpoint Transformation Matrix TI)');
        f.write(s);


