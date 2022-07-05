import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import copy
# Local 
from welib.tools.clean_exceptions import *
from welib.weio import FASTInputFile

from welib.fast.fast_mesh import PointMesh
from welib.yams.rotations import Rodriguez_A


class Morison:

    def __init__(self, graph, File, WtrDpth, MSL2SWL):
        """ 
        graph: graph generated from hydrodyn file
        File : File content of hydrodyn input file
        """
        self.graph        = copy.deepcopy(graph)
        self.File         = File

        # Internal
        self.p={}
        self.m={}
        self.p['WtrDpth'] = WtrDpth
        self.p['MSL2SWL'] = MSL2SWL


        # --- Important init 
        # Reindex Nodes from 1-N to match what HydroDyn use
        self.graph.reindexNodes(offset=0)
        # Create "Simulation Nodes"
        self._divisionDone=False
        self.NodesBeforeSwap = self.generateSimulationNodes()

    def __repr__(self):
        s='<{} object>:\n'.format(type(self).__name__)
        s+='|properties:\n'
        s+='|- File: (input file data)\n'
        s+='|parameters (p):\n'
        s+='| - WtrDpth: {}\n'.format(p.WtrDpth)
        s+='| - MSL2SWL: {}\n'.format(p.MSL2SWL)
        s+='|methods:\n'
        s+='|- init\n'
        return s

    # --------------------------------------------------------------------------------}
    # --- Functions 
    # --------------------------------------------------------------------------------{
    def init(self, initData):
        """
        Initialize HydroDyn model 

        Gravity: position of transition point
        """
        f = self.File
        graph = self.graph
        self._setupMembers()
        #------------------------ set up joint (or joint-node) properties --
        #       ! Redundant work (these are already assigned to the member data arrays, 
        #       ! but is needed on the joint data because we report the tMG, and MGDensity at each Joint node in the Summary File
        #       call SetNodeMG( InitInp%NMGDepths, InitInp%MGDepths, InitInp%Nodes(i), InitInp%MSL2SWL, InitInp%Nodes(i)%tMG, InitInp%Nodes(i)%MGDensity )
        NJoints = len(graph.Nodes)
        NNodes  = len(self.NodesBeforeSwap)
        self.m['nodeInWater']   = np.zeros(NNodes)
        self.m['vrel']          = np.zeros((    3, NNodes   ))
        self.m['FV']            = np.zeros((    3, NNodes   ))
        self.m['FA']            = np.zeros((    3, NNodes   ))
        self.m['FDynP']         = np.zeros(       NNodes     )
        self.m['F_I_End']       = np.zeros((    3, NJoints))
        self.m['F_BF_End']      = np.zeros((    6, NJoints))
        self.m['F_A_End']       = np.zeros((    3, NJoints))
        self.m['F_D_End']       = np.zeros((    3, NJoints))
        self.m['F_B_End']       = np.zeros((    6, NJoints))
        self.m['F_IMG_End']     = np.zeros((    6, NJoints))
        self.p['An_End']        = np.zeros((    3, NJoints))
        self.p['DragConst_End'] = np.zeros(        NJoints  )
        self.p['I_MG_End']      = np.zeros(( 3, 3, NJoints))
        self.p['F_WMG_End']     = np.zeros((    3, NJoints))
        self.p['Mass_MG_End']   = np.zeros(        NJoints )
        self.p['AM_End']        = np.zeros(( 3, 3, NJoints))
        self.p['DP_Const_End']  = np.zeros((    3, NJoints))
        self.p['nodeInWater']   = initData['nodeInWater']  # At each time step
        self.p['WaveVel']       = initData['WaveVel']
        self.p['WaveAcc']       = initData['WaveAcc']
        self.p['WaveDynP']      = initData['WaveDynP']
        self.p['WaveTime']      = initData['WaveTime']
        self.p['Gravity']       = initData['Gravity']
        self.p['WtrDens']       = initData['WtrDens']

        # --- Connectivity
        Connectivity=[]
        for e in graph.Elements:
            n1 = e.nodes[0]
            n2 = e.nodes[1]
            i1 = graph.Nodes.index(n1)
            i2 = graph.Nodes.index(n2)
            mem      = e.MorisonData
            idx      = mem['nodeIDs']
            if len(idx)>2:
                for ii, _ in enumerate(idx[:-1]):
                    i1=idx[ii]
                    i2=idx[ii+1]
                    Connectivity.append([i1,i2])
            else:
                Connectivity.append([i1,i2])




   
        # ---- Input Mesh
        u = dict()
        u['Mesh']= PointMesh(NNodes, RefPoint=np.array([0,0,0]), Connectivity=Connectivity) # TODO TODO RefPoint will be used for rigid body rotations, might need TPRef here
        for i,pos in enumerate(self.NodesBeforeSwap):
            # Here positions are relative to MSL not SWL
            pos[2] = pos[2] + self.p['MSL2SWL']
            u['Mesh'].Position[i,:] = pos
        # --- Output Mesh is duplicate
        y = dict()
        y['Mesh'] = copy.deepcopy(u['Mesh'])
        # Define initial system states here:
        #m%LastIndWave              = 1
        
        p = self.p
        # --- loop through joints to calculate joint quantities (the joints are the first NJoints nodes)
        for i,n in enumerate(graph.Nodes):
            An        = 0.0
            Vn        = 0.0
            I_n       = 0.0
            MGdens    = 0.0
            tMG       = -999.0
            An_drag   = 0.0
            #print('--------------------------- JOINT',i+1)
            if n.point[2] >= -p['WtrDpth']:
                # loop through each member attached to the joint, getting the radius of its appropriate end

                elem = graph.node2Elements(n)
                # --- Loop on elements connected to joint
                for e in elem:
                    member = e.MorisonData
                    if n.ID == e.nodes[0].ID:
                        MemberEndIndx = member['NElements']
                    else:
                        MemberEndIndx = 0
                    #print('Member',e.ID, 'Indx',MemberEndIndx+1)
                    # Compute the signed area*outward facing normal of this member
                    sgn = 1.0
                    if MemberEndIndx == 0:
                        sgn = -1.0                                # Local coord sys points into member at starting node, so flip sign of local z vector
                    else:
                        sgn = 1.0                                 # Local coord sys points out of member at ending node, so leave sign of local z vector
                    # Account for reordering of what the original node for the end was -- This affects the sign of the An term which can pose a problem for members crossing the waterline
                    if e.flipped:
                        sgn = -1.0 * sgn
                    # Compute the signed quantities for this member end (for drag regardless of PropPot value), and add them to the joint values
                    An_drag = An_drag + sgn* member['k']*np.pi*(member['RMG'][MemberEndIndx])**2     # area-weighted normal vector
                    # For the following quantities, the attached member cannot be modeled using WAMIT if we're to count it
                    if not e.data['Pot']:
                        # Compute the signed quantities for this member end, and add them to the joint values
                        An = An + sgn* member['k']*np.pi*(member['RMG'][MemberEndIndx])**2     # area-weighted normal vector
                        Vn = Vn + sgn* member['k']*      (member['RMG'][MemberEndIndx])**3     # r^3-weighted normal vector used for mass
                        I_n=I_n + sgn* member['k']*np.pi*(member['RMG'][MemberEndIndx])**4     # r^4-weighted normal vector used for moments of inertia
                        if tMG == -999.0:
                            # All member nodes at this joint will have the same MG thickness and density, so only do this once
                            tMG =    member['tMG'][MemberEndIndx      ]
                            MGdens = member['MGdensity'][MemberEndIndx]
                #Vn *= 2*np.pi/3   # Semisphere volume is Vn = 2/3 pi \sum (r_MG^3 k)   # TODO TODO TODO uncomment when OpenFAST fixed
                p['An_End'][:,i] = An_drag 
                Amag_drag =  An_drag.dot(An_drag)
                Amag      =  An.dot(An)
                if Amag_drag == 0:
                    p['DragConst_End'][i] =  0
                else:
                    p['DragConst_End'][i] = n.data['JAxCd']*p['WtrDens'] / ( 4.* Amag_drag )
                # magnitudes of normal-weighted values
                Amag = np.sqrt(Amag)
                Vmag = np.linalg.norm(Vn)
                Imag = np.linalg.norm(I_n)
                # Constant part of the external hydrodynamic added mass term
                if Vmag > 0.0:
                    v2D = Vn.reshape((3,1))
                    p['AM_End'][:,:,i] = (n.data['JAxCa']*p['WtrDens']/Vmag)*  (v2D).dot(v2D.T) 
                # Constant part of the external hydrodynamic dynamic pressure force
                if Amag > 0.0:
                    p['DP_Const_End'][:,i] = -n.data['JAxCp']*An 
                # marine growth mass/inertia magnitudes
                p['Mass_MG_End'][i] = MGdens * tMG * Amag
                p['F_WMG_End'][2,i] =        -MGdens * tMG * Amag * p['Gravity']  # Z component of the directional force due to marine growth mass at joint
                Ir_MG_end   =  0.25 * MGdens * tMG * Imag  # radial moment of inertia magnitude
                Il_MG_end   =  0.5  * MGdens * tMG * Imag  # axial moment of inertia magnitude
                ## get rotation matrix for moment of inertia orientations
                R_I = Rodriguez_A(I_n)
                ## globally-oreinted moment of inertia matrix for joint
                Irl_mat = np.zeros((3,3))
                Irl_mat[0,0] = Ir_MG_end
                Irl_mat[1,1] = Ir_MG_end
                Irl_mat[2,2] = Il_MG_end
                #p['I_MG_End'][:,:,i] = MatMul( MatMul(R_I, Irl_mat), Transpose(R_I) ) # final moment of inertia matrix for node
                #print('An    ',An)
                #print('Vn    ',Vn)
                #print('I_n   ',I_n)
                #print('An_Drag  ', An_drag)
                #print('An_Drag_0', sgn* member['k']*np.pi*(member['RMG'][MemberEndIndx])**2  )
                #print('Amag_d', Amag_drag)
                #print('Amag ', Amag)
                #print('AM_end       ',p['AM_End']      [0,:,i])
                #print('AM_end       ',p['AM_End']      [1,:,i])
                #print('AM_end       ',p['AM_End']      [2,:,i])
                #print('DP_Const_End ',p['DP_Const_End'][:,i])
                #print('MassMGEnd    ',p['Mass_MG_End'] [i])
                #print('F_WMG_End    ',p['F_WMG_End']   [:,i])
                #print('I_MG_End     ',p['I_MG_End']    [0,:,i])
                #print('I_MG_End     ',p['I_MG_End']    [1,:,i])
                #print('I_MG_End     ',p['I_MG_End']    [2,:,i])
            # ----END IF  # InitInp%InpJoints(i)%Position(3) >= -p%WtrDpth
        # ----END Loop on Joints to compute End Parameters

        # ---  WriteOutputs TODO
        # IF ( p%OutSwtch > 0) then  #@mhall: moved this "if" to after allocations
        #    CALL MrsnOUT_Init( InitInp, y, p, InitOut, errStat, errMsg )
        #       # Determine if we need to perform output file handling
        #    IF ( p%OutSwtch == 1 .OR. p%OutSwtch == 3 ) THEN  
        #       CALL MrsnOUT_OpenOutput( Morison_ProgDesc%Name, TRIM(InitInp%OutRootName)//'.HD', p, InitOut, errStat, errMsg )
        #    END IF
        # END IF  

        # We call CalcOutput to compute the loads for the initial reference position
        # Then we can use the computed load components in the Summary File
        # NOTE: Morison module has no states, otherwise we could not do this.
        self.calcOutput(t=0, u=u, y=y)

        #       ! Write Summary information now that everything has been initialized. 
        #    CALL WriteSummaryFile( InitInp%UnSum, InitInp%Gravity, InitInp%MSL2SWL, InitInp%WtrDpth, InitInp%NJoints, InitInp%NNodes, InitInp%Nodes, p%NMembers, p%Members, &
        #                           p%NumOuts, p%OutParam, p%NMOutputs, p%MOutLst,  p%NJOutputs, p%JOutLst, u%Mesh, y%Mesh, &
        #                           p, m, errStat, errMsg )
        # For convenience
        self._y = y
        self._u = u
        return u, y

    # --------------------------------------------------------------------------------}
    # ---  
    # --------------------------------------------------------------------------------{
    def _MorisonPositions(self):
        """ Return nodal positions used for Morison, WaveKin and Current.
        NOTE: MSL2SWL is removed from "z" """
        if not self._divisionDone:
            raise Exception('Call generateSimulationNodes first')
        graph   = self.graph
        MSL2SWL = self.p['MSL2SWL']
        MorisonPos=[]
        for n in graph.Nodes:
            MorisonPos.append([n.x,n.y,n.z]) 
        for e in graph.Elements:
            if hasattr(e, 'MorisonData'):
                for p in e.MorisonData['SubNodesPositions']:
                    MorisonPos.append([p[0], p[1], p[2]])
        MorisonPos=np.asarray(MorisonPos)
        return MorisonPos

    def generateSimulationNodes(self):
        """ Morison_GenerateSimulationNodes in Morison.f90 """
        if not self._divisionDone:
            graph = self.graph
            # NOTE: we change the reference for the z-coordinates to SWL for all nodes
            nodeCount=0
            for n in graph.Nodes:
                n.z -= self.p['MSL2SWL']
                nodeCount+=1
            for e in graph.Elements:
                if e.data['Pot'] is False:
                    numDiv = np.ceil(e.length/e.data['DivSize']).astype(int)
                    MorisonData = {}
                    MorisonData['NElements'] = numDiv
                    MorisonData['dl']        = e.length/numDiv
                    MorisonData['refLength'] = e.length
                    n1, n2 = e.nodes
                    SubNodesPositions = np.zeros((numDiv-1,3))
                    MorisonData['nodeIDs']     = np.zeros(numDiv+1).astype(int)
                    MorisonData['nodeIDs'][0]  = e.nodeIDs[0]
                    MorisonData['nodeIDs'][-1] = e.nodeIDs[1]
                    for j in range(numDiv-1):
                        s = (j+1)/numDiv
                        SubNodesPositions[j,:] = n1.point * (1-s) + n2.point * s
                        MorisonData['nodeIDs'][j+1] = nodeCount
                        nodeCount+=1
                    MorisonData['SubNodesPositions'] = SubNodesPositions
                    e.MorisonData=MorisonData
        self._divisionDone=True
        Nodes=self._MorisonPositions()
        return Nodes

    # --------------------------------------------------------------------------------}
    # --- Setup of Members 
    # --------------------------------------------------------------------------------{
    def _allocateMemberDataArrays(self, member):
        NElements = member['NElements']
        for k in ['dRdl_mg','dRdl_in','floodstatus','alpha','alpha_fb','alpha_fb_star',
                'm_fb_l','m_fb_u','h_cfb_l','h_cfb_u','I_lfb_l','I_lfb_u','I_rfb_l','I_rfb_u',
                'm_mg_l','m_mg_u','h_cmg_l','h_cmg_u','I_lmg_l','I_lmg_u','I_rmg_l','I_rmg_u',
                'Cfl_fb', 'Cfr_fb','CM0_fb']:
            member[k] = np.zeros(NElements)
        for k in ['R','RMG','Rin','tMG','MGdensity''Cd','Ca','Cp','AxCd','AxCa','AxCp']:
            member[k] = np.zeros(NElements+1)
        member['Loads']=dict()
        for k in ['F_D','F_A','F_B','F_BF','F_I','F_If','F_WMG','F_IMG']:
            member['Loads'][k] = np.zeros((6,NElements+1))


    def _setMemberProperties(self, member, m):
        """  see SetMemberProperties in Morison.f90 """
        WtrDepth = self.p['WtrDpth']                 # TODO which one is it
        N        = member['NElements']
        dl       = member['dl']
        vec      = m.nodes[1].point-m.nodes[0].point
        # calculate reference orientation information.  Note: members are straight to start
        memLength = member['refLength']
        member['k'] = (vec/memLength)  # vector along member from start to end point
        tk = np.array([[member['k'][0], member['k'][1], member['k'][2]]])
        member['kkt']    = tk.T.dot(tk)
        member['Ak']     =  np.eye(3) - member['kkt']
        phi = np.arccos(vec[2]/memLength)  # incline angle   
        sinPhi = np.sin(phi)
        cosPhi = np.cos(phi)  
        member['cosPhi_ref'] = cosPhi
        # --- MG TODO
        # These are all per node and not done here, yet
        #  do i = 1, member%NElements+1
        #     call SetNodeMG( InitInp%NMGDepths, InitInp%MGDepths, InitInp%Nodes(member%NodeIndx(i)), InitInp%MSL2SWL, member%tMG(i), member%MGDensity(i) )
        member['MGdensity']=np.zeros(N+1) # TODO
        member['tMG'] = np.zeros(N+1)   # TODO
        prop1 = m.nodeProps[0]  # NOTE: t&D are not stored in nodes since they are member dependent
        prop2 = m.nodeProps[1] 
        t             = np.linspace(prop1.data['t']  , prop2.data['t']  , N+1)
        member['R']   = np.linspace(prop1.data['D']/2, prop2.data['D']/2, N+1)
        member['RMG'] = member['R']+member['tMG']
        member['Rin'] = member['R']-t

        # --- see SetExternalHydroCoefs from Morison.f90
        C = self.graph.NodePropertySets['SimpleCoefs'][0] # we only have one of these sets
        if m.data['CoefMod'] == 'SimpleCoefs':
            # Simple model : all nodes receive the same coefficients
            member['Cd']  = np.array([C['CdMG']   if t>0 else C['Cd']   for t in member['tMG']])
            member['Ca']  = np.array([C['CaMG']   if t>0 else C['Ca']   for t in member['tMG']])
            member['Cp']  = np.array([C['CpMG']   if t>0 else C['Cp']   for t in member['tMG']])
            member['AxCd']= np.array([C['AxCdMG'] if t>0 else C['AxCd'] for t in member['tMG']])
            member['AxCa']= np.array([C['AxCaMG'] if t>0 else C['AxCa'] for t in member['tMG']])
            member['AxCp']= np.array([C['AxCpMG'] if t>0 else C['AxCp'] for t in member['tMG']])
        elif m.data['CoefMod'] == 'DepthCoefs':
            # Depth-based model: coefficients are set using depth-based table data
            raise NotImplementedError()
            #do i = 1, member%NElements + 1
            #CALL SetDepthBasedCoefs( nodes(member%NodeIndx(i))%Position(3)+MSL2SWL,  member%tMG(i), NCoefDpth, CoefDpths, member%Cd(i), member%Ca(i), member%Cp(i), member%AxCd(i), member%AxCa(i), member%AxCp(i) )
        elif m.data['CoefMod'] == 'MemberCoefs':
            # Member-based model: coefficients set using member-specific coefficient tables
            # do i = 1, member%NElements + 1
            #  ! Pull member  end-node data from the tables and then linearly interpolate it onto the interior member nodes    
            #  s = (real(i,ReKi)-1.0) / real(member%NElements,ReKi)
            #  if ( member%tMG(i) > 0.0_ReKi ) then
            #     member%Cd    (i) = CoefMembers(MmbrCoefIDIndx)%MemberCdMG1*(1-s) + CoefMembers(MmbrCoefIDIndx)%MemberCdMG2*s
            #     member%Ca    (i) = CoefMembers(MmbrCoefIDIndx)%MemberCaMG1*(1-s) + CoefMembers(MmbrCoefIDIndx)%MemberCaMG2*s
            #     member%Cp    (i) = CoefMembers(MmbrCoefIDIndx)%MemberCpMG1*(1-s) + CoefMembers(MmbrCoefIDIndx)%MemberCpMG2*s 
            #     member%AxCd  (i) = CoefMembers(MmbrCoefIDIndx)%MemberAxCaMG1*(1-s) + CoefMembers(MmbrCoefIDIndx)%MemberAxCdMG2*s
            #     member%AxCa  (i) = CoefMembers(MmbrCoefIDIndx)%MemberAxCaMG1*(1-s) + CoefMembers(MmbrCoefIDIndx)%MemberAxCaMG2*s
            #     member%AxCp  (i) = CoefMembers(MmbrCoefIDIndx)%MemberAxCpMG1*(1-s) + CoefMembers(MmbrCoefIDIndx)%MemberAxCpMG2*s
            #  else
            #     member%Cd    (i) = CoefMembers(MmbrCoefIDIndx)%MemberCd1 *(1-s) + CoefMembers(MmbrCoefIDIndx)%MemberCd2 *s
            #     member%Ca    (i) = CoefMembers(MmbrCoefIDIndx)%MemberCa1 *(1-s) + CoefMembers(MmbrCoefIDIndx)%MemberCa2 *s
            #     member%Cp    (i) = CoefMembers(MmbrCoefIDIndx)%MemberCp1 *(1-s) + CoefMembers(MmbrCoefIDIndx)%MemberCp2 *s
            #     member%AxCd  (i) = CoefMembers(MmbrCoefIDIndx)%MemberAxCd1  *(1-s) + CoefMembers(MmbrCoefIDIndx)%MemberAxCd2  *s
            #     member%AxCa  (i) = CoefMembers(MmbrCoefIDIndx)%MemberAxCa1  *(1-s) + CoefMembers(MmbrCoefIDIndx)%MemberAxCa2  *s
            #     member%AxCp  (i) = CoefMembers(MmbrCoefIDIndx)%MemberAxCp1  *(1-s) + CoefMembers(MmbrCoefIDIndx)%MemberAxCp2  *s
            #  end if
            raise NotImplementedError()
        else:
            raise NotImplementedError()
        # calculate reference incline angle and heading, and related trig values.  Note: members are straight to start
        Za = m.nodes[0].point[2]
        Zb = m.nodes[1].point[2]
        # find fill location of member (previously in SetElementFillProps)
        # TODO Filled members
        if 'FillGroups' in self.graph.NodePropertySets.keys():
            # Find index 
            raise NotImplementedError()
        ID = -1
        if ID>0:
            raise NotImplementedError()
        #             if ( MmbrFilledIDIndx > 0 ) then    
        #                member%FillDens     =  InitInp%FilledGroups(MmbrFilledIDIndx)%FillDens
        #                member%FillFSLoc    =  InitInp%FilledGroups(MmbrFilledIDIndx)%FillFSLoc - InitInp%MSL2SWL
        #                 if (member%FillFSLoc >= Zb) then
        #                   member%z_overfill = member%FillFSLoc - Zb
        #                   member%l_fill = member%RefLength
        #                   member%memfloodstatus = 1  ! fully flooded   
        #                 elseif (Za >= member%FillFSLoc) then
        #                    ! No ballast
        #                   member%memfloodstatus = 0  
        #                   member%z_overfill = 0.0_ReKi
        #                   member%l_fill = 0.0_ReKi
        #                else
        #                   member%z_overfill =0
        #                   if ( Zb <= -InitInp%WtrDpth ) then
        #                      member%memfloodstatus = 0  ! member fully buried in seabed
        #                      member%l_fill = 0
        #                   else
        #                      member%memfloodstatus = 2  ! partially flooded member
        #                      member%l_fill = (member%FillFSLoc - Za)/cosPhi
        else:
            member['FillDens']       = 0.0
            member['FillFSLoc']      = 0.0  # Future calculations for ballasting MUST verify that MbrFilledIDIndx > 0 for any ballasting calcs or this value will cause errors
            member['z_overfill']     = 0
            member['l_fill']         = 0
            member['memfloodstatus'] = 0
        # Check the member does not exhibit any of the following conditions
        if not m.data['Pot']:
            if abs(Zb) < abs(member['RMG'][-1]*sinPhi):
                raise Exception('The upper end-plate of a member must not cross the water plane.  This is not true for Member ID ', m)
            if abs(Za) < abs(member['RMG'][0]*sinPhi):
                raise Exception('The lower end-plate of a member must not cross the water plane.  This is not true for Member ID ', m)
             #if ( ( Za < -WtrDepth and Zb >= -WtrDepth ) and ( phi > 10.0*np.pi/180  or  abs((member['Rmg'][-1] - member['Rmg'][i))/member%RefLength)>0.1 ) ) then
             #   call SetErrStat(ErrID_Fatal, 'A member which crosses the seabed must not be inclined more than 10 degrees from vertical or have a taper larger than 0.1.  This is not true for Member ID '//trim(num2lstr(member%MemberID)), errStat, errMsg, 'SetMemberProperties' )   
        # calculate h_floor if seabed-piercing
        member['h_floor'] = 0
        member['i_floor'] = member['NElements']+1  # Default to entire member is below the seabed
        member['doEndBuoyancy'] = False
        if Za < -WtrDepth:
            raise NotImplementedError()
            # TODO PYTHON MAKE SURE each time i_floor is used, the proper "python indexing is used" not Fortran
            # do i= 2, member%NElements+1
            #                   Za = InitInp%Nodes(member%NodeIndx(i))%Position(3)
            #                   if (Za > -WtrDepth) then            ! find the lowest node above the seabed
            #                      
            #                      if (cosPhi < 0.173648178 ) then ! phi > 80 degrees and member is seabed crossing
            #                         call SetErrStat(ErrID_Fatal, 'A seabed crossing member must have an inclination angle of <= 80 degrees from vertical.  This is not true for Member ID '//trim(num2lstr(member%MemberID)), errStat, errMsg, 'SetMemberProperties' )
            #                      end if
            #                      
            #                      member%h_floor = (-WtrDepth-Za)/cosPhi  ! get the distance from the node to the seabed along the member axis (negative value)
            #                      member%i_floor = i-1                    ! record the number of the element that pierces the seabed
            #                      member%doEndBuoyancy = .true.
            #                      exit
            #                   else if ( EqualRealNos(Za, -WtrDepth ) ) then
            #                      member%doEndBuoyancy = .true.
            #                   end if
            #                end do
        else:
            member['i_floor'] = -1 # lower end is at or above the seabed
        #print('Member',Za, member['doEndBuoyancy'], member['i_floor'])
        # calculate element-level values
        member['dRdl_mg'] = np.diff(member['RMG'])/dl
        member['dRdl_in'] = np.diff(member['Rin'])/dl

        def GetAlpha(R1,R2):
            # calculates relative center of volume location for a (tapered) cylindrical element
            # R1: interior radius of element at node point
            # R2: interior radius of other end of part-element
            return (R1*R1 + 2.0*R1*R2 + 3.0*R2*R2)/4.0/(R1*R1 + R1*R2 + R2*R2)
        member['alpha']    = [GetAlpha(member['RMG'][i], member['RMG'][i+1]) for i in range(N)]
        member['alpha_fb'] = [GetAlpha(member['Rin'][i], member['Rin'][i+1]) for i in range(N)]
        member['Vinner']   = 0 # Total  volume of member without marine growth
        member['Vouter']   = 0 # Total outer volume of member including marine growth
        member['Vballast'] = 0 # Total ballasted volume of member
        member['Vsubmerged'] = 0 



        # --- force-related constants for each element
        from welib.hydro.tools import tapered_cylinder_geom, tapered_cylinder_prop_MG
        Positions = np.vstack( (m.nodes[0].point, member['SubNodesPositions'], m.nodes[1].point) )
        for i in range(N):
            Za = Positions[i  ,2]  # z location of node i
            Zb = Positions[i+1,2]  # z location of node i+1
            # ------------------ marine growth weight and inertia ------------------------------------------------
            Vinner_l   = 0.0
            Vouter_l   = 0.0
            Vinner_U   = 0.0
            Vouter_U   = 0.0
            if i > member['i_floor']:
                # full marine growth: get the properties for each half-element lumped to the appropriate node
                Rmid   = 0.5*(member['R']  [i] + member['R']  [i+1])  # radius at middle of segment, where division occurs
                RmidMG = 0.5*(member['RMG'][i] + member['RMG'][i+1])  # radius with marine growth at middle of segment, where division occurs
                Lmid   = 0.5*dl   # = 0.5*(R2-R1)/m  half-length of segment
                Vinner_l, Vouter_l, member['m_mg_l'][i], member['h_cmg_l'][i], member['I_lmg_l'][i], member['I_rmg_l'][i] = tapered_cylinder_prop_MG(member['R'][i  ], Rmid, member['RMG'][i  ], RmidMG, Lmid, member['MGdensity'][i])
                Vinner_u, Vouter_u, member['m_mg_u'][i], member['h_cmg_u'][i], member['I_lmg_u'][i], member['I_rmg_u'][i] = tapered_cylinder_prop_MG(member['R'][i+1], Rmid, member['RMG'][i+1], RmidMG,-Lmid, member['MGdensity'][i])
            elif i == member['i_floor']:
                raise NotImplementedError()
                #  ! crossing seabed: get the properties for part-element above the seabed and lump to the upper node      
                #  Rmid   = (-member%h_floor*member%R(  i) +(dl+member%h_floor)*member%R(  i+1))/dl
                #  RmidMG = (-member%h_floor*member%RMG(i) +(dl+member%h_floor)*member%RMG(i+1))/dl
                #  Lmid   = -member%h_floor
                #  CALL MarineGrowthPartSegment(member%R(i+1), Rmid, member%RMG(i+1),RmidMG, -Lmid, member%MGDensity(i),  Vinner_u, Vouter_u, member%m_mg_u(i), member%h_cmg_u(i), member%I_lmg_u(i), member%I_rmg_u(i))   ! get precomputed quantities for upper half-segment
                #  Vinner_l   = 0.0
                #  Vouter_l   = 0.0
            # ------------------ flooded ballast inertia ---------------------------------------------------------
            Vballast_l = 0.0
            Vballast_U = 0.0
            if (member['memfloodstatus'] > 0 and (member['FillFSLoc'] > Za)):
                raise NotImplementedError()
                #  # Fully filled element, so split in middle
                #  if ((i > member%i_floor) .and. (member%FillFSLoc >= Zb)) then
                #     ! get the properties for each half-element lumped to the appropriate node
                #     Rmidin = 0.5*(member%Rin(i)+member%Rin(i+1))  ! radius of member interior at middle of segment, where division occurs
                #     Lmid   = 0.5*dl   ! = 0.5*(R2-R1)/m  half-length of segment
                #     CALL FloodedBallastPartSegment(member%Rin(i  ), Rmidin,  Lmid, member%FillDens, Vballast_l, member%m_fb_l(i), member%h_cfb_l(i), member%I_lfb_l(i), member%I_rfb_l(i))   ! get precomputed quantities for lower half-segment
                #     CALL FloodedBallastPartSegment(member%Rin(i+1), Rmidin, -Lmid, member%FillDens, Vballast_u, member%m_fb_u(i), member%h_cfb_u(i), member%I_lfb_u(i), member%I_rfb_u(i))   ! get precomputed quantities for upper half-segment
                #  # partially filled element, so split at FillFSLoc
                #  else if ((i > member%i_floor)  .AND. (member%FillFSLoc < Zb)) then
                #     ! get the properties for each partial-element lumped to the appropriate node
                #     Lmid   = member%FillFSLoc - Za 
                #     Rmidin = member%Rin(i)+(Lmid/(Zb-Za))*(member%Rin(i+1)-member%Rin(i))  ! radius of member interior at middle of segment, where division occurs
                #     CALL FloodedBallastPartSegment(member%Rin(i  ), Rmidin,  Lmid, member%FillDens, Vballast_l, member%m_fb_l(i), member%h_cfb_l(i), member%I_lfb_l(i), member%I_rfb_l(i))   ! get precomputed quantities for lower half-segment
                #     CALL FloodedBallastPartSegment(member%Rin(i+1), Rmidin, -Lmid, 0.0, Vballast_u, member%m_fb_u(i), member%h_cfb_u(i), member%I_lfb_u(i), member%I_rfb_u(i))   ! get precomputed quantities for upper half-segment
                #  else if (i == member%i_floor) then     ! Hopefully we don't have a partially filled element crossing the seabed.
                #     ! crossing seabed: get the properties for part-element above the seabed and lump to the upper node
                #     RmidMG = (-member%h_floor*member%RMG(i) +(dl+member%h_floor)*member%RMG(i+1))/dl
                #     Rmidin = (-member%h_floor*member%Rin(i) +(dl+member%h_floor)*member%Rin(i+1))/dl
                #     Lmid   = -member%h_floor
                #     CALL FloodedBallastPartSegment(member%Rin(i+1), Rmidin, -Lmid, member%FillDens,  Vballast_u, member%m_fb_u(i), member%h_cfb_u(i), member%I_lfb_u(i), member%I_rfb_u(i))   ! get precomputed quantities for upper half-segment
                #     Vballast_l = 0.0
            else:  # Either no ballast flooding in member, or this particular element isn't flooded at all
                Vballast_u           = 0.0
                Vballast_l           = 0.0
                member['m_fb_u'] [i] = 0.0
                member['h_cfb_u'][i] = 0.0
                member['I_lfb_u'][i] = 0.0
                member['I_rfb_u'][i] = 0.0
            # Determine volumes to add to Non-WAMIT modeled members, etc.
            if not m.data['Pot']:
                if Zb < -WtrDepth:
                    pass # fully buried element, do not add these volume contributions to totals
                elif Zb<0.0:
                    # fully submerged elements.  
                    # NOTE: For an element which is fractionaly in the seabed, the entire element volume is added to totals
                    member['Vinner']     += Vinner_l + Vinner_u
                    member['Vouter']     += Vouter_l + Vouter_u
                    member['Vsubmerged'] += Vouter_l + Vouter_u
                elif (0.0 > Za) and (0.0 <= Zb):
                    if (i == 0):
                        raise Exception('The lowest element of a member must not cross the free surface.  This is true for MemberID ',m)
                    # partially submerged element
                    member['Vinner'] += Vinner_l + Vinner_u
                    member['Vouter'] += Vouter_l + Vouter_u
                    # compute volume portion which is submerged
                    Lmid = -Za/cosPhi 
                    Vouter_l, h_c = tapered_cylinder_geom( member['RMG'][i], member['RMG'][i]+Lmid*member['dRdl_mg'][i], Lmid)
                    member['Vsubmerged'] += Vouter_l 
                else: # fully above the water
                    member['Vinner'] += Vinner_l + Vinner_u
                    member['Vouter'] += Vouter_l + Vouter_u
            # ------------------ flooded ballast weight (done) --------------------
            # NOTE: this section of code is somewhat redundant with "flooded ballast inertia" section above
            li = dl*(i-1)
            if Zb < -WtrDepth:
                # fully buried element
                member['floodstatus'][i] = 0
            elif member['memfloodstatus'] > 0 and member['FillFSLoc'] > Zb:
                # fully filled elements 
                raise NotImplementedError()
#                     member%floodstatus(i) = 1
#                     member%Vballast = member%Vballast + Vballast_l + Vballast_u
#                     ! depth-adjusted force distribution constant
#                     member%alpha_fb_star(i) = member%alpha_fb(i)*( Zb - member%FillFSLoc )**3 / ( ( (1-member%alpha_fb(i))*(Za - member%FillFSLoc))**3 + member%alpha_fb(i)*(Zb - member%FillFSLoc)**3 )
#                     ! force and moment magnitude constants
#                     member%Cfl_fb(i) = TwoPi * member%dRdl_in(i) * member%FillDens * gravity * dl *( (li - member%l_fill)*member%Rin(i) + 0.5*((li - member%l_fill)* member%dRdl_in(i) + member%Rin(i))*dl + 1.0/3.0* member%dRdl_in(i)*dl**2 )
#                     member%Cfr_fb(i) =    Pi *                     member%FillDens * gravity * dl *( member%Rin(i)**2 +  member%dRdl_in(i)*member%Rin(i)*dl +1.0/3.0 * member%dRdl_in(i)**2 *dl**2 )
#                     member%CM0_fb(i) = TwoPi *                     member%FillDens * gravity * dl *( 0.25*dl**3* member%dRdl_in(i)**4 + 0.25*dl**3* member%dRdl_in(i)**2 + dl**2* member%dRdl_in(i)**3*member%Rin(i) + 2.0/3.0*dl**2* member%dRdl_in(i)*member%Rin(i) + 1.5*dl* member%dRdl_in(i)**2*member%Rin(i)**2 + 0.5*dl*member%Rin(i)**2 +  member%dRdl_in(i)*member%Rin(i)**3 )
            elif (member['memfloodstatus'] > 0) and (member['FillFSLoc'] > Za) and (member['FillFSLoc'] < Zb): 
                # partially filled element
                # Need to enforce the modeling requirement that the first/bottom-most element of a member be fully flooded
                if (i == 0):
                    raise Exception('The modeling of partially flooded/ballested members requires that the first/bottom-most element of a member must be fully flooded. This is not true for MemberID ',m)
                # Need to enforce the modeling requirement that a partially flooded member must not be close to horizontal
                #if (InitInp%Nodes(member%NodeIndx(N+1))%Position(3) - member%Rin(N+1)*sinPhi) < member%FillFSLoc :
                #    raise Exception('The modeling of partially flooded/ballested members requires the the member not be near horizontal.  This is not true for MemberID ', m)
                member['floodstatus'][i] = 2
                # length along axis from node i to fill level
                member['h_fill'] = member['l_fill'] - i*dl
                raise NotImplementedError()
                #Since this element is only partially flooded/ballasted, compute the Volume fraction which is filled
                #   call TaperCalc( member%Rin(i), member%Rin(i)+member%h_fill*member%dRdl_in(i), member%h_fill, Vballast_l, h_c)
                #   Vballast_u = 0.0
                #   member%Vballast = member%Vballast + Vballast_l + Vballast_u ! Note: Vballast_l will match calculations above
                #   ! depth-adjusted force distribution constant
                #   member%alpha_fb_star(i) = (1 - member%alpha_fb(i))*( Za - member%FillFSLoc )**3 / ( ( (1-member%alpha_fb(i))*(Za - member%FillFSLoc))**3 - member%alpha_fb(i)*(Zb - member%FillFSLoc)**3 )
                #   ! force and moment magnitude constants
                #   member%Cfl_fb(i) = TwoPi * member%dRdl_in(i) * member%FillDens * gravity * member%h_fill *( (li - member%l_fill)*member%Rin(i) + 0.5*((li - member%l_fill)*member%dRdl_in(i) + member%Rin(i))*member%h_fill + 1.0/3.0*member%dRdl_in(i)*member%h_fill**2 )
                #   member%Cfr_fb(i) =    Pi * member%FillDens * gravity * member%h_fill *( member%Rin(i)**2 + member%dRdl_in(i)*member%Rin(i)*member%h_fill +1.0/3.0 *member%dRdl_in(i)**2 *member%h_fill**2 )
                #   member%CM0_fb(i) = TwoPi * member%FillDens * gravity * member%h_fill *( 0.25*member%h_fill**3*member%dRdl_in(i)**4 + 0.25*member%h_fill**3*member%dRdl_in(i)**2 + member%h_fill**2*member%dRdl_in(i)**3*member%Rin(i) + 2.0/3.0*member%h_fill**2*member%dRdl_in(i)*member%Rin(i)  &
                #                                                                           + 1.5*member%h_fill*member%dRdl_in(i)**2*member%Rin(i)**2 + 0.5*member%h_fill*member%Rin(i)**2 + member%dRdl_in(i)*member%Rin(i)**3 ) &
                #                              -0.25 * member%FillDens * gravity * Pi * (  member%Rin(i) + member%h_fill*member%dRdl_in(i))**4
                # 
            else: # unflooded element
                member['floodstatus'][i] = 0


    def _setupMembers(self):
        """ see SetupMembers in Morison.f90 
          - swap nodes if needed 
          - compute member properties
        """
        graph = self.graph
        # --- Swap Nodes if needed
        for m in graph.Elements:
            # See FlipMemberNodeData in Morison .f90
            doSwap=False # Y1>Y2
            pos1 = m.nodes[0].point
            pos2 = m.nodes[1].point
            if pos1[2] == pos2[2]:  # Z1=Z2
                if pos1[0] == pos2[0]: # X1=X2
                    if pos1[1] > pos2[1]: 
                        doSwap=True # Y1>Y2
                elif pos1[0] > pos2[0]:
                    doSwap=True # X1>X2
            elif pos1[2] > pos2[2]:
                doSwap=True # Z1>Z2

            if doSwap:
                #print('>>> Swapping mID',m.ID, 'NodeIDs:', m.nodeIDs)
                m.swapNodes()
                m.MorisonData['SubNodesPositions'] =  m.MorisonData['SubNodesPositions'][-1::-1,:]
                m.MorisonData['nodeIDs']           =  m.MorisonData['nodeIDs'][-1::-1]
            m.flipped=doSwap
        # Update connectivity
        graph.updateConnectivity()
        # --- Set Member properties
        for m in graph.Elements:
            # see SetMemberProperties in Morison.f90
            member = m.MorisonData
            self._allocateMemberDataArrays(member)
            self._setMemberProperties(member, m)


    def calcOutput(self, t, x=None, xd=None, xo=None, u=None, y=None, opts=None):
        """ NOTE: Morison has no state
        u:


        """


        graph = self.graph
        p = self.p
        m = self.m
        umesh = u['Mesh']
        ymesh = y['Mesh']
        g = p['Gravity']

        # Calculation options
        defaultOpts ={'verbose':False, 'MG':True, 'Buoyancy':True, 'Ballast':True, 'HydroD':True, 'HydroA':True, 'HydroI':True, 'End':True}
        if opts is not None:
            defaultOpts.update(opts)
        bVerbose = defaultOpts['verbose'] 
        bMG      = defaultOpts['MG']
        bBuoy    = defaultOpts['Buoyancy']
        bBallast = defaultOpts['Ballast']
        bHydroD  = defaultOpts['HydroD']
        bHydroA  = defaultOpts['HydroA']
        bHydroI  = defaultOpts['HydroI']
        bEnd     = defaultOpts['End']

        def myprint(*args, **kwargs):
            if bVerbose:
                print(*args, **kwargs)
            else:
                pass


        myprint()
        if bVerbose:
            print('---------------------------------------------------------------------------')
            print('--------------------------- CALCOUTPUT {:10.4f}'.format(t))
            print('---------------------------------------------------------------------------')
            umesh.printDebug()


        #===============================================================================================
        # Calculate the fluid kinematics at all mesh nodes and store for use in the equations below
        #    InterpolationSlope = GetInterpolationSlope(Time, p, m, IntWrapIndx)
        iTime = 0 # TODO TODO interpolate wave kinematics at time
        myprint('---------------------------FLUID KINEMATICS')
        for j, pos in enumerate(self.NodesBeforeSwap):    
            m['nodeInWater'][j] = p['nodeInWater'][iTime,j]
            m['FDynP'][j]       = p['WaveDynP'   ][iTime,j]   # TODO InterpolateWithSlope(InterpolationSlope, m%LastIndWave, p%WaveDynP(:,j))
            m['FA'][0,j]        = p['WaveAcc'    ][iTime,j,0] # TODO InterpolateWithSlope(InterpolationSlope, m%LastIndWave, p%WaveAcc(:,j,i)) 
            m['FA'][1,j]        = p['WaveAcc'    ][iTime,j,1] 
            m['FA'][2,j]        = p['WaveAcc'    ][iTime,j,2] 
            m['FV'][0,j]        = p['WaveVel'    ][iTime,j,0] # TODO InterpolateWithSlope(InterpolationSlope, m%LastIndWave, p%WaveVel(:,j,i)) 
            m['FV'][1,j]        = p['WaveVel'    ][iTime,j,1] 
            m['FV'][2,j]        = p['WaveVel'    ][iTime,j,2] 
            m['vrel'][:,j]      = m['FV'][:,j] - umesh.TranslationVel[j,:]
            #m['vrel'][0,j]      = 1.5
            #m['vrel'][2,j]      = 0.5
            #m['FA'][0,j]       = 0.5
            #m['FA'][2,j]       = 0.1
            #print('Node',j, m['nodeInWater'][j], m['FDynP'][j], m['FV'][:,j])
        # ==============================================================================================
        # Calculate instantaneous loads on each member except for the hydrodynamic loads on member ends.
        # This covers aspects of the load calculations previously in CreateDistributedMesh.  

        from welib.FEM.utils import DCM
        # Zero out previous time-steps loads (these are loads which are computed at the member-level and summed onto a node, 
        # so they need to be zeroed out before the summations happen)
        m['F_BF_End'] *= 0
        m['F_B_End']  *= 0
        ymesh.Force  *= 0
        ymesh.Moment *= 0
        # Very important transfer disp/rot
        umesh.transferMotion2IdenticalMesh(ymesh)

        for im, e in enumerate(graph.Elements):
            myprint('---------------------------MEMBER LOADS ',im+1)
            mem      = e.MorisonData
            N        = mem['NElements']
            memLoads = mem['Loads']
            idx      = mem['nodeIDs']
            #zero member loads
            for k in ['F_D','F_A','F_B','F_BF','F_I','F_If','F_WMG','F_IMG']:
                mem['Loads'][k] *= 0
            # --- Internal Hydrodynamic loads
            # Loop through member elements
            for i in range(N):
                myprint('---------------------------ELEMENT LOADS ',im+1, i+1)
                #print('idx',np.array(idx)+1)
                # calculate instantaneous incline angle and heading, and related trig values
                pos1    = umesh.TranslationDisp[idx[i]  ,:] + umesh.Position[idx[i],:]
                pos1[2] = pos1[2] - p['MSL2SWL']
                pos2    = umesh.TranslationDisp[idx[i+1],:] + umesh.Position[idx[i+1],:]
                pos2[2] = pos2[2] - p['MSL2SWL']
                phi, sinPhi, cosPhi, tanPhi, sinBeta, cosBeta, k_hat= getOrientationAngles(pos1, pos2)
                CMatrix = Morison_DirCosMtrx(pos1, pos2) # TODO harmony with DCM from FEM
                #print('pos1',pos1)
                #print('pos2',pos2)
                #print('phi',phi, sinBeta, k_hat)
                #print('CMatrix',CMatrix[0,:])
                #print('CMatrix',CMatrix[1,:])
                #print('CMatrix',CMatrix[2,:])
                CTrans  = CMatrix.T
                # save some commonly used variables   
                dl      = mem['dl']
                z1      = pos1[2]           # get node z locations from input mesh
                z2      = pos2[2]
                r1      = mem['RMG'][i ]   # outer radius element nodes including marine growth
                r2      = mem['RMG'][i+1]
                dRdl_mg = mem['dRdl_mg'][i] # Taper of element including marine growth
                a_s1    = umesh.TranslationAcc[idx[i  ], :]
                alpha_s1= umesh.RotationAcc   [idx[i  ], :]
                omega_s1= umesh.RotationVel   [idx[i  ], :]
                a_s2    = umesh.TranslationAcc[idx[i+1], :]
                alpha_s2= umesh.RotationAcc   [idx[i+1], :]
                omega_s2= umesh.RotationVel   [idx[i+1], :]
                if bMG:
                    if not e.data['Pot']: # Member is NOT modeled with Potential Flow Theory
                        # should i_floor theshold be applied to below calculations to avoid wasting time on computing zero-valued things? <<<<<
                        # should lumped half-element coefficients get combined at initialization? <<<
                        # --------------------------------------------------------------------------------}
                        # --- Marine growth 
                        # --------------------------------------------------------------------------------{
                        # Sides: Section 4.1.2
                        F_WMG = np.zeros(6)
                        # lower node
                        F_WMG[2] = - mem['m_mg_l'][i]*g # weight force  : Note: this is a constant
                        F_WMG[3] = - mem['m_mg_l'][i]*g * mem['h_cmg_l'][i]* sinPhi * sinBeta# weight force
                        F_WMG[4] =   mem['m_mg_l'][i]*g * mem['h_cmg_l'][i]* sinPhi * cosBeta# weight force
                        memLoads['F_WMG'][:,i] += F_WMG
                        ymesh.Force [idx[i],:] += F_WMG[:3]
                        ymesh.Moment[idx[i],:] += F_WMG[3:]
                        # upper node
                        F_WMG[2] = - mem['m_mg_u'][i]*g # weight force  : Note: this is a constant 
                        F_WMG[3] = - mem['m_mg_u'][i]*g * mem['h_cmg_u'][i]* sinPhi * sinBeta# weight force
                        F_WMG[4] =   mem['m_mg_u'][i]*g * mem['h_cmg_u'][i]* sinPhi * cosBeta# weight force
                        memLoads['F_WMG'][:,i+1] += F_WMG  
                        ymesh.Force [idx[i+1],:] += F_WMG[:3]
                        ymesh.Moment[idx[i+1],:] += F_WMG[3:]
                        F_IMG = np.zeros(6)
                        Imat  = np.zeros((3,3))
                        # lower node 
                        Ioffset   = mem['h_cmg_l'][i]*mem['h_cmg_l'][i]*mem['m_mg_l'][i]
                        Imat[0,0] = mem['I_rmg_l'][i] - Ioffset
                        Imat[1,1] = mem['I_rmg_l'][i] - Ioffset
                        Imat[2,2] = mem['I_lmg_l'][i] - Ioffset
                        Imat      =  CMatrix.dot(Imat).dot(CTrans)
                        iArm = mem['h_cmg_l'][i] * k_hat
                        iTerm     = ( -a_s1 - np.cross(omega_s1, np.cross(omega_s1,iArm )) - np.cross(alpha_s1,iArm) ) * mem['m_mg_l'][i]
                        F_IMG[:3] = iTerm
                        F_IMG[3:] = - np.cross(a_s1 * mem['m_mg_l'][i], mem['h_cmg_l'][i] * k_hat) + Imat.dot(alpha_s1)  - np.cross(omega_s1,Imat.dot(omega_s1))
                        memLoads['F_IMG'][:,i] +=F_IMG
                        ymesh.Force [idx[i], :] +=F_IMG[:3]
                        ymesh.Moment[idx[i], :] +=F_IMG[3:]
                        # upper node
                        Ioffset   = mem['h_cmg_u'][i]*mem['h_cmg_u'][i]*mem['m_mg_u'][i]
                        Imat[0,0] = mem['I_rmg_u'][i] - Ioffset
                        Imat[1,1] = mem['I_rmg_u'][i] - Ioffset
                        Imat[2,2] = mem['I_lmg_u'][i] - Ioffset
                        Imat      = CMatrix.dot(Imat).dot(CTrans)
                        iArm = mem['h_cmg_u'][i] * k_hat
                        iTerm     = ( -a_s2 - np.cross(omega_s2, np.cross(omega_s2,iArm )) - np.cross(alpha_s2,iArm) ) * mem['m_mg_u'][i]
                        F_IMG[:3] = iTerm
                        F_IMG[3:] = - np.cross(a_s2 * mem['m_mg_u'][i], mem['h_cmg_u'][i] * k_hat) + Imat.dot(alpha_s2) - np.cross(omega_s2,Imat.dot(omega_s2))
                        memLoads['F_IMG'][:,i+1] += F_IMG
                        ymesh.Force [idx[i+1], :] += F_IMG[:3]
                        ymesh.Moment[idx[i+1], :] += F_IMG[3:]
        
                if bBuoy:
                    if not e.data['Pot']: # Member is NOT modeled with Potential Flow Theory
                        # --------------------------------------------------------------------------------}
                        # ---Buoyancy loads
                        # --------------------------------------------------------------------------------{
                        # sides: Sections 3.1 and 3.2 ------------------------
                        if z1 < 0:  # if segment is at least partially submerged ...
                            if z1*z2 <= 0: # special calculation if the slice is partially submerged
                                # Check that this is not the 1st element of the member
                                if i==0:
                                    raise Exception('The lowest element of a Morison member has become partially submerged!  This is not allowed.  Please review your model and create a discretization such that even with displacements, the lowest element of a member does not become partially submerged.')
                                h0 = -z1/cosPhi             # distances along element centerline from point 1 to the waterplane
                                if abs(dRdl_mg)< 0.0001:      # untapered cylinder case
                                    Vs =    np.pi*r1*r1*h0   # volume of total submerged portion
                                    if Vs ==0:
                                        cx = 0.0  # Avoid singularity, but continue to provide the correct solution
                                    else:
                                        cr = 0.25*r1*r1*tanPhi/h0
                                        cl = 0.5*h0 + 0.125*r1*r1*tanPhi*tanPhi/h0
                                        cx = cr*cosPhi + cl*sinPhi
                                else: # inclined tapered cylinder case (note I've renamed r0 to rh here##)
                                    #===================
                                    #Per plan equations
                                    # NOTE:  Variable changes of Plan     vs       Code
                                    #---------------------------------------------------
                                    #                             V                 Vs
                                    #                             a_h               a0
                                    #                             b_h               b0
                                    #                             x_c               cx
                                    #                             h                 h0
                                    #                             r1                r_MG,i
                                    #                             r_c               cr
                                    #                             h_c               cl
                                    # NOTE: a0 and b0 always appear as a0b0, never separately.
                                    rh   = r1 + h0*dRdl_mg    # radius of element at point where its centerline crosses the waterplane
                                    C_1  = 1.0 - dRdl_mg**2 * tanPhi**2
                                    # waterplane ellipse shape
                                    b0   = rh/np.sqrt(C_1)
                                    a0   = rh/((C_1)*cosPhi)             # simplified from what's in ConicalCalcs.ipynb
                                    a0b0 = a0*b0
                                    C_2  = a0b0*rh*cosPhi - r1**3
                                    cl   = -(-0.75*a0b0*rh**2*cosPhi + 0.75*r1**4*C_1 + r1*C_1*C_2) / (dRdl_mg*C_1*C_2)
                                    cr   = (0.75*a0b0*dRdl_mg*rh**2*sinPhi)/(C_1*C_2)
                                    cx   = cr*cosPhi + cl*sinPhi 
                                    Vs   = np.pi*(a0b0*rh*cosPhi - r1**3)/(3.0*dRdl_mg)       
                                    # End per plan equations
                                    #===================
                                pwr = 3
                                alpha    = (1.0-mem['alpha'][i])*z1**pwr/(-mem['alpha'][i]*z2**pwr + (1.0-mem['alpha'][i])*z1**pwr)
                                Fb  = Vs*p['WtrDens']*g       #buoyant force
                                Fr  = -Fb*sinPhi     #radial component of buoyant force
                                Fl  = Fb*cosPhi      #axial component of buoyant force
                                Moment = -Fb*cx      #This was matt's code        #moment induced about the center of the cylinder's bottom face
                                # calculate (imaginary) bottom plate forces/moment to subtract from displacement-based values
                                Fl  = Fl  + p['WtrDens']*g*z1* np.pi *r1*r1        
                                Moment  = Moment  + p['WtrDens']*g* sinPhi * np.pi/4.0*r1**4       
                                # reduce taper-based moment to remove (not double count) radial force distribution to each node 
                                Moment  = Moment + Fr*(1.0-alpha)*dl
                                F_B1, F_B2 = DistributeElementLoads(Fl, Fr, Moment, sinPhi, cosPhi, sinBeta, cosBeta, alpha)
                                #print('Case 1')
                                #print('>>> FB_1',F_B1[:3])
                                #print('>>> FB_2',F_B2[:3])
                                memLoads['F_B'][:, i]   += F_B1  # alpha
                                memLoads['F_B'][:, i-1] += F_B2  # 1-alpha
                                ymesh.Force [idx[i  ], :] += F_B1[:3]
                                ymesh.Moment[idx[i  ], :] += F_B1[3:]
                                ymesh.Force [idx[i-1], :] += F_B2[:3]
                                ymesh.Moment[idx[i-1], :] += F_B2[3:]
                            else: # normal, fully submerged case
                                Fl = -2.0*np.pi*dRdl_mg*p['WtrDens']*g*dl*( z1*r1 + 0.5*(z1*dRdl_mg + r1*cosPhi)*dl + 1.0/3.0*(dRdl_mg*cosPhi*dl*dl) )   # from CylinderCalculationsR1.ipynb
                                Fr = -np.pi*p['WtrDens']*g*dl*(r1*r1 + dRdl_mg*r1*dl + (dRdl_mg**2*dl**2)/3.0)*sinPhi                          # from CylinderCalculationsR1.ipynb
                                Moment = -np.pi*dl*g*p['WtrDens']*(3.0*dl**3*dRdl_mg**4 + 3.0*dl**3*dRdl_mg**2 + 12.0*dl**2*dRdl_mg**3*r1 + 8.0*dl**2*dRdl_mg*r1 + 18.0*dl*dRdl_mg**2*r1*r1 + 6.0*dl*r1*r1 + 12.0*dRdl_mg*r1**3)*sinPhi/12.0   # latest from CylinderCalculationsR1.ipynb
             
                                # precomputed as mem['alpha[i] ... alpha0 = (r1*r1 + 2*r1*r2 + 3*r2**2)/4/(r1*r1 + r1*r2 + r2**2)
                                z1d = -min(0.0,z1)
                                z2d = -min(0.0,z2)
                                pwr = 3
                                alpha = mem['alpha'][i]*z2d**pwr/(mem['alpha'][i]*z2d**pwr+(1-mem['alpha'][i])*z1d**pwr)
                                # reduce moment to remove (not double count) radial force distribution to each node
                                Moment = Moment - Fr*alpha*dl
                                F_B1, F_B2 = DistributeElementLoads(Fl, Fr, Moment, sinPhi, cosPhi, sinBeta, cosBeta, alpha)
                                #print('Case 2')
                                #print('>>> FB_1',F_B1[:3])
                                #print('>>> FB_2',F_B2[:3])
                                memLoads['F_B'][:,i+1] += F_B1  # alpha
                                memLoads['F_B'][:, i]  += F_B2  # 1-alpha
                                ymesh.Force [idx[i  ], :] += F_B2[:3]
                                ymesh.Moment[idx[i  ], :] += F_B2[3:]
                                ymesh.Force [idx[i+1], :] += F_B1[:3]
                                ymesh.Moment[idx[i+1], :] += F_B1[3:]
                    # --- End If not potential element
                # --------------------------------------------------------------------------------}
                # --- Flooded ballast (for Pot or not Pot)
                # --------------------------------------------------------------------------------{
                if bBallast:
                    # --- Inertia Section 6.1.1 
                    Imat = np.zeros((3,3))
                    F_If = np.zeros(6)
                    # lower node
                    Ioffset   = mem['h_cfb_l'][i]*mem['h_cfb_l'][i]*mem['m_fb_l'][i]
                    Imat[0,0] = mem['I_rfb_l'][i] - Ioffset
                    Imat[1,1] = mem['I_rfb_l'][i] - Ioffset
                    Imat[2,2] = mem['I_lfb_l'][i] - Ioffset
                    iArm = mem['h_cfb_l'][i] * k_hat
                    iTerm     = ( -a_s1  - np.cross(omega_s1, np.cross(omega_s1,iArm ))  -  np.cross(alpha_s1,iArm) ) * mem['m_fb_l'][i]
                    F_If[:3] =  iTerm
                    F_If[3:] =  - np.cross(a_s1 * mem['m_fb_l'][i], mem['h_cfb_l'][i] * k_hat) + Imat.dot(alpha_s1) - np.cross(omega_s1,Imat.dot(omega_s1)) 
                    memLoads['F_If'][:,i] += F_If
                    ymesh.Force [idx[i],:] += F_If[:3]
                    ymesh.Moment[idx[i],:] += F_If[3:6]
                    # upper node
                    Ioffset   = mem['h_cfb_u'][i]*mem['h_cfb_u'][i]*mem['m_fb_u'][i]
                    Imat[0,0] = mem['I_rfb_u'][i] - Ioffset
                    Imat[1,1] = mem['I_rfb_u'][i] - Ioffset
                    Imat[2,2] = mem['I_lfb_u'][i] - Ioffset
                    iArm = mem['h_cfb_u'][i] * k_hat
                    iTerm     = ( -a_s2  - np.cross(omega_s2, np.cross(omega_s2,iArm ))  -  np.cross(alpha_s2,iArm) ) * mem['m_fb_u'][i]
                    F_If[:3] = iTerm
                    F_If[3:] = - np.cross(a_s2 * mem['m_fb_u'][i], mem['h_cfb_u'][i] * k_hat) + Imat.dot(alpha_s2) - np.cross(omega_s2,Imat.dot(omega_s2))
                    memLoads['F_If'][:,i+1]  += F_If
                    ymesh.Force [idx[i+1],:] += F_If[:3]
                    ymesh.Moment[idx[i+1],:] += F_If[3:]  

                    # -- flooded ballast weight : sides : Section 5.1.2 & 5.2.2  : Always compute regardless of PropPot setting
                    # NOTE: For memfloodstatus and floodstatus: 0 = fully buried or not ballasted, 1 = fully flooded, 2 = partially flooded
                    # fully filled elements
                    if mem['floodstatus'][i]==1:
                        # Compute lstar
                        if mem['memfloodstatus'] == 2:
                            # partially flooded MEMBER
                            lstar = dl*(i-1) - mem['l_fill']
                        elif (cosPhi >= 0.0 ):
                            lstar = dl*(i-N-1) 
                        else:
                            lstar = dl*(i-1)
                        Fl =2*np.pi * mem['dRdl_in'][i] * mem['FillDens'] * p['gravity'] * dl *( -( mem['Rin'][i] + 0.5* mem['dRdl_in'][i]*dl )*mem['z_overfill'] + ( lstar*mem['Rin'][i] + 0.5*(lstar*mem['dRdl_in'][i] + mem['Rin'][i] )*dl + mem['dRdl_in'][i]*dl**2/3.0 )*cosphi )
                        # forces and moment in tilted coordinates about node i
                        Fr = mem['Cfr_fb'][i]*sinPhi     
                        Moment  = mem['CM0_fb'][i]*sinPhi - Fr*mem['alpha_fb_star'][i]*dl
                        # calculate full vector and distribute to nodes
                        F_B1, F_B2 = DistributeElementLoads(Fl, Fr, Moment, sinPhi, cosPhi, sinBeta, cosBeta, (1-mem['alpha_fb_star'][i]))
                        memLoads['F_BF'][:, i]   += F_B2 # 1-alpha
                        memLoads['F_BF'][:, i+1] += F_B1 # alpha
                        ymesh.Force [idx[i  ], :] += F_B2[:3]
                        ymesh.Moment[idx[i  ], :] += F_B2[3:]
                        ymesh.Force [idx[i+1], :] += F_B1[:3]
                        ymesh.Moment[idx[i+1], :] += F_B1[3:]
                    elif mem['floodstatus'][i] == 2:
                        # partially filled element
                        # forces and moment in tilted coordinates about node i
                        Fl     = mem['Cfl_fb'][i]*cosPhi
                        Fr     = mem['Cfr_fb'][i]*sinPhi
                        Moment = mem['CM0_fb'][i]*sinPhi + Fr*(1 - mem['alpha_fb_star'][i])*dl
                        # calculate full vector and distribute to nodes
                        F_B1, F_B2 = DistributeElementLoads(Fl, Fr, Moment, sinPhi, cosPhi, sinBeta, cosBeta, mem['alpha_fb_star'][i])
                        memLoads['F_BF'][:, i]   += F_B1 # alpha
                        memLoads['F_BF'][:, i-1] += F_B2 # 1- alpha
                        ymesh.Force [idx[i  ], :] += F_B1[:3]
                        ymesh.Moment[idx[i  ], :] += F_B1[3:]
                        ymesh.Force [idx[i-1], :] += F_B2[:3]
                        ymesh.Moment[idx[i-1], :] += F_B2[3:]    
                    else:
                        pass
                    # no load for unflooded element or element fully below seabed
            # --- End loop on Member Elements i=1,N

            # --------------------------------------------------------------------------------}
            # --- Hydrodynamic loads on sides 
            # --------------------------------------------------------------------------------{
            # NOTE: All geometry-related calculations are based on the undisplaced configuration of the structure
            for i in range(N+1): # Loop through member nodes
                myprint('---------------------------NODAL LOADS ',im+1, i+1)
                # We need to subtract the MSL2SWL offset to place this in the SWL reference system
                z1 = umesh.Position[idx[i],2] - p['MSL2SWL']
                if i > mem['i_floor'] and z1 <= 0.0:  # node is above (or at? TODO: check) seabed and below or at free-surface)
                    # TODO: Note that for computational efficiency, we could precompute h_c and deltal for each element when we are NOT using wave stretching
                    # We would still need to test at time marching for nodes just below the free surface because that uses the current locations not the reference locations
                    # see table in Section 7.1.1
                    if i == 0:
                        deltal = mem['dl']/2.0
                        h_c    = mem['dl']/4.0
                    elif (i == N):
                        deltal =  mem['dl']/2.0
                        h_c    = -mem['dl']/4.0
                    elif mem['i_floor'] == i+1 : # This node is the upper node of an element which crosses the seabed
                        deltal = mem['dl']/2.0 - mem['h_floor']  # TODO: h_floor is negative valued, should we be subrtracting it from dl/2? GJH
                        h_c    = 0.5*(mem['dl']/2.0 + mem['h_floor'])
                    else:
                        # We need to subtract the MSL2SWL offset to place this  in the SWL reference system
                        pos1    = umesh.Position[idx[i], :]
                        pos1[2] = pos1[2] - p['MSL2SWL']
                        pos2    = umesh.Position[idx[i+1],:]
                        pos2[2] = pos2[2] - p['MSL2SWL']
                        if pos1[2] <= 0.0 and 0.0 < pos2[2]: # This node is just below the free surface #TODO: Needs to be augmented for wave stretching
                            # We need to subtract the MSL2SWL offset to place this  in the SWL reference system
                            #TODO: Fix this one
                            pos1    = umesh.Position[idx[i],:] # use reference position for following equation
                            pos1[2] = pos1[2] - p['MSL2SWL']
                            h       = ( pos1[2] ) / mem['cosPhi_ref'] #TODO: Needs to be augmented for wave stretching
                            deltal  = mem['dl']/2.0 + h
                            h_c     = 0.5*(h-mem['dl']/2.0)
                        else:
                            # This node is a fully submerged interior node
                            deltal = mem['dl']
                            h_c    = 0.0
                    if i == 0:
                        dRdl_p  = abs(mem['dRdl_mg'][i])
                        dRdl_pp = mem['dRdl_mg'][i]   
                    elif i > 0 and i<N:
                        dRdl_p  = 0.5*( abs(mem['dRdl_mg'][i-1]) + abs(mem['dRdl_mg'][i]) )
                        dRdl_pp = 0.5*( mem['dRdl_mg'][i-1] + mem['dRdl_mg'][i] )
                    else:
                        dRdl_p  = abs(mem['dRdl_mg'][-1])
                        dRdl_pp = mem['dRdl_mg'][-1]
                    # ------------------- hydrodynamic drag loads: sides: Section 7.1.2 ------------------------ 
                    if bHydroD:
                        vec = mem['Ak'].dot(m['vrel'][:,idx[i]] )
                        dotp = mem['k'].dot( m['vrel'][:,idx[i]].flatten())
                        KV   = dotp*mem['kkt'] # 3x3 matrix
                        vec2 = KV .dot( m['vrel'][:,idx[i]])
                        f_hydro = mem['Cd'][i]*p['WtrDens']*mem['RMG'][i]* np.linalg.norm(vec)*vec  + 0.5*mem['AxCd'][i]*p['WtrDens']*np.pi*mem['RMG'][i]*dRdl_p * vec2
                        myprint('f_hydro_d {:16.4f}{:16.4f}{:16.4f}'.format(*f_hydro))
                        #print('t1 ', mem['Cd'][i]*p['WtrDens']*mem['RMG'][i]* np.linalg.norm(vec)*vec )
                        #print('t10', mem['Cd'][i]*p['WtrDens']*mem['RMG'][i]                   )
                        #print('t12', np.linalg.norm(vec)*vec )
                        #print('t2 ', 0.5*mem['AxCd'][i]*p['WtrDens']*np.pi*mem['RMG'][i]*dRdl_p * vec2)
                        #print('t20', 0.5*mem['AxCd'][i]*p['WtrDens']*np.pi*mem['RMG'][i]*dRdl_p)
                        #print('t21', vec2)
                        memLoads['F_D'][:, i]  = LumpDistrHydroLoads( f_hydro, mem['k'], deltal, h_c )
                        ymesh.Force [idx[i], :] += memLoads['F_D'][:3, i]
                        ymesh.Moment[idx[i], :] += memLoads['F_D'][3:, i]
                    if not e.data['Pot']:
                        # ------------------- hydrodynamic added mass loads: sides: Section 7.1.3 ------------------------
                        if bHydroA:
                            Am = mem['Ca'][i]*p['WtrDens']*np.pi*mem['RMG'][i]*mem['RMG'][i]*mem['Ak'] + 2.0*mem['AxCa'][i]*p['WtrDens']*np.pi*mem['RMG'][i]*mem['RMG'][i]*dRdl_p*mem['kkt']
                            f_hydro = - Am.dot( umesh.TranslationAcc[idx[i], :] )
                            myprint('f_hydro_a {:16.4f}{:16.4f}{:16.4f}'.format(*f_hydro))
                            memLoads['F_A'][:, i]  = LumpDistrHydroLoads( f_hydro, mem['k'], deltal, h_c )
                            ymesh.Force [idx[i], :] += memLoads['F_A'][:3, i]
                            ymesh.Moment[idx[i], :] += memLoads['F_A'][3:, i]
                        # ------------------- hydrodynamic inertia loads: sides: Section 7.1.4 ------------------------
                        if bHydroI:
                            t1 = (mem['Ca'][i]+mem['Cp'][i])*p['WtrDens']*np.pi*mem['RMG'][i]*mem['RMG'][i]*          mem['Ak'] .dot( m['FA'][:,idx[i]] )
                            t2 =          2.0*mem['AxCa'][i]*p['WtrDens']*np.pi*mem['RMG'][i]*mem['RMG'][i]*dRdl_p *  mem['kkt'].dot( m['FA'][:,idx[i]] )
                            t3 =          2.0*m['FDynP'][idx[i]]*mem['AxCp'][i]*np.pi*mem['RMG'][i]*dRdl_pp*mem['k'] 
                            f_hydro = t1 + t2 +t3
                            myprint('f_hydro_i {:16.4f}{:16.4f}{:16.4f}'.format(*f_hydro))
                            memLoads['F_I'][:, i] = LumpDistrHydroLoads( f_hydro, mem['k'], deltal, h_c)
                            ymesh.Force [idx[i], :] += memLoads['F_I'][:3, i]
                            ymesh.Moment[idx[i], :] += memLoads['F_I'][3:, i]
                # --- End loop through nodes

            # --------------------------------------------------------------------------------}
            # --- End plate loads, per member
            # --------------------------------------------------------------------------------{
            if bEnd:
                # Any end plate loads that are modeled on a per-member basis
                # reassign convenience variables to correspond to member ends
                # We need to subtract the MSL2SWL offset to place this  in the SWL reference system
                pos1    = umesh.TranslationDisp[idx[0], :] + umesh.Position[idx[0], :] 
                pos2    = umesh.TranslationDisp[idx[1], :] + umesh.Position[idx[1], :] 
                pos1[2] -= p['MSL2SWL']
                pos2[2] -= p['MSL2SWL']
                z1 = pos1[2]
                phi1, sinPhi1, cosPhi1, tanPhi, sinBeta1, cosBeta1, k_hat1 = getOrientationAngles(pos1, pos2)
                if N == 1:       # Only one element in member
                    sinPhi2  = sinPhi1
                    cosPhi2  = cosPhi1
                    sinBeta2 = sinBeta1
                    cosBeta2 = cosBeta1
                else:
                    pos1    = umesh.TranslationDisp[idx[-2], :] + umesh.Position[idx[-2], :]
                    pos2    = umesh.TranslationDisp[idx[-1], :] + umesh.Position[idx[-1], :]
                    pos1[2] -= p['MSL2SWL']
                    pos2[2] -= p['MSL2SWL']
                    phi2, sinPhi2, cosPhi2, tanPhi, sinBeta2, cosBeta2, k_hat2 = getOrientationAngles(pos1, pos2)
                pos2    = umesh.TranslationDisp[idx[-1], :] + umesh.Position[idx[-1], :]
                pos2[2] -= p['MSL2SWL']
                z2 = pos2[2]
                # Check the member does not exhibit any of the following conditions
                if not e.data['Pot']:
                    if abs(z2) < abs(mem['RMG'][-1]*sinPhi2):
                        raise Exception('The upper end-plate of a member must not cross the water plane.  This is not true for Member ID ',e)
                    if abs(z1) < abs(mem['RMG'][0]*sinPhi1):
                       raise Exception('The lower end-plate of a member must not cross the water plane.  This is not true for Member ID ',e)
                # TODO: Do the equations below still work if z1 > z2 ?
                #TODO, should not have to test seabed crossing in time-marching loop
                if mem['i_floor'] == -1:   # both ends are above seabed
                   #--- Water ballast buoyancy ---
                   # if member is fully flooded
                   if mem['memfloodstatus'] == 1:
                       Fl      = -mem['FillDens'] * g * np.pi *     mem['Rin'][  0]**2* (mem['z_overfill'] + max(z2-z1, 0.0))
                       Moment  =  mem['FillDens'] * g * np.pi *0.25*mem['Rin'][  0]**4*sinPhi
                       m['F_BF_End'][:, idx[0]] += endLoad(Fl, Moment, sinPhi1, cosPhi1, sinBeta1, cosBeta1)
                       Fl      =   mem['FillDens'] * g * np.pi *     mem['Rin'][-1]**2* (mem['z_overfill'] + max(z1-z2, 0.0))
                       Moment  =  -mem['FillDens'] * g * np.pi *0.25*mem['Rin'][-1]**4*sinPhi            
                       m['F_BF_End'][:, idx[-1]] += endLoad(Fl, Moment, sinPhi2, cosPhi2, sinBeta2, cosBeta2)
                   elif mem['l_fill'] > 0:  # if member is partially flooded
                       Fl      = -mem['FillDens'] * g * np.pi      *mem['Rin'][0]**2*mem['l_fill']*cosPhi
                       Moment  =  mem['FillDens'] * g * np.pi *0.25*mem['Rin'][0]**4*sinPhi
                       m['F_BF_End'][:, idx[0]] += endLoad(Fl, Moment, sinPhi1, cosPhi1, sinBeta1, cosBeta1)
                   else:
                       pass # no load if member is not flooded at all
                elif  mem['i_floor'] < mem['NElements']:  # upper node is still above the seabed, but lower node is below seabed
                    if (mem['memfloodstatus'] == 1):
                        Fl      =   mem['FillDens'] * g * np.pi *     mem['Rin'][-1]**2* (mem['z_overfill'] + max(z1-z2, 0.0))
                        Moment  =  -mem['FillDens'] * g * np.pi *0.25*mem['Rin'][-1]**4*sinPhi            
                        m['F_BF_End'][:, idx[-1]] += endLoad(Fl, Moment, sinPhi2, cosPhi2, sinBeta2, cosBeta2)
                else :
                    pass  # no loads because both end nodes are below seabed
                # --- no inertia loads from water ballast modeled on ends
                # --- external buoyancy loads: ends ---
                if not e.data['Pot']:
                    pos1    = umesh.TranslationDisp[idx[0 ], :] + umesh.Position[idx[0 ], :]
                    pos2    = umesh.TranslationDisp[idx[-1], :] + umesh.Position[idx[-1], :]
                    pos1[2] -= p['MSL2SWL']
                    pos2[2] -= p['MSL2SWL']
                    z1 = pos1[2]
                    z2 = pos2[2]
                    if mem['i_floor'] == -1:  # both ends above or at seabed
                        if z2<= 0.0:
                            # Compute loads on both ends
                            Fl      = -p['WtrDens'] * g * np.pi *     mem['RMG'][0]**2*z1
                            Moment  = -p['WtrDens'] * g * np.pi *0.25*mem['RMG'][0]**4*sinPhi
                            m['F_B_End'][:, idx[0]] += endLoad(Fl, Moment, sinPhi1, cosPhi1, sinBeta1, cosBeta1) 
                            Fl      = p['WtrDens'] * g * np.pi *     mem['RMG'][-1]**2*z2
                            Moment  = p['WtrDens'] * g * np.pi *0.25*mem['RMG'][-1]**4*sinPhi
                            m['F_B_End'][:, idx[-1]] += endLoad(Fl, Moment, sinPhi2, cosPhi2, sinBeta2, cosBeta2)
                        elif z1< 0.0:
                            # Compute loads only on lower end
                            Fl      = -p['WtrDens'] * g * np.pi *     mem['RMG'][0]**2*z1
                            Moment  = -p['WtrDens'] * g * np.pi *0.25*mem['RMG'][0]**4*sinPhi
                            m['F_B_End'][:, idx[0]] += endLoad(Fl, Moment, sinPhi1, cosPhi1, sinBeta1, cosBeta1)
                        else:
                            pass # Entire member is above the still water line
                    elif mem['doEndBuoyancy'] and z2<= 0.0: # The member crosses the seabed line so only the upper end could have bouyancy effects, if at or below free surface
                        # Only compute the buoyancy contribution from the upper end
                        Fl      = p['WtrDens'] * g * np.pi *     mem['RMG'][-1]**2*z2
                        Moment  = p['WtrDens'] * g * np.pi *0.25*mem['RMG'][-1]**4*sinPhi
                        m['F_B_End'][:, idx[-1]] += endLoad(Fl, Moment, sinPhi2, cosPhi2, sinBeta2, cosBeta2)
                    else:
                        print('>>> Skipping')
                        pass # entire member is buried below the seabed         
                # end if  PropPot
            # --- end do loop through members

        if bEnd:
            # --- Hydrodynamic drag loads: joints
            # NOTE:  All wave kinematics have already been zeroed out above the SWL or instantaneous wave height (for WaveStMod > 0), so loads derived from the kinematics will be correct
            #        without the use of a nodeInWater value, but other loads need to be multiplied by nodeInWater to zero them out above the SWL or instantaneous wave height.
            NJoints = len(graph.Nodes)
            for j in range(NJoints):
                # Obtain the node index because WaveVel, WaveAcc, and WaveDynP are defined in the node indexing scheme, not the markers
                # Compute the dot product of the relative velocity vector with the directional Area of the Joint
                vmag =  m['nodeInWater'][j] * ( m['vrel'][0,j]*p['An_End'][0,j] + m['vrel'][1,j]*p['An_End'][1,j] + m['vrel'][2,j]*p['An_End'][2,j] )
                #NOTE: The PropPot values are only for members, and when the p['AM_End, p['DP_Const_End, p['Mass_MG_End, and p['I_MG_End are computed at init,
                #      contributions to these values are added only if the member connecting to the joint is NOT modeled with potential flow theory
                #      However, the p['An_End term used data from ALL members attached to a node, regardless of the PropPot setting.
                # Lumped added mass loads
                qdotdot_t = umesh.TranslationAcc[j,:]
                qdotdot_r = umesh.RotationAcc   [j,:]
                m['F_A_End'][:,j] = m['nodeInWater'][j] * p['AM_End'][:,:,j].dot(- qdotdot_t) 
                m['F_I_End'][:,j] = p['DP_Const_End'][:,j] * m['FDynP'][j] + p['AM_End'][:,:,j].dot(m['FA'][:,j])
                # Marine growth inertia: ends: Section 4.2.2  
                m['F_IMG_End'][:3,j] = -m['nodeInWater'][j] * p['Mass_MG_End'][j]*qdotdot_t
                m['F_IMG_End'][3:,j] = -m['nodeInWater'][j] * ( p['I_MG_End'][:,:,j].dot(qdotdot_r) - np.cross(umesh.RotationVel[j,:], p['I_MG_End'][:,:,j].dot(umesh.RotationVel[j,:]) ))

                F_end = np.zeros(6)
                for i in range(6):
                    # We are now combining the dynamic pressure term into the inertia term
                    if i < 3:
                        m['F_D_End'][i,j] =  p['An_End'][i,j]*p['DragConst_End'][j]*abs(vmag)*vmag  # Note: vmag is zero if node is not in the water
                        F_end[i] = m['F_D_End'][i,j] + m['F_I_End'][i,j] + p['F_WMG_End'][i,j] + m['F_B_End'][i,j] + m['F_BF_End'][i,j] + m['F_A_End'][i,j] + m['F_IMG_End'][i,j]
                    else:
                        F_end[i] = m['F_B_End'][i,j] + m['F_BF_End'][i,j]  + m['F_IMG_End'][i,j]
                ymesh.Force [j, :] += F_end[:3]
                ymesh.Moment[j, :] += F_end[3:]
                #if j==3:
                myprint('---------------------------JOINTS LOADS ',j+1)
                #    print('vmag',vmag)
                #    print('An_End',p['An_End'][:,j])
                #    print('DragConst', p['DragConst_End'][j])
                #    print('F_D_End'  ,m['F_D_End']  [:,j])
                #    print('F_I_End'  ,m['F_I_End']  [:,j])
                #    print('F_WMG_End',p['F_WMG_End'][:,j])
                #    print('F_B_End'  ,m['F_B_End']  [:,j])
                #    print('F_BF_End' ,m['F_BF_End'] [:,j])
                #    print('F_A_End'  ,m['F_A_End']  [:,j])
                #    print('F_IMG_End',m['F_IMG_End'][:,j])

                if (bVerbose and sum(abs(F_end))>1e-6):
                    print('F_end {:16.4f}{:16.4f}{:16.4f}'.format(*F_end[:3]))
                    print('M_end {:16.4f}{:16.4f}{:16.4f}'.format(*F_end[3:]))
        if bVerbose:
#             print('SUPER HACK')
#             print('SUPER HACK')
#             ymesh.Force *=0
#             ymesh.Moment *=0
#             ymesh.Force[2,0]=1
            print('---------------------------YMESH FORCE')
            for j in range(len(self.NodesBeforeSwap)):
                print('F {:5d} {:16.3f} {:16.3f} {:16.3f}'.format(j+1,*ymesh.Force [j,:]))
            print('---------------------------YMESH MOMENT')
            for j in range(len(self.NodesBeforeSwap)):
                print('M {:5d} {:16.3f} {:16.3f} {:16.3f}'.format(j+1,*ymesh.Moment [j,:]))
        # --- Write Outputs TODO
#          # OutSwtch determines whether or not to actually output results via the WriteOutput array
#          # 1 = Morison will generate an output file of its own.  2 = the caller will handle the outputs, but
#          # Morison needs to provide them.  3 = Both 1 and 2, 0 = No one needs the Morison outputs provided
#          # via the WriteOutput array.
#       IF ( p['OutSwtch > 0 ) THEN
#             # Map calculated results into the AllOuts Array
#          CALL MrsnOut_MapOutputs(Time, y, p, u, m, AllOuts, errStat, errMsg)
#             # Put the output data in the WriteOutput array
#          DO I = 1,p['NumOuts
#             y%WriteOutput[i] = p['OutParam[i]%SignM * AllOuts( p['OutParam[i]%Indx )
#          END DO
#             # Generate output into the output file
#          IF ( p['OutSwtch == 1 .OR. p['OutSwtch == 3 ) THEN
#             CALL MrsnOut_WriteOutputs( p['UnOutFile, Time, y, p, errStat, errMsg )         
#          END IF
#       END IF
        return y


    # --------------------------------------------------------------------------------}
    # --- Useful properties
    # --------------------------------------------------------------------------------{
    @property
    def VolumeStructure(self): return np.sum([e.MorisonData['Vinner'] for e in self.graph.Elements])

    @property
    def VolumeSubmerged(self): return np.sum([e.MorisonData['Vsubmerged'] for e in self.graph.Elements])

    @property
    def VolumeMG(self): return np.sum([e.MorisonData['Vouter']-e.MorisonData['Vinner'] for e in self.graph.Elements])

    @property
    def VolumeBallast(self): return np.sum([e.MorisonData['Vballast'] for e in self.graph.Elements])

    @property
    def MassMG(self): return np.sum([np.sum(e.MorisonData['m_mg_l'])+np.sum(e.MorisonData['m_mg_l']) for e in self.graph.Elements])

    # --------------------------------------------------------------------------------}
    # --- IO/Converters
    # --------------------------------------------------------------------------------{
    def writeSummary(self, filename=None, fid=None):
        morisonToSum(self, filename, fid)

    def toYAMSData(self):
        """ 
        Convert to Data needed by YAMS
        """
        p=dict()
        return p

def LumpDistrHydroLoads(f_hydro, k_hat, dl, h_c):
    """ """
    lumpedLoad = np.zeros(6)
    lumpedLoad[:3] = f_hydro*dl
    lumpedLoad[3:] = np.cross(k_hat*h_c, f_hydro)*dl
    return lumpedLoad

def endLoad(Fl, M, sinPhi, cosPhi, sinBeta, cosBeta):
    """ Takes loads on end node i and converts to 6DOF loads
    INPUTS:
    - Fl        : (N)   axial load about node i
    - M         : (N-m) radial moment about node i, positive in direction of tilt angle
    - sinPhi    : trig functions of  tilt angle 
    - cosPhi   
    - sinBeta   : trig functions of heading of tilt
    - cosBeta  
    OUTPUT:
    - Fi(6)     : (N, Nm) force/moment vector for end node i
    """
    Fi=np.zeros(6)
    Fi[0] =  + Fl*sinPhi*cosBeta
    Fi[1] =  + Fl*sinPhi*sinBeta
    Fi[2] =  + Fl*cosPhi
    Fi[3] =  - M*sinBeta
    Fi[4] =  + M*cosBeta
    return Fi


def morisonToSum(mor, filename=None, fid=None, more=False):
    """ 
    Write a summary file, similar to HydroDyn
    """
    graph  = mor.graph
    MSL2SWL = mor.p['MSL2SWL']

    ExtBuoyancy   = 0.0
    totalFillMass = 0.0
    totalDisplVol = 0.0
    totalVol      = 0.0
    totalMGVol    = 0.0
    totalFillVol  = 0.0
    totalMGMass   = 0.0
    COB           = 0.0
    NNodes = len(mor.NodesBeforeSwap)
    F_B           = np.zeros((6,NNodes))
    F_BF          = np.zeros((6,NNodes))
    F_WMG         = np.zeros((6,NNodes))
    for e in graph.Elements:
         m = e.MorisonData
         totalVol      += m['Vouter']
         totalMGVol    += m['Vouter'] - m['Vinner']
         totalDisplVol += m['Vsubmerged']
         totalFillVol  += m['Vballast']
         totalMGMass += np.sum(m['m_mg_l'])
         totalMGMass += np.sum(m['m_mg_u'])
         # TODO TODO
         idx      = m['nodeIDs']
         memLoads = m['Loads']
         for i in range(m['NElements']+1):
            F_B  [:,idx[i]] += memLoads['F_B']  [:,i]
            F_BF [:,idx[i]] += memLoads['F_BF'] [:,i]
            F_WMG[:,idx[i]] += memLoads['F_WMG'][:,i]

    # --- Helper functions
    s=''
    s+='Strip-Theory Volume Calculations(m^3)\n'
    s+='-------------------------------------\n'
    s+='  Structure Volume     :{:15.5e}\n'.format(totalVol).replace('e+','E+').replace('e-','E-')
    s+='  Submerged Volume     :{:15.5e}\n'.format(totalDisplVol).replace('e+','E+').replace('e-','E-')
    s+='  Marine Growth Volume :{:15.5e}\n'.format(totalMGVol  ).replace('e+','E+').replace('e-','E-')
    s+='  Ballasted Volume     :{:15.5e}\n'.format(totalFillVol).replace('e+','E+').replace('e-','E-')
    s+='              NOTE: Structure, Submerged and Marine Growth volumes are based on members not modelled with WAMIT\n'
    s+='                      Ballasted volume is computed from all members which are marked as filled in the HydroDyn input file, regardless of PropPot flag\n'
    s+='\n'
    s+='\n'

    # Attach the external distributed buoyancy loads to the distributed mesh so they can be transferred to the WRP
    # Because of wave stretching and user-supplied waves, we may have loads above the still water line (SWL) which will be used
    # in the hydrodynamics for conditions where the wave height is > SWL.  So we now need to check that the vertical position
    # is <= SWL for this summary file calculation.
    ymesh = copy.deepcopy(mor._y['Mesh'])
    numJoints = len(mor.graph.Nodes)
    m=mor.m
    ptLoad = np.zeros(6)
    for j in range(ymesh.nNodes):
        if (ymesh.Position[j,2] <= MSL2SWL): # need to check relative to MSL2SWL offset because the Mesh Positons are relative to MSL
            if j < numJoints:
                ptLoad = F_B[:,j] + m['F_B_End'][:,j]
            else:
                ptLoad = F_B[:,j]
            ymesh.Force [j,:]  = ptLoad[:3]
            ymesh.Moment[j,:]  = ptLoad[3:]
        else:
            ymesh.Force [j,:]  = 0
            ymesh.Moment[j,:]  = 0
    ExtBuoyancy = np.zeros(6)
    ExtBuoyancy[:3], ExtBuoyancy[3:] = ymesh.mapLoadsToPoint((0,0,0))

    for j in range(ymesh.nNodes):
        if j < numJoints:
            ptLoad = F_BF[:,j] + m['F_BF_End'][:,j]
        else:
            ptLoad = F_BF[:,j]
        ymesh.Force [j,:]  = ptLoad[:3]
        ymesh.Moment[j,:]  = ptLoad[3:]
    IntBuoyancy = np.zeros(6)
    IntBuoyancy[:3], IntBuoyancy[3:] = ymesh.mapLoadsToPoint((0,0,0))
    TotBuoyancy = ExtBuoyancy + IntBuoyancy
    #print('ExtBuoyancy',np.around(ExtBuoyancy,4))

    s+='\n'
    s+='Total Buoyancy loads summed about ( 0.0, 0.0, 0.0 )\n'
    s+='---------------------------------------------------\n'
    s+='                                BuoyFxi               BuoyFyi               BuoyFzi               BuoyMxi               BuoyMyi               BuoyMzi \n'
    s+='                                  (N)                   (N)                   (N)                  (N-m)                 (N-m)                 (N-m)  \n'
    s+=' External:        {:22.6e}{:22.6e}{:22.6e}{:22.6e}{:22.6e}{:22.6e}\n'.format(*ExtBuoyancy).replace('e+','E+').replace('e-','E-')
    s+=' Internal:        {:22.6e}{:22.6e}{:22.6e}{:22.6e}{:22.6e}{:22.6e}\n'.format(*IntBuoyancy).replace('e+','E+').replace('e-','E-')
    s+=' Total   :        {:22.6e}{:22.6e}{:22.6e}{:22.6e}{:22.6e}{:22.6e}\n'.format(*TotBuoyancy).replace('e+','E+').replace('e-','E-')
    s+='              NOTE: External buoyancy is based on members not modelled with WAMIT\n'
    s+='                      Internal buoyancy is computed from all members which are marked as filled in the HydroDyn input file, regardless of PropPot flag\n'
    s+='                    Total buoyancy does not include WAMIT-modelled buoyancy contribution\n'
    s+='\n'
    s+='\n'
#       !   ! Now compute marine growth weight at the WRP
#       DO J = 1, yMesh%Nnodes
#          if (J <= numJoints) then
#             yMesh%Force(:,J)   = F_WMG(1:3,J) + p%F_WMG_End(:,J)
#          else
#             yMesh%Force(:,J)   = F_WMG(1:3,J)
#          end if
#          yMesh%Moment(:,J)  = F_WMG(4:6,J)
#       END DO ! DO J
#       MG_Wt = 0.0
#       CALL Transfer_Point_to_Point( yMesh, WRP_Mesh, M_P_2_P, errStat, errMsg, uMesh, WRP_Mesh_position )
#       MG_Wt(1:3) = WRP_Mesh%Force(:,1)
#       MG_Wt(4:6) = WRP_Mesh%Moment(:,1)
#        CALL MeshMapDestroy( M_P_2_P, errStat, errMsg ); IF ( errStat /= ErrID_None ) CALL WrScr(TRIM(errMsg))
    s+='\n'
    s+='Weight loads about ( 0.0, 0.0, 0.0 )\n'
    s+='------------------------------------\n'
    s+='                                 MGFxi                 MGFyi                 MGFzi                 MGMxi                 MGMyi                 MGMzi  \n'
    s+='                                  (N)                   (N)                   (N)                  (N-m)                 (N-m)                 (N-m)  \n'
    s+=' Marine Growth:   {:22.6e}{:22.6e}{:22.6e}{:22.6e}{:22.6e}{:22.6e}\n'.format(0,0,0,0,0,0).replace('e+','E+').replace('e-','E-')
    s+='\n'
    s+='\n'
    s+='\n'
    s+='Nodes (first [{:4d} ] are joints, remainder are internal nodes)\n'.format(len(graph.Nodes))
    s+='\n'
    s+='\n'
    s+='   i      MbrIndx      Nxi         Nyi         Nzi           R          t          tMG        MGDens     PropPot    FilledFlag  FilledMass      Cd          Ca          Cp         AxCd        AxCa        AxCp        JAxCd       JAxCa       JAxCp  \n'
    s+='  (-)       (-)        (m)         (m)         (m)          (m)        (m)         (m)       (kg/m^3)      (-)         (-)        (kg)          (-)         (-)         (-)         (-)         (-)         (-)         (-)         (-)         (-)   \n'
    # Write the node data
    for i,n in enumerate(graph.Nodes):
        # need to add MSL2SWL offset from this because the Positons are relative to SWL, but we should report them relative to MSL here
        d = n.data
        pos = n.point
        pos[2] += MSL2SWL
        s0= ' {:5d}      -       {:10.4f}  {:10.4f}  {:10.4f}      -           -       {:10.3e}  {:10.3e}'.format(i+1,pos[0],pos[1],pos[2],0, 0)# TODO TODO nodes(i)%tMG,  nodes(i)%MGdensity
        s0+= '      -           -           -           -           -           -           -           -           -      '
        s0+= ' {:10.3e}  {:10.3e}  {:10.3e}\n'.format(d['JAxCd'], d['JAxCa'], d['JAxCp'])
        s0 = s0.replace('e+','E+').replace('e-','E-')
        s+=s0
    c = len(graph.Nodes)
    for j,e in enumerate(graph.Elements):
        m = e.MorisonData
        for i in np.arange(1,m['NElements']):
            c += 1
            fillFlag = m['l_fill'] - m['dl']*i > 0.0
            # need to add MSL2SWL offset from this because the Positons are relative to SWL, but we should report them relative to MSL here
            if e.flipped:
                ii = m['NElements']-i
            else:
                ii=i
#             ii=i
            pos = m['SubNodesPositions'][ii-1,:]
            pos[2] += MSL2SWL
            s0= ' {:5d}  {:10d}  {:10.4f}  {:10.4f}  {:10.4f}  {:10.3e}  {:10.3e}  {:10.3e}  {:10.3e}'.format(c,e.ID, pos[0], pos[1], pos[2], m['R'][ii], m['R'][ii]-m['Rin'][ii], m['tMG'][ii], m['MGdensity'][ii])
            spropot='T' if e.data['Pot'] else 'F'
            sfill  ='T' if fillFlag      else 'F'
            s0+= '           {:s}           {:s}'.format(spropot, sfill)
            s0+= '  {:10.3e}  {:10.3e}  {:10.3e}  {:10.3e}  {:10.3e}  {:10.3e}  {:10.3e}'.format(m['m_fb_u'][ii]+m['m_fb_l'][ii],  m['Cd'][ii], m['Ca'][ii], m['Cp'][ii], m['AxCd'][ii], m['AxCa'][ii], m['AxCp'][ii])
            s0+= '         -           -           -  \n'
            s0 = s0.replace('e+','E+').replace('e-','E-')
            s+=s0
    s+='\n'
    s+='\n'
    s+='\n'
    s+=' Members\n'
    s+='\n'
    s+='\n'
    s+=' MemberID  joint1  joint2      Length       NElem         Volume       MGVolume          R1           t1             R2           t2          PropPot      FilledFlag   FillDensity    FillFSLoc     FillMass         Cd1            Ca1          Cp1          AxCd1         AxCa1         AxCp1        JAxCd1        JAxCa1       JAxCp1           Cd2           Ca2           Cp2          AxCd2         AxCa2         AxCp2        JAxCd2        JAxCa2        JAxCp2   \n'
    s+='   (-)      (-)     (-)         (m)          (-)          (m^3)         (m^3)            (m)          (m)            (m)          (m)           (-)           (-)        (kg/m^3)         (-)          (kg)           (-)            (-)          (-)           (-)           (-)           (-)           (-)           (-)          (-)            (-)           (-)           (-)           (-)           (-)           (-)           (-)           (-)           (-)    \n'
    for i, e in enumerate(graph.Elements):
        d = e.data
        m = e.MorisonData
        n1 = e.nodes[0]
        n2 = e.nodes[1]
        if d['Pot']:
            MGvolume  = 0
            memberVol = 0
        else:
            memberVol   = m['Vouter']
            MGvolume    = m['Vouter'] - m['Vinner']
        if m['l_fill'] > 0.0:
            filledFlag = True
            mass_fill   = m['FillDens']*m['Vballast']
        else:
            filledFlag = False
            mass_fill  = 0.0
        spropot='T' if e.data['Pot'] else 'F'
        sfill  ='T' if filledFlag    else 'F'
        #print('----------------------------------------------------------------')
        #print('>>> Member ',i+1)
        #print('IDs',m['nodeIDs'])
        #print('RMG',m['RMG'])
        #print('Rin',m['Rin'])

        s0=  ' {:8d}  {:6d}  {:6d}  '.format(e.ID, e.nodeIDs[0]+1, e.nodeIDs[1]+1)
        s0+= '{:12.5e}  {:12d} '.format(m['refLength'], m['NElements'])
        s0+= ' {:12.5e}  {:12.5e}  {:12.5e}  {:12.5e}  {:12.5e}  {:12.5e}'.format(memberVol, MGvolume, m['RMG'][0], m['RMG'][0]-m['Rin'][0], m['RMG'][-1], m['RMG'][-1]-m['Rin'][-1]) 
        s0+= '             {:s}             {:s}'.format(spropot, sfill)
        s0+= '  {:12.5e}  {:12.5e}  {:12.5e}'.format(m['FillDens'], m['FillFSLoc'], mass_fill)
        s0+= '  {:12.5e}  {:12.5e}  {:12.5e}  {:12.5e}  {:12.5e}  {:12.5e}  {:12.5e}  {:12.5e}  {:12.5e}'.format(m['Cd'][0], m['Ca'][0], m['Cp'][0], m['AxCd'][0], m['AxCa'][0], m['AxCp'][0], n1.data['JAxCd'], n1.data['JAxCa'], n1.data['JAxCp'])
        s0+= '  {:12.5e}  {:12.5e}  {:12.5e}  {:12.5e}  {:12.5e}  {:12.5e}  {:12.5e}  {:12.5e}  {:12.5e}'.format(m['Cd'][-1], m['Ca'][-1], m['Cp'][-1], m['AxCd'][-1], m['AxCa'][-1], m['AxCp'][-1], n2.data['JAxCd'], n2.data['JAxCa'], n2.data['JAxCp'])
        s0+= '\n'
        s0 = s0.replace('e+','E+').replace('e-','E-')
        s+=s0
    s+='\n'
    s+='\n'
    s+='\n'
    s+='Requested Member Outputs\n'
    s+='\n'
    s+='\n'
    s+='   Label         Xi          Yi          Zi       MemberID    StartXi     StartYi     StartZi      EndXi       EndYi       EndZi        Loc    \n'
    s+='    (-)          (m)         (m)         (m)        (-)         (m)         (m)         (m)         (m)         (m)         (m)         (-)    \n'
    # TODO
    s+='\n'
    s+='\n'
    s+='\n'
    s+=' Requested Joint Outputs\n'
    s+='\n'
    s+='\n'
    s+='   Label         Xi          Yi          Zi      InpJointID\n'
    s+='    (-)          (m)         (m)         (m)        (-)    \n'
    # TODO
    s+='\n'

    if fid is None:
        with open(filename, 'w') as fid:
            fid.write(s);
    else:
        fid.write(s);



def getOrientationAngles(p1, p2):
    """ 
    p1: Point 1
    p2: Point 2
    """
    # calculate instantaneous incline angle and heading, and related trig values
    # the first and last NodeIndx values point to the corresponding Joint nodes idices which are at the start of the Mesh
    vec      = np.asarray(p2) - np.asarray(p1)
    vecLen   = np.sqrt(vec.dot(vec))
    vecLen2D = np.sqrt(vec[0]**2+vec[1]**2)
    if vecLen < 0.000001:
        raise Exception('An element of the Morison structure has co-located endpoints!  This should never occur.  Please review your model.')
    k_hat = vec / vecLen 
    phi   = np.arctan2(vecLen2D, vec[2])  # incline angle   
    if phi==0:
        beta = 0
    else:
        beta = np.arctan2(vec[1], vec[0])# heading of incline     
    sinPhi  = np.sin(phi)
    cosPhi  = np.cos(phi)  
    tanPhi  = np.tan(phi)     
    sinBeta = np.sin(beta)
    cosBeta = np.cos(beta)
    return phi, sinPhi, cosPhi, tanPhi, sinBeta, cosBeta, k_hat


def Morison_DirCosMtrx(pos0, pos1):
    """  Compute the direction cosine matrix given two points along the axis of a cylinder
    TODO TODO Merge with DCM of FEM.
    Main axis is x????
    """
    x0 = pos0[0]
    y0 = pos0[1]
    z0 = pos0[2]
    x1 = pos1[0]
    y1 = pos1[1]
    z1 = pos1[2]
    xz  = np.sqrt((x0-x1)*(x0-x1)+(z0-z1)*(z0-z1))
    xyz = np.sqrt((x0-x1)*(x0-x1)+(y0-y1)*(y0-y1)+(z0-z1)*(z0-z1))

    if xz==0:
        if y1<y0:
            DirCos = np.array([[1,0,0], [0,0,-1], [0,1,0]])
        else:
            DirCos = np.array([[1,0,0], [0,0, 1], [0,-1,0]])
    else:
        DirCos = np.zeros((3,3))
        DirCos[0, 0] = (z1-z0)/xz
        DirCos[0, 1] = -(x1-x0)*(y1-y0)/(xz*xyz)
        DirCos[0, 2] = (x1-x0)/xyz
        DirCos[1, 0] = 0.0
        DirCos[1, 1] = xz/xyz
        DirCos[1, 2] = (y1-y0)/xyz
        DirCos[2, 0] = -(x1-x0)/xz
        DirCos[2, 1] = -(y1-y0)*(z1-z0)/(xz*xyz)
        DirCos[2, 2] = (z1-z0)/xyz
    return DirCos

def DistributeElementLoads(Fl, Fr, M, sinPhi, cosPhi, sinBeta, cosBeta, alpha):
    """
    Takes loads on node i in element tilted frame and converts to 6DOF loads at node i and adjacent node
    INPUTS:
    - Fl     : (N)   axial load about node i
    - Fr     : (N)   radial load about node i in direction of tilt
    - M      : (N-m) radial moment about node i, positive in direction of tilt angle
    - sinPhi : trig functions of  tilt angle 
    - cosPhi :
    - sinBeta: trig functions of heading of tilt
    - cosBeta:
    - alpha  : fraction of load staying with node i (1-alpha goes to other node)  
    OUTPUS:
    - F1(6): (N, Nm) force/moment vector for node i
    - F2(6): (N, Nm) force/moment vector for the other node (whether i+1, or i-1)
    """
    F1=np.zeros(6)
    F1[0] =  cosBeta*(Fl*sinPhi + Fr*cosPhi)*alpha
    F1[1] =  sinBeta*(Fl*sinPhi + Fr*cosPhi)*alpha
    F1[2] =         (Fl*cosPhi - Fr*sinPhi)*alpha
    F1[3] =  -sinBeta * M                    *alpha
    F1[4] =  cosBeta * M                    *alpha
    F1[5] =  0.0
       
    F2=np.zeros(6)
    F2[0] =  cosBeta*(Fl*sinPhi + Fr*cosPhi)*(1-alpha)
    F2[1] =  sinBeta*(Fl*sinPhi + Fr*cosPhi)*(1-alpha)
    F2[2] =          (Fl*cosPhi - Fr*sinPhi)*(1-alpha)
    F2[3] =  -sinBeta * M                    *(1-alpha)
    F2[4] =  cosBeta * M                    *(1-alpha)
    F2[5] = 0.0
    return F1, F2



if __name__ == '__main__':
    import sys
    if len(sys.argv)>=1:
        filename=sys.argv[1]
    else:
        filename='_SparNoRNA_HD_RefH.dat'
