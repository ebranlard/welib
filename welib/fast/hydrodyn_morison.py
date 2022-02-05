import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# Local 
from welib.tools.clean_exceptions import *
from welib.weio import FASTInputFile


class Morison:

    def __init__(self, graph, File, WtrDpth, MSL2SWL):
        """ 
        graph: graph generated from hydrodyn file
        File : File content of hydrodyn input file
        """
        import copy
        self.graph        = copy.deepcopy(graph)
        self.File         = File

        # Internal
        self.p={}
        self.m={}
        self.p['WtrDpth'] = WtrDpth
        self.p['MSL2SWL'] = MSL2SWL

        # ---
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
    def init(self, initData, gravity = 9.81, WtrDens=1025):
        """
        Initialize HydroDyn model 

        gravity: position of transition point
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
   
        # ---- Input Mesth
        # u%Mesh 
        # Here positions are relative to MSL not SWL
        #    pos(3) = pos(3) + InitInp%MSL2SWL
        # y%Mesh is duplicate
        # Define initial system states here:
        #m%LastIndWave              = 1
        
        # --- loop through joints to calculate joint quantities (the joints are the first NJoints nodes)
        for n in graph.Nodes:
            An        = 0.0
            Vn        = 0.0
            I_n       = 0.0
            MGdens    = 0.0
            tMG       = -999.0
            An_drag   = 0.0
            #       IF ( InitInp%InpJoints(i)%Position(3) >= -p%WtrDpth ) THEN
            #          ! loop through each member attached to the joint, getting the radius of its appropriate end
            #          DO J = 1, InitInp%InpJoints(I)%NConnections
            #             ! identify attached member and which end to us
            #             IF (InitInp%InpJoints(I)%ConnectionList(J) > 0) THEN         ! set up for end node 1
            #                !TODO: Should not perform a copy here?  A pointer to data would be better?
            #                member = p%Members(InitInp%InpJoints(I)%ConnectionList(J))   
            #                MemberEndIndx = 1
            #             ELSE     
            #                ! set up for end node N+1
            #                ! NOTE:  %ConnectionList(J) is negative valued if InitInp%Morison%InpMembers(I)%MJointID2 == InitInp%Morison%InpJoints(J)%JointID.  See HydroDynInput_ProcessInitData, members section
            #                member = p%Members(-InitInp%InpJoints(I)%ConnectionList(J))
            #                MemberEndIndx = member%NElements + 1
            #             END IF
            #             ! Compute the signed area*outward facing normal of this member
            #             sgn = 1.0
            #             IF ( MemberEndIndx == 1 ) THEN
            #                sgn = -1.0                                ! Local coord sys points into member at starting node, so flip sign of local z vector
            #             ELSE
            #                sgn = 1.0                                 ! Local coord sys points out of member at ending node, so leave sign of local z vector
            #             END IF
            #             ! Account for reordering of what the original node for the end was -- This affects the sign of the An term which can pose a problem for members crossing the waterline
            #             if (member%Flipped)   sgn = -1.0 * sgn
            #             ! Compute the signed quantities for this member end (for drag regardless of PropPot value), and add them to the joint values
            #             An_drag = An_drag + sgn* member%k*Pi*(member%RMG(MemberEndIndx))**2     ! area-weighted normal vector
            #             ! For the following quantities, the attached member cannot be modeled using WAMIT if we're to count it
            #             IF  (.NOT. member%PropPot) THEN
            #                ! Compute the signed quantities for this member end, and add them to the joint values
            #                An = An + sgn* member%k*Pi*(member%RMG(MemberEndIndx))**2     ! area-weighted normal vector
            #                Vn = Vn + sgn* member%k*   (member%RMG(MemberEndIndx))**3     ! r^3-weighted normal vector used for mass
            #                I_n=I_n + sgn* member%k*Pi*(member%RMG(MemberEndIndx))**4     ! r^4-weighted normal vector used for moments of inertia
            #                if (tMG == -999.0) then
            #                   ! All member nodes at this joint will have the same MG thickness and density, so only do this once
            #                   tMG = member%tMG(MemberEndIndx)
            #                   MGdens = member%MGdensity(MemberEndIndx) 
            #                end if
            #             END IF
            #          END DO   !J = 1, InitInp%InpJoints(I)%NConnections
            #          p%An_End(:,i) = An_drag 
            #          Amag_drag = Dot_Product(An_drag ,An_drag)
            #          Amag = Dot_Product(An ,An)
            #          IF (EqualRealNos(Amag_drag, 0.0_ReKi)) THEN
            #             p%DragConst_End(i) =  0.0
            #          ELSE
            #             p%DragConst_End(i) = InitInp%Nodes(i)%JAxCd*p%WtrDens / ( 4.0_ReKi * Amag_drag )
            #          END IF
            #          ! magnitudes of normal-weighted values
            #          Amag = sqrt(Amag)
            #          Vmag = sqrt(Dot_Product(Vn ,Vn))
            #          Imag = sqrt(Dot_Product(I_n,I_n))
            #          ! Constant part of the external hydrodynamic added mass term
            #          if ( Vmag > 0.0 ) then
            #             v2D(:,1) = Vn        
            #             p%AM_End(:,:,i) = (InitInp%Nodes(I)%JAxCa*InitInp%WtrDens/ Vmag)*matmul(transpose(v2D), v2D) 
            #          end if
            #          ! Constant part of the external hydrodynamic dynamic pressure force
            #          if ( Amag > 0.0 ) then
            #             p%DP_Const_End(:,i) = -InitInp%Nodes(i)%JAxCp*An 
            #          endif
            #          ! marine growth mass/inertia magnitudes
            #          p%Mass_MG_End(i) = MGdens * tMG * Amag
            #          p%F_WMG_End(3,i) =        -MGdens * tMG * Amag * InitInp%Gravity  ! Z component of the directional force due to marine growth mass at joint
            #          Ir_MG_end   =  0.25 * MGdens * tMG * Imag  ! radial moment of inertia magnitude
            #          Il_MG_end   =  0.5  * MGdens * tMG * Imag  ! axial moment of inertia magnitude
            #          ! get rotation matrix for moment of inertia orientations
            #          call RodrigMat(I_n, R_I, errStat, errMsg)
            #          IF ( errStat >= AbortErrLev ) RETURN
            #          ! globally-oreinted moment of inertia matrix for joint
            #          Irl_mat = 0.0
            #          Irl_mat(1,1) = Ir_MG_end
            #          Irl_mat(2,2) = Ir_MG_end
            #          Irl_mat(3,3) = Il_MG_end
            #          p%I_MG_End(:,:,i) = MatMul( MatMul(R_I, Irl_mat), Transpose(R_I) ) ! final moment of inertia matrix for node
            #       END IF  ! InitInp%InpJoints(i)%Position(3) >= -p%WtrDpth
            #          ! Initialize the outputs      
            #    IF ( p%OutSwtch > 0) then  !@mhall: moved this "if" to after allocations
            #       CALL MrsnOUT_Init( InitInp, y, p, InitOut, errStat, errMsg )
            #          ! Determine if we need to perform output file handling
            #       IF ( p%OutSwtch == 1 .OR. p%OutSwtch == 3 ) THEN  
            #          CALL MrsnOUT_OpenOutput( Morison_ProgDesc%Name, TRIM(InitInp%OutRootName)//'.HD', p, InitOut, errStat, errMsg )
            #       END IF
            #    END IF  
            #    ! We will call CalcOutput to compute the loads for the initial reference position
            #    ! Then we can use the computed load components in the Summary File
            #    ! NOTE: Morison module has no states, otherwise we could no do this. GJH
            #    call Morison_CalcOutput(0.0_DbKi, u, p, x, xd, z, OtherState, y, m, errStat, errMsg )
            #       ! Write Summary information now that everything has been initialized. 
            #    CALL WriteSummaryFile( InitInp%UnSum, InitInp%Gravity, InitInp%MSL2SWL, InitInp%WtrDpth, InitInp%NJoints, InitInp%NNodes, InitInp%Nodes, p%NMembers, p%Members, &
            #                           p%NumOuts, p%OutParam, p%NMOutputs, p%MOutLst,  p%NJOutputs, p%JOutLst, u%Mesh, y%Mesh, &
            #                           p, m, errStat, errMsg )

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
            for n in graph.Nodes:
                n.z -= self.p['MSL2SWL']
            for e in graph.Elements:
                if e.data['Pot'] is False:
                    numDiv = np.ceil(e.length/e.data['DivSize']).astype(int)
                    MorisonData = {}
                    MorisonData['NElements'] = numDiv
                    MorisonData['dl']        = e.length/numDiv
                    MorisonData['refLength'] = e.length
                    n1, n2 = e.nodes
                    SubNodesPositions = np.zeros((numDiv-1,3))
                    for j in range(numDiv-1):
                        s = (j+1)/numDiv
                        SubNodesPositions[j,:] = n1.point * (1-s) + n2.point * s
                    MorisonData['SubNodesPositions'] = SubNodesPositions
                    e.MorisonData=MorisonData
        self._divisionDone=True
        Nodes=self._MorisonPositions()
        return Nodes

    # --------------------------------------------------------------------------------}
    # ---  
    # --------------------------------------------------------------------------------{
    def _setupMembers(self):
        """ see SetupMembers in Morison.f90 
          - swap nodes if needed 
          - compute member properties
        """
        WtrDepth = self.p['WtrDpth'] # TODO which one is it
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
                print('>>> Swapping mID',m.ID, 'NodeIDs:', m.nodeIDs)
                print('TODO TODO TODO SWAPPING IS INCOMPLETE/BUGGY')
                m.swapNodes()
                m.MorisonData['SubNodesPositions'] =  m.MorisonData['SubNodesPositions'][-1::-1,:]
            m.flipped=doSwap
        # Update connectivity
        graph.updateConnectivity()
        # --- Set Member properties
        for m in graph.Elements:
            # see SetMemberProperties in Morison.f90
            member = m.MorisonData
            N  = member['NElements']
            dl = member['dl']
            vec     = m.nodes[1].point-m.nodes[0].point
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
#             do i = 1, member%NElements+1
#                call SetNodeMG( InitInp%NMGDepths, InitInp%MGDepths, InitInp%Nodes(member%NodeIndx(i)), InitInp%MSL2SWL, member%tMG(i), member%MGDensity(i) )
            member['MGdensity']=np.zeros(N) # TODO
            member['tMG'] = np.zeros(N+1)   # TODO
            t             = np.linspace(m.nodes[0].data['t']  , m.nodes[1].data['t']  , N+1)
            member['R']   = np.linspace(m.nodes[0].data['D']/2, m.nodes[1].data['D']/2, N+1)
            member['RMG'] = member['R']+member['tMG']
            member['Rin'] = member['R']-t

            # --- see SetExternalHydroCoefs from Morison.f90
            C = graph.NodePropertySets['SimpleCoefs'][0] # we only have one of these sets
            if m.data['CoefMod'] == 'SimpleCoefs':
                # Simple model : all nodes receive the same coefficients
                member['Cd']=np.zeros(N+1)
                member['Cd']=np.zeros(N+1)
                member['Cd']=np.zeros(N+1)
                member['Cd']=np.zeros(N+1)
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
            if 'FillGroups' in graph.NodePropertySets.keys():
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
            # calculate element-level values
            member['dRdl_mg'] = np.diff(member['RMG'][1:])/dl
            member['dRdl_in'] = np.diff(member['Rin'][1:])/dl

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
            member['m_fb_u']      = np.zeros(N) # per element
            member['m_fb_l']      = np.zeros(N) # per element
            member['h_cfb_u']     = np.zeros(N)
            member['I_lfb_u']     = np.zeros(N)
            member['I_rfb_u']     = np.zeros(N)
            member['I_lmg_u']     = np.zeros(N)
            member['I_lmg_l']     = np.zeros(N)
            member['I_rmg_u']     = np.zeros(N)
            member['I_rmg_l']     = np.zeros(N)
            member['m_mg_u']      = np.zeros(N)
            member['m_mg_l']      = np.zeros(N)
            member['h_cmg_u']     = np.zeros(N)
            member['h_cmg_l']     = np.zeros(N)
            member['floodstatus'] = np.zeros(N)
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
                    if (InitInp%Nodes(member%NodeIndx(N+1))%Position(3) - member%Rin(N+1)*sinPhi) < member%FillFSLoc :
                        raise Exception('The modeling of partially flooded/ballested members requires the the member not be near horizontal.  This is not true for MemberID ', m)
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





def morisonToSum(mor, filename=None, fid=None, more=False):
    """ 
    Write a summary file, similar to HydroDyn
    """
    graph  = mor.graph
    MSL2SWL = mor.File['MSL2SWL']

    ExtBuoyancy   = 0.0
    totalFillMass = 0.0
    totalDisplVol = 0.0
    totalVol      = 0.0
    totalMGVol    = 0.0
    totalFillVol  = 0.0
    totalMGMass   = 0.0
    COB           = 0.0
    F_B           = 0.0
    F_BF          = 0.0
    F_WMG         = 0.0
    for e in graph.Elements:
         m = e.MorisonData
         n1 = e.nodes[0]
         n2 = e.nodes[-1]
         totalVol      += m['Vouter']
         totalMGVol    += m['Vouter'] - m['Vinner']
         totalDisplVol += m['Vsubmerged']
         totalFillVol  += m['Vballast']
         totalMGMass += np.sum(m['m_mg_l'])
         totalMGMass += np.sum(m['m_mg_u'])
     # TODO TODO
#          do i = 1, mem%NElements+1 
#             F_B  (:,mem%NodeIndx(i)) = F_B  (:,mem%NodeIndx(i)) + m%memberLoads(j)%F_B  (:,i) # TODO memberLoads
#             F_BF (:,mem%NodeIndx(i)) = F_BF (:,mem%NodeIndx(i)) + m%memberLoads(j)%F_BF (:,i)
#             F_WMG(:,mem%NodeIndx(i)) = F_WMG(:,mem%NodeIndx(i)) + m%memberLoads(j)%F_WMG(:,i)
      



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
#          ! Attach the external distributed buoyancy loads to the distributed mesh so they can be transferred to the WRP
#          ! Because of wave stretching and user-supplied waves, we may have loads above the still water line (SWL) which will be used
#          ! in the hydrodynamics for conditions where the wave height is > SWL.  So we now need to check that the vertical position
#          ! is <= SWL for this summary file calculation.
#       DO J = 1, yMesh%Nnodes
#          if ( yMesh%Position(3,J) <= MSL2SWL ) then  ! need to check relative to MSL2SWL offset because the Mesh Positons are relative to MSL
#             if (J <= numJoints) then
#                ptLoad = F_B(:,J) + m%F_B_end(:,J)
#             else
#                ptLoad = F_B(:,J)
#             end if
#             yMesh%Force(:,J)   = ptLoad(1:3)
#             yMesh%Moment(:,J)  = ptLoad(4:6)
#          else
#             yMesh%Force(:,J)   = 0.0
#             yMesh%Moment(:,J)  = 0.0
#          end if              ! <= still water line check
#       END DO ! DO J
#          ! Transfer the loads from the distributed mesh to the (0,0,0) point mesh
#       CALL MeshMapCreate           ( yMesh, WRP_Mesh, M_P_2_P, errStat, errMsg                )
#         !CALL CheckError( errStat, 'Message from MeshMapCreate HD_M_L_2_ED_P: '//NewLine//errMsg )
#       CALL Transfer_Point_to_Point( yMesh, WRP_Mesh, M_P_2_P, errStat, errMsg, uMesh, WRP_Mesh_position )
#       ExtBuoyancy(1:3) = WRP_Mesh%Force (:,1)
#       ExtBuoyancy(4:6) = WRP_Mesh%Moment(:,1)
#          ! Now compute internal Buoyancy
#       DO J = 1, yMesh%Nnodes
#          if (J <= numJoints) then
#             ptLoad = F_BF(:,J) + m%F_BF_end(:,J)
#          else
#             ptLoad = F_BF(:,J)
#          end if
#          yMesh%Force(:,J)   = ptLoad(1:3)
#          yMesh%Moment(:,J)  = ptLoad(4:6)
#       END DO ! DO J
#       IntBuoyancy = 0.0
#       CALL Transfer_Point_to_Point( yMesh, WRP_Mesh, M_P_2_P, errStat, errMsg, uMesh, WRP_Mesh_position )
#       IntBuoyancy(1:3) = WRP_Mesh%Force(:,1)
#       IntBuoyancy(4:6) = WRP_Mesh%Moment(:,1)
    s+='\n'
    s+='Total Buoyancy loads summed about ( 0.0, 0.0, 0.0 )\n'
    s+='---------------------------------------------------\n'
    s+='                                BuoyFxi               BuoyFyi               BuoyFzi               BuoyMxi               BuoyMyi               BuoyMzi \n'
    s+='                                  (N)                   (N)                   (N)                  (N-m)                 (N-m)                 (N-m)  \n'
    s+=' External:        {:22.6e}{:22.6e}{:22.6e}{:22.6e}{:22.6e}{:22.6e}\n'.format(0,0,0,0,0,0).replace('e+','E+').replace('e-','E-')
    s+=' Internal:        {:22.6e}{:22.6e}{:22.6e}{:22.6e}{:22.6e}{:22.6e}\n'.format(0,0,0,0,0,0).replace('e+','E+').replace('e-','E-')
    s+=' Total   :        {:22.6e}{:22.6e}{:22.6e}{:22.6e}{:22.6e}{:22.6e}\n'.format(0,0,0,0,0,0).replace('e+','E+').replace('e-','E-')
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
            pos = m['SubNodesPositions'][ii-1, :]
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
        sfill  ='T' if filledFlag      else 'F'

        s0=  ' {:8d}  {:6d}  {:6d}  '.format(e.ID, e.nodeIDs[0], e.nodeIDs[1])
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


if __name__ == '__main__':
    filename='../../data/Monopile/MT100_HD.dat'
    filename='../../data/Monopile/TetraSpar_HydroDyn_v2.dat'
    #filename='../../data/SparNoRNA/SparNoRNA_HD_RefH.dat'
    #sumfile ='../../data/SparNoRNA/Main.HD_python.sum'
    filename='_SparNoRNA_HD_RefH.dat'
    sumfile ='_Main.HD_python.sum'


#     hd = weio.FASTInputFile(filename)
#     hd.write('Out.dat')
#     Graph = hd.toGraph()
    hd = HydroDyn(filename)
    hd.init()
    hd.MorisonPositions
    hd.writeSummary(sumfile)


#     Graph.divideElements(3)
#     print(Graph)
#     import numpy as np
#     import matplotlib.pyplot as plt
#     from matplotlib import collections  as mc
#     from mpl_toolkits.mplot3d import Axes3D
#     fig = plt.figure()
#     ax = fig.add_subplot(1,2,1,projection='3d')
#     lines=Graph.toLines(output='coord')
#     for l in lines:
#     #     ax.add_line(l)
#         ax.plot(l[:,0],l[:,1],l[:,2])

#     ax.autoscale()
    # ax.set_xlim([-40,40])
    # ax.set_ylim([-40,40])
    # ax.set_zlim([-40,40])
    # ax.margins(0.1)
#     plt.show()

# 
# if __name__ == '__main__':
#     pass
