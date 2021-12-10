""" 
Tools to perform direct elimination, similar to SubDyn implementation

"""

import numpy as np


def nDOF_c(self):
    """ Returns number of DOF after constraint reduction (via the matrix T)
    TODO
    """
    nDOF_c = 0
    # Rigid assemblies contribution
    nDOF_c = nDOF_c + 6*size(RA)
    # Contribution from all the other joints
    for iNode, node = in enumerate(self.Nodes):
        #m = Init%NodesConnE(iNode,1) ! Col1: number of elements connected to this joint
        NodeType = Init%Nodes(iNode,iJointType)
        if   NodeType == idJointPin  :
            nDOF_c += 5 + 1*m
            print('Node',iNode, ' is a pin joint, number of members involved: ',m)
        elif NodeType == idJointUniversal :
            nDOF_c = nDOF_c + 4 + 2*m
            print('Node',iNode, ' is an universal joint, number of members involved: ',m)
        elif NodeType == idJointBall :
            nDOF_c = nDOF_c + 3 + 3*m
            print('Node',iNode, ' is a ball joint, number of members involved: ',m)
        elif NodeType == idJointCantilever:
            if ( NodeHasRigidElem(iNode, Init, p, er)):
                # This joint is involved in a rigid link assembly, we skip it (accounted for above)
                print'(4X,A,I5,A,I5)','Node',iNode, ' is involved in a Rigid assembly'
            else:
                # That's a regular Cantilever joint
                nDOF_c = nDOF_c + 6

def rigidLinkAssemblies(self):
    """ Setup a list of rigid link assemblies (RA)
    !! Variables created by this routine:
    !! - RA(ia)= [e1,..,en]  list of elements forming each rigid link assembly "ia".
    !!                       Needed for BuildTMatrix
    !! - RAm1(e)=(RA^-1(e)= a) : for a given element give the index of a rigid assembly. 
    !!                       Needed for BuildTMatrix
    type(IList), dimension(:)   :: RA   !< RA(a) = [e1,..,en]  list of elements forming a rigid link assembly
    integer(IntKi), dimension(:):: RAm1 !< RA^-1(e) = a , for a given element give the index of a rigid assembly
    """
#         allocate(RAm1(1:Init%NElem)) ! NOTE: do not deallocate, this is an "output" of this function
#         RAm1(1:Init%NElem) = -1
# 
#         ! --- Establish a list of rigid link elements
#         Er = RigidLinkElements(Init, p, ErrStat2, ErrMsg2)
#         nRA=0
#         do while (len(Er)>0)
#         nRA=nRA+1
#         ! Creating List Ea of elements of a given assembly
#         call init_list(Ea, 0, 0, ErrStat2, ErrMsg2);
#         e0 = pop(Er, ErrStat2, ErrMsg2);
#         call append(Ea, e0, ErrStat2, ErrMsg2);
#         call AddNeighbors(e0, Er, Ea)
#         if (DEV_VERSION) then
#         call print_list(Ea,'Rigid assembly (loop 1) element list')
#         endif
#         do ie = 1, len(Ea)
#         e0 = get(Ea, ie, ErrStat2, ErrMsg2)
#         RAm1(e0) = nRA ! Index of rigid assembly that this element belongs to
#         enddo
#         call destroy_list(Ea, ErrStat2, ErrMsg2)
#         enddo
#         call destroy_list(Er, ErrStat2, ErrMsg2)
# 
#         ! --- Creating RA, array of lists of assembly elements.
#         ! Note: exactly the same as all the Ea created above, but we didn't know the total number of RA
#         allocate(RA(1:nRA)) ! NOTE: do not deallocate, this is an "output" of this function
#         do ia = 1, nRA
#         call init_list(RA(ia), 0, 0, ErrStat2, ErrMsg2)
#         enddo
#         do ie = 1, Init%NElem
#         ia = RAm1(ie) ! Index of the assembly the element belongs to: RA^{-1}(ie) = ia
#         if (ia>0) then
#         call append(RA(ia), ie, ErrStat2, ErrMsg2)
#         endif
#         enddo
#         if (DEV_VERSION) then
#         do ia = 1, nRA
#         call print_list(RA(ia),'Rigid assembly (loop 2) element list')
#         enddo
#         endif
# CONTAINS
#    !> The neighbor-elements of element e0 (that are found within the list Er) are added to the list Ea  
#    RECURSIVE SUBROUTINE AddNeighbors(e0, Er, Ea) 
#       integer(IntKi), intent(in) :: e0  !< Index of an element
#       type(IList), intent(inout) :: Er  !< List of rigid elements
#       type(IList), intent(inout) :: Ea  !< List of elements in a rigid assembly
#       type(IList)     :: En             !< List of neighbors of e0
#       integer (IntKi) :: ik
#       integer (IntKi) :: ek, ek2
#       integer (IntKi) :: iWhichNode_e0, iWhichNode_ek
#       call init_list(En, 0, 0, ErrStat2, ErrMsg2)
#       ! Loop through all elements, setup list of e0-neighbors, add them to Ea, remove them from Er
#       ik=0
#       do while (ik< len(Er))
#          ik=ik+1
#          ek = Er%List(ik)
#          if (ElementsConnected(p, e0, ek, iWhichNode_e0, iWhichNode_ek)) then
#             if (DEV_VERSION) then
#                print*,'Element ',ek,'is connected to element',e0,'via its node',iWhichNode_ek
#             endif
#             ! Remove element from Er (a rigid element can belong to only one assembly)
#             ek2 =  pop(Er, ik,  ErrStat2, ErrMsg2) ! same as ek before
#             ik=ik-1
#             if (ek/=ek2) then
#                print*,'Problem in popping',ek,ek2
#                STOP
#             endif
#             call append(En, ek, ErrStat2, ErrMsg2)
#             call append(Ea, ek, ErrStat2, ErrMsg2)
#          endif
#       enddo
#       ! Loop through neighbors and recursively add neighbors of neighbors
#       do ik = 1, len(En)
#          ek = En%List(ik)
#          call AddNeighbors(ek, Er, Ea)
#       enddo
#       call destroy_list(En, ErrStat2, ErrMsg2)

def buildTMatrix(self):
    """ 
    !------------------------------------------------------------------------------------------------------
    !> Build transformation matrix T, such that x= T.x~ where x~ is the reduced vector of DOF
    !! Variables set by this routine
    !! - p%NodesDOFred(iNode)=[list of DOF]: Created for each node, the list of DOF of this node in the 
    !!         reduced system. 
    !!         NOTE: follower nodes in rigid assembly have no DOFred (convention)
    !! - p%nDOF_red: number of DOF in reduced system (<= nDOF)
    !! - p%reduced: true if a reduction is needed, i.e. a T matrix is needed, and nDOF_red<nDOF
    !!
    !! Variables returned:
    !! - T_red: retuction matrix such that x= T_red.x~ where x~ is the reduced vector of DOF
    """
    pass
#    use IntegerList, only: init_list, find, pop, destroy_list, len
#    use IntegerList, only: print_list
#    TYPE(SD_InitType),            INTENT(IN   ) :: Init
#    TYPE(SD_ParameterType),target,INTENT(INOUT) :: p
#    type(IList), dimension(:),    INTENT(IN   ) :: RA   !< RA(a) = [e1,..,en]  list of elements forming a rigid link assembly
#    integer(IntKi), dimension(:), INTENT(IN   ) :: RAm1 !< RA^-1(e) = a , for a given element give the index of a rigid assembly
#    INTEGER(IntKi),               INTENT(  OUT) :: ErrStat     ! Error status of the operation
#    CHARACTER(*),                 INTENT(  OUT) :: ErrMsg      ! Error message if ErrStat /= ErrID_None
#    real(FEKi), dimension(:,:), allocatable :: T_red !< Transformation matrix for DOF elimination
#    ! Local  
#    real(ReKi), dimension(:,:), allocatable   :: Tc
#    integer(IntKi), dimension(:), allocatable :: INodesID !< List of unique nodes involved in Elements
#    integer(IntKi), dimension(:), allocatable :: IDOFOld !< 
#    integer(IntKi), dimension(:), pointer :: IDOFNew !< 
#    real(ReKi), dimension(6,6) :: I6       !< Identity matrix of size 6
#    integer(IntKi) :: iPrev
#    type(IList) :: IRA !< list of rigid assembly indices to process
#    integer(IntKi) :: aID, ia ! assembly ID, and index in IRA
#    integer(IntKi) :: iNode, iNodeSel, iNodeRemaining, iiNodeRemaining
#    integer(IntKi) :: er !< Index of one rigid element belong to a rigid assembly
#    integer(IntKi) :: JType
#    integer(IntKi) :: I
#    integer(IntKi) :: nc !< Number of DOF after constraints applied
#    integer(IntKi) :: nj
#    real(ReKi)  :: phat(3) !< Directional vector of the joint
#    type(IList), dimension(:), allocatable :: RA_DOFred ! DOF indices for each rigid assembly, in reduced system
#    INTEGER(IntKi)       :: ErrStat2
#    CHARACTER(ErrMsgLen) :: ErrMsg2
#    ErrStat = ErrID_None
#    ErrMsg  = ""
# 
#    ! --- Misc inits
#    nullify(IDOFNew)
#    I6(1:6,1:6)=0; do i = 1,6 ; I6(i,i)=1_ReKi; enddo ! I6 =  eye(6)
#    allocate(p%NodesDOFred(1:p%nNodes), stat=ErrStat2); if(Failed()) return; ! Indices of DOF for each joint, in reduced system
#    allocate(RA_DOFred(1:size(RA)), stat=ErrStat2); if(Failed()) return; ! Indices of DOF for each rigid assmbly, in reduced system
# 
#    p%nDOF_red = nDOF_ConstraintReduced()
#    p%reduced  = reductionNeeded()      ! True if reduction needed, allow for optimization if not needed
# 
#    if (DEV_VERSION) then
#       print*,'nDOF constraint elim', p%nDOF_red , '/' , p%nDOF
#    endif
#    CALL AllocAry( T_red, p%nDOF, p%nDOF_red, 'p%T_red',  ErrStat2, ErrMsg2); if(Failed()) return; ! system stiffness matrix 
#    T_red=0.0_FeKi
#    call init_list(IRA, size(RA), 0, ErrStat2, ErrMsg2); if(Failed()) return;
#    IRA%List(1:size(RA)) = (/(ia , ia = 1,size(RA))/)
#    if (DEV_VERSION) then
#       call print_list(IRA, 'List of RA indices')
#    endif
# 
#    ! --- For each node:
#    !  - create list of indices I      in the assembled vector of DOF
#    !  - create list of indices Itilde in the reduced vector of DOF
#    !  - increment iPrev by the number of DOF of Itilde
#    iPrev =0 
#    do iNode = 1, p%nNodes
#       iNodeSel = iNode ! Unless changed by Rigid assembly, using this index
#       if (allocated(Tc))      deallocate(Tc)
#       if (allocated(IDOFOld)) deallocate(IDOFOld)
#       JType = int(Init%Nodes(iNodeSel,iJointType))
#       if(JType == idJointCantilever ) then
#          if ( NodeHasRigidElem(iNodeSel, Init, p, er)) then ! return True and element index "er" if element is rigid
#             ! --- The joint is involved in a rigid link assembly
#             aID = RAm1(er) ! ID of rigid assembly for element er
#             if (aID<0) then
#                call Fatal('No rigid assembly attributed to node'//trim(Num2LStr(iNodeSel))//'. RAm1 wrong'); return
#             endif
#             ia  = find(IRA, aID, ErrStat2, ErrMsg2); if(Failed()) return ! We "pop" IRA, so index and ID are different
#             if (DEV_VERSION) then
#                print'(4X,A,I5,A,I5,A,I5)','Node',iNodeSel, ' is involved in RA:', aID, '. Current index in list of RA', ia
#             endif
#             if ( ia <= 0) then
#                ! --- This rigid assembly has already been processed, simple triggers below
#                ! OLD: The DOF list is taken from the stored RA DOF list
#                ! call init_list(p%NodesDOFred(iNodeSel), RA_DOFred(aID)%List, ErrStat2, ErrMsg2)
#                ! NEW: this node has no DOFs, so we set an empty list of DOFred for this node
#                !call init_list(p%NodesDOFred(iNodeSel), 0, 0, ErrStat2, ErrMsg2)
#                if (DEV_VERSION) then
#                   print*,'   The RA',aID,', has already been processed!'! The following node has no reduced DOF'
#                   !print*,'   but based on its RA, we can list its Itilde DOF:'
#                   !print*,'   N',iNodeSel,'I ',p%NodesDOF(iNodeSel)%List(1:6)
#                   !print*,'   N',iNodeSel,'It',RA_DOFred(aID)%List
#                endif
#                cycle ! We pass to the next joint, important so that:
#                !     - we don't increase iPrev
#                !     - we don't set Tc
#                !     - p%NodesDOFred is not set (assuming it has already been done)
#             else
#                ! --- Proceeding the rigid assembly
#                ! Returns TC and INodesID, do not change other variables
#                call RAElimination( RA(aID)%List, Tc, INodesID, Init, p, ErrStat2, ErrMsg2); if(Failed()) return;
#                aID = pop(IRA, ia, ErrStat2, ErrMsg2) ! this assembly has been processed, remove it from IRA list
#                nj = size(INodesID) ! Number of nodes in this rigid assembly
#                allocate(IDOFOld(1:6*nj))
#                do I=1, nj
#                   IDOFOld( (I-1)*6+1 : I*6 ) = p%NodesDOF(INodesID(I))%List(1:6)
#                enddo
# 
#                ! Storing DOF list for this RA (Note: same as NodesDOFred below, only for debug)
#                nc=size(Tc,2) ! Should be 6 
#                call init_list(RA_DOFred(aID), (/ (iprev + i, i=1,nc) /), ErrStat2, ErrMsg2);
# 
#                ! --- Processing trigger for leader/follower Nodes
#                iNodeSel = INodesID(1)  ! The first index returned is the leader of the assembly, we use this from now on
#                do iiNodeRemaining=2,size(INodesID) ! start at 2 because 1 is always the leader
#                   iNodeRemaining = INodesID(iiNodeRemaining)
#                   ! OLD: The DOF list is taken from the stored RA DOF list
#                   ! call init_list(p%NodesDOFred(iNode), RA_DOFred(aID)%List, ErrStat2, ErrMsg2)
#                   ! NEW: this node has no DOFs, so we set an empty list of DOFred for this node
#                   call init_list(p%NodesDOFred(iNodeRemaining), 0, 0, ErrStat2, ErrMsg2)
#                   if (DEV_VERSION) then
#                      print'(4X,A,I5,A,I5,I5)','Node',iNodeRemaining,' has no reduced DOF since its the follower of leader node ',INodesID(1),iNodeSel
#                   endif
#                enddo
#             endif
#          else
#             ! --- Regular cantilever joint
#             ! TODO/NOTE: We could apply fixed constraint/BC here, returning Tc as a 6xn matrix with n<6
#             !            Extreme case would be Tc: 6*0, in which case NodesDOFred would be empty ([])
#             allocate(Tc(1:6,1:6))
#             allocate(IDOFOld(1:6))
#             Tc=I6
#             IDOFOld = p%NodesDOF(iNodeSel)%List(1:6)
#          endif
#       else
#          ! --- Ball/Pin/Universal joint
#          allocate(IDOFOld(1:len(p%NodesDOF(iNodeSel))))
#          IDOFOld(:) = p%NodesDOF(iNodeSel)%List(:)
#          phat = Init%Nodes(iNodeSel, iJointDir:iJointDir+2)
#          ! Return Tc, do not change other variable
#          call JointElimination(Init%NodesConnE(iNodeSel,:), JType, phat, p, Tc, ErrStat2, ErrMsg2); if(Failed()) return
#       endif ! Cantilever or Special Joint
#       nc=size(Tc,2) 
#       call init_list(p%NodesDOFred(iNodeSel), nc, 0, ErrStat2, ErrMsg2)
#       p%NodesDOFred(iNodeSel)%List(1:nc) = (/ (iprev + i, i=1,nc) /)
#       IDOFNew => p%NodesDOFred(iNodeSel)%List(1:nc) ! alias to shorten notations
#       if (DEV_VERSION) then
#          ! KEEP ME, VERY USEFUL
#          print*,'N',iNodeSel,'I ',IDOFOld
#          print*,'N',iNodeSel,'It',IDOFNew
#       endif
#       T_red(IDOFOld, IDOFNew) = Tc
#       iPrev = iPrev + nc
#    enddo
#    if (DEV_VERSION) then
#       print'(A)','--- End of BuildTMatrix'
#       print*,'   - T_red set'
#       print*,'   - p%nDOF_red', p%nDOF_red
#       print*,'   - p%reduced ', p%reduced
#       print*,'   - p%NodesDOFred: (list of reduced DOF indices per node) '
#       do iNode = 1, p%nNodes
#          print*,'N',iNode, 'It', p%NodesDOFred(iNode)%List(:)
#       enddo
#    endif
# 
#    ! --- Safety checks
#    if (len(IRA)>0) then 
#       call Fatal('Not all rigid assemblies were processed'); return
#    endif
#    if (iPrev /= p%nDOF_red) then 
#       call Fatal('Inconsistency in number of reduced DOF'); return
#    endif
# # 
# 
