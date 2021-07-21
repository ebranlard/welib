!>  Takes nPart particle and nCPs control points
subroutine fui_part_c(CPs,ioff,iCP_beg, nCPs_loc, nCPs, iComputeGrad, Uind_out, Grad_out, &
        nPartTot, Part, iPartCurPos, PartSmoothModel, PartSmoothParam, PartIntensities) 
    use PrecisionMod, only: C_INT, C_DOUBLE
    implicit none
    ! Arguments declarations 
    integer(C_INT), intent(in)                        :: iComputeGrad    !< 
    integer(C_INT), intent(in)                        :: ioff
    integer(C_INT), intent(in)                        :: iCP_beg
    integer(C_INT), intent(in)                        :: nCPs_loc
    integer(C_INT), intent(in)                        :: nCPs
    real(C_DOUBLE), dimension(3,nCPs), intent(in)     :: CPs             !< 
    real(C_DOUBLE), dimension(3,nCPs), intent(inout)    :: Uind_out        !< Induced velocity vector - No Side effects!!!
    real(C_DOUBLE), dimension(9,nCPs), intent(inout)    :: Grad_out        !< Gradient vector - No Side effects
    integer(C_INT), intent(in)                        :: nPartTot
    real(C_DOUBLE), dimension(3,nPartTot), intent(in) :: Part
    integer(C_INT) ,intent(in)                        :: iPartCurPos
    integer(C_INT), intent(in)                        :: PartSmoothModel
    real(C_DOUBLE), dimension(nPartTot), intent(in)   :: PartSmoothParam
    real(C_DOUBLE), dimension(3,nPartTot), intent(in) :: PartIntensities
    integer :: icp,ip

    ! 
    real(C_DOUBLE),parameter :: MINNORM=1e-4_C_DOUBLE
    real(C_DOUBLE), parameter :: fourpi_inv    = 0.079577471545947667884441881686257181017229822870228224373833_C_DOUBLE
    ! --------------------------------------------------------------------------------
    ! ---  Variables for particle 11
    ! --------------------------------------------------------------------------------
    ! Kindof arguments
    real(C_DOUBLE), dimension(9) :: Grad       !< no side effects
    real(C_DOUBLE), dimension(3) :: Uind       !< no side effects
    real(C_DOUBLE), dimension(3) :: DeltaP     !< CP-PP "control point - particle point"
    real(C_DOUBLE), dimension(3) :: Omega      !< 
    real(C_DOUBLE)               :: Eps        !< 
    ! 
    real(C_DOUBLE),dimension(3)  :: C          !< Cross product of Omega and r
    real(C_DOUBLE)               :: E          !< Exponential poart for the mollifider
    real(C_DOUBLE)               :: Ebar       !< 1-E
    real(C_DOUBLE)               :: r3_inv     !< 
    real(C_DOUBLE)               :: n5_inv     !< 
    real(C_DOUBLE)               :: n2_inv     !< TODO, remove this one or r3..but watch out, look for the Grad trick of model 2
    real(C_DOUBLE)               :: rDeltaP    !< norm , distance between point and particle
    real(C_DOUBLE)               :: ScalarPart !< the part containing the inverse of the distance, but not 4pi, Mollifier

    
    !print*,'Part fortran',(iCP_beg+1),(iCP_beg+1)+nCPS_loc-1,nCPs
    !$OMP PARALLEL default(shared)
    !$OMP do private(&
    !$OMP& icp,ip,Uind,Grad,DeltaP,Omega,Eps,C,E,Ebar,r3_inv,n5_inv,n2_inv,rDeltaP,ScalarPart&
    !$OMP& ) schedule(runtime)
    ! loop on CPs 
    do icp=(iCP_beg+1),(iCP_beg+1)+nCPS_loc-1
        ! loop on particles
        do ip=1,iPartCurPos
            ! Setting up, kindof arguments for particle 11
            Uind(1:3) = 0.0_C_DOUBLE
            Grad(1:9) = 0.0_C_DOUBLE
            DeltaP = CPs(:,icp)-Part(1:3,ip)
            Eps    = PartSmoothParam(ip)
            Omega  = PartIntensities(1:3,ip)
            ! --------------------------------------------------------------------------------
            ! --- Inlining Particle 11 
            ! --------------------------------------------------------------------------------
            rDeltaP=sqrt(DeltaP(1)**2+ DeltaP(2)**2+ DeltaP(3)**2)! norm
            ! TODO, escape should depend on the 
            if (rDeltaP<MINNORM) then! to avoid singularity, save a bit of computation and use limit values
                !--------------------------------------------------------------------------------
                !--- Exactly on the Singularity 
                !--------------------------------------------------------------------------------
                !return

                ! NEW setting the gradient to 0, to avoid self influence (artifact of the low order Smoothing Kernel)
                ! => Commenting what is below
                !             if( SmoothModel==0) then 
                !                 !return
                !             end if
                !             if(iComputeGrad) then
                !                 if( SmoothModel==2) then 
                !                     r3_inv=1_C_DOUBLE/sqrt(Eps) ! Epsilon^3
                !                 else
                !                     r3_inv=1_C_DOUBLE/Eps ! Epsilon ^3
                !                 endif
                !                 ! For model 1 and 2 is the same limit
                !                 Grad(2) = fourpi_inv*( -Omega(3) *r3_inv) ! 12
                !                 Grad(3) = fourpi_inv*( +Omega(2) *r3_inv) ! 13
                !                 Grad(4) = fourpi_inv*( +Omega(3) *r3_inv) ! 21
                !                 Grad(6) = fourpi_inv*( -Omega(1) *r3_inv) ! 23
                !                 Grad(7) = fourpi_inv*( -Omega(2) *r3_inv) ! 31
                !                 Grad(8) = fourpi_inv*( +Omega(1) *r3_inv) ! 33
                !             endif
                ! TODO far distance
            else
                !--------------------------------------------------------------------------------
                !--- Normal Procedure 
                !--------------------------------------------------------------------------------


                C(1) = Omega(2) * DeltaP(3) - Omega(3) * DeltaP(2)
                C(2) = Omega(3) * DeltaP(1) - Omega(1) * DeltaP(3)
                C(3) = Omega(1) * DeltaP(2) - Omega(2) * DeltaP(1)
                ! influence of all control points 
                E    = 0.0_C_DOUBLE
                Ebar = 1.0_C_DOUBLE
                select case (PartSmoothModel) !
                case (0) ! No mollification, exept on zero
                    r3_inv=1._C_DOUBLE/(rDeltaP**3)
                    ScalarPart=r3_inv
                case (2) !
                    !Exponential mollifier
                    r3_inv=1._C_DOUBLE/(rDeltaP**3)
                    E=exp(-rDeltaP**3/Eps**3)
                    Ebar=1._C_DOUBLE-E
                    ScalarPart=Ebar*r3_inv
                case (3) !
                    !Exponential mollifier to the power 4
                    r3_inv=1._C_DOUBLE/(rDeltaP**3)
                    E=exp(-rDeltaP**4/Eps**4)
                    Ebar=1._C_DOUBLE-E
                    ScalarPart=Ebar*r3_inv
                case (22) !
                    ! Finite support
                    r3_inv= 1._C_DOUBLE/sqrt(Eps**6+rDeltaP**6)
                    ScalarPart= r3_inv
                case default 
                    print*,'SmoothModel',PartSmoothModel
                    print*,'Blob: wrong viscous model'
                    STOP
                end select 

                Uind(1:3)=C*ScalarPart*fourpi_inv

                !          write(*,*) 'DeltaP', DeltaP
                !          write(*,*) 'Omega', Omega
                !          write(*,*) 'Smoothparam', Eps
                !          write(*,*) 'SmoothModel', SmoothModel
                !          write(*,*) 'C', C
                !          write(*,*) 'r3_inv', r3_inv
                !          write(*,*) 'ScalarPart', ScalarPart
                !          write(*,*) 'UIout', UIout
                !          STOP
                if (iComputeGrad==1) then 
                    n5_inv=1._C_DOUBLE/rDeltaP**5
                    select case ( PartSmoothModel) !
                    ! Omega is my alpha
                    ! C = (alpha x r)  The Omega term comes from grad(alpha x r)
                    case (0) !
                        Grad(1) = fourpi_inv*(                   -3.*DeltaP(1) * C(1) *n5_inv) ! 11
                        Grad(2) = fourpi_inv*( -Omega(3) *r3_inv -3.*DeltaP(2) * C(1) *n5_inv) ! 12
                        Grad(3) = fourpi_inv*( +Omega(2) *r3_inv -3.*DeltaP(3) * C(1) *n5_inv) ! 13

                        Grad(4) = fourpi_inv*( +Omega(3) *r3_inv -3.*DeltaP(1) * C(2) *n5_inv) ! 21
                        Grad(5) = fourpi_inv*(                   -3.*DeltaP(2) * C(2) *n5_inv) ! 22
                        Grad(6) = fourpi_inv*( -Omega(1) *r3_inv -3.*DeltaP(3) * C(2) *n5_inv) ! 23

                        Grad(7) = fourpi_inv*( -Omega(2) *r3_inv -3.*DeltaP(1) * C(3) *n5_inv) ! 31
                        Grad(8) = fourpi_inv*( +Omega(1) *r3_inv -3.*DeltaP(2) * C(3) *n5_inv) ! 33
                        Grad(9) = fourpi_inv*(                   -3.*DeltaP(3) * C(3) *n5_inv) ! 32

                    case (2) !
                        ! Exponential mollifier 
                        n2_inv=1._C_DOUBLE/rDeltaP**2
                        ! Almost the same as above 
                        ! Remember I chose Eps=epsilon**3, so don't cube it twice!!
                        Grad(1) = fourpi_inv*(                       + (-3.*n5_inv*Ebar + 3./Eps**3*n2_inv*E)*DeltaP(1)*C(1))!11
                        Grad(2) = fourpi_inv*( -Omega(3)*r3_inv*Ebar + (-3.*n5_inv*Ebar + 3./Eps**3*n2_inv*E)*DeltaP(2)*C(1))!12
                        Grad(3) = fourpi_inv*( +Omega(2)*r3_inv*Ebar + (-3.*n5_inv*Ebar + 3./Eps**3*n2_inv*E)*DeltaP(3)*C(1))!13

                        Grad(4) = fourpi_inv*( +Omega(3)*r3_inv*Ebar + (-3.*n5_inv*Ebar + 3./Eps**3*n2_inv*E)*DeltaP(1)*C(2))!21
                        Grad(5) = fourpi_inv*(                       + (-3.*n5_inv*Ebar + 3./Eps**3*n2_inv*E)*DeltaP(2)*C(2))!22
                        Grad(6) = fourpi_inv*( -Omega(1)*r3_inv*Ebar + (-3.*n5_inv*Ebar + 3./Eps**3*n2_inv*E)*DeltaP(3)*C(2))!23

                        Grad(7) = fourpi_inv*( -Omega(2)*r3_inv*Ebar + (-3.*n5_inv*Ebar + 3./Eps**3*n2_inv*E)*DeltaP(1)*C(3))!31
                        Grad(8) = fourpi_inv*( +Omega(1)*r3_inv*Ebar + (-3.*n5_inv*Ebar + 3./Eps**3*n2_inv*E)*DeltaP(2)*C(3))!33
                        Grad(9) = fourpi_inv*(                       + (-3.*n5_inv*Ebar + 3./Eps**3*n2_inv*E)*DeltaP(3)*C(3))!32

                    case (3) !
                        ! Exponential mollifier 
                        n2_inv=1._C_DOUBLE/rDeltaP !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!FAKE
                        ! Almost the same as above 
                        ! Remember I chose Eps=epsilon**3, so don't cube it twice!!
                        Grad(1) = fourpi_inv*(                       + (-3.*n5_inv*Ebar + 4./Eps**4*n2_inv*E)*DeltaP(1)*C(1))!11
                        Grad(2) = fourpi_inv*( -Omega(3)*r3_inv*Ebar + (-3.*n5_inv*Ebar + 4./Eps**4*n2_inv*E)*DeltaP(2)*C(1))!12
                        Grad(3) = fourpi_inv*( +Omega(2)*r3_inv*Ebar + (-3.*n5_inv*Ebar + 4./Eps**4*n2_inv*E)*DeltaP(3)*C(1))!13

                        Grad(4) = fourpi_inv*( +Omega(3)*r3_inv*Ebar + (-3.*n5_inv*Ebar + 4./Eps**4*n2_inv*E)*DeltaP(1)*C(2))!21
                        Grad(5) = fourpi_inv*(                       + (-3.*n5_inv*Ebar + 4./Eps**4*n2_inv*E)*DeltaP(2)*C(2))!22
                        Grad(6) = fourpi_inv*( -Omega(1)*r3_inv*Ebar + (-3.*n5_inv*Ebar + 4./Eps**4*n2_inv*E)*DeltaP(3)*C(2))!23

                        Grad(7) = fourpi_inv*( -Omega(2)*r3_inv*Ebar + (-3.*n5_inv*Ebar + 4./Eps**4*n2_inv*E)*DeltaP(1)*C(3))!31
                        Grad(8) = fourpi_inv*( +Omega(1)*r3_inv*Ebar + (-3.*n5_inv*Ebar + 4./Eps**4*n2_inv*E)*DeltaP(2)*C(3))!33
                        Grad(9) = fourpi_inv*(                       + (-3.*n5_inv*Ebar + 4./Eps**4*n2_inv*E)*DeltaP(3)*C(3))!32

                    case (22) !
                        ! Finite support with Epsilon 6 
                        n5_inv= r3_inv**3 * rDeltaP**4 ! Fake n5_inv
                        ! EXACTLY the same as case 0 due to fake variables n3 and n5 
                        Grad(1) = fourpi_inv*(                   -3.*DeltaP(1) * C(1) *n5_inv) ! 11
                        Grad(2) = fourpi_inv*( -Omega(3) *r3_inv -3.*DeltaP(2) * C(1) *n5_inv) ! 12
                        Grad(3) = fourpi_inv*( +Omega(2) *r3_inv -3.*DeltaP(3) * C(1) *n5_inv) ! 13

                        Grad(4) = fourpi_inv*( +Omega(3) *r3_inv -3.*DeltaP(1) * C(2) *n5_inv) ! 21
                        Grad(5) = fourpi_inv*(                   -3.*DeltaP(2) * C(2) *n5_inv) ! 22
                        Grad(6) = fourpi_inv*( -Omega(1) *r3_inv -3.*DeltaP(3) * C(2) *n5_inv) ! 23

                        Grad(7) = fourpi_inv*( -Omega(2) *r3_inv -3.*DeltaP(1) * C(3) *n5_inv) ! 31
                        Grad(8) = fourpi_inv*( +Omega(1) *r3_inv -3.*DeltaP(2) * C(3) *n5_inv) ! 33
                        Grad(9) = fourpi_inv*(                   -3.*DeltaP(3) * C(3) *n5_inv) ! 32
                    end select  ! Smooth Model
                else
                    Grad(1:9)=0.0_C_DOUBLE
                end if  ! iComputeGrad
                !print*,'Uind',Uind
            end if  ! Norm Big enough
            ! --------------------------------------------------------------------------------
            ! --- END of particle 11 inlining 
            ! --------------------------------------------------------------------------------
            !icp-icp_beg+1
            Uind_out(1:3,ioff+icp-icp_beg)=Uind_out(1:3,ioff+icp-icp_beg)+Uind(1:3)
            if(iComputeGrad==1) then
                Grad_out(1:9,ioff+icp-icp_beg)=Grad_out(1:9,ioff+icp-icp_beg)+Grad(1:9)
            endif
        enddo! loop particles
    enddo ! loop Cps
    !$OMP END DO 
    !$OMP END PARALLEL
end subroutine

