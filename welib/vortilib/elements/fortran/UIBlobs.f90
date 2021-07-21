!> Module for blobs induced velocities
module UIBlobs
    use PrecisionMod, only: MK
    implicit none

    real(MK),parameter :: MINNORM=1e-4

contains

    !>  Takes nPart particle and nCPs control points
    subroutine fUi_Particle(CPs,Part, Omega, SmoothModel, SmoothParam, bComputeGrad, Uiout, Gradout,nCPs,nPart)
        ! Arguments declarations 
        integer, intent(in) :: nCPs
        integer, intent(in) :: nPart
        logical, intent(in) :: bComputeGrad !< 
        real(MK), dimension(3,nCPs), intent(in) :: CPs   !< side effects
        real(MK), dimension(3,nCPs), intent(inout) :: UiOut   !< side effects
        real(MK), dimension(9,nCPs), intent(inout)  :: Gradout !< side effects
        real(MK), dimension(3,nPart), intent(in)   :: Part      !< 
        real(MK), dimension(3,nPart), intent(in)   :: Omega      !< 
        integer, intent(in)    :: SmoothModel !< 
        real(MK), dimension(nPart), intent(in)   :: SmoothParam !< 
        ! Variable declarations 
        real(MK), dimension(3) :: UItmp   !< 
        real(MK), dimension(9) :: Gradtmp !< 
        real(MK), dimension(3) :: DP      !< 
        integer :: icp,ip

        ! loop on CPs 
        do icp=1,nCPs
            ! loop on particles
            do ip=1,nPart
                UItmp=0.0_MK
                Gradtmp=0.0_MK 
                DP=CPs(:,icp)-Part(1:3,ip)
                call  fUi_Particle_11(DP, Omega(1:3,ip), SmoothModel ,&
                    SmoothParam(ip), bComputeGrad, UItmp,Gradtmp)
                UIout(1:3,icp)=UIout(1:3,icp)+UItmp(1:3)
                if(bComputeGrad) then
                    Gradout(1:9,icp)=Gradout(1:9,icp)+Gradtmp(1:9)
                endif
            enddo! loop particles
        enddo ! loop Cps
    end subroutine

    !>  Takes one particle and one point
    subroutine fUi_Particle_11(DeltaP, Omega, SmoothModel, SmoothParam, bComputeGrad, Ui, Grad)
        use MathConstants, only: fourpi_inv
        ! Arguments declarations 
        logical, intent(in) :: bComputeGrad  !<  
        real(MK), dimension(9), intent(out) :: Grad  !< no side effects
        real(MK), dimension(3), intent(out) :: Ui  !<  no side effects
        real(MK), dimension(3), intent(in) :: DeltaP  !<  CP-PP "control point - particle point"
        real(MK), dimension(3), intent(in) :: Omega  !<  
        integer, intent(in) :: SmoothModel  !<  
        real(MK), intent(in) :: SmoothParam  !<  
        ! Variable declarations 
        real(MK),dimension(3) :: C  !<  Cross product of Omega and r
        real(MK) :: E  !<   Exponential poart for the mollifider
        real(MK) :: Ebar  !<  1-E
        real(MK) :: r3_inv  !<  
        real(MK) :: n5_inv  !<  
        real(MK) :: n2_inv  !<   TODO, remove this one or r3..but watch out, look for the Grad trick of model 2
        real(MK) :: rDeltaP  !< norm , distance between point and particle 
        real(MK) :: ScalarPart  !< the part containing the inverse of the distance, but not 4pi, Mollifier

        Ui(1:3)   = 0.0_MK
        Grad(1:9) = 0.0_MK
        rDeltaP=sqrt(DeltaP(1)**2+ DeltaP(2)**2+ DeltaP(3)**2)! norm
        ! TODO, escape should depend on the 
        if (rDeltaP<MINNORM) then! to avoid singularity, save a bit of computation and use limit values
            !--------------------------------------------------------------------------------
            !--- Exactly on the Singularity 
            !--------------------------------------------------------------------------------
            return
        
            ! NEW setting the gradient to 0, to avoid self influence (artifact of the low order Smoothing Kernel)
            ! => Commenting what is below
            !             if( SmoothModel==0) then 
            !                 return
            !             end if
            !             if(bComputeGrad) then
            !                 if( SmoothModel==2) then 
            !                     r3_inv=1_MK/sqrt(SmoothParam) ! Epsilon^3
            !                 else
            !                     r3_inv=1_MK/SmoothParam ! Epsilon ^3
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
            select case (SmoothModel) !
            case (0) ! No mollification, exept on zero
                r3_inv=1._MK/(rDeltaP**3)
                ScalarPart=r3_inv
            case (2) !
                !Exponential mollifier
                r3_inv=1._MK/(rDeltaP**3)
                E=exp(-rDeltaP**3/SmoothParam**3)
                Ebar=1._MK-E
                ScalarPart=Ebar*r3_inv
            case (3) !
                !Exponential mollifier to the power 4
                r3_inv=1._MK/(rDeltaP**3)
                E=exp(-rDeltaP**4/SmoothParam**4)
                Ebar=1._MK-E
                ScalarPart=Ebar*r3_inv
            case (22) !
                ! Finite support
                r3_inv= 1._MK/sqrt(SmoothParam**6+rDeltaP**6)
                ScalarPart= r3_inv
            case default 
                print*,'SmoothModel',SmoothModel
                print*,'Blob: wrong viscous model'
            end select 

            Ui(1:3)=C*ScalarPart*fourpi_inv

        !          write(*,*) 'DeltaP', DeltaP
        !          write(*,*) 'Omega', Omega
        !          write(*,*) 'Smoothparam', SmoothParam
        !          write(*,*) 'SmoothModel', SmoothModel
        !          write(*,*) 'C', C
        !          write(*,*) 'r3_inv', r3_inv
        !          write(*,*) 'ScalarPart', ScalarPart
        !          write(*,*) 'UIout', UIout
        !          STOP
            if (bComputeGrad) then 
                n5_inv=1._MK/rDeltaP**5
                select case ( SmoothModel) !
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
                    n2_inv=1._MK/rDeltaP**2
                    ! Almost the same as above 
                    ! Remember I chose SmoothParam=epsilon**3, so don't cube it twice!!
                    Grad(1) = fourpi_inv*(                       + (-3.*n5_inv*Ebar + 3./SmoothParam**3*n2_inv*E)*DeltaP(1)*C(1))!11
                    Grad(2) = fourpi_inv*( -Omega(3)*r3_inv*Ebar + (-3.*n5_inv*Ebar + 3./SmoothParam**3*n2_inv*E)*DeltaP(2)*C(1))!12
                    Grad(3) = fourpi_inv*( +Omega(2)*r3_inv*Ebar + (-3.*n5_inv*Ebar + 3./SmoothParam**3*n2_inv*E)*DeltaP(3)*C(1))!13

                    Grad(4) = fourpi_inv*( +Omega(3)*r3_inv*Ebar + (-3.*n5_inv*Ebar + 3./SmoothParam**3*n2_inv*E)*DeltaP(1)*C(2))!21
                    Grad(5) = fourpi_inv*(                       + (-3.*n5_inv*Ebar + 3./SmoothParam**3*n2_inv*E)*DeltaP(2)*C(2))!22
                    Grad(6) = fourpi_inv*( -Omega(1)*r3_inv*Ebar + (-3.*n5_inv*Ebar + 3./SmoothParam**3*n2_inv*E)*DeltaP(3)*C(2))!23

                    Grad(7) = fourpi_inv*( -Omega(2)*r3_inv*Ebar + (-3.*n5_inv*Ebar + 3./SmoothParam**3*n2_inv*E)*DeltaP(1)*C(3))!31
                    Grad(8) = fourpi_inv*( +Omega(1)*r3_inv*Ebar + (-3.*n5_inv*Ebar + 3./SmoothParam**3*n2_inv*E)*DeltaP(2)*C(3))!33
                    Grad(9) = fourpi_inv*(                       + (-3.*n5_inv*Ebar + 3./SmoothParam**3*n2_inv*E)*DeltaP(3)*C(3))!32

                case (3) !
                    ! Exponential mollifier 
                    n2_inv=1._MK/rDeltaP !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!FAKE
                    ! Almost the same as above 
                    ! Remember I chose SmoothParam=epsilon**3, so don't cube it twice!!
                    Grad(1) = fourpi_inv*(                       + (-3.*n5_inv*Ebar + 4./SmoothParam**4*n2_inv*E)*DeltaP(1)*C(1))!11
                    Grad(2) = fourpi_inv*( -Omega(3)*r3_inv*Ebar + (-3.*n5_inv*Ebar + 4./SmoothParam**4*n2_inv*E)*DeltaP(2)*C(1))!12
                    Grad(3) = fourpi_inv*( +Omega(2)*r3_inv*Ebar + (-3.*n5_inv*Ebar + 4./SmoothParam**4*n2_inv*E)*DeltaP(3)*C(1))!13

                    Grad(4) = fourpi_inv*( +Omega(3)*r3_inv*Ebar + (-3.*n5_inv*Ebar + 4./SmoothParam**4*n2_inv*E)*DeltaP(1)*C(2))!21
                    Grad(5) = fourpi_inv*(                       + (-3.*n5_inv*Ebar + 4./SmoothParam**4*n2_inv*E)*DeltaP(2)*C(2))!22
                    Grad(6) = fourpi_inv*( -Omega(1)*r3_inv*Ebar + (-3.*n5_inv*Ebar + 4./SmoothParam**4*n2_inv*E)*DeltaP(3)*C(2))!23

                    Grad(7) = fourpi_inv*( -Omega(2)*r3_inv*Ebar + (-3.*n5_inv*Ebar + 4./SmoothParam**4*n2_inv*E)*DeltaP(1)*C(3))!31
                    Grad(8) = fourpi_inv*( +Omega(1)*r3_inv*Ebar + (-3.*n5_inv*Ebar + 4./SmoothParam**4*n2_inv*E)*DeltaP(2)*C(3))!33
                    Grad(9) = fourpi_inv*(                       + (-3.*n5_inv*Ebar + 4./SmoothParam**4*n2_inv*E)*DeltaP(3)*C(3))!32

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
                Grad(1:9)=0.0_MK
            end if  ! bComputeGrad
        end if  ! Norm Big enough
    end subroutine fUi_Particle_11
end module UIBlobs

