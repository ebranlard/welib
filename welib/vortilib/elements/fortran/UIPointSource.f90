!> 
module UIPointSource
    use PrecisionMod, only: MK
    implicit none

    real(MK),parameter :: MINNORM=1e-4
contains

    subroutine fUi_PointSource(CPs,SrcP,Intensities, SmoothModel,SmoothParam,bComputeGrad,Uind, Grad,nCPs,nSrc)
        ! Arguments declarations 
        integer, intent(in) :: nCPs
        integer, intent(in) :: nSrc
        logical, intent(in) :: bComputeGrad !< 
        real(MK), dimension(3,nCPs), intent(in)    :: CPs         !< side effects
        real(MK), dimension(3,nCPs), intent(inout) :: Uind        !< side effects
        real(MK), dimension(9,nCPs), intent(inout) :: Grad        !< side effects
        real(MK), dimension(3,nSrc+1), intent(in)  :: SrcP        !< 
        real(MK), dimension(nSrc), intent(in)      :: Intensities !< 
        integer, intent(in)                        :: SmoothModel !< 
        real(MK), dimension(nSrc), intent(in)      :: SmoothParam !< 
        ! Variable declarations 
        integer :: icp,is
        real(MK), dimension(3) :: Uindtmp !< 
        real(MK), dimension(9) :: Gradtmp !< 
        real(MK), dimension(3) :: DP      !< 
        ! loop on control points 
        do icp=1,nCPs
            ! loop on Sources
            do is=1,nSrc
                Uindtmp=0.0_MK;
                Gradtmp=0.0_MK;
                DP=CPs(1:3,icp)-SrcP(1:3,is)
                call fUi_PointSource_11(DP,Intensities(is),SmoothModel,SmoothParam(is),bComputeGrad,&
                    Uindtmp,Gradtmp)
                Uind(1:3,icp) = Uind(1:3,icp)+Uindtmp(1:3)
                Grad(1:9,icp) = Grad(1:9,icp)+Gradtmp(1:9)
            enddo 
        enddo
    end subroutine

    subroutine fUi_PointSource_11(DeltaP,Intensity, SmoothModel,SmoothParam,bComputeGrad,Uind, Grad)
        use MathConstants ,only: fourpi_inv
        ! Input/output arguments 
        real(MK), dimension(3), intent(in) :: DeltaP       !< 3 x 1   Pcp-Ps
        real(MK), intent(in)               :: Intensity    !< 
        integer, intent(in)                :: SmoothModel  !< 
        real(MK), intent(in)               :: SmoothParam  !< 
        logical, intent(in)                :: bComputeGrad !< 
        real(MK), dimension(3),intent(out) :: Uind         !< No Side effects
        real(MK), dimension(9),intent(out) :: Grad         !< No Side effects
        ! Variables declaration 
        real(MK) :: r      !< 
        real(MK) :: r3_inv !< 
        real(MK) :: r2_inv !< 
        real(MK) :: E      !< 
        real(MK) :: Ebar   !< 
        real(MK) :: factor   !< 
        real(MK) :: Kv   !< Mollifier

        r=sqrt(DeltaP(1)**2+ DeltaP(2)**2+ DeltaP(3)**2)

        Uind(1:3)=0.0_MK
        Grad(1:9)=0.0_MK
        ! TODO: escape should depend on physics
        if (r<MINNORM) then! to avoid singularity, save a bit of computation and use limit values
            return
            !--------------------------------------------------------------------------------
            !--- Exactly on the Singularity 
            !--------------------------------------------------------------------------------
            !if( SmoothModel==0) then 
            !    return
            !end if
            ! NOTE: We use to return the value obtained analytically for the gradient on the singularity
            ! Now, we just return 0
            !
            !if(bComputeGrad) then
            !    !if( SmoothModel==2) then 
            !    !    r3_inv=1._MK/SmoothParam**3 ! Epsilon^3
            !    !else
            !    !    r3_inv=1._MK/SmoothParam**3 ! Epsilon ^3
            !    !endif
            !    ! For model 1 and 2 is the same limit
            !    !Grad(1) = fourpi_inv*Intensity*r3_inv ! 11
            !    !Grad(5) = fourpi_inv*Intensity*r3_inv ! 22
            !    !Grad(9) = fourpi_inv*Intensity*r3_inv ! 33
            !endif
            ! TODO far distance
        else
            !--------------------------------------------------------------------------------
            !--- Normal Procedure 
            !--------------------------------------------------------------------------------
            ! Velocity 
            select case(SmoothModel)
            case(0)
                r3_inv=1._MK/(r**3);
                Kv=1.0_MK
            case (1) !
                !Exponential mollifier SmoothParam=Epsilon^3
                r3_inv=1._MK/(r**3)
                E=exp(-r**3/SmoothParam**3)
                Ebar=1._MK-E
                Kv=Ebar
            case (2) !
                ! Finite support with SmoothParam=Epsilon^6  !!!
                r3_inv= 1._MK/sqrt(SmoothParam**6+r**6)
                Kv= 1.0_MK
            case default
                print*,'Wrong smooth model for point source'
            end select

            Uind(1:3)=Kv*Intensity*fourpi_inv*(DeltaP)*r3_inv 

            ! Gradient 
            if (bComputeGrad) then 
                select case ( SmoothModel) !
                case (0) !
                    r2_inv=1._MK/r**2
                    Grad(1) = fourpi_inv*r3_inv*Intensity*(1 -3*DeltaP(1)*DeltaP(1)*r2_inv) ! 11
                    Grad(2) = fourpi_inv*r3_inv*Intensity*(  -3*DeltaP(2)*DeltaP(1)*r2_inv) ! 12
                    Grad(3) = fourpi_inv*r3_inv*Intensity*(  -3*DeltaP(3)*DeltaP(1)*r2_inv) ! 13

                    Grad(4) = fourpi_inv*r3_inv*Intensity*(  -3*DeltaP(1)*DeltaP(2)*r2_inv) ! 21
                    Grad(5) = fourpi_inv*r3_inv*Intensity*(1 -3*DeltaP(2)*DeltaP(2)*r2_inv) ! 22
                    Grad(6) = fourpi_inv*r3_inv*Intensity*(  -3*DeltaP(3)*DeltaP(2)*r2_inv) ! 23

                    Grad(7) = fourpi_inv*r3_inv*Intensity*(  -3*DeltaP(1)*DeltaP(3)*r2_inv) ! 31
                    Grad(8) = fourpi_inv*r3_inv*Intensity*(  -3*DeltaP(2)*DeltaP(3)*r2_inv) ! 32
                    Grad(9) = fourpi_inv*r3_inv*Intensity*(1 -3*DeltaP(3)*DeltaP(3)*r2_inv) ! 33

                case (1) !
                    ! Exponential mollifier 
                    r2_inv=1._MK/r**2
                    ! Almost the same as above 
                    ! Remember I chose SmoothParam=epsilon**3, so don't cube it twice!!
                    factor=-3*Ebar*r2_inv+3*r*E/SmoothParam**3
                    Grad(1) = fourpi_inv*r3_inv*Intensity*( Ebar  + factor*DeltaP(1)*DeltaP(1))!11
                    Grad(2) = fourpi_inv*r3_inv*Intensity*(         factor*DeltaP(2)*DeltaP(1))!12
                    Grad(3) = fourpi_inv*r3_inv*Intensity*(         factor*DeltaP(3)*DeltaP(1))!13

                    Grad(4) = fourpi_inv*r3_inv*Intensity*(         factor*DeltaP(1)*DeltaP(2))!21
                    Grad(5) = fourpi_inv*r3_inv*Intensity*( Ebar  + factor*DeltaP(2)*DeltaP(2))!22
                    Grad(6) = fourpi_inv*r3_inv*Intensity*(         factor*DeltaP(3)*DeltaP(2))!23

                    Grad(7) = fourpi_inv*r3_inv*Intensity*(         factor*DeltaP(1)*DeltaP(3))!31
                    Grad(8) = fourpi_inv*r3_inv*Intensity*(         factor*DeltaP(2)*DeltaP(3))!32
                    Grad(9) = fourpi_inv*r3_inv*Intensity*( Ebar  + factor*DeltaP(3)*DeltaP(3))!33

                case (2) !
                    ! Finite support with Epsilon 6 
                    r2_inv= r**4/(SmoothParam**6+r**6) ! Fake r2_inv
                    ! EXACTLY the same as case 0 due to fake variables r3 and r5 
                    Grad(1) = fourpi_inv*r3_inv*Intensity*(1 -3*DeltaP(1)*DeltaP(1)*r2_inv) ! 11
                    Grad(2) = fourpi_inv*r3_inv*Intensity*(  -3*DeltaP(2)*DeltaP(1)*r2_inv) ! 12
                    Grad(3) = fourpi_inv*r3_inv*Intensity*(  -3*DeltaP(3)*DeltaP(1)*r2_inv) ! 13

                    Grad(4) = fourpi_inv*r3_inv*Intensity*(  -3*DeltaP(1)*DeltaP(2)*r2_inv) ! 21
                    Grad(5) = fourpi_inv*r3_inv*Intensity*(1 -3*DeltaP(2)*DeltaP(2)*r2_inv) ! 22
                    Grad(6) = fourpi_inv*r3_inv*Intensity*(  -3*DeltaP(3)*DeltaP(2)*r2_inv) ! 23

                    Grad(7) = fourpi_inv*r3_inv*Intensity*(  -3*DeltaP(1)*DeltaP(3)*r2_inv) ! 31
                    Grad(8) = fourpi_inv*r3_inv*Intensity*(  -3*DeltaP(2)*DeltaP(3)*r2_inv) ! 32
                    Grad(9) = fourpi_inv*r3_inv*Intensity*(1 -3*DeltaP(3)*DeltaP(3)*r2_inv) ! 33
                end select  ! Smooth Model
            end if  ! bComputeGrad
        end if  ! Norm Big enough
    end subroutine
end module
