!> 
module VortexPoint2D
    use PrecisionMod,   only: MK
    implicit none
    real(MK),parameter :: MINNORM2=1e-8
contains
    subroutine fUi_VortexPoint2D_11(DeltaP,Intensity,SmoothModel,SmoothParam,bComputeGrad,Uind, Grad)
        use MathConstants , only: twopi_inv
        implicit none 
        ! Input/output arguments 
        real(MK), dimension(3), intent(in) :: DeltaP       !< 3 x 1   Pcp-Pv
        real(MK), intent(in)               :: Intensity    !<  Circulation Gamma
        integer, intent(in) ,optional      :: SmoothModel  !< 
        real(MK), intent(in),optional      :: SmoothParam  !< 
        logical, intent(in)                :: bComputeGrad !< 
        real(MK), dimension(3),intent(out) :: Uind         !< No Side effects. Note that
        real(MK), dimension(9),intent(out) :: Grad         !< No Side effects
        !         ! Variables declaration 
        real(MK) :: r2     !< 
        real(MK) :: factor !< 
        real(MK) :: E      !< 
        real(MK) :: rho2   !< 
        real(MK) :: Qm     !< 
        ! 
        ! No need of gradient in 2D
        if(bComputeGrad) then
            Grad(1:9)=0.0_MK !
        else
            Grad(1:9)=0.0_MK
        endif
        Uind=0.0_MK

        r2  = DeltaP(1)**2+ DeltaP(2)**2
        if (r2<MINNORM2) then! to avoid strong singularity
            Uind(3)=0.0_MK
            Grad(1:9)=0.0_MK 
!         elseif (r2>100) then 
!             Uind(3)=0.0_MK
!             Grad(1:9)=0.0_MK 
        else
            select case (SmoothModel) !
            case (0) ! No mollification, exept on zero
                factor=1._MK/r2
            case (22) !
                ! Finite support 
                factor=1._MK/sqrt(SmoothParam**4+r2**2)
                ! --------------------------------------------------------------------------------
                ! --- Majda Kernels
                ! --------------------------------------------------------------------------------
            case (2) ! Second order - Gaussian - Majada
                rho2   = r2/SmoothParam**2;
                E      = exp(-rho2);
                Qm     = 1._MK;
                factor = (1._MK-Qm*E)/r2;
            case(4)  ! Fourth Order - Gaussian - Majda
                rho2   = r2/SmoothParam**2;
                E      = exp(-rho2);
                Qm     = 1._MK-rho2;
                factor = (1._MK-Qm*E)/r2;
            case(6)  ! Sixth Order - Gaussian - Majda
                rho2   = r2/SmoothParam**2;
                E      = exp(-rho2);
                Qm     = (1._MK-2._MK*rho2+rho2**2/2._MK);
                factor = (1._MK-Qm*E)/r2;
            case(8)  ! Eigth Order - Gaussian - Majda
                rho2   = r2/SmoothParam**2;
                E      = exp(-rho2);
                Qm     = (1._MK-3_MK*rho2+3._MK*rho2**2/2._MK-rho2**3/6._MK);
                factor = (1._MK-Qm*E)/r2;
            case(10)  ! Tenth Order - Gaussian - Mads
                rho2   = r2/SmoothParam**2;
                E      = exp(-rho2/2._MK);  !!!!!!!! Note the 2 in the Gaussian
                Qm     = (1._MK-2._MK*rho2+3./4.*rho2**2-1./12.*rho2**3+1./384.*rho2**4);
                factor = (1._MK-Qm*E)/r2;
                ! --------------------------------------------------------------------------------
                ! --- Mads Kernel (they are the same, it's just the fact that there is a 2 in the exponential that changes things)
                ! --------------------------------------------------------------------------------
            case(102) ! Second order Gaussian - Mads
                rho2   = r2/SmoothParam**2;
                E      = exp(-rho2/2._MK);  !Note the 2 in the Gaussian
                Qm     = 1._MK;
                factor = (1._MK-Qm*E)/r2;
            case(104) ! Fourth order Gaussian - Mads
                rho2   = r2/SmoothParam**2;
                E      = exp(-rho2/2._MK);  !Note the 2 in the Gaussian
                Qm     = 1._MK-1._MK/2._MK*rho2;
                factor = (1._MK-Qm*E)/r2;
            case(106)  ! Sixth Order - Gaussian - Mads
                rho2   = r2/SmoothParam**2;
                E      = exp(-rho2/2._MK);    !Note the 2 in the Gaussian
                Qm     = (1._MK-rho2+1._MK/8._MK*rho2**2);
                factor = (1._MK-Qm*E)/r2;
            case(108)  ! Eigth Order - Gaussian - Mads
                rho2   = r2/SmoothParam**2;
                E      = exp(-rho2/2._MK);                    !Note the 2 in the Gaussian
                Qm     = (1._MK-3._MK/2._MK*rho2+3._MK/8._MK*rho2**2-1._MK/48._MK*rho2**3);
                factor = (1._MK-Qm*E)/r2;
            case(110)  ! Tenth Order - Gaussian - Mads
                rho2   = r2/SmoothParam**2;
                E      = exp(-rho2/2._MK);                               !Note the 2 in the Gaussian
                Qm     = (1._MK-2._MK*rho2+3./4.*rho2**2-1./12.*rho2**3+1./384.*rho2**4);
                factor = (1._MK-Qm*E)/r2;
            case default 
                write(*,*)'2D vortex point: wrong smooth model'
                STOP -1
            end select 
            factor = factor * Intensity * twopi_inv;
            Uind(1)=-DeltaP(2)*factor
            Uind(2)= DeltaP(1)*factor
            Uind(3)=0.0_MK

        endif

    end subroutine

end module VortexPoint2D
