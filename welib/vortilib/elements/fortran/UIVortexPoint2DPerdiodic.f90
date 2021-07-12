!> 
module VortexPoint2DPeriodic
    use PrecisionMod,   only: MK
    implicit none
    real(MK),parameter :: MINNORM2=1e-8
contains
    subroutine fUi_VortexPoint2DPeriodicX_11(DeltaP,Intensity,Lambda,SmoothModel,SmoothParam,bComputeGrad,Uind, Grad)
        use MathConstants , only: twopi
        implicit none 
        ! Input/output arguments 
        real(MK), dimension(3), intent(in) :: DeltaP       !< 3 x 1   Pcp-Pv
        real(MK), intent(in)               :: Intensity    !<  Circulation Gamma
        real(MK), intent(in)               :: Lambda       !<  Period in the X direction
        integer, intent(in) ,optional      :: SmoothModel  !< 
        real(MK), intent(inout),optional      :: SmoothParam  !< 
        logical, intent(in)                :: bComputeGrad !< 
        real(MK), dimension(3),intent(out) :: Uind         !< No Side effects. Note that
        real(MK), dimension(9),intent(out) :: Grad         !< No Side effects
        !         ! Variables declaration 
        real(MK) :: x     !< 
        real(MK) :: y     !< 
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
        factor = 1._MK
        if(SmoothModel>0 .and. SmoothParam>0) then
            rho2   = r2/SmoothParam**2;
            E      = exp(-rho2);
            Qm     = 1._MK;
            factor = (1._MK-Qm*E);
            !
        else
            SmoothParam=0.0_MK
        endif
        Uind=0.0_MK

        r2  = DeltaP(1)**2+ DeltaP(2)**2
        if (r2<MINNORM2) then! to avoid strong singularity
            Uind(3)=0.0_MK
            Grad(1:9)=0.0_MK 
        else
            x=DeltaP(1)*twopi/Lambda
            y=DeltaP(2)*twopi/Lambda
            Uind(1)= -(Intensity*Sinh(y))/(2.*(-Cos(x) + Cosh(y)+SmoothParam**2))
            Uind(2)=  (Intensity* Sin(x))/(2.*(-Cos(x) + Cosh(y)+SmoothParam**2))
            Uind(3)=0.0_MK

        endif

    end subroutine
            
end module VortexPoint2DPeriodic
