!> 
module UIVortexCylinders
    use PrecisionMod, only: MK
    implicit none

    real(MK),parameter :: MINNORM=1e-4
contains

    !> Semi Infinite vortex cylinder with tangential vorticiy
    subroutine fUi_VortexCylinderTangSemiInf_11(CP,PCenter,R,gamma_t,N, SmoothModel,SmoothParam,bComputeGrad,Uind,Grad)
        ! Arguments declarations 
        real(MK), dimension(3), intent(in)    :: CP           !< 
        real(MK), dimension(3), intent(in)    :: PCenter      !< PCenter : Position of ring center
        real(MK), intent(in)                  :: R            !< Cylinder radius
        real(MK), intent(in)                  :: gamma_t      !< Tangential vorticity Intensity
        real(MK), dimension(3), intent(in)    :: N            !< Vector norm is gamma_t, direction is normal to the
        integer, intent(in)                   :: SmoothModel  !< 
        real(MK), intent(in)                  :: SmoothParam  !< 
        logical, intent(in)                   :: bComputeGrad !< 
        real(MK), dimension(3), intent(inout) :: Uind         !< side effects
        real(MK), dimension(9), intent(inout) :: Grad         !< side effects
        ! Variable declarations 
        real(MK), dimension(3) :: ex  !<  
        real(MK), dimension(3) :: ey  !< 
        real(MK), dimension(3) :: ez  !<  
        real(MK), dimension(3,3) :: T_r2g
        real(MK), dimension(3,3) :: T_g2r
        integer :: imin
        real(MK), dimension(3) :: CP_loc  !< 
        real(MK) :: R_cp
        real(MK) :: Z_cp
        real(MK) :: phi
        real(MK), dimension(3) :: Uindtmp         !< 
        real(MK), dimension(9) :: Gradtmp         !<
        real(MK) :: urad
        real(MK) :: upsi
        !         ! Note: theories return only vr (index 1) and vz (index 3), there is some cartesian handling though, but does it work ? Is it consistent with the rest of the codes. Does it make use of the inputed vpsi to reduce the error in the angle estimate. 
        ! 
        !! Setting up an orthonornal cartesian system normal to the ring
        ez(1:3)=N(1:3)/norm2(N) ! just in case N is not unitary
        imin=minloc(abs(ez),1)
        ey=ez
        ey(imin)=0
        ey=ey([ 2,3,1 ] )
        ey(imin)=-ey(imin)
        ey=ey/norm2(ey)
        ex=cross(ey,ez)
!         print*,ex
!         print*,ey
!         print*,ez
        !! Now setting up rotation matrix 
        T_r2g(1:3,1)=ex(1:3)
        T_r2g(1:3,2)=ey(1:3)
        T_r2g(1:3,3)=ez(1:3)
        T_g2r=transpose(T_r2g)

        !! Coordinate transformation to a cartesian frame normal to the ring
        CP_loc(1)=T_g2r(1,1)*(CP(1)-PCenter(1)) + T_g2r(1,2)*(CP(2)-PCenter(2)) + T_g2r(1,3)*(CP(3)-PCenter(3))
        CP_loc(2)=T_g2r(2,1)*(CP(1)-PCenter(1)) + T_g2r(2,2)*(CP(2)-PCenter(2)) + T_g2r(2,3)*(CP(3)-PCenter(3))
        CP_loc(3)=T_g2r(3,1)*(CP(1)-PCenter(1)) + T_g2r(3,2)*(CP(2)-PCenter(2)) + T_g2r(3,3)*(CP(3)-PCenter(3))

        !! Coordinate transform to a canonical polar coordinates
        Z_cp=CP_loc(3) ! Since we are already in the ring coordinate system 
        R_cp=sqrt(CP_loc(1)**2+CP_loc(2)**2)
        !print*,'CP_loc',CP_loc

        !! Computing velocity in polar coordinates
        Uindtmp=0.0_MK
        Gradtmp=0.0_MK
        !print*,'gamma_t',gamma_t
        ! We use the Smooth Param of the cylinder no matter the SmoothModel
        if(SmoothModel==0) then
            call fUi_CylinderTangSemiInf11_loc(R_cp, Z_cp, gamma_t, R, SmoothParam, Uindtmp)
            !print*,'Uind pol',Uindtmp
        else
            call fUi_CylinderTangSemiInf11_loc(R_cp, Z_cp, gamma_t, R, SmoothParam, Uindtmp)
        endif

        !! Converting to cartesian 
        phi=atan2(CP_loc(2),CP_loc(1)) 
        !print*,'phi',phi*180/3.141593
        urad=Uindtmp(1)
        upsi=0.0_MK
        Uindtmp(1)=+urad*cos(phi)
        Uindtmp(2)=+urad*sin(phi)
        !print*,'Uind cart',Uindtmp

        !! Going back to global frame 
        Uind(1)=Uind(1) + T_r2g(1,1)*Uindtmp(1) + T_r2g(1,2)*Uindtmp(2) + T_r2g(1,3)*Uindtmp(3)
        Uind(2)=Uind(2) + T_r2g(2,1)*Uindtmp(1) + T_r2g(2,2)*Uindtmp(2) + T_r2g(2,3)*Uindtmp(3)
        Uind(3)=Uind(3) + T_r2g(3,1)*Uindtmp(1) + T_r2g(3,2)*Uindtmp(2) + T_r2g(3,3)*Uindtmp(3)

        !print*,'Uind glob',Uind
        !! TODO gradient 
        if (bComputeGrad) then 
            Grad = 0.0D0  
        else
            Grad = 0.0D0 
        end if 
    contains
        !> Cross product of two vectors ! See crossprod in UTILS...
        function cross(a, b)
            real(MK), dimension(3) :: cross
            real(MK), dimension(3), intent(in) :: a, b
            cross(1) = a(2) * b(3) - a(3) * b(2)
            cross(2) = a(3) * b(1) - a(1) * b(3)
            cross(3) = a(1) * b(2) - a(2) * b(1)
        end function cross
    end subroutine

    !> Induction from a Cylinder at a control point, within the referential of a Cylinder: coordinates (r,theta, z), the ring is at z=0
    subroutine fUi_CylinderTangSemiInf11_loc(Rcp, Zcp, gamma_t, R, SmoothParam, Uind)
        use EllipticIntegrals, only: elliptic_ke, elliptic_pi
        use MathConstants, only:  pi
        ! Arguments declarations 
        real(MK), intent(in) :: Rcp   !< Control point radius
        real(MK), intent(in) :: Zcp   !< Control point z location
        real(MK), intent(in) :: gamma_t !< Cylinder intensity
        real(MK), intent(in) :: R     !< Cylinder radius
        real(MK), intent(in) :: SmoothParam  !< 
        real(MK), dimension(3), intent(inout) :: Uind  !< side effects, in coord (r, theta, z)
        ! Variable declarations 
        real(MK), parameter :: MIN_DIST=1e-7 ! threshold for using axis formula
        real(MK) :: KK  !<  
        real(MK) :: EE  !<  
        real(MK) :: PPI  !<  
        real(MK) :: k  !<  
        real(MK) :: k_2  !<  
        real(MK) :: k0_2  !<  

        if (rcp<MIN_DIST*R) then
            !! Using Axis formula : v_z=-Gamma/(2R) *1 / (1+(z/R)^2)^(3/2) See Naca 1184 - Watsle, or Glauert Contraction of slipstream of an airscrew 1926 
            Uind(1)=0
            Uind(3)=gamma_t/(2_MK)*(1_MK + Zcp/sqrt(Zcp**2+R**2))
        else
            ! Eliptic integrals 
            k_2 = 4_MK*rcp*R/((R+rcp)**2+zcp**2+SmoothParam**2)
            k = sqrt(k_2)
            k0_2 = 4_MK*rcp*R/((R+rcp)**2+SmoothParam**2)
            call elliptic_ke( k_2 , KK,EE) 
            call elliptic_pi(k0_2,90._MK,k_2,PPI)
            ! ur 
            Uind(1)=-gamma_t/(2_MK*pi)*sqrt(R/rcp)*( (2_MK-k_2)/k*KK - 2_MK/k*EE )
            ! ur 
            Uind(3)= gamma_t/2_MK*( 1._MK/2._MK*(1+(R-rcp)/((rcp+R)*sqrt(1-k0_2))) &
                +zcp*k/(2_MK*pi*sqrt(rcp*R))*(KK + (R-rcp)/(R+rcp)*PPI ) )
        endif
    end subroutine fUi_CylinderTangSemiInf11_loc

    subroutine fui_cylinder11_loc_c(Rcp, Zcp, gamma_t, R,SmoothParam, Uind) bind(C,name='fui_cylinder11_loc_c')
        use PrecisionMod, only: C_DOUBLE
        !DEC$ ATTRIBUTES DLLEXPORT :: fui_cylinder11_loc_c
        !GCC$ ATTRIBUTES DLLEXPORT :: fui_cylinder11_loc_c
        ! Arguments declarations 
        real(C_DOUBLE), intent(in) :: Rcp   !< Control point radius
        real(C_DOUBLE), intent(in) :: Zcp   !< Control point z location
        real(C_DOUBLE), intent(in) :: gamma_t !< cylinder intensity
        real(C_DOUBLE), intent(in) :: R     !< cylinder radius
        real(C_DOUBLE), intent(in) :: SmoothParam     !< cylinder radius
        real(C_DOUBLE), dimension(3), intent(inout) :: Uind  !< side effects, in coord (r, theta, z)
        real(MK), dimension(3) :: Uind_MK 

        Uind_MK=0.0_MK
        call fUi_CylinderTangSemiInf11_loc(Rcp, Zcp, gamma_t, R, SmoothParam, Uind_MK)
        Uind=real(Uind_MK,C_DOUBLE)

    end subroutine 

end module
