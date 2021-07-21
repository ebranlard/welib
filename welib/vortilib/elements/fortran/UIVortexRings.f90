!> 
module UIVortexRings
    use PrecisionMod, only: MK
    implicit none

    real(MK),parameter :: MINNORM=1e-4
contains

    !     subroutine fUi_VortexRings(CPs,PRingCenter,R,GammaOriented, SmoothModel,SmoothParam,bComputeGrad,Uind,Grad)
    !         use MathConstants, only: fourpi_inv
    !         ! Arguments declarations 
    !         integer, intent(in) :: nCPs
    !         integer, intent(in) :: nRings
    !         logical, intent(in) :: bComputeGrad !< 
    !         real(MK), dimension(3,nCPs), intent(in)    :: CPs         !< side effects
    !         real(MK), dimension(3,nCPs), intent(inout) :: Uind        !< side effects
    !         real(MK), dimension(9,nCPs), intent(inout) :: Grad        !< side effects
    !         real(MK), dimension(3,nRings), intent(in)  :: PRingCenter !< 
    !         real(MK), dimension(nSrc), intent(in)      :: Intensities !< 
    !         integer, intent(in)                        :: SmoothModel !< 
    !         real(MK), dimension(nSrc), intent(in)      :: SmoothParam !< 
    ! !         ! Variable declarations 
    !     end subroutine
    !
    subroutine fUi_VortexRing_11(CP,PRingCenter,R,GammaOriented, SmoothModel,SmoothParam,bComputeGrad,Uind,Grad)
        ! Arguments declarations 
        real(MK), dimension(3), intent(in)    :: CP           !< 
        real(MK), dimension(3), intent(in)    :: PRingCenter  !< PRingCenter : Position of ring center 
        real(MK), intent(in)                  :: R            !< Ring radius
        real(MK), dimension(3), intent(in)    :: GammaOriented!< Vector norm is ring intensity, direction is normal to the ring
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
        real(MK) :: Gamma  !<  
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
        Gamma=norm2(GammaOriented)
        ez(1:3)=GammaOriented(1:3)/Gamma
        imin=minloc(abs(ez),1)
        ey=ez
        ey(imin)=0
        ey=ey([ 2,3,1 ] )
        ey(imin)=-ey(imin)
        ey=ey/norm2(ey)
        ex=cross(ey,ez)
        !! Now setting up rotation matrix 
        T_r2g(1:3,1)=ex(1:3)
        T_r2g(1:3,2)=ey(1:3)
        T_r2g(1:3,3)=ez(1:3)
        T_g2r=transpose(T_r2g)

        !! Coordinate transformation to a cartesian frame normal to the ring
        CP_loc(1)=T_g2r(1,1)*(CP(1)-PRingCenter(1)) + T_g2r(1,2)*(CP(2)-PRingCenter(2)) + T_g2r(1,3)*(CP(3)-PRingCenter(3))
        CP_loc(2)=T_g2r(2,1)*(CP(1)-PRingCenter(1)) + T_g2r(2,2)*(CP(2)-PRingCenter(2)) + T_g2r(2,3)*(CP(3)-PRingCenter(3))
        CP_loc(3)=T_g2r(3,1)*(CP(1)-PRingCenter(1)) + T_g2r(3,2)*(CP(2)-PRingCenter(2)) + T_g2r(3,3)*(CP(3)-PRingCenter(3))

        !! Coordinate transform to a canonical polar coordinates
        Z_cp=CP_loc(3) ! Since we are already in the ring coordinate system 
        R_cp=sqrt(CP_loc(1)**2+CP_loc(2)**2)
        !print*,'CP_loc',CP_loc

        !! Computing velocity in polar coordinates
        Uindtmp=0.0_MK
        Gradtmp=0.0_MK
        if(SmoothModel==0) then
            call fUi_Ring11_loc(R_cp, Z_cp, Gamma, R, Uindtmp)
            !print*,'Uind pol',Uindtmp
        else
            print*,'Smooth param',SmoothParam ! for unused var compiler warning
            write(*,*)'Smooth Ring not done'
            STOP -1
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

    subroutine test_ring()
        ! Arguments declarations 
        real(MK):: Gamma !< Ring intensity
        real(MK):: R     !< Ring radius
        real(MK), dimension(3) :: Uind  !< side effects, in coord (r, theta, z)
        real(MK), dimension(3) :: CP           !< 
        real(MK), dimension(3):: PRingCenter  !< PRingCenter : Position of ring center 
        real(MK), dimension(3):: GammaOriented!< Vector norm is ring intensity, direction is normal to the ring
        integer:: SmoothModel  !< 
        real(MK):: SmoothParam  !< 
        logical:: bComputeGrad !< 
        real(MK), dimension(9):: Grad         !< side effects
       
        R=3._MK
        Gamma=2._MK
        Uind=0_MK
        call fUi_Ring11_loc(1._MK, 2._MK, Gamma, R, Uind)
        GammaOriented=[2,0,0];
        PRingCenter=[3,2,1];
        CP=[5,2,1]
        CP=CP+[0.,cos(0.1),sin(0.1)];

        SmoothModel=0
        SmoothParam=0.0
        bComputeGrad=.false.
        Uind=[0,0,0]
        Grad(1:9)=0
        call fUi_VortexRing_11(CP,PRingCenter,R,GammaOriented, SmoothModel,SmoothParam,bComputeGrad,Uind,Grad)
    end subroutine


    !> Induction from a ring at a control point, within the referential of a ring: coordinates (r,theta, z), the ring is at z=0
    subroutine fUi_Ring11_loc(Rcp, Zcp, Gamma, R, Uind)
        use EllipticIntegrals, only: elliptic_ke
        use MathConstants, only: fourpi_inv
        ! Arguments declarations 
        real(MK), intent(in) :: Rcp   !< Control point radius
        real(MK), intent(in) :: Zcp   !< Control point z location
        real(MK), intent(in) :: Gamma !< Ring intensity
        real(MK), intent(in) :: R     !< Ring radius
        real(MK), dimension(3), intent(inout) :: Uind  !< side effects, in coord (r, theta, z)
        ! Variable declarations 
        real(MK), parameter :: MIN_DIST=1e-7 ! threshold for using axis formula
        real(MK) :: a  !<  
        real(MK) :: AA !<  
        real(MK) :: B  !<  
        real(MK) :: E  !<  
        real(MK) :: I1  !< 
        real(MK) :: I2  !<
        real(MK) :: K  !<  
        real(MK) :: m  !<  
        !real(MK) :: aL,x,rr,cons,phi,a2  !<  Lewis notations

        if (Rcp<MIN_DIST*R) then
            !! Using Axis formula : v_z=-Gamma/(2R) *1 / (1+(z/R)^2)^(3/2) See Naca 1184 - Watsle, or Glauert Contraction of slipstream of an airscrew 1926 
            Uind(1)=0
            Uind(3)=Gamma/(2_MK*R)*(1_MK/((1_MK+(Zcp/R)**2)**(1.5_MK))) ! no minus sign... it's a matter of convention
        else
            ! Formulation uses Formula from Yoon 2004.. it's not my favorite though.. 
            !Rcp=sqrt(Xcp.**2+Ycp.**2)
            !Zcp=Zcp-z ! !!!! In this case, rememeber to remove z to be in Ring coordinate system
            a=sqrt((Rcp+R)**2+(Zcp)**2)
            m=4_MK*Rcp*R/(a**2)
            AA=(Zcp)**2+Rcp**2+R**2
            B=-2_MK*Rcp*R
            call elliptic_ke(m, K, E) 
            I1=4_MK/a*K
            I2=4_MK/a**3*E/(1_MK-m)
            Uind(1)=( (Zcp)/B )*(I1-AA*I2)      ! vr  !! Watchout GammaR/4pi removed , norm2(r_i) for yoon is R
            Uind(3)=((R+Rcp*AA/B)*I2-Rcp/B*I1 ) ! vz 
            ! General scaling 
            Uind(1:3)=Uind(1:3)*Gamma*fourpi_inv*R

            ! Lewis
!             x =Zcp/R;
!             rr=Rcp/R;
!             aL    = x**2+(rr-1._MK)**2
!             a2    = x**2+(rr+1._MK)**2
!             cons = 0.5_MK/(pi*R*sqrt(x**2+(rr+1._MK)**2))*Gamma
!             m    = 4._MK*rr/a2
!             phi=atan(sqrt(4._MK*rr/aL))
!             if(m>0.999) then
!                 K=log(4./cos(phi))
!                 E=1.+0.5*(K-1.0/1.2) * cos(phi)**2
!             else
!                 call elliptic_ke(m, K, E) 
!             endif
!             Uind(3)= -cons*     (K-(1.+2.*(rr-1.)/aL)*E);
!             Uind(1)=  cons*x/rr*(K-(1.+2.*     rr/aL)*E);
            ! Oye with Grad
!                 zc=Zcp_prime(i)
!                 rc=Rcp_prime(i)
!                 a=sqrt((rc+R)**2+(zc)**2)
!                 m=4*rc*R/a**2
!                 A=(zc)**2+rc**2+R**2
!                 B=-2*rc*R
!                 call ellipke( m , K,E) 
!                 I1=4/a*K
!                 I2=4/a**3*E/(1-m)
!                 uip(i,2)=0
!                 uip(i,1)=Gamma/(4*pi)*R*( (zc)/B )*(I1-A*I2) ! vr hopefully
!                 uip(i,3)=Gamma/(4*pi)*R*((R+rc*A/B)*I2-rc/B*I1 ) ! vz hopefully
!                 if (bComputeGrad) then 
!                     a=(R + rc)**2 + zc**2
!                     b=(R + rc)**2 + zc**2
!                     ! X and Z are the same 
!                     ! Y and R are the same 
!                     Uzz=(1/2).*a.**(-3/2).*b.**(-2).*pi.**(-1).*zc.*(((-7).*R.**4+6.*R.**2.*rc.**2+rc.**4+2.*((-3).*R.**2+rc.**2).*zc.**2+zc.**4).*E+(-1).*b.*((-1).*R.**2+rc.**2+zc.**2).*K)
! 
!                     Uzr=(1/2).*a.**(-3/2).*b.**(-2).*pi.**(-1).*rc.**(-1).*(((R.**2+(-1).*rc.**2).**2.*(R.**2+rc.**2)+2.*(R.**4+(-6).*R.**2.*rc.**2+rc.**4).*zc.**2+(R.**2+rc.**2).*zc.**4).*E+(-1).*b.*((R.**2+(-1).*rc.**2).**2+(R.**2+rc.**2).*zc.**2).*K)
! 
!                     Urr= (-1/2).*a.**(-3/2).*b.**(-2).*pi.**(-1).*rc.**(-2).*zc.*((2.*rc.**6+5.* rc.**4.*(R.**2+zc.**2)+4.*rc.**2.*((-2).*R.**2+zc.**2).*(R.**2+zc.**2)+( R.**2+zc.**2).**3).*E+(-1).*b.*(2.*rc.**4+3.*rc.**2.*((-1).* R+zc).*(R+zc)+(R.**2+zc.**2).**2).*K)
! 
!                     Urz=(1/2).*a.**(-3/2).*b.**(-2).*pi.**(-1).*rc.**(-1).*(((R.**2+(-1).* rc.**2).**2.*(R.**2+rc.**2)+2.*(R.**4+(-6).*R.**2.*rc.**2+rc.**4).*zc.**2+( R.**2+rc.**2).*zc.**4).*E+(-1).*b.*((R.**2+(-1).*rc.**2).**2+ (R.**2+rc.**2).*zc.**2).*K)
! 
!                     grad(i,1)=Urr*Gamma
!                     grad(i,3)=Urz*Gamma
!                     grad(i,7)=Uzr*Gamma
!                     grad(i,9)=Uzz*Gamma
!                 end if 




        endif
    end subroutine fUi_Ring11_loc

    subroutine fui_ring11_loc_c(Rcp, Zcp, Gamma, R, Uind) bind(C,name='fui_ring11_loc_c')
        use PrecisionMod, only: C_DOUBLE
        !DEC$ ATTRIBUTES DLLEXPORT :: fui_ring11_loc_c
        !GCC$ ATTRIBUTES DLLEXPORT :: fui_ring11_loc_c
        ! Arguments declarations 
        real(C_DOUBLE), intent(in) :: Rcp   !< Control point radius
        real(C_DOUBLE), intent(in) :: Zcp   !< Control point z location
        real(C_DOUBLE), intent(in) :: Gamma !< Ring intensity
        real(C_DOUBLE), intent(in) :: R     !< Ring radius
        real(C_DOUBLE), dimension(3), intent(inout) :: Uind  !< side effects, in coord (r, theta, z)
        real(MK), dimension(3) :: Uind_MK 

        Uind_MK=0.0_MK
        call fUi_Ring11_loc(Rcp, Zcp, Gamma, R, Uind_MK)
        Uind=real(Uind_MK,C_DOUBLE)

    end subroutine 

end module
