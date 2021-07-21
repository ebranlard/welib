!> 
module HillsVortex
    use  PrecisionMod, only: MK
    implicit none
            
contains
            
    !> Velocity induced by Hills vortex (hence in the frame with no free stream at infinity). 
    ! Induced velocity returned in Cartesian coordinates
    !
    ! Cylindrical coordinate system (rho,phi,z)
    subroutine fUi_HillsVortex_1(CP,a,us,z0,Uind,Grad,Defm,Vort,Strm)
        ! Arguments
        real(MK),dimension(3), intent(in) :: CP !< Control Point
        real(MK), intent(in) :: a  !< hill vortex radius
        real(MK), intent(in) :: us !< hills self induced velocity us=-5/2 U0  (ie convects along -ez)
        real(MK), intent(in) :: z0 !< Position on the z axis
        real(MK),dimension(3),intent(out) :: Uind !< Induced velocity
        real(MK),dimension(9),intent(out) :: Grad !< Gradient
        real(MK),dimension(3),intent(out) :: Defm !< Deformation
        real(MK),dimension(3),intent(out) :: Vort !< Vorticity
        real(MK),dimension(3),intent(out),optional :: Strm !< Stremfunction (optional)
        ! Variables
        real(MK)              :: r3d,rho
        real(MK)              :: x,y,z
        real(MK)              :: urho
        real(MK)              :: uz
        real(MK),dimension(3) :: e_phi  ! Tangential coordinate
        real(MK),dimension(3) :: e_rho  ! Cylindrical radial
        real(MK)              :: om_phi
        real(MK)              :: defm_phi
        real(MK)              :: U0 !< I don't want to think : U0 =-2/5 us
        real(MK)              :: Stokes 
        real(MK)              :: Stokes0
        !
        U0=-2._MK/5._MK*us
        !
        !
        x=CP(1)
        y=CP(2)
        z=CP(3)-z0  ! Note the allowed offset

        r3d=sqrt(x**2 + y**2 + z**2)
        rho=sqrt(x**2 + y**2)
        ! Tangential coordinate (tangent to (y,z))
        if (rho/a<1.e-12) then
            ! This should save us from singularities
            e_phi(1:3)=0.0_MK
            e_rho(1:3)=0.0_MK
        else
            e_phi(1)=-y/rho
            e_phi(2)= x/rho
            e_phi(3)=0
            e_rho(1)= x/rho
            e_rho(2)= y/rho
            e_rho(3)=0
        endif
        if(r3d<a) then 
            ! --- Inside the sphere 
            ! Velocity
            uz     =3._MK/5._MK*us*(1._MK-(2._MK*rho**2+z**2)/(a**2))+us*2._MK/5._MK ;
            urho   =3._MK/5._MK*us*(rho*z)/(a**2);
            ! Vorticity (along e_phi)
            om_phi= 3._MK*us*rho/a**2  != -15/2  * (u0 rho)/a^2 = 3 us rho/a**2
            ! Deformation
            defm_phi= 9._MK/5._MK*us**2/a**2*rho*z ! =45/4 * uo^2/a^2 * rho z
            ! Gradient
            Grad(1:9)=0.0_MK !TODO
            ! Stokes Stream function
            Stokes=-3.0_MK/4.0_MK*U0*rho**2*(1-(rho**2+z**2)/(a**2))
            Stokes0=-rho**2*U0/2.0_MK !< Due to the addition of a free stream term in uz_in
        else 
            ! --- Outside of the sphere 
            ! Velocity
            uz    =2._MK/5._MK*us* (((a**2)/(z**2+rho**2))**(5._MK/2._MK))*(2._MK*z**2-rho**2)/(2._MK*a**2);
            urho  =3._MK/5._MK*us*rho*z/(a**2)*(((a**2)/(z**2+rho**2))**(5._MK/2._MK));
            ! Vorticity
            om_phi   = 0.0_MK
            ! Deformation
            defm_phi = 0.0_MK
            ! Gradient 
            Grad(1:9)=0.0_MK !TODO
            ! Stokes Streamfunction
            Stokes=1.0_MK/2.0_MK*U0*rho**2*(1-a**3/((rho**2+z**2)**1.5_MK))
            Stokes0=-rho**2*U0/2.0_MK !< Due to the addition of a free stream term in uz_in
        endif
        !
        Uind(1) = urho * e_rho(1)
        Uind(2) = urho * e_rho(2)
        Uind(3) = uz
        !
        Defm(1) =  defm_phi * e_phi(1)
        Defm(2) =  defm_phi * e_phi(2)
        Defm(3) = 0.0_MK
        !
        Vort(1) = om_phi * e_phi(1)
        Vort(2) = om_phi * e_phi(2)
        Vort(3) = 0.0_MK
        ! Streamfunction
        if(present(Strm)) then
            if (rho/a>1.e-12) then
                Strm(1) = (Stokes0+Stokes)/rho * e_phi(1)
                Strm(2) = (Stokes0+Stokes)/rho * e_phi(2)
                Strm(3) = 0.0_MK
            else
                Strm=0.0_MK
            endif
        endif

    end subroutine 

end module HillsVortex
