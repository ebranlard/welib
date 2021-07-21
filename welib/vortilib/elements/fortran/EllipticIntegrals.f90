!> 
module EllipticIntegrals
    use PrecisionMod, only: MK
    implicit none

    ! For Elliptic pi
    real(MK), dimension(10),parameter :: T=(/.9931285991850949_MK,.9639719272779138_MK,&
              .9122344282513259_MK,.8391169718222188_MK,&
              .7463319064601508_MK,.6360536807265150_MK,&
              .5108670019508271_MK,.3737060887154195_MK,&
              .2277858511416451_MK,.07652652113349734_MK/)
    real(MK), dimension(10),parameter :: W=(/.01761400713915212_MK,.04060142980038694_MK,&
              .06267204833410907_MK,.08327674157670475_MK,&
              .1019301198172404_MK,.1181945319615184_MK,&
              .1316886384491766_MK,.1420961093183820_MK,&
              .1491729864726037_MK,.1527533871307258_MK/)
contains


    !> Complete elliptic integrals of first and second kind
    ! Just like matlab function ellipke
    subroutine elliptic_ke(m,K,E)
        use PrecisionMod, only: precision_equal
        !  m=k^2          K              E
        ! ---------------------------------
        !  .00      1.570796      1.570796
        !  .50      1.854074      1.350643
        ! 1.00       ì            1.000000
        real(MK), intent(in) :: m  !< m=k^2 
        real(MK), intent(out) :: K, E
        ! Variables
        real(MK) :: pk ,ak, bk, ae, be
        pk=1.0_MK-m

        !if(m>0.99999_MK) then
        !    ! Using asymptotic expression (cant, phi is not the one from m)
        !    K=log(4./cos(phi))
        !    E=1.+0.5*(K-1.0/1.2) * cos(phi)**2
        !endif
        if (precision_equal(m,1._MK)) then
            k=1.0_MK+300
            e=1.0_MK
        else
            !
            ak=(((.01451196212_MK*pk+.03742563713_MK)*pk +.03590092383_MK)*pk+.09666344259_MK)*pk+ 1.38629436112_MK
            bk=(((.00441787012_MK*pk+.03328355346_MK)*pk +.06880248576_MK)*pk+.12498593597_MK)*pk+.5_MK
            k=ak-bk*log(pk)

            ae=(((.01736506451_MK*pk+.04757383546_MK)*pk+.0626060122_MK)*pk+.44325141463_MK)*pk+1.0_MK
            be=(((.00526449639_MK*pk+.04069697526_MK)*pk+.09200180037_MK)*pk+.2499836831_MK)*pk
            e=ae-be*log(pk)
        endif
    end subroutine 

    !> Compute the elliptic integral of the third kind using Gauss-Legendre quadrature
    ! Just like matlab function ellipticPi(n,phi,m) BUT PHI IN DEGREES
    ! Use phi=90 for complete elliptic integral Pi
    subroutine elliptic_pi(n,phi,m,EL3)
        use PrecisionMod, only: precision_equal
        real(MK),intent(in) :: n !< parameter [0 1]
        real(MK),intent(in) :: phi !< argument in degrees
        real(MK),intent(in) :: m !< modulus [0 1] 
        real(MK),intent(out) :: EL3   !< Elliptic integral value Pi
        integer :: i
        real(MK) :: c0,c1,c2,t1,t2,f1,f2
        real(MK) :: k,c
        logical :: lb1,lb2
        k   =sqrt(m)
        c=n

        lb1= precision_equal(k,1.0_MK)  .and. abs(phi-90.0).le.1.0d-8
        lb2= precision_equal(c,1.0_MK)  .and. abs(phi-90.0).le.1.0d-8
        if (lb1.or.lb2) then
            el3=1.0_MK+300
            return
        endif
        c1=0.87266462599716d-2*phi
        c2=c1
        el3=0.0_MK
        do i=1,10
           c0=c2*T(i)
           t1=c1+c0
           t2=c1-c0
           f1=1.0_MK/((1.0_MK-c*sin(t1)*sin(t1))*sqrt(1.0_MK-k*k*sin(t1)*sin(t1)))
           f2=1.0_MK/((1.0_MK-c*sin(t2)*sin(t2))*sqrt(1.0_MK-k*k*sin(t2)*sin(t2)))
           el3=el3+W(i)*(f1+f2)
        enddo
        el3=c1*el3
    end subroutine

    subroutine elliptic_ke_c(m,k,e) bind(C,name='elliptic_ke')
        use PrecisionMod, only: C_DOUBLE
        !DEC$ ATTRIBUTES DLLEXPORT :: elliptic_ke_c
        !GCC$ ATTRIBUTES DLLEXPORT :: elliptic_ke_c
        real(C_DOUBLE),intent(in)   :: m !< 
        real(C_DOUBLE), intent(out) :: k !< 
        real(C_DOUBLE), intent(out) :: e !< 
        call elliptic_ke(m,k,e)
    end subroutine 

    subroutine elliptic_pi_c(n,phi,m,pi) bind(C,name='elliptic_pi')
        use PrecisionMod, only: C_DOUBLE
		!DEC$ ATTRIBUTES DLLEXPORT :: elliptic_pi_c
	    !GCC$ ATTRIBUTES DLLEXPORT :: elliptic_pi_c
        real(C_DOUBLE),intent(in)    :: n     !< 
        real(C_DOUBLE),intent(in)    :: phi     !< 
        real(C_DOUBLE),intent(in)    :: m     !< 
        real(C_DOUBLE), intent(out)    :: pi      !< 
        call elliptic_pi(n,phi,m,pi)
    end subroutine 
end module EllipticIntegrals
