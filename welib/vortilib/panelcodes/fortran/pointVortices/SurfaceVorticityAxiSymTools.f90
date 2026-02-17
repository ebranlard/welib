!> 
module SurfaceVorticityAxiSymTools
    use PrecisionMod, only: CDK, CIK, MK
    implicit none
            
contains
            
    subroutine svas_data_preparation(m,XP,RP,XCP,RCP,ds,sine,cosine,slope,curve,te)&
            bind(C,name='svas_data_preparation')
        ! Arguments
        integer(CIK), intent(in) :: m
        real(CDK), dimension(*), intent(in) :: XP,RP   !< Panel coordinates (X and R), size m+1!!!
        real(CDK), dimension(*), intent(out)  :: XCP,RCP !< Control Point coordinatse (X and R)
        real(CDK), dimension(*), intent(out) :: ds       !< Panel Length
        real(CDK), dimension(*), intent(out) :: sine     !
        real(CDK), dimension(*), intent(out) :: cosine
        real(CDK), dimension(*), intent(out) :: slope
        real(CDK), dimension(*), intent(out) :: curve
        integer(CIK), intent(in) :: te !< number of subdivisions
        ! Variables
        real(CDK) :: x1,x2,r1,r2,abscos,t
        integer(CIK) :: i
        real(CDK), parameter :: pi=acos(-1._CDK)
        real(CDK), parameter :: ex=0.000001 ! TODO, precision


        do i=1,m
            x1 = XP(i)
            r1 = RP(i)
            x2 = XP(i+1)
            r2 = RP(i+1)
            ! Control Points Between panels
            XCP(i) = (x1+x2)/2
            RCP(i) = (r1+r2)/2
            !
            ds(i)     = sqrt((x2-x1)**2+(r2-r1)**2)
            sine(i)   = (r2-r1)/ds(i);
            cosine(i) = (x2-x1)/ds(i);
            abscos    = abs(cosine(i))
            ! Slope
            if(abscos>ex) then
                t=atan(sine(i)/cosine(i))
                if(cosine(i)>ex) then
                    slope(i) =t
                else ! cosine<-ex
                    ! <<<< Here small difference compared to sv2d
                    if(i>te) then
                        slope(i) =t-pi
                    else
                        slope(i) =t+pi
                    endif
                endif
            else
                slope(i)=sine(i) / abs(sine(i))*pi/2
            endif
        enddo
        ! <<<< Here Computing curve ompared to sv2d
        do i=2,m-1
            curve(i) = (slope(i+1)-slope(i-1)) / (8*pi)
        enddo
        curve(1)   = (slope(2)-slope(m)  -2*pi) / (8*pi)
        curve(m)   = (slope(1)-slope(m-1)-2*pi) / (8*pi)
        curve(te)  = 0
        curve(te+1)= 0

    end subroutine




    !> Coupling coefficients computing in a fashion similar to Lewis: 
    subroutine svas_coupling_coefficients(m,lda,XCP,RCP,ds,sine,cosine,curve,coup)&
            bind(C,name='svas_coupling_coefficients')
        ! Arguments
        integer(CIK), intent(in) :: m
        integer(CIK), intent(in) :: lda
        real(CDK), dimension(*), intent(in)    :: XCP,RCP    !< Control Point coordinatse
        real(CDK), dimension(*), intent(in)    :: ds         !< Panel Length
        real(CDK), dimension(*), intent(in)    :: sine     !
        real(CDK), dimension(*), intent(in)    :: cosine
        real(CDK), dimension(*), intent(in)    :: curve
        real(CDK), dimension(lda,*), intent(out) :: coup
        ! Variables
        real(CDK) :: ux,ur
        real(CDK), parameter :: pi      = acos(-1._CDK)
        integer(CIK) :: i,j
        real(CDK) :: cons !

        ! --- Self induced coefficients (diagonal) using curvatue correction from Lewis
        do i=1,m
            cons=4*pi*RCP(i)/ds(i)
            coup(i,i) = -0.5_CDK-(log(2*cons) - 0.25_CDK)/cons*cosine(i)-curve(i) ! Lewis equation 4.22
        enddo
        ! --- Coupling coefficients off diagonal
        do i=1,m
            do j=1,m
                if(j/=i) then
                    call svas_ui_n1(XCP(j),RCP(j),XCP(i),RCP(i),ux,ur)
                    ! Velocity (unit point vortex from ring i at CP j)
                    coup(j,i) = (ux*cosine(j) + ur*sine(j)) * ds(i)
                endif
            enddo
        enddo

    end subroutine

    !> cf routine "velocities" of Lewis p 522
    subroutine svas_ui_n1(x1,r1,xm,rm,ux,ur) bind(C,name='svas_ui_n1')
        use UIVortexRings, only: fUi_Ring11_loc
        !use EllipticIntegrals, only: elliptic_ke
        !use MathConstants, only: pi
        !
        real(CDK), intent(in) :: x1, r1 !< Control points
        real(CDK), intent(in) :: xm, rm !< Ring parameters
        real(CDK), intent(out) :: ux,ur
        !
        real(MK), dimension(3) :: Uind
        !real(MK) :: x, r,ux0,ur0
        !real(MK) :: aL, cons, fi1, K1,E1, k,phi,m,a2,K1as,E1as

        ! Induced velocity due to a unitary ring of radius 1 
        call fUi_Ring11_loc(r1, x1-xm, 1.0_MK, rm, Uind)
        !
        ux= - real(Uind(3),CDK); ! NOTE MYSTERIOUS MINUS SIGN
        ur= - real(Uind(1),CDK);
        !
        ! Lewis:
        ! Scaling as funciton of radius (radius becomes 1)
        !x=real((x1-xm)/rm,MK)
        !r=real(     r1/rm,MK)
        !aL    = x**2+(r-1._MK)**2
        !a2    = x**2+(r+1._MK)**2
        !cons = 0.5_MK/ (pi*rm*sqrt(x**2+(r+1._MK)**2))
        !m    = 4._MK*r/a2
        !phi=atan(sqrt(4._MK*r/aL))

        !if(m>1) then
        !    K1=log(4./cos(phi))
        !    E1=1.+0.5*(K1-1.0/1.2) * cos(phi)**2
        !else
        !    call elliptic_ke(m, K1, E1) 
        !endif
        !ux=real(-cons*    (K1-(1.+2.*(r-1.)/aL)*E1),CDK);
        !ur=real( cons*x/r*(K1-(1.+2.*r/aL)     *E1),CDK);

        !print'(A,2F12.5)','a',ux,ur
        !print'(A,2F12.5)','b',ux0,ur0
    end subroutine

end module SurfaceVorticityAxiSymTools
