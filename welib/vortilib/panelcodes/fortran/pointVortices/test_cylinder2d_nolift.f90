!-  Coordinate system and panelling when starting from LE: (see figure 1.7 p21 of Lewis)
! .>phi    ^ y    <.   + theta 
! i        |        j 
!          |
!LE        |----TE---> x
!-  Coordinate system and panelling when starting from TE (and using theta negative)
!          ^ y    <.   + theta 
!          |        j 
!          |
!LE        |----TE---> x
!
!
!
!
! This performs the same as Program 1.1 of Lewis (p 499):
!  - Computes the surface vorticity on a 2D circular cylinder
!
! For coordinate system used by Lewis, see OmnivorDoc and Lewis
!
!
program test_cylinder2d_nolift
    use FortLibTest
    use SurfaceVorticity2DTools
    implicit none
    !
!     integer, parameter :: m = 19
    integer, parameter :: m = 5
    real(CDK), dimension(m,m) :: MCoeff
    real(CDK), dimension(m+1) :: XP,YP
    real(CDK), dimension(m) :: XCP,YCP
    real(CDK), dimension(m) :: sine, cosine, ds, rhs,  slope
    real(CDK), dimension(m) :: gammas
    real(CDK), dimension(m) :: gammas_exact
    real(CDK), dimension(m) :: vphi
    real(CDK), dimension(m) :: delx,dely
    real(CDK) :: R, phi, dphi, Ux
    integer :: pCoordSystem
    integer :: i
    integer :: ndiv, lda
    ! Test
    testname='Jerboa Cylind.2D No Lift, surf.vort'
    bTestPlot=.false.

    ! --------------------------------------------------------------------------------
    ! --- Parameters
    ! --------------------------------------------------------------------------------
    Ux=1;
    R = 1
    pCoordSystem=1 !< 1= Lewis, 2= Mine
    ndiv=1
    !pi=acos(-1._CDK)

    ! --------------------------------------------------------------------------------
    ! --- Panel points around the cylinder
    ! --------------------------------------------------------------------------------
    do i=1,m+1
        dphi=2*pi/m
        if(pCoordSystem==1) then
            ! Starts at LE, go suction side, and include the LE twice
            phi   = (i-1)*dphi
            XP(i) = R*(1- cos(phi));
            YP(i) = R*(   sin(phi));
        elseif(pCoordSystem==3) then
            ! Starts at LE, go suction side, and include the TE twice
            phi =-(i-1)*dphi
            XP(i) = R*(1-cos(phi));
            YP(i) = R*  sin(phi);
        else
            ! Starts at TE, go pressure side, and include the TE twice
            phi =-(i-1)*dphi
            XP(i) = R*cos(phi);
            YP(i) = R*sin(phi);
        endif
    enddo

    call sv2d_data_preparation(m,XP,YP,XCP,YCP,ds,sine,cosine,slope,ndiv,delx,dely)

    print*,'XP',XP
    print*,'YP',YP
    print*,'ds   ',ds
    print*,'sin  ',sine
    print*,'cos  ',cosine
    print*,'slope',slope
    print*,'slope',maxval(slope)
    print*,'slope',minval(slope)
    print*,'xcp  ',xcp
    print*,'ycp  ',ycp

!     sine=-sine
!     slope=-slope

    ! RHS
    do i=1,m
        !print'(A,4F8.2)','P,CP',XP(i),YP(i),XCP(i),YCP(i)
        if(pCoordSystem==1) then
            rhs(i)          = -Ux*cosine(i)
            phi             =  (i-0.5)*dphi
            gammas_exact(i) =  2.0*Ux*sin(phi)
        elseif(pCoordSystem==3) then
            rhs(i)          = -Ux*cosine(i)
            phi             = -(i-0.5)*dphi
            gammas_exact(i) = -2.0*Ux*sin(phi)
        else
            rhs(i)          =  Ux*cosine(i)
            phi             = -(i-0.5)*dphi
            gammas_exact(i) = -2.0*Ux*sin(phi)
        endif
        vphi(i)=phi
    enddo
    lda=size(MCoeff,1)
    call sv2d_coupling_coefficients(m,lda,XCP,YCP,ds,sine,cosine,slope,MCoeff,ndiv,delx,dely)
    print*,'rhs  ',rhs
    print*,'ge   ',gammas_exact
    print*,''
    print*,'M1',MCoeff(1,:)
    print*,'M2',MCoeff(2,:)
    print*,'M3',MCoeff(3,:)
    print*,'M4',MCoeff(4,:)
    print*,''
    !
    call invert_by_pivot(m,MCoeff)

    call AxV(m,MCoeff,RHS,gammas)


    !if(bTestPlot) then
    !    call gnuplot_2curves(vphi,gammas,vphi,gammas_exact,1,'phi','gamma','num','exact')
    !endif
    print*,'ge   ',gammas_exact
    print*,'gn   ',gammas
    print*,''

    call  test_almost_equal('Gammas',gammas,gammas_exact,1.d-2,bTestStop,bTestPrint)
    print*,''
end program
