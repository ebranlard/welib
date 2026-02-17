! This performs the same as Program 1.3 of Lewis (p 502):
!  - Computes the surface vorticity on a 2D ellipse, with or without subdivisions
!
! For coordinate system used by Lewis, see OmnivorDoc and Lewis
program test_ellipse_2d_no_lift
    use FortLibTest
    use SurfaceVorticity2DTools
    implicit none
    !
    integer, parameter :: MM = 20 ! big storage, just for lazy allocation
    real(CDK), dimension(MM,MM) :: MCoeff
    real(CDK), dimension(MM+1) :: XP,YP
    real(CDK), dimension(MM) :: XCP,YCP
    real(CDK), dimension(MM) :: sine, cosine, ds, rhs,  slope
    real(CDK), dimension(MM) :: gammas
    real(CDK), dimension(MM) :: gammas_exact
    real(CDK), dimension(MM) :: vphi
    real(CDK), dimension(MM) :: delx,dely
    real(CDK) :: phi, dphi, Ux,Uy,U0, theta
    real(CDK) :: major ! major axis
    real(CDK) :: ratio ! minor/major
    real(CDK) :: alpha ! angle of attack
    real(CDK) :: f,abyr  ! ellipse param
    integer :: ndiv
    integer :: i,iCase
    integer :: m, lda
    ! Test
    bTestPlot=.false.

    !do iCase =1,2
    do iCase =1,1
        ! --------------------------------------------------------------------------------
        ! --- Parameters, for each case
        ! --------------------------------------------------------------------------------
        U0=1;
        major=1
        if(iCase==1) then
            ! Case 1 
            testname='Jerboa Ellipse2D No Lift, No Subdiv'
            !m=20     ! m=20 works well but that's just lucky
            !m=5     ! m=20 works well but that's just lucky
            m=MM     ! m=20 works well but that's just lucky
            ratio=0.2
            alpha=30.*pi/180.
            ndiv=0
        else
            ! Case 2 
            testname='Jerboa Ellipse2D No Lift,    Subdiv'
            !m=30     
            m=MM
            ratio=0.05
            alpha=10.*pi/180.
            ndiv=2 !<<<<<<<<<<<<<< 2 subdivisions are enough
        endif
        ! --------------------------------------------------------------------------------
        ! ---  
        ! --------------------------------------------------------------------------------
        ! Velocity components
        Ux= U0*cos(alpha)
        Uy= U0*sin(alpha)

        ! --- Panel points around the ellipse (see coordinate system up the script)
        dphi=2*pi/m
        do i=1,m+1
            phi   =(i-1)*dphi
            XP(i) = 0.5*major       *(1.-cos(phi));
            YP(i) = 0.5*major*ratio*  sin(phi);
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

        ! --- RHS (tangential velocity to the panels)
        do i=1,m
            ! Projection of Uinf on the tangent of the panel
            rhs(i)  = -Ux*cosine(i)-Uy*sine(i)
        enddo
        ! --- Building and solving
        lda=size(MCoeff,1)
        call sv2d_coupling_coefficients(m,lda,XCP,YCP,ds,sine,cosine,slope,MCoeff,ndiv,delx,dely)
        print*,'rhs  ',rhs
        !print*,'ge   ',gammas_exact
        print*,''
        print*,'M1',MCoeff(1,:)
        print*,'M2',MCoeff(2,:)
        print*,'M3',MCoeff(3,:)
        print*,'M4',MCoeff(4,:)
        print*,'M5',MCoeff(5,:)
        print*,''
        !
        !call sv2d_back_diagonal_correction(m,ds,MCoeff)
        !
        call invert_by_pivot(m,MCoeff)
        ! 
        call AxV(m,MCoeff,RHS,gammas)

        ! 

        ! --- Analytical velocity for an ellipse
        abyr= sqrt((1.-ratio)/(1.+ratio)); ! see Lewis p 50
        do i=1,m
            phi     = (i-0.5)*dphi
            vphi(i) = phi
            theta   = pi-phi
            f       = sqrt(1.0+ (abyr**4) - 2.*(abyr**2)*cos(2.*theta) )
            gammas_exact(i) =  2.0*U0*sin(theta-alpha)/f  ! See Lewis p50
        enddo
        print*,'ge   ',gammas_exact
        print*,'gn   ',gammas
        print*,''
        if(bTestPlot) then
            !call gnuplot_2curves(vphi,gammas,vphi,gammas_exact+1,1,'phi','gamma','num','exact')
            !call gnuplot_2curves(XCP(1:1),YCP(1:1),XCP(3:m),YCP(3:m),1,'x','y','First','Rest')
            ! CP-like
            !call gnuplot_2curves(XCP,gammas_exact,XCP,gammas,1,'x','gamma','num','')
            ! phi-gammas
            !call gnuplot_2curves(vphi,gammas,vphi,gammas_exact,1,'phi','gamma','num','exact')
        endif
        ! 
        ! 
        call  test_almost_equal('Gammas',gammas(1:m),gammas_exact(1:m),1.0d-1,bTestStop,bTestPrint)
    enddo
end program
