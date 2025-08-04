! This performs the same as Program 2.2 of Lewis (p 510):
!  - Computes the surface vorticity of a lifting 2D ellipse with prescribed circulation
!   (Comparison with table  2.3 or figure 2.6 of Lewis)
!
! For coordinate system used by Lewis, see OmnivorDoc and Lewis
program test_ellipse2d_lift
    use FortLibTest
    use SurfaceVorticity2DTools
    implicit none
    !
    integer, parameter :: MM = 30 ! big storage, just for lazy allocation
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
    real(CDK) :: ratio ! minor/major [lambda in Lewis notations)
    real(CDK) :: alpha ! angle of attack
    real(CDK) :: f,abyr  ! ellipse param
    real(CDK) :: Gamma_th   ! Circulation around the ellipse
    real(CDK) :: Gamma_num  ! Circulation around the ellipse
    integer :: ndiv
    integer :: i,iCase
    integer :: m, te, m_new, lda
    ! Test
    bTestPlot=.false.

    !do iCase =1,2
    do iCase =2,2
        ! --------------------------------------------------------------------------------
        ! --- Parameters, for each case
        ! --------------------------------------------------------------------------------
        U0=1;
        major=1
        !
        if(iCase==1) then
            ! Case 1 (Comparison with table  2.3 or figure 2.6 of Lewis)
            testname='Jerboa Ellipse2D    Lift, Prescrib.'
            m=MM
            ratio=0.05
            alpha=10.*pi/180.
            ndiv=1
        else
            testname='Jerboa Ellipse2D    Lift, KuttaCond'
            !m=30     
            m=MM
            ratio=0.05
            alpha=10.*pi/180.
            ndiv=1 
            te=m/2 !< Index of the trailing edge
            !print*,'Gamma',Gamma
        endif
        Gamma_th=pi*(1+ratio)*major*sin(alpha)
        Gamma_num=0.0
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
            if(iCase==1) then
                rhs(i)  = -Ux*cosine(i)-Uy*sine(i)+Gamma_th !<<<<<<<<<<<<<<<< Note Prescribing gamma
            else
                rhs(i)  = -Ux*cosine(i)-Uy*sine(i)
            endif
        enddo
        ! --- Building and solving
        lda=size(MCoeff,1)
        call sv2d_coupling_coefficients(m,lda,XCP,YCP,ds,sine,cosine,slope,MCoeff,ndiv,delx,dely)
        print*,''
        print*,'rhs  ',rhs
        print*,''
        print*,'M1',MCoeff(1,:)
        print*,'M2',MCoeff(2,:)
        print*,'M3',MCoeff(3,:)
        print*,'M4',MCoeff(4,:)
        print*,'M5',MCoeff(5,:)
        print*,'M6',MCoeff(6,:)
        print*,''
        !
        print*,'>>> Back correction <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'
        call sv2d_back_diagonal_correction(m,lda,ds,MCoeff)
        !
        print*,''
        print*,'M1',MCoeff(1,:)
        print*,'M2',MCoeff(2,:)
        print*,'M3',MCoeff(3,:)
        print*,'M4',MCoeff(4,:)
        print*,'M5',MCoeff(5,:)
        print*,'M6',MCoeff(6,:)
        print*,''
        if(iCase==1) then
            call sv2d_bound_vortex_correction(m,lda,ds,MCoeff)
            m_new=m
        else
            ! One equation is removed
            call sv2d_kutta_condition(m,lda,te,MCoeff,rhs,m_new)
        endif
        
        print*,''
        print*,'M1',MCoeff(1,:)
        print*,'M2',MCoeff(2,:)
        print*,'M3',MCoeff(3,:)
        print*,'M4',MCoeff(4,:)
        print*,'M5',MCoeff(5,:)
        print*,'M6',MCoeff(6,:)
        print*,''
        print*,'rhs  ',rhs
        print*,''
        !
        call invert_by_pivot(m_new,MCoeff)
        ! 
        call AxV(m_new,MCoeff,RHS,gammas)
        ! 
        if(iCase==2) then
           print*,''
           print*,'gn (bef)  ',gammas
            ! --- Reverting to original matrix size
            call sv2d_revert_kutta(m,te,gammas)
           print*,'gn (aft)  ',gammas
           print*,''
        endif

        ! --- Computing full circulation
        do i=1,m
            Gamma_num=Gamma_num+gammas(i)*ds(i)
        enddo
        ! --- Analytical velocity for an ellipse
        abyr= sqrt((1.-ratio)/(1.+ratio)); ! see Lewis p 50
        do i=1,m
            phi     = (i-0.5)*dphi
            vphi(i) = phi
            theta   = pi-phi
            f       = sqrt(1.0+ (abyr**4) - 2.*(abyr**2)*cos(2.*theta) )
            gammas_exact(i) =  (2.0*U0*sin(theta-alpha)+2*sin(alpha) )/f  ! See Lewis p50
        enddo
        if(bTestPlot) then
            !call gnuplot_2curves(vphi,gammas,vphi,gammas_exact+1,1,'phi','gamma','num','exact')
            !call gnuplot_2curves(XCP(1:1),YCP(1:1),XCP(3:m),YCP(3:m),1,'x','y','First','Rest')
            ! CP-like
            !call gnuplot_2curves(XCP,gammas_exact,XCP,gammas,1,'x','gamma','theory','num')
            ! phi-gammas
            !call gnuplot_2curves(vphi,gammas,vphi,gammas_exact,1,'phi','gamma','num','exact')
        endif
        print*,'ge   ',gammas_exact
        print*,'gn   ',gammas
        print*,''
        ! 
         do i=1,m
             print*,i,gammas_exact(i),gammas(i),gammas(i)-gammas_exact(i)
         enddo
         print*,'Gamma',Gamma_th,Gamma_num
        ! 
        ! NOTE: first and last point is a bit off (seems to be consistent with Fig2.6 of Lewis)
        call  test_almost_equal('Gammas',gammas(2:m-1),gammas_exact(2:m-1),1.2d-1,bTestStop,bTestPrint)
    enddo
end program
