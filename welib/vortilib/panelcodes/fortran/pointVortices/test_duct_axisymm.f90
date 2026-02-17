! This performs the same as Program 4.3 of Lewis (p 524):
!  - Computes the surface vorticity on an axisymmetrix duct
!
! For coordinate system used by Lewis, see OmnivorDoc and Lewis
!
!
program test_duct_axisymm
    use FortLib
    use JerboaLib
    implicit none
    !
    integer, parameter :: MM = 100 ! big storage, just for lazy allocation
    real(CDK), dimension(MM,MM) :: MCoeff
    real(CDK), dimension(MM+1) :: XP,RP
    real(CDK), dimension(MM) :: XCP,RCP
    real(CDK), dimension(MM) :: sine, cosine, ds, rhs,  slope, curve
    real(CDK), dimension(MM) :: gammas
    real(CDK), dimension(MM) :: CP
    !real(CDK), dimension(MM) :: gammas_exact
    real(CDK) :: U0
    real(CDK) :: Gamma_num
    integer :: i,m,m_new,lda
    integer :: te
    testname='Jerboa Duct Axisymmetric,         ,'
    bTestPlot=.true.

    ! --------------------------------------------------------------------------------
    ! --- Parameters
    ! --------------------------------------------------------------------------------
    U0=1;

    ! --------------------------------------------------------------------------------
    ! ---  Panel points, and control points
    ! --------------------------------------------------------------------------------
    call get_naca662015()
    m=40
    te=m/2
    call svas_data_preparation(m,XP,RP,XCP,RCP,ds,sine,cosine,slope,curve,te)
    !call gnuplot_plot(XP,RP,1,'x','r')
    !call gnuplot_2curves(XP,RP,XCP,RCP,1,'x','r','Panels','CPs')

    ! --- RHS (tangential velocity to the panels)
    do i=1,m
        ! Projection of Uinf on the tangent of the panel
        rhs(i)  = -cosine(i)
    enddo
    ! --- Building and solving
    lda=size(MCoeff,1)
    call svas_coupling_coefficients(m,lda,XCP,RCP,ds,sine,cosine,curve,MCoeff)
    !
    call sv2d_back_diagonal_correction(m,lda,ds,MCoeff)
    !
    call sv2d_kutta_condition(m,lda,te,MCoeff,rhs,m_new)
    !
    call invert_by_pivot(m_new,MCoeff)
    ! 
    call AxV(m_new,MCoeff,RHS,gammas)
    ! 
    ! --- Reverting to original matrix size
    call sv2d_revert_kutta(m,te,gammas)

    ! --- Computing full circulation
    Gamma_num=0.0_MK
    do i=1,m
        Gamma_num=Gamma_num+gammas(i)*ds(i)
        Cp(i)    = 1-(gammas(i)/U0)**2
    enddo
    if(bTestPlot) then
        ! CP-like
        !call gnuplot_plot(XCP(1:m),gammas(1:m),1,'x','gamma')
        ! Cp
        call gnuplot_plot(XCP(1:m),-CP(1:m),1,'x','-Cp')
        !call gnuplot_plot(XCP(1:m),gammas(1:m),1,'x','gamma')
        ! phi-gammas
        !call gnuplot_2curves(vphi,gammas,vphi,gammas_exact,1,'phi','gamma','num','exact')
       ! do i=1,m
       !     print*,i,gammas_exact(i),gammas(i),gammas(i)-gammas_exact(i)
       ! enddo
        print*,'Gamma',Gamma_num
    endif
!     ! 
!     ! 
!     ! NOTE: first and last point is a bit off (seems to be consistent with Fig2.6 of Lewis)
!     call  test_almost_equal('Gammas',gammas(2:m-1),gammas_exact(2:m-1),1.2d-1,bTestStop,bTestPrint)
! enddo


contains

    subroutine get_naca662015()
        XP(1)  = 0.000000    ;  RP(1)  =  0.835000
        XP(2)  = 0.006156    ;  RP(2)  =  0.847523
        XP(3)  = 0.024472    ;  RP(3)  =  0.857191
        XP(4)  = 0.054497    ;  RP(4)  =  0.867361
        XP(5)  = 0.095492    ;  RP(5)  =  0.877616
        XP(6)  = 0.146447    ;  RP(6)  =  0.887300
        XP(7)  = 0.206107    ;  RP(7)  =  0.895706
        XP(8)  = 0.273005    ;  RP(8)  =  0.902498
        XP(9)  = 0.345492    ;  RP(9)  =  0.907284
        XP(10) = 0.421783    ;  RP(10) =  0.909725
        XP(11) = 0.500000    ;  RP(11) =  0.909500
        XP(12) = 0.578217    ;  RP(12) =  0.905888
        XP(13) = 0.654509    ;  RP(13) =  0.898088
        XP(14) = 0.726995    ;  RP(14) =  0.885847
        XP(15) = 0.793893    ;  RP(15) =  0.872291
        XP(16) = 0.853553    ;  RP(16) =  0.859551
        XP(17) = 0.904509    ;  RP(17) =  0.849009
        XP(18) = 0.945503    ;  RP(18) =  0.841442
        XP(19) = 0.975528    ;  RP(19) =  0.837324
        XP(20) = 0.993844    ;  RP(20) =  0.835504
        XP(21) = 1.000000    ;  RP(21) =  0.835000
        XP(22) = 0.993844    ;  RP(22) =  0.834496
        XP(23) = 0.975528    ;  RP(23) =  0.832676
        XP(24) = 0.945503    ;  RP(24) =  0.828558
        XP(25) = 0.904509    ;  RP(25) =  0.820991
        XP(26) = 0.853553    ;  RP(26) =  0.810449
        XP(27) = 0.793893    ;  RP(27) =  0.797709
        XP(28) = 0.726995    ;  RP(28) =  0.784153
        XP(29) = 0.654509    ;  RP(29) =  0.771912
        XP(30) = 0.578217    ;  RP(30) =  0.764112
        XP(31) = 0.500000    ;  RP(31) =  0.760500
        XP(32) = 0.421783    ;  RP(32) =  0.760275
        XP(33) = 0.345492    ;  RP(33) =  0.762716
        XP(34) = 0.273005    ;  RP(34) =  0.767502
        XP(35) = 0.206107    ;  RP(35) =  0.774294
        XP(36) = 0.146447    ;  RP(36) =  0.782700
        XP(37) = 0.095492    ;  RP(37) =  0.792384
        XP(38) = 0.054497    ;  RP(38) =  0.802639
        XP(39) = 0.024472    ;  RP(39) =  0.812809
        XP(40) = 0.006156    ;  RP(40) =  0.822477
        XP(41) = 0.000000    ;  RP(41) =  0.835000
    end subroutine
!
! Lewis Figure 4.7
!
!    x         r
!0.000000     0.835000
!0.006156     0.847523
!0.024472     0.857191
!0.054497     0.867361
!0.095492     0.877616
!0.146447     0.887300
!0.206107     0.895706
!0.273005     0.902498
!0.345492     0.907284
!0.421783     0.909725
!0.500000     0.909500
!0.578217     0.905888
!0.654509     0.898088
!0.726995     0.885847
!0.793893     0.872291
!0.853553     0.859551
!0.904509     0.849009
!0.945503     0.841442
!0.975528     0.837324
!0.993844     0.835504
!1.000000     0.835000
!0.993844     0.834496
!0.975528     0.832676
!0.945503     0.828558
!0.904509     0.820991
!0.853553     0.810449
!0.793893     0.797709
!0.726995     0.784153
!0.654509     0.771912
!0.578217     0.764112
!0.500000     0.760500
!0.421783     0.760275
!0.345492     0.762716
!0.273005     0.767502
!0.206107     0.774294
!0.146447     0.782700
!0.095492     0.792384
!0.054497     0.802639
!0.024472     0.812809
!0.006156     0.822477
!0.000000     0.835000



! Lewis table 4.2
!   x         r
!0.000000  0.000000
!0.003074  0.031214
!0.012179  0.061229
!0.026965  0.088891
!0.046863  0.113137
!0.07109   0.133035
!0.098771  0.147821
!0.128786  0.156926
!0.160000  0.160000
!0.189048  0.160000
!0.218095  0.160000
!0.247143  0.160000
!0.276190  0.160000
!0.305238  0.160000
!0.334286  0.160000
!0.363333  0.160000
!0.392381  0.160000
!0.421429  0.160000
!0.450476  0.160000
!0.479524  0.160000
!0.508571  0.160000
!0.537619  0.160000
!0.566667  0.160000
!0.595714  0.160000
!0.624762  0.160000
!0.653809  0.160000 
!0.682857  0.160000
!0.711905  0.160000
!0.740952  0.160000
!0.770000  0.160000
!0.798435  0.152381
!0.826869  0.144762
!0.855304  0.137143
!0.883739  0.129524
!0.912173  0.121905
!0.940608  0.114286
!0.969043  0.106667
!0.997477  0.099048
!1.025912  0.091429
!1.054347  0.083810
!1.082781  0.076190
!1.111216  0.068571
!1.139651  0.060952
!1.168085  0.053333
!1.196520  0.045714
!1.224954  0.038095
!1.253389  0.030476
!1.281824  0.022857
!1.310258  0.015238
!1.338693  0.007619
!1.367128  0.000000








end program
