!> 
module SurfaceVorticity2DTools
    use PrecisionMod

contains

    subroutine invert_by_pivot(m,Mat)
        integer(CIK), intent(in) :: m
        real(CDK), dimension(:,:), intent(inout) :: Mat
        ! Variables
        real(CDK), dimension(:),allocatable :: pivot
        real(CDK) :: a,b
        integer :: i,j,k
        allocate(pivot(1:m))

        do i=1,m
            a=Mat(i,i)
            Mat(i,i)=1.0
            do j=1,m
                pivot(j) = Mat(j,i)/a;
                Mat(j,i) = pivot(j)
            enddo
            do j=1,m
                if(i/=j) then
                    b   = Mat(i,j);
                    Mat(i,j) = 0.0_MK
                    do k=1,m
                        Mat(k,j) = Mat(k,j)-b*pivot(k)
                    enddo
                endif
            enddo
        enddo
    end subroutine

    subroutine AxV(m,Mat,RHS,Res)
        !
        integer(CIK), intent(in) :: m
        real(CDK), dimension(:,:), intent(in) :: Mat
        real(CDK), dimension(:), intent(in) :: RHS
        real(CDK), dimension(:), intent(out) :: Res
        !
        integer :: i,j

        do i=1,m
            Res(i)=0._CDK
            do j=1,m
                Res(i)=Res(i)+Mat(i,j)*RHS(j)
            enddo
        enddo

    end subroutine

    subroutine sv2d_data_preparation(m,XP,YP,XCP,YCP,ds,sine,cosine,slope,ndiv,delx,dely)&
            bind(C,name='sv2d_data_preparation')
        ! Arguments
        integer(CIK), intent(in) :: m
        real(CDK), dimension(*), intent(in) :: XP,YP     !< Panel coordinates
        real(CDK), dimension(*), intent(out)  :: XCP,YCP !< Control Point coordinatse
        real(CDK), dimension(*), intent(out) :: ds       !< Panel Length
        real(CDK), dimension(*), intent(out) :: sine     !
        real(CDK), dimension(*), intent(out) :: cosine
        real(CDK), dimension(*), intent(out) :: slope
        integer(CIK), intent(in) :: ndiv !< number of subdivisions
        real(CDK), dimension(*), intent(out) :: delx !< delta in the x direction for one subdivision
        real(CDK), dimension(*), intent(out) :: dely !< 
        ! Variables
        real(CDK) :: x1,x2,y1,y2,abscos,t
        integer(CIK) :: i
        real(CDK), parameter :: pi=acos(-1._CDK)
        real(CDK), parameter :: ex=0.000001 ! TODO, precision


        do i=1,m
            x1 = XP(i)
            y1 = YP(i)
            x2 = XP(i+1)
            y2 = YP(i+1)
            ! Control Points Between panels
            XCP(i) = (x1+x2)/2
            YCP(i) = (y1+y2)/2
            ! Subdivisions?
            if(ndiv<=0) then
                delx(i)=(x2-x1)/1
                dely(i)=(y2-y1)/1
            else
                delx(i)=(x2-x1)/ndiv
                dely(i)=(y2-y1)/ndiv
            endif
            !
            ds(i)     = sqrt((x2-x1)**2+(y2-y1)**2)
            sine(i)   = (y2-y1)/ds(i);
            cosine(i) = (x2-x1)/ds(i);
            abscos    = abs(cosine(i))
            ! Slope
            if(abscos>ex) then
                t=atan(sine(i)/cosine(i))
                if(cosine(i)>ex) then
                    slope(i) =t
                else ! cosine<-ex
                    slope(i) =t-pi
                endif
            else
                slope(i)=sine(i) / abs(sine(i))*pi/2
            endif
        enddo

    end subroutine


    !> Coupling coefficients computing in a fashion similar to Lewis: 
    ! - A panel of constant surface vorticity is represented by one or more 2D vortex points
    !   If one point is used the vorticity is concentrated at the middle of the panel, where the control point actually is.
    !   That way, the influence of i on j is strongly linked to the influence of j on i
    ! - The curvature correction of Lewis is also used for the self induced velocity (using the slope)
    subroutine sv2d_coupling_coefficients(m,lda,XCP,YCP,ds,sine,cosine,slope,coup,ndiv,delx,dely)&
            bind(C,name='sv2d_coupling_coefficients')
        ! Arguments
        integer(CIK), intent(in) :: m
        integer(CIK), intent(in) :: lda
        real(CDK), dimension(*), intent(in)    :: XCP,YCP    !< Control Point coordinatse
        real(CDK), dimension(*), intent(in)    :: ds         !< Panel Length
        real(CDK), dimension(*), intent(in)    :: sine     !
        real(CDK), dimension(*), intent(in)    :: cosine
        real(CDK), dimension(*), intent(in)    :: slope
        real(CDK), dimension(lda,*), intent(out) :: coup
        integer(CIK), intent(in) :: ndiv
        real(CDK), dimension(*), intent(in) :: delx
        real(CDK), dimension(*), intent(in) :: dely
        ! Variables
        real(CDK) :: r2,u,v
        real(CDK), parameter :: pi      = acos(-1._CDK)
        real(CDK), parameter :: twopi   = 2*pi
        real(CDK), parameter :: eightpi = 8*pi
        integer(CIK) :: i,j,k
        real(CDK) :: xv,yv ! Vortex point location (when using subdivision)

        ! --- Self induced coefficients (diagonal) using curvatue correction from Lewis

        print*,'curv', (slope(2)-slope(m  )-twopi)
        do i=2,m-1
            print*,'curv', (slope(i+1)-slope(i-1))
        enddo
        print*,'curv', (slope(1)-slope(m-1)-twopi)

        coup(1,1) = -0.5_CDK-(slope(2)-slope(m  )-twopi)/(eightpi)
        coup(m,m) = -0.5_CDK-(slope(1)-slope(m-1)-twopi)/(eightpi)
        do i=2,m-1
            coup(i,i) = -0.5_CDK-(slope(i+1)-slope(i-1))/(eightpi)
        enddo
        ! --- Coupling coefficients, off diagonal
        if(ndiv<=1) then ! We use the standard procedure

            ! --- Coupling coefficients off diagonal
            do i=1,m
                xv=XCP(i)
                yv=YCP(i)
                !do j=i,m
                do j=1,m
                    if(j/=i) then
                        ! Distance between two CONTROL POINTS (ie assumes the point vortex is 
                        r2         = (XCP(j)-XCP(i))**2 +(YCP(j)-YCP(i))**2 
                        ! Velocity (unit point vortex)
                        u         =  (YCP(j)-YCP(i))/(twopi*r2)
                        v         = -(XCP(j)-XCP(i))/(twopi*r2)
!                         print*,'i,j,u,v',i,j,u,v, xv,yv
                        ! NOTE: if j=i,j:
                        ! Coupling, making use the fact that control points and vortex points are at the same location
                        !coup(j,i) = (u*cosine(j) + v*sine(j)) * ds(i)
                        !coup(i,j) =-(u*cosine(i) + v*sine(i)) * ds(j)
                        ! From ndiv formula:
                        coup(j,i) = (u*cosine(j) + v*sine(j)) * ds(i)
                    endif
                enddo
            enddo
        else ! We divide each panel into ndiv points

            do i=1,m
                do j=1,m ! can't make use of symmetry here 
                    if(j/=i) then
                        u=0
                        v=0
                        do k=1,ndiv
                            ! 
                            xv=XCP(i) + (k-0.5_CDK*(1+ndiv))*delx(i)
                            yv=YCP(i) + (k-0.5_CDK*(1+ndiv))*dely(i)
                            ! Distance between two CONTROL POINTS (ie assumes the point vortex is 
                            r2         = (XCP(j)-xv)**2 +(YCP(j)-yv)**2 
                            ! Velocity (unit point vortex)
                            u         = u +(YCP(j)-yv)/(r2) ! two pi removed
                            v         = v -(XCP(j)-xv)/(r2)
                           !print*,'i,j,u,v',i,j,u,v, xv,yv
                        enddo
                        u=u/(twopi*ndiv)
                        v=v/(twopi*ndiv)
                        coup(j,i) = (u*cosine(j) + v*sine(j)) * ds(i)
                        !coup(i,j) =-(u*cosine(i) + v*sine(i)) * ds(j)
                    endif
                enddo
            enddo
        endif

    end subroutine

    subroutine sv2d_back_diagonal_correction(m,lda,ds,coup) bind(C,name='sv2d_back_diagonal_correction')
        ! Arguments
        integer(CIK), intent(in) :: m
        integer(CIK), intent(in) :: lda
        real(CDK), dimension(*), intent(in)      :: ds   !< Panel Length
        real(CDK), dimension(lda,*), intent(inout) :: coup !< Coupling coefficient matrix
        ! Variables
        real(CDK) :: cum_sum
        integer(CIK) :: i,j

        do i=1,m
            cum_sum=0.0_CDK
            do j=1,m
                if(j/=m+1-i) then
                    print*,'>>>', j, coup(j,i)*ds(j)
                    cum_sum=cum_sum+coup(j,i)*ds(j)
                endif
            enddo
            print*,'CumSum i',i-1,m+1-i-1, cum_sum, coup(m+1-i,i), -cum_sum / ds(m+1-i)
            coup(m+1-i,i) = -cum_sum / ds(m+1-i)
        enddo

    end subroutine


    subroutine sv2d_bound_vortex_correction(m,lda,ds,coup) bind(C,name='sv2d_bound_vortex_correction')
        ! Arguments
        integer(CIK), intent(in) :: m
        integer(CIK), intent(in) :: lda
        real(CDK), dimension(*), intent(in)    :: ds         !< Panel Length
        real(CDK), dimension(lda,*), intent(inout) :: coup !< Coupling coefficient matrix
        ! Variables
        integer(CIK) :: i,j

        do j=1,m
            do i=1,m
                coup(i,j) = coup(i,j)+ds(j)
            enddo
        enddo

    end subroutine

    !> Enfore the kutta-condition in the matrix. The matrix size is reduced by 1:
    ! Column te+1 is subtracted from clumn te
    subroutine sv2d_kutta_condition(m,lda,te,coup,rhs,m_new) bind(C,name='sv2d_kutta_condition')
        ! Arguments
        integer(CIK), intent(in) :: m
        integer(CIK), intent(in) :: lda
        integer(CIK), intent(in) :: te
        real(CDK), dimension(lda,*), intent(inout) :: coup
        real(CDK), dimension(*), intent(inout) :: rhs
        integer(CIK), intent(out) :: m_new
        ! Variables
        integer(CIK) :: i,j

        ! Substracting column te+1 from column te
        do j=te,m
            do i=1,m
                if(j>te) then
                    ! moving data from next columns to current column
                    coup(i,j)=coup(i,j+1)
                else ! j==te
                    ! Substracting column te+1 from te
                    coup(i,j)=coup(i,j)-coup(i,j+1)
                endif
            enddo
        enddo

        ! Substracting row te+1 from row te
        do i=te,m
            do j=1,m
                if(i>te) then
                    ! moving data from next row to current row
                    coup(i,j)=coup(i+1,j)
                else ! i==te
                    ! Substracting row te+1 from te
                    coup(i,j)=coup(i,j)-coup(i+1,j)
                endif
            enddo
        enddo

        ! Performing substraction of rows for rhs also

        do i=te+1,m
            if(i>te+1) then
                ! moving data from next row to current row
                rhs(i-1)=rhs(i)
            else ! i==te
                ! Substracting row te+1 from te
                rhs(i-1)=rhs(i-1)-rhs(i)
            endif
        enddo

        ! SAFETY
        coup(m,1:m)=-9999999
        coup(:,m)=-9999999
        rhs (m)  =-9999999

        m_new = m-1

    end subroutine


    subroutine sv2d_revert_kutta(m,te,Sol)
        integer(CIK), intent(in) :: m
        integer(CIK), intent(in) :: te
        real(CDK), dimension(:), intent(inout) :: Sol
        !
        integer :: i
        ! Shift lower surface values along the array by one place
        do i=te,m-2
            Sol(m-i+te)=Sol(m-i+te-1)
        enddo
        ! Replate lower trailing edge value by minus te value
        Sol(te+1)=-Sol(te)
    end subroutine


    ! Velocity due to m panels at XPv and YPv 
    subroutine sv2d_ui_n1(X,Y,m,XPv,YPv,Gammas,U,V) bind(C,name='sv2d_ui_n1')
        ! Arguments
        integer(CIK), intent(in) :: m
        real(CDK), intent(in)    :: X,Y    !< Control Point coordinatse
        real(CDK), dimension(*), intent(in)    :: XPv,YPv !< the XCP and YCP of the other subroutines above
        real(CDK), dimension(*), intent(in)    :: Gammas ! = gammas*ds [m*m/s]
        real(CDK), intent(out)    :: U,V    !< Velocity at control point
        ! Variables
        real(CDK) :: r2,dx,dy, contrib
        real(CDK), parameter :: pi = acos(-1._CDK)
        real(CDK), parameter :: twopi =2*pi
        integer(CIK) :: i
        u=0
        v=0
        ! --- Loop on vortex elements
        do i=1,m
            ! Distance between two CONTROL POINTS (ie assumes the point vortex is 
            dx  = (X-XPv(i))
            dy  = (Y-YPv(i))
            r2  = (dx)**2 +(dy)**2 
            contrib=Gammas(i) /(twopi*r2)
            ! Cumulative Velocity 
            u         = u+dy*contrib
            v         = v-dx*contrib
        enddo
    end subroutine



end module SurfaceVorticity2DTools
