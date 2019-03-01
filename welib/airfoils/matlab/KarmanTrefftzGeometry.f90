program mapsl
    !  Based on: Justin E. Kerwin December, 2000
    !------------------------------------------------------------------------------­

    !	Declare the variables----------------------------------------------------­
    implicit none
    character(len=1) yesno
    character(len=20) :: buffer
    character(len=20) :: params
    real, parameter :: pi=3.1415927e00, rad=1.7453293e-02
    real :: x,y,theta,rsq,u,v,r,xc,yc,alpha_deg,alpha,rc,beta,rcsq,g,gamma, &
        lambda,tau_deg,xi,eta,x_circle,y_circle,fl2,u_foil,v_foil, &
        theta1,theta2,dtheta,rmax,dr,cp,xle,xte,chord,cl
    real, allocatable, dimension(:) :: xi_foil,cp_foil,eta_foil
    complex :: z,zp1,zm1,zeta,w_circle,dzeta_dz,w_foil
    integer :: n,m,ntheta,nradial

    !input the variables defining foil geometry and angle of attack ----------­
    write(*,'(a)') ' enter xc,yc,trailing edge angle ' ! position of center of circle--------­
    read(*,*) xc,yc,tau_deg ! note: xc must be <= 0 ------------­
    lambda=2.0-tau_deg/180.0 ! the exponent in the k-t mapping function -­
    fl2=4.0*lambda**2 ! pre-comute a mapping function constant----­
    rcsq=(1.0-xc)**2+yc**2 ! radius of circle passing through (1,0) ---­

    rc=sqrt(rcsq) 
    beta=atan(yc/(1.0-xc)) 

    alpha_deg=0
    theta1=0
    theta2=360
    dtheta=1
    ntheta=nint((theta2+theta1)/dtheta)+1 ! number of radial grid lines 
    !	alocate the arrays and open output files --------------------------------­
    allocate(xi_foil(ntheta),eta_foil(ntheta),cp_foil(ntheta)) ! arrays for foil cp -----­

    buffer=leading_zero(xc)
    params='_'//trim(buffer)
    buffer=leading_zero(yc)
    params=trim(params)//'_'//trim(buffer)
    write(buffer,'(I0)')int(tau_deg)
    params=trim(params)//'_'//trim(buffer)
!     write(buffer,'(I0)')int(alpha_deg)
!     params=trim(params)//'_'//trim(buffer)

!     open(2,file='cp-karman-trefftz'//trim(params)//'.dat',status='unknown',form='formatted')
    open(3,file='geom-karman-trefftz'//trim(params)//'.dat',status='unknown',form='formatted')

    !-----generate velocity and presure field -------------------------------------­
    m=1
    r=rc+real(m-1)*dr
    rsq=r**2
    ! grid radial coordinate in the z (circle) plane--­
    do n=1,ntheta
        theta=rad*(theta1-real(n-1)*dtheta) ! grid angular coordinate--­
        x=xc+r*cos(theta) ! convert to cartesian-----­
        y=yc+r*sin(theta)
        !-----------compute the velocity field around the circle-----------------------­
        u=cos(alpha)-rcsq*cos(2.0*theta-alpha)/rsq-g*sin(theta)/r
        v=sin(alpha)-rcsq*sin(2.0*theta-alpha)/rsq+g*cos(theta)/r

        !-----------express the field point position and velocity in complex form------­
        z=cmplx(x,y) ! form complex number (x+iy)-------­
        w_circle=cmplx(u,-v)
        ! complex velocity in z plane------­
        !-----------use the karman-trefftz transformation to map points to zeta plane--­
        zp1=(z+1.0)**lambda
        zm1=(z-1.0)**lambda
        zeta=lambda*(zp1+zm1)/(zp1-zm1)
        xi=real(zeta) ! transformed x coordinate---------­
        eta=aimag(zeta) ! transformed y coordinate --------­
        !-----------compute the derivative of the mapping function; transform velocities
        dzeta_dz=fl2*((z-1.0)**(lambda-1.0)*(z+1.0)**(lambda-1.0))/ (zp1-zm1)**2

        if(abs(dzeta_dz)>0.0001) then
            w_foil=w_circle/dzeta_dz ! complex velocity in foil plane------­

        else
            w_foil=0.0
        end if
        u_foil=real(w_foil)
        v_foil=-aimag(w_foil)

        ! avoids divide by zero at trailing edge

        !-----------compute the presure coefficient, cp, and output to plotting file---­

        cp=abs(w_foil)**2-1.0

        if(m==1) then
            ! save presure on foil surface for later use--­
            xi_foil(n)=xi
            eta_foil(n)=eta
            cp_foil(n)=cp
        end if
    end do
    !-----scale the chordwise coordinate to (0,1) and output foil pressure dist.---­
    xle=10.0
    xte=-10.0
    do n=1,ntheta
        xle=min(xle,xi_foil(n))
        xte=max(xte,xi_foil(n))
    end do
    ! Coordinate scalling
    chord=xte-xle
    do n=1,ntheta
        xi_foil(n)=(xi_foil(n)-xle)/chord
        eta_foil(n)=eta_foil(n)/chord
    end do
!     write(2,'(2f10.5)') (xi_foil(n),cp_foil(n),n=1,ntheta) 
!     close(2) 

    write(3,'(2f10.5)') (xi_foil(n),eta_foil(n),n=1,ntheta) 
    close(3) 

contains

    function leading_zero(x) result(buffer)
        real :: x
        character(len=20) :: buffer
        integer :: ix
        if (x.lt.0) then
            ix = ceiling(x)
        else
            ix = floor(x)
        endif
        if (x.lt.0) then
            write(buffer,1) abs(ix), abs(x)-floor(abs(x))
            1   format('-',i0,f0.1)    
        else
            write(buffer,2) floor(x), abs(x)-floor(abs(x))
            2   format(i0,f0.1)
        endif

    end function




end program mapsl
