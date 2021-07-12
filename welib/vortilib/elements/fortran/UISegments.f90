!> Segments Induced velocity functions.
module UISegments
    use PrecisionMod, only: MK,PRECISION_UI
    implicit none

    ! Norm for the simple test exactly on singulartiy
    real(MK),parameter :: MINNORMSIMP=1e-6_MK
    ! Min norm for escaping if denominator or cross prod too small ! Has quite some influence!!!
    ! real(MK),parameter :: MINNORM=1e-08_MK
    ! real(MK),parameter :: MINNORM2=1e-16_MK
    ! Added to denominator for floating point
    real(MK),parameter :: MINDENOM=1e-15_MK
    real(MK),parameter :: MIN_EXP_VALUE=-10.0_MK
    !integer:: iunitSeg
    !logical :: bOutputSeg=.false.
contains
    !>
    !> Takes n segments and n control points
    subroutine fUi_SegmentCst()
        implicit none
    end subroutine

    !> Takes n segments one ofter the other and n control points
    subroutine fUi_ContinousLine(CPs,SgmtP,Intensities, SmoothModel,SmoothParam,bComputeGrad,Uind, Grad,nCPs,nSgmt)
        implicit none
        ! Arguments declarations 
        integer, intent(in) :: nCPs
        integer, intent(in) :: nSgmt
        logical, intent(in) :: bComputeGrad !< 
        real(MK), dimension(3,nCPs), intent(in)    :: CPs          !< side effects
        real(MK), dimension(3,nCPs), intent(inout) :: Uind         !< side effects
        real(MK), dimension(9,nCPs), intent(inout) :: Grad         !< side effects
        real(MK), dimension(3,nSgmt+1), intent(in) :: SgmtP        !< 
        real(MK), dimension(nSgmt), intent(in)   :: Intensities  !< 
        integer, intent(in)                        :: SmoothModel  !< 
        real(MK), dimension(nSgmt), intent(in)     :: SmoothParam  !< 
        ! Variable declarations 
        real(MK), dimension(3) :: Uindtmp      !<  
        real(MK), dimension(9) :: Gradtmp      !<  
        real(MK), dimension(3) :: DP1, DP2
        integer :: icp,is
        ! loop on control points 
        do icp=1,nCPs
            ! loop on segments
            do is=1,nSgmt
                Uindtmp=0.0_MK;
                Gradtmp=0.0_MK;
                DP1=CPs(:,icp)-SgmtP(:,is)
                DP2=CPs(:,icp)-SgmtP(:,is+1)
                call fUi_SegmentCst_11(DP1, DP2, Intensities(is), SmoothModel,&
                    SmoothParam(is), bComputeGrad, Uindtmp, Gradtmp)
                Uind(1:3,icp) = Uind(1:3,icp)+Uindtmp(1:3)
                Grad(1:9,icp) = Grad(1:9,icp)+Gradtmp(1:9)
            enddo 
        enddo
    end subroutine


    !> Used to be fUi_VortexLine, takes one segment and one control points
    subroutine fUi_SegmentCst_11_b(DeltaPa, DeltaPb, norm2_r0, Intensity, visc_model , visc_param, bComputeGrad, Ui, Grad)
        use MathConstants 
        implicit none 
        ! Input/output arguments 
        real(MK), dimension(3), intent(in) :: DeltaPa      !< 3 x 1   Pcp-P1 !!!!!!!!!!
        real(MK), dimension(3), intent(in) :: DeltaPb      !< 3 x 1   Pcp-P2
        real(MK), intent(in)               :: Intensity    !< 
        real(MK), intent(in)               :: norm2_r0      !<  norm squared of segment
        integer, intent(in)                :: visc_model   !< 
        real(MK), intent(in)               :: visc_param   !< 
        logical, intent(in)                :: bComputeGrad !< 
        real(MK), dimension(3),intent(out) :: Ui           !<  No Side effects
        real(MK), dimension(9),intent(out) :: Grad         !<  No Side effects
        ! Variables declaration 
        real(MK),dimension(3) :: crossprod      !< 
        !real(MK),dimension(3) :: D              !< 
        real(MK)              :: denom          !< 
        real(MK)              :: denominator    !< 
        real(MK)              :: h2             !< Square of h
!         real(MK)              :: h2_inv         !< Square of h
        real(MK)              :: h              !< Only used by one model
        real(MK)              :: Kv             !< 
        !real(MK)              :: Kv2            !< Gradient kernel only for Lamb-Oseen
!         real(MK)              :: Kv_tmp         !< 
        real(MK)              :: norm_a         !< 
        real(MK)              :: norm_b         !< 
        real(MK)              :: norm2_crossprod !< 
        !real(MK)              :: norm2_r0tmp !< 
        real(MK)              :: xa             !< 
        real(MK)              :: ya             !< 
        real(MK)              :: za             !< 
        real(MK)              :: xb             !< 
        real(MK)              :: yb             !< 
        real(MK)              :: zb             !< 
        real(MK)              :: visc_param2    !< Square of viscous param
        real(MK)              :: exp_value      !< 
        ! 
        xa=DeltaPa(1)
        ya=DeltaPa(2)
        za=DeltaPa(3)
        xb=DeltaPb(1)
        yb=DeltaPb(2)
        zb=DeltaPb(3)
        !--------------------------------------------------------------------------------
        !--- Simple test if on an extremity point
        !--------------------------------------------------------------------------------
        if(abs(xa)+abs(ya)+abs(za)<MINNORMSIMP) then
            Ui(1:3)=0.0_MK
            Grad(1:9)=0.0_MK
        elseif(abs(xb)+abs(yb)+abs(zb)<MINNORMSIMP) then
            Ui(1:3)=0.0_MK
            Grad(1:9)=0.0_MK
        else
            norm_b = sqrt(xb*xb + yb*yb + zb*zb)
            norm_a = sqrt(xa**2 + ya**2 + za**2)
            !denominator = norm_a*norm_b*(norm_a*norm_b + xa*xb+ya*yb+za*zb)
            denom       =               (norm_a*norm_b + xa*xb+ya*yb+za*zb)
            denominator = norm_a*norm_b*denom
            crossprod(1) = ya*zb-za*yb
            crossprod(2) = za*xb-xa*zb
            crossprod(3) = xa*yb-ya*xb
            norm2_crossprod=crossprod(1)**2+ crossprod(2)**2+ crossprod(3)**2
            !print "(A,6F21.16)","xa xb:",xa,ya,za,xb,yb,zb
            ! check for singularity  
            ! TODO this check is viscous model dependent(e.g. crossprod)... have to go towards different funcitons soon..
            ! to avoid singularity, save a bit of computation and use limit values
            if (denominator<PRECISION_UI .or. norm2_crossprod<PRECISION_UI) then
                !--------------------------------------------------------------------------------
                !--- Exactly on the Singularity, velocity is zero, but gradient may be constant (is is wanted?)
                !--------------------------------------------------------------------------------
                !                 if(denominator<MINNORM) then
                !                     print*,'because of denom'
                !                 endif
                !                 if(norm2_crossprod<MINNORM2) then
                !                     print*,'because of crossprod'
                !                 endif
                Ui(1:3)=0.0_MK
                Grad(1:9)=0.0_MK
                if( visc_model==0) then 
                    return
                end if
                if(bComputeGrad) then
                    ! denom=(norm_a*norm_b + xa*xb+ya*yb+za*zb)+MINDENOM
                    ! if(norm_a <MINNORM) norm_a=MINNORM
                    ! if(norm_b <MINNORM) norm_b=MINNORM
                    ! if(denom <MINNORM) denom=MINNORM
                    ! ! For now only one limit for all models...
                    ! h2 = MINNORM2
                    ! Kv = 1.0_MK-exp(-1.25643_MK*h2/visc_param2)
                    ! Grad(2) = Kv**2*Intensity*fourpi_inv*1/denom*(1/norm_a+1/norm_b)*(zb-za) ! 12
                    ! Grad(3) = Kv**2*Intensity*fourpi_inv*1/denom*(1/norm_a+1/norm_b)*(ya-yb) ! 13
                    ! Grad(4) = Kv**2*Intensity*fourpi_inv*1/denom*(1/norm_a+1/norm_b)*(za-zb) ! 21
                    ! Grad(6) = Kv**2*Intensity*fourpi_inv*1/denom*(1/norm_a+1/norm_b)*(xb-xa) ! 23
                    ! Grad(7) = Kv**2*Intensity*fourpi_inv*1/denom*(1/norm_a+1/norm_b)*(yb-ya) ! 31
                    ! Grad(8) = Kv**2*Intensity*fourpi_inv*1/denom*(1/norm_a+1/norm_b)*(xa-xb) ! 32
                    ! print*,Grad(4)
                endif
            else
                !norm2_r0tmp = (xa-xb)*(xa-xb) + (ya-yb)*(ya-yb) +(za-zb)*(za-zb)
                ! if(abs(norm2_r0-norm2_r0tmp)>1e-3) then
                !     print*,'Norm mismatch',norm2_r0tmp,norm2_r0
                !     STOP
                ! endif
                denominator=denominator+MINDENOM
                ! TODO far distance
                !--------------------------------------------------------------------------------
                !--- Normal Procedure 
                !--------------------------------------------------------------------------------
                ! viscous model  
                select case ((visc_model)) !
                case ( 0 ) !! No vortex core model
                    Kv=1.0_MK
                case ( 1 ) !!Rankine - t<=>rc
                    ! orthogonal distance r1xr2/r0 
                    h2 = norm2_crossprod/ norm2_r0 
                    visc_param2=visc_param**2
                    if (h2<visc_param2) then 
                        Kv=h2/visc_param2
                    else
                        Kv=1.0_MK 
                    end if 
                case ( 2 ) !!Lamb-Oseen - vsic_param<=>rc 
                    ! BUG WITH OPENMP ??
                    ! orthogonal distance r1xr2/r0 
                    h2 = norm2_crossprod/ norm2_r0 
                    visc_param2=visc_param**2
                    exp_value = -1.25643_MK*(h2)/(visc_param2)
                    if(exp_value<MIN_EXP_VALUE) then
                        Kv = 1.0_MK
                    else
                        Kv = 1.0_MK-exp(exp_value)
                    endif
                case ( 3 ) !!Vatistas n=2 - visc_param<=>rc 
                    ! orthogonal distance r1xr2/r0 
                    h2 = norm2_crossprod/ norm2_r0 
                    visc_param2=visc_param**2
                    ! h = (norm_a+norm_b)/2; 
                    Kv = h2/sqrt(visc_param2**2+h2**2)
                case ( 4 ) !!Cut-off radius as function of length: delta^2*norm(r0)^2  
                    Kv=1.0_MK
                    denominator=denominator+visc_param*visc_param*norm2_r0
                case ( 5 ) !!Cut-off radius dimension fixed: (delta l_0)^2  
                    Kv=1.0_MK
                    denominator=denominator+visc_param*visc_param
                case ( 33 ) !!Vatistas n=2 - visc_paramt<=>rc  
                    ! See Hoydonck 2012. and Script MainRingSensitivityN
                    ! Not matture enough
                    h = (norm_a+norm_b)/2.0_MK
                    visc_param2=visc_param**2
                    if(h<2*sqrt(norm2_r0)) then
                        ! use Vatistas n=2, normal one
                        h2 = norm2_crossprod/norm2_r0
                    end if
                    Kv = (h2)/sqrt(visc_param2**2+h2**2)
                case DEFAULT
                    Kv=1.0_MK
                end select 

                Kv=Intensity*fourpi_inv*Kv*(norm_a+norm_b)/denominator
                Ui(1:3) = Kv*crossprod(1:3)
                !print "(A,3F21.16)","ui   :",Ui(1:3)

                if (bComputeGrad) then 
                    Grad(1:9)=0.0_MK
                    !denom=(norm_a*norm_b + xa*xb+ya*yb+za*zb)+MINDENOM
!                     denom=denom+MINDENOM
!                     D=-(DeltaPa/norm_a**3+DeltaPb/norm_b**3)&
!                         -(1./norm_a+1./norm_b)*1/denom*(DeltaPa/norm_a*norm_b+ DeltaPb/norm_b*norm_a+ DeltaPa+DeltaPb)
! 
!                     Grad(1) = crossprod(1)*D(1)                              ! 11
!                     Grad(2) = crossprod(1)*D(2) +(1./norm_a+1./norm_b)*(zb-za) ! 12
!                     Grad(3) = crossprod(1)*D(3) +(1./norm_a+1./norm_b)*(ya-yb) ! 13
! 
!                     Grad(4) = crossprod(2)*D(1) +(1./norm_a+1./norm_b)*(za-zb) ! 21
!                     Grad(5) = crossprod(2)*D(2)                              ! 22
!                     Grad(6) = crossprod(2)*D(3) +(1./norm_a+1./norm_b)*(xb-xa) ! 23
! 
!                     Grad(7) = crossprod(3)*D(1) +(1./norm_a+1./norm_b)*(yb-ya) ! 31
!                     Grad(8) = crossprod(3)*D(2) +(1./norm_a+1./norm_b)*(xa-xb) ! 32
!                     Grad(9) = crossprod(3)*D(3)                              ! 33
! 
!                     Grad=Intensity*fourpi_inv*1./denom*Grad;
!                     if(visc_model>0) then
!                         ! h2 and visc_param2 are the same for everybody
!                         ! We use Lamb-Oseen model (it's the easiest)
!                         exp_value = -1.25643_MK*(h2)/(visc_param2)
!                         if(exp_value<MIN_EXP_VALUE) then
!                             ! Far away, we then use inviscid Gradient from above.
!                         else
!                             ! The Lamb-Oseen kernel (may be different than above Kv)
!                             Kv2 = 1.0_MK-exp(exp_value)
!                             ! Computing the Grad of Kv2
!                             D=0._MK
!                             ! Grad(raxrb) ra x rb 
!                             ! !!!!! Watcth out, need the segment extremity coordinates, but we don't have it
!                             ! It's the opposite formula hence, since xa-xb=x1-x2
!                             D(1)=(zb-za)*crossprod(2)+(ya-yb)*crossprod(3);
!                             D(2)=(za-zb)*crossprod(1)+(xb-xa)*crossprod(3);
!                             D(3)=(yb-ya)*crossprod(1)+(xa-xb)*crossprod(2);
!                             D=-2._MK*1.25643_MK/(norm2_r0 *visc_param2)*exp(exp_value)*D;
!                             !                         print*,'Kv',Kv
!                             !                         print*,'Grad',Grad
! 
!                             ! Grad_viscous=Kv2*Grad + u_invisc * D
!                             Grad=Kv2*Grad ! Intensity already in Grad and Ui
                            ! Since Ui contains already Kv, we have to divide by it. Alternative: store Ui inviscid
!                            print*,'Kv',Kv,visc_param,exp_value,Kv2! this print is here to remind that c code has no grad
                           ! also a singulartiy can occur with the grid.
!                             Grad(1) = Grad(1) + Ui(1)*D(1)/Kv
!                             Grad(2) = Grad(2) + Ui(1)*D(2)/Kv
!                             Grad(3) = Grad(3) + Ui(1)*D(3)/Kv
! 
!                             Grad(4) = Grad(4) + Ui(2)*D(1)/Kv
!                             !                         print*, 'Term3', Ui(2)*D(1)/Kv
!                             Grad(5) = Grad(5) + Ui(2)*D(2)/Kv
!                             Grad(6) = Grad(6) + Ui(2)*D(3)/Kv
! 
!                             Grad(7) = Grad(7) + Ui(3)*D(1)/Kv
!                             Grad(8) = Grad(8) + Ui(3)*D(2)/Kv
!                             Grad(9) = Grad(9) + Ui(3)*D(3)/Kv
!                         endif
!                     end if  ! Visc model >0
                else
                    Grad(1:9)=0.0_MK
                end if  ! bComputeGrad
            end if ! denominator size or distances too small
        end if ! 
    end subroutine fUi_SegmentCst_11_b

    !> Used to be fUi_VortexLine, takes one segment and one control points
    subroutine fUi_SegmentCst_11(DeltaPa, DeltaPb, Intensity, visc_model , visc_param, bComputeGrad, Ui, Grad)
        use MathConstants 
        implicit none 
        ! Input/output arguments 
        real(MK), dimension(3), intent(in) :: DeltaPa      !< 3 x 1   Pcp-P1 !!!!!!!!!!
        real(MK), dimension(3), intent(in) :: DeltaPb      !< 3 x 1   Pcp-P2
        real(MK), intent(in)               :: Intensity    !< 
        integer, intent(in)                :: visc_model   !< 
        real(MK), intent(in)               :: visc_param   !< 
        logical, intent(in)                :: bComputeGrad !< 
        real(MK), dimension(3),intent(out) :: Ui           !<  No Side effects
        real(MK), dimension(9),intent(out) :: Grad         !<  No Side effects
        ! Variables declaration 
        real(MK),dimension(3) :: crossprod      !< 
        real(MK),dimension(3) :: D              !< 
        real(MK)              :: denom          !< 
        real(MK)              :: denominator    !< 
        real(MK)              :: h2             !< Square of h
        real(MK)              :: h              !< Only used by one model
        real(MK)              :: Kv             !< 
        real(MK)              :: norm_a         !< 
        real(MK)              :: norm_b         !< 
        real(MK)              :: norm2_r0        !< 
        real(MK)              :: norm2_crossprod !< 
        real(MK)              :: xa             !< 
        real(MK)              :: ya             !< 
        real(MK)              :: za             !< 
        real(MK)              :: xb             !< 
        real(MK)              :: yb             !< 
        real(MK)              :: zb             !< 
        real(MK)              :: visc_param2    !< Square of viscous param
        real(MK)              :: exp_value      !< 
        ! 
        xa=DeltaPa(1)
        ya=DeltaPa(2)
        za=DeltaPa(3)
        xb=DeltaPb(1)
        yb=DeltaPb(2)
        zb=DeltaPb(3)
        !--------------------------------------------------------------------------------
        !--- Simple test if on an extremity point
        !--------------------------------------------------------------------------------
        if(abs(xa)+abs(ya)+abs(za)<MINNORMSIMP) then
            Ui(1:3)=0.0_MK
            Grad(1:9)=0.0_MK
            !print*,'Exit1'
        elseif(abs(xb)+abs(yb)+abs(zb)<MINNORMSIMP) then
            Ui(1:3)=0.0_MK
            Grad(1:9)=0.0_MK
            !print*,'Exit2'
        else
            norm_b = sqrt(xb*xb + yb*yb + zb*zb)
            norm_a = sqrt(xa**2 + ya**2 + za**2)
            norm2_r0 = (xa-xb)*(xa-xb) + (ya-yb)*(ya-yb) +(za-zb)*(za-zb)
            !denominator = norm_a*norm_b*(norm_a*norm_b + xa*xb+ya*yb+za*zb)
            denom       =               (norm_a*norm_b + xa*xb+ya*yb+za*zb)
            denominator = norm_a*norm_b*denom
            crossprod(1) = ya*zb-za*yb
            crossprod(2) = za*xb-xa*zb
            crossprod(3) = xa*yb-ya*xb
            norm2_crossprod=crossprod(1)**2+ crossprod(2)**2+ crossprod(3)**2
            !if(abs(norm2_r0-norm2_seg)>1e-20_MK) then
            !    print*,'Norm mismatch'
            !endif
            ! check for singularity  
            !if(bOutputSeg) then
            !    write(iunitSeg,'(a,3EN14.2E2)') 'a', xa,ya,za
            !    write(iunitSeg,'(a,3EN14.2E2)') 'b', xb,yb,zb
            !    write(iunitSeg,'(a,3EN14.2E2)') 'sum', norm_a+norm_b, sqrt(norm2_r0),norm_a+norm_b-sqrt(norm2_r0)
            !    write(iunitSeg,'(3EN14.2E2)') norm2_r0 , norm_a, norm_b
            !    write(iunitSeg,'(a, 2EN14.2E2)') 'den', denominator, MINNORM*norm2_r0 
            !    write(iunitSeg,'(a, 3EN14.2E2)') 'cor', norm2_crossprod, MINNORM2*norm2_r0, PRECISION_UI
            !endif
            ! TODO this check is viscous model dependent(e.g. crossprod)... have to go towards different funcitons soon..
            ! to avoid singularity, save a bit of computation and use limit values
            if (norm2_r0<PRECISION_UI .or. denominator<PRECISION_UI .or. norm2_crossprod<PRECISION_UI) then
                !--------------------------------------------------------------------------------
                !--- Exactly on the Singularity, velocity is zero, but gradient may be constant (is is wanted?)
                !--------------------------------------------------------------------------------
                if(norm2_r0<PRECISION_UI) then
                    !print*,'Exit3, because of segment too small'
                else
                    if(denominator<PRECISION_UI) then
                        !print*,'Exit3, because of denom'
                    endif
                    if(norm2_crossprod<PRECISION_UI) then
                        !print*,'Exit3, because of crossprod'
                    endif
                endif
                !if(bOutputSeg) then
                !    write(iunitSeg,'(a)') 'exitc'
                !endif
                Ui(1:3)=0.0_MK
                Grad(1:9)=0.0_MK
                if( visc_model==0) then 
                    return
                end if
                if(bComputeGrad) then
                    ! denom=(norm_a*norm_b + xa*xb+ya*yb+za*zb)+MINDENOM
                    ! if(norm_a <MINNORM) norm_a=MINNORM
                    ! if(norm_b <MINNORM) norm_b=MINNORM
                    ! if(denom <MINNORM) denom=MINNORM
                    ! ! For now only one limit for all models...
                    ! h2 = MINNORM2
                    ! Kv = 1.0_MK-exp(-1.25643_MK*h2/visc_param2)
                    ! Grad(2) = Kv**2*Intensity*fourpi_inv*1/denom*(1/norm_a+1/norm_b)*(zb-za) ! 12
                    ! Grad(3) = Kv**2*Intensity*fourpi_inv*1/denom*(1/norm_a+1/norm_b)*(ya-yb) ! 13
                    ! Grad(4) = Kv**2*Intensity*fourpi_inv*1/denom*(1/norm_a+1/norm_b)*(za-zb) ! 21
                    ! Grad(6) = Kv**2*Intensity*fourpi_inv*1/denom*(1/norm_a+1/norm_b)*(xb-xa) ! 23
                    ! Grad(7) = Kv**2*Intensity*fourpi_inv*1/denom*(1/norm_a+1/norm_b)*(yb-ya) ! 31
                    ! Grad(8) = Kv**2*Intensity*fourpi_inv*1/denom*(1/norm_a+1/norm_b)*(xa-xb) ! 32
                    ! print*,Grad(4)
                endif
            else
                !denominator=denominator+MINDENOM
                denominator=denominator
                ! TODO far distance
                !--------------------------------------------------------------------------------
                !--- Normal Procedure 
                !--------------------------------------------------------------------------------
                ! viscous model  
                select case ((visc_model)) !
                case ( 0 ) !! No vortex core model
                    Kv=1.0_MK
                case ( 1 ) !!Rankine - t<=>rc
                    ! orthogonal distance r1xr2/r0 
                    h2 = norm2_crossprod/ norm2_r0 
                    visc_param2=visc_param**2
                    if (h2<visc_param2) then 
                        Kv=h2/visc_param2
                    else
                        Kv=1.0_MK 
                    end if 
                case ( 2 ) !!Lamb-Oseen - vsic_param<=>rc 
                    ! BUG WITH OPENMP ??
                    ! orthogonal distance r1xr2/r0 
                    h2 = norm2_crossprod/ norm2_r0 
                    visc_param2=visc_param**2
                    exp_value = -1.25643_MK*(h2)/(visc_param2)
                    if(exp_value<MIN_EXP_VALUE) then
                        Kv = 1.0_MK
                    else
                        Kv = 1.0_MK-exp(exp_value)
                    endif
                case ( 3 ) !!Vatistas n=2 - visc_param<=>rc 
                    ! orthogonal distance r1xr2/r0 
                    h2 = norm2_crossprod/ norm2_r0 
                    visc_param2=visc_param**2
                    ! h = (norm_a+norm_b)/2; 
                    Kv = h2/sqrt(visc_param2**2+h2**2)
                case ( 4 ) !!Cut-off radius no dimensions - visc_param<=>delta**2 
                    h2 = norm2_crossprod/ norm2_r0 
                    Kv=1.0_MK
                    ! delta*norm(r0)^2  
                    denominator=denominator+visc_param*norm2_r0
                    visc_param2=visc_param**2 ! TODO
                case ( 5 ) !!Cut-off radius dimension - visc_param<=>(delta l)**2 
                    h2 = norm2_crossprod/ norm2_r0 
                    Kv=1.0_MK
                    ! (delta l_0)^2  
                    denominator=denominator+visc_param
                    visc_param2=visc_param**2 ! TODO
                case ( 33 ) !!Vatistas n=2 - visc_paramt<=>rc  
                    ! See Hoydonck 2012. and Script MainRingSensitivityN
                    ! Not matture enough
                    h = (norm_a+norm_b)/2.0_MK
                    visc_param2=visc_param**2
                    if(h<2*sqrt(norm2_r0)) then
                        ! use Vatistas n=2, normal one
                        h2 = norm2_crossprod/norm2_r0
                    else
                        h2 = 1._MK ! TODO
                    endif
                    Kv = (h2)/sqrt(visc_param2**2+h2**2)
                case DEFAULT
                    Kv=1.0_MK
                end select 

                Kv=Intensity*fourpi_inv*Kv*(norm_a+norm_b)/denominator
                Ui(1:3) = Kv*crossprod(1:3)
                !if(bOutputSeg) then
                !    write(iunitSeg,'(A,E14.2E2)') 'Kv', Kv
                !    write(iunitSeg,'(A,3E14.2E2)') 'crossprod', crossprod
                !    write(iunitSeg,'(A,E14.2E2)') 'Intensity', Intensity
                !    write(iunitSeg,'(A,3E14.2E2)') 'Ui', Ui
                !endif

                if (bComputeGrad) then 
                    denom=(norm_a*norm_b + xa*xb+ya*yb+za*zb)+MINDENOM
                    D=-(DeltaPa/norm_a**3+DeltaPb/norm_b**3)&
                        -(1/norm_a+1/norm_b)*1/denom*(DeltaPa/norm_a*norm_b+ DeltaPb/norm_b*norm_a+ DeltaPa+DeltaPb)

                    Grad(1) = crossprod(1)*D(1)                              ! 11
                    Grad(2) = crossprod(1)*D(2) +(1/norm_a+1/norm_b)*(zb-za) ! 12
                    Grad(3) = crossprod(1)*D(3) +(1/norm_a+1/norm_b)*(ya-yb) ! 13

                    Grad(4) = crossprod(2)*D(1) +(1/norm_a+1/norm_b)*(za-zb) ! 21
                    Grad(5) = crossprod(2)*D(2)                              ! 22
                    Grad(6) = crossprod(2)*D(3) +(1/norm_a+1/norm_b)*(xb-xa) ! 23

                    Grad(7) = crossprod(3)*D(1) +(1/norm_a+1/norm_b)*(yb-ya) ! 31
                    Grad(8) = crossprod(3)*D(2) +(1/norm_a+1/norm_b)*(xa-xb) ! 32
                    Grad(9) = crossprod(3)*D(3)                              ! 33

                    Grad=Intensity*fourpi_inv*1/denom*Grad;
                    if(visc_model>0) then
                        ! h2 and visc_param2 are the same for everybody
                        ! We use Lamb Oseen model !!! TODO, expensive if far away!!!
                        exp_value = -1.25643_MK*(h2)/(visc_param2)
                        if(exp_value<MIN_EXP_VALUE) then
                            ! Far away, not need to compute smooth gradient
                            Kv = 1.0_MK
                        else
                            Kv = 1.0_MK-exp(exp_value)
                        !                     end select
                            ! Computing the Grad of Kv
                            D=0._MK
                            ! Grad(raxrb) ra x rb 
                            ! !!!!! Watcth out, need the segment extremity coordinates, but we don't have it
                            ! It's the opposite formula hence, since xa-xb=x1-x2
                            D(1)=(zb-za)*crossprod(2)+(ya-yb)*crossprod(3);
                            D(2)=(za-zb)*crossprod(1)+(xb-xa)*crossprod(3);
                            D(3)=(yb-ya)*crossprod(1)+(xa-xb)*crossprod(2);
                            D=-2*1.25643/(norm2_r0 *visc_param2)*exp(exp_value)*D;
                            !                         print*,'Kv',Kv
                            !                         print*,'Grad',Grad

                            Grad=Kv*Grad ! Intensity already in Grad and Ui
                            ! Since Ui contains already Kv, we have to divide by it..
                            Grad(1) = Grad(1) + Ui(1)*D(1)/Kv
                            Grad(2) = Grad(2) + Ui(1)*D(2)/Kv
                            Grad(3) = Grad(3) + Ui(1)*D(3)/Kv

                            Grad(4) = Grad(4) + Ui(2)*D(1)/Kv
                            !                         print*, 'Term3', Ui(2)*D(1)/Kv
                            Grad(5) = Grad(5) + Ui(2)*D(2)/Kv
                            Grad(6) = Grad(6) + Ui(2)*D(3)/Kv

                            Grad(7) = Grad(7) + Ui(3)*D(1)/Kv
                            Grad(8) = Grad(8) + Ui(3)*D(2)/Kv
                            Grad(9) = Grad(9) + Ui(3)*D(3)/Kv
                        endif
                    end if  ! Visc model >0
                else
                    Grad(1:9)=0.0_MK
                end if  ! bComputeGrad
            end if ! denominator size or distances too small
        end if ! 
        ! printf("!4.3f !4.3f !4.3f !4.3f !4.3fNewLine",Uout(1),Uout(2),Uout(3),Kv,denominator);  
    end subroutine fUi_SegmentCst_11


end module UISegments


