
subroutine fui_sgmt_c(CPs, ioff, iCP_beg, nCPs_loc, nCPs, iComputeGrad, Uind_out, Grad_out, &
        nSgmtTot, nSgmtPTot, SgmtP, Sgmt, iSgmtCurPos, iSgmtStart, &
        SgmtNorm2, SgmtSmoothModel, SgmtSmoothParam, SgmtIntensities)
    use MathConstants 
    use PrecisionMod, only: C_DOUBLE, C_INT, PRECISION_UI
    implicit none
    integer(C_INT), intent(in)                        :: iComputeGrad    !< 
    integer(C_INT), intent(in)                        :: ioff
    integer(C_INT), intent(in)                        :: iCP_beg
    integer(C_INT), intent(in)                        :: nCPs_loc
    integer(C_INT), intent(in)                        :: nCPs
    real(C_DOUBLE), dimension(3,nCPs), intent(in)     :: CPs             !< 
    real(C_DOUBLE), dimension(3,nCPs), intent(inout)  :: Uind_out      !< Induced velocity vector - Side effects!!!
    real(C_DOUBLE), dimension(9,nCPs), intent(inout)  :: Grad_out      !< Gradient vector - Side effects
    integer(C_INT), intent(in)                        :: nSgmtTot
    integer(C_INT), intent(in)                        :: nSgmtPTot
    integer(C_INT), dimension(2,nSgmtTot), intent(in)  :: Sgmt
    real(C_DOUBLE), dimension(3,nSgmtPTot), intent(in) :: SgmtP
    integer(C_INT),intent(in)                         :: iSgmtCurPos
    integer(C_INT),intent(in)                         :: iSgmtStart
    real(C_DOUBLE), dimension(nSgmtTot), intent(in)   :: SgmtNorm2
    integer(C_INT), intent(in)                        :: SgmtSmoothModel
    real(C_DOUBLE), dimension(nSgmtTot), intent(in)   :: SgmtSmoothParam
    real(C_DOUBLE), dimension(nSgmtTot), intent(in)   :: SgmtIntensities
    ! Variables
    integer(C_INT) :: icp,is
    ! --------------------------------------------------------------------------------
    ! --- Variables for inlined C function 
    ! --------------------------------------------------------------------------------
    ! Kindof arguments
    real(C_DOUBLE), dimension(3) :: Uind          !<
    real(C_DOUBLE), dimension(9) :: Grad          !< 
    real(C_DOUBLE), dimension(3) :: P1            !< Point1 of a given segment
    real(C_DOUBLE), dimension(3) :: P2            !< Point2 of a given segment
    real(C_DOUBLE), dimension(3) :: DP1, DP2
    real(C_DOUBLE)               :: norm2_r0
    real(C_DOUBLE)               :: Intensity
    integer(C_INT)               :: visc_model   !< 
    real(C_DOUBLE)               :: visc_param   !< 
    ! Variables declaration 
    real(C_DOUBLE),dimension(3) :: crossprod      !< 
    real(C_DOUBLE)              :: denom          !< 
    real(C_DOUBLE)              :: denominator    !< 
    real(C_DOUBLE)              :: h2             !< Square of h
    real(C_DOUBLE)              :: h              !< Only used by one model
    real(C_DOUBLE)              :: Kv             !< 
    real(C_DOUBLE)              :: norm_a         !< 
    real(C_DOUBLE)              :: norm_b         !< 
    real(C_DOUBLE)              :: norm2_crossprod !< 
    real(C_DOUBLE)              :: xa             !< 
    real(C_DOUBLE)              :: ya             !< 
    real(C_DOUBLE)              :: za             !< 
    real(C_DOUBLE)              :: xb             !< 
    real(C_DOUBLE)              :: yb             !< 
    real(C_DOUBLE)              :: zb             !< 
    real(C_DOUBLE)              :: visc_param2    !< Square of viscous param
    real(C_DOUBLE)              :: exp_value      !< 
    !
    real(MK),parameter :: MIN_EXP_VALUE=-10.0_C_DOUBLE
    real(MK),parameter :: MINDENOM=1e-15_C_DOUBLE
    real(MK),parameter :: MINNORMSIMP=1e-6_C_DOUBLE

    !print*,'Sgmt fortran',(iCP_beg+1),(iCP_beg+1)+nCPS_loc-1,nCPs
    visc_model=SgmtSmoothModel

    !$OMP PARALLEL default(shared)
    !$OMP do private(&
    !$OMP& icp,is,Uind,Grad,P1,P2,DP1,DP2,norm2_r0,Intensity,visc_param,&
    !$OMP& crossprod,denom,denominator,h2,h,Kv,norm_a,norm_b,norm2_crossprod,xa,ya,za,xb,yb,zb,visc_param2,exp_value&
    !$OMP& ) schedule(runtime)
    ! loop on CPs 
    do icp=(iCP_beg+1),(iCP_beg+1)+nCPS_loc-1
        do is=iSgmtStart,iSgmtCurPos ! loop on selected segments 
            ! Kind of arguments for Segment 11
            Uind=0.0_C_DOUBLE
            Grad=0.0_C_DOUBLE 
            P1         = SgmtP(1:3, Sgmt(1,is)) ! Segment extremity points
            P2         = SgmtP(1:3, Sgmt(2,is))
            DP1        = CPs(:,icp)-P1;
            DP2        = CPs(:,icp)-P2
            norm2_r0   = SgmtNorm2(is)
            Intensity  = SgmtIntensities(is)
            visc_param = SgmtSmoothParam(is)
            ! --------------------------------------------------------------------------------
            ! --- inlining of Segment 11 
            ! --------------------------------------------------------------------------------
            ! 
            xa=DP1(1)
            ya=DP1(2)
            za=DP1(3)
            xb=DP2(1)
            yb=DP2(2)
            zb=DP2(3)
            !--------------------------------------------------------------------------------
            !--- Simple test if on an extremity point
            !--------------------------------------------------------------------------------
            if(abs(xa)+abs(ya)+abs(za)<MINNORMSIMP) then
                !Uind(1:3)=0.0_C_DOUBLE
                !Grad(1:9)=0.0_C_DOUBLE
            elseif(abs(xb)+abs(yb)+abs(zb)<MINNORMSIMP) then
                !Uind(1:3)=0.0_C_DOUBLE
                !Grad(1:9)=0.0_C_DOUBLE
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
                    if( visc_model==0) then 
                        !return
                    end if
                    if(iComputeGrad==1) then
                        ! denom=(norm_a*norm_b + xa*xb+ya*yb+za*zb)+MINDENOM
                        ! if(norm_a <MINNORM) norm_a=MINNORM
                        ! if(norm_b <MINNORM) norm_b=MINNORM
                        ! if(denom <MINNORM) denom=MINNORM
                        ! ! For now only one limit for all models...
                        ! h2 = MINNORM2
                        ! Kv = 1.0_C_DOUBLE-exp(-1.25643_C_DOUBLE*h2/visc_param2)
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
                    h2 = 0.0_MK
                    select case ((visc_model)) !
                    case ( 0 ) !! No vortex core model
                        Kv=1.0_C_DOUBLE
                    case ( 1 ) !!Rankine - t<=>rc
                        ! orthogonal distance r1xr2/r0 
                        h2 = norm2_crossprod/ norm2_r0 
                        visc_param2=visc_param**2
                        if (h2<visc_param2) then 
                            Kv=h2/visc_param2
                        else
                            Kv=1.0_C_DOUBLE 
                        end if 
                    case ( 2 ) !!Lamb-Oseen - vsic_param<=>rc 
                        ! BUG WITH OPENMP ??
                        ! orthogonal distance r1xr2/r0 
                        h2 = norm2_crossprod/ norm2_r0 
                        visc_param2=visc_param**2
                        exp_value = -1.25643_C_DOUBLE*(h2)/(visc_param2)
                        if(exp_value<MIN_EXP_VALUE) then
                            Kv = 1.0_C_DOUBLE
                        else
                            Kv = 1.0_C_DOUBLE-exp(exp_value)
                        endif
                    case ( 3 ) !!Vatistas n=2 - visc_param<=>rc 
                        ! orthogonal distance r1xr2/r0 
                        h2 = norm2_crossprod/ norm2_r0 
                        visc_param2=visc_param**2
                        ! h = (norm_a+norm_b)/2; 
                        Kv = h2/sqrt(visc_param2**2+h2**2)
                    case ( 4 ) !!Cut-off radius as function of length: delta^2*norm(r0)^2  
                        Kv=1.0_C_DOUBLE
                        denominator=denominator+visc_param*visc_param*norm2_r0
                    case ( 5 ) !!Cut-off radius dimension fixed: (delta l_0)^2  
                        Kv=1.0_C_DOUBLE
                        denominator=denominator+visc_param*visc_param
                    case ( 33 ) !!Vatistas n=2 - visc_paramt<=>rc  
                        ! See Hoydonck 2012. and Script MainRingSensitivityN
                        ! Not matture enough
                        h = (norm_a+norm_b)/2.0_C_DOUBLE
                        visc_param2=visc_param**2
                        if(h<2*sqrt(norm2_r0)) then
                            ! use Vatistas n=2, normal one
                            h2 = norm2_crossprod/norm2_r0
                        end if
                        Kv = (h2)/sqrt(visc_param2**2+h2**2)
                    case DEFAULT
                        Kv=1.0_C_DOUBLE
                    end select 

                    Kv=Intensity*fourpi_inv*Kv*(norm_a+norm_b)/denominator
                    Uind(1:3) = Kv*crossprod(1:3)
                    !print "(A,3F21.16)","ui   :",Uind(1:3)

                    if (iComputeGrad==1) then 
                        Grad(1:9)=0.0_C_DOUBLE
                        !denom=(norm_a*norm_b + xa*xb+ya*yb+za*zb)+MINDENOM
                        !                     denom=denom+MINDENOM
                        !                     D=-(DP1/norm_a**3+DP2/norm_b**3)&
                        !                         -(1./norm_a+1./norm_b)*1/denom*(DP1/norm_a*norm_b+ DP2/norm_b*norm_a+ DP1+DP2)
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
                        !                         exp_value = -1.25643_C_DOUBLE*(h2)/(visc_param2)
                        !                         if(exp_value<MIN_EXP_VALUE) then
                        !                             ! Far away, we then use inviscid Gradient from above.
                        !                         else
                        !                             ! The Lamb-Oseen kernel (may be different than above Kv)
                        !                             Kv2 = 1.0_C_DOUBLE-exp(exp_value)
                        !                             ! Computing the Grad of Kv2
                        !                             D=0._C_DOUBLE
                        !                             ! Grad(raxrb) ra x rb 
                        !                             ! !!!!! Watcth out, need the segment extremity coordinates, but we don't have it
                        !                             ! It's the opposite formula hence, since xa-xb=x1-x2
                        !                             D(1)=(zb-za)*crossprod(2)+(ya-yb)*crossprod(3);
                        !                             D(2)=(za-zb)*crossprod(1)+(xb-xa)*crossprod(3);
                        !                             D(3)=(yb-ya)*crossprod(1)+(xa-xb)*crossprod(2);
                        !                             D=-2._C_DOUBLE*1.25643_C_DOUBLE/(norm2_r0 *visc_param2)*exp(exp_value)*D;
                        !                             !                         print*,'Kv',Kv
                        !                             !                         print*,'Grad',Grad
                        ! 
                        !                             ! Grad_viscous=Kv2*Grad + u_invisc * D
                        !                             Grad=Kv2*Grad ! Intensity already in Grad and Ui
                        ! Since Ui contains already Kv, we have to divide by it. Alternative: store Ui inviscid
                        !                            print*,'Kv',Kv,visc_param,exp_value,Kv2! this print is here to remind that c code has no grad
                        ! also a singulartiy can occur with the grid.
                        !                             Grad(1) = Grad(1) + Uind(1)*D(1)/Kv
                        !                             Grad(2) = Grad(2) + Uind(1)*D(2)/Kv
                        !                             Grad(3) = Grad(3) + Uind(1)*D(3)/Kv
                        ! 
                        !                             Grad(4) = Grad(4) + Uind(2)*D(1)/Kv
                        !                             !                         print*, 'Term3', Uind(2)*D(1)/Kv
                        !                             Grad(5) = Grad(5) + Uind(2)*D(2)/Kv
                        !                             Grad(6) = Grad(6) + Uind(2)*D(3)/Kv
                        ! 
                        !                             Grad(7) = Grad(7) + Uind(3)*D(1)/Kv
                        !                             Grad(8) = Grad(8) + Uind(3)*D(2)/Kv
                        !                             Grad(9) = Grad(9) + Uind(3)*D(3)/Kv
                        !                         endif
                        !                     end if  ! Visc model >0
                    else
                        Grad(1:9)=0.0_C_DOUBLE
                    end if  ! iComputeGrad
                end if ! denominator size or distances too small
            end if ! 
            ! --------------------------------------------------------------------------------
            ! --- END inlining of Segment 11 
            ! --------------------------------------------------------------------------------
            !icp-icp_beg+1
            Uind_out(1:3,ioff+icp-icp_beg) = Uind_out(1:3,ioff+icp-icp_beg)+Uind(1:3)
            if(iComputeGrad==1) then
                Grad_out(1:9,ioff+icp-icp_beg) = Grad_out(1:9,ioff+icp-icp_beg)+Grad(1:9)
            endif
        end do
    enddo
    !$OMP END DO 
    !$OMP END PARALLEL


end subroutine 

