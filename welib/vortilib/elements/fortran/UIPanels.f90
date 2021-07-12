!>
module UIPanels
    use PrecisionMod,only:MK
    implicit none
    integer, save:: nWarning=0

contains
    !>
    subroutine fUi_QuadDoubletCst_1n(P1,P2,P3,P4,CPs,Intensity, SmoothModel, SmoothParam,bComputeGrad,bTrailedOnly,UI,Grad,nCPs)
        !use MathConstants,only : pi,fourpi
        use UISegments, only : fUi_SegmentCst_11
        ! Arguments declarations
        integer, intent(in)                        :: nCPs        ! >
        logical, intent(in)                        :: bComputeGrad                     !< 
        logical, intent(in)                        :: bTrailedOnly
        real(MK), dimension(3), intent(in)         :: P1,P2,P3,P4                      !< 3
        real(MK), dimension(3,nCPs), intent(in)    :: CPs                              !< 3 x n
        real(MK), dimension(3,nCPs), intent(inout) :: UI                               !< 3 x n
        real(MK), dimension(9,nCPs), intent(inout) :: Grad                             !< 9 x n
        real(MK), intent(in)                       :: Intensity                        !< 
        integer, intent(in)                        :: SmoothModel        ! >
        real(MK), intent(in)                       :: SmoothParam                  ! >
        ! Variables declarations
        real(MK), dimension(3)                     :: UItmp
        real(MK), dimension(9)                     :: Gradtmp
        real(MK), dimension(3) :: DP1          !< 
        real(MK), dimension(3) :: DP2          !< 
        integer :: i
        ! UItmp, is incremented at each call
        !          write(*,*) 'P1', P1
        !          write(*,*) 'P2', P2
        !          write(*,*) 'P3', P3
        !          write(*,*) 'P4', P4
        do i=1,nCPs
            !              write(*,*) 'CPs', CPs(:,i)
            ! --------------------------------------------------------------------------------
            ! --- "Shed" segments 
            ! --------------------------------------------------------------------------------
            if(.not.bTrailedOnly) then
                ! 1-2 segment
                DP1=CPs(:,i)-P1; DP2=CPs(:,i)-P2;
                call fUi_SegmentCst_11(DP1, DP2,Intensity,SmoothModel,SmoothParam,bComputeGrad,UItmp,Gradtmp)
                UI(:,i)=UI(:,i)+UItmp
                Grad(:,i)=Grad(:,i)+Gradtmp
                !              write(*,*) '1', UItmp, UI(:,i)
                ! 3-4 segment
                DP1=CPs(:,i)-P3; DP2=CPs(:,i)-P4 ;
                call  fUi_SegmentCst_11(DP1,DP2,Intensity,SmoothModel,SmoothParam,bComputeGrad,UItmp,Gradtmp)
                UI(:,i)=UI(:,i)+UItmp
                Grad(:,i)=Grad(:,i)+Gradtmp
                !              write(*,*) '3', UItmp, UI(:,i)
            endif
            ! --------------------------------------------------------------------------------
            ! --- "Trailed" segments 
            ! --------------------------------------------------------------------------------
            ! 2-3 segment
            DP1=CPs(:,i)-P2; DP2=CPs(:,i)-P3;
            call  fUi_SegmentCst_11(DP1, DP2, Intensity, SmoothModel , SmoothParam, bComputeGrad, UItmp,Gradtmp)
            UI(:,i)=UI(:,i)+UItmp
            Grad(:,i)=Grad(:,i)+Gradtmp
            !              write(*,*) '2', UItmp, UI(:,i)
            ! 4-1 segment
            DP1=CPs(:,i)-P4; DP2=CPs(:,i)-P1;
            call  fUi_SegmentCst_11(DP1, DP2, Intensity, SmoothModel , SmoothParam, bComputeGrad, UItmp,Gradtmp)
            UI(:,i)=UI(:,i)+UItmp
            Grad(:,i)=Grad(:,i)+Gradtmp
            !              write(*,*) '4', UItmp, UI(:,i)
            !              UI(:,i)=UItmp
        enddo
    end subroutine  fUi_QuadDoubletCst_1n


    !> 
    subroutine fUi_QuadSourcePlaneCst_11(CP,Sigma,xi,eta,RefPoint,TransfoMat,Viscous,bComputeGrad,UI,Grad)
        use MathConstants,only : pi,fourpi
        use PrecisionMod, only : precision_equal
        real(MK),parameter :: eps_quadsource=1e-6_MK !!!!!!!!!!!!!!!!!! !< Used if z coordinate close to zero
        ! Arguments declarations 
        logical, intent(in)                  :: bComputeGrad        !< 
        real(MK), dimension(3), intent(in)   :: CP                  !< 3 x n
        real(MK), dimension(3), intent(out)  :: UI                  !< 3 x n
        real(MK), dimension(9), intent(out)  :: Grad                !< 9 x n
        real(MK), intent(in)                 :: Sigma              !< 
        real(MK), dimension(3), intent(in)   :: RefPoint            !< 3 x n !coordinate of panel origin in ref coordinates
        real(MK), dimension(4), intent(in)   :: xi!1,xi2,xi3,xi4     !< 4 x n
        real(MK), dimension(4), intent(in)   :: eta!1,eta2,eta3,eta4 !< 4 x n
        real(MK), dimension(3,3), intent(in) :: TransfoMat          !< 3 x 3 x n
        real(MK), intent(in)                 :: Viscous ! TODO
        ! Variable declarations 
        real(MK), dimension(3,3) :: A                          !< 
        real(MK), dimension(3,3) :: tA                         !< 
        real(MK)                 :: d12, d23, d34, d41         !< 
        real(MK)                 :: m12, m23, m34, m41         !< 
        real(MK)                 :: xi1,  xi2,  xi3,  xi4      !< 
        real(MK)                 :: eta1,  eta2,  eta3,  eta4  !< 
        real(MK)                 :: e1,  e2,  e3,  e4          !< 
        real(MK)                 :: h1,  h2,  h3,  h4          !< 
        real(MK)                 :: r1,  r2,  r3,  r4          !< 
        real(MK)                 :: RJ12, RJ23, RJ34, RJ41     !< 
        real(MK)                 :: TAN12, TAN23, TAN34, TAN41 !< 
        real(MK), dimension(3) :: V     !< 
        real(MK), dimension(3) :: P    !< 
        real(MK), dimension(3) :: Pe   !< 
        xi1=xi(1)
        xi2=xi(2)
        xi3=xi(3)
        xi4=xi(4)
        eta1=eta(1)
        eta2=eta(2)
        eta3=eta(3)
        eta4=eta(4)
        ! TODO viscous is not used yet 

        Grad=0.0_MK
        V = 0.0_MK 
        !param that are constant for each panel, distances and slopes - The slopes can be if divided by zero NaN => security required 
        d12=sqrt((xi2-xi1)**2+(eta2-eta1)**2)
        d23=sqrt((xi3-xi2)**2+(eta3-eta2)**2)
        d34=sqrt((xi4-xi3)**2+(eta4-eta3)**2)
        d41=sqrt((xi1-xi4)**2+(eta1-eta4)**2)

        ! transform control points in panel coordinate system using matrix - vectorial way 
        A=TransfoMat(1:3,1:3)
        tA=transpose(A)

        P(1:3)=CP(1:3)-RefPoint(1:3)
        Pe = matmul(A,P) ! transfo in element coordinate system, noted x,y,z, but in fact xi eta zeta
        ! scalars
        r1=sqrt((Pe(1)-xi1)**2 + (Pe(2)-eta1)**2 + Pe(3)**2)
        r2=sqrt((Pe(1)-xi2)**2 + (Pe(2)-eta2)**2 + Pe(3)**2)
        r3=sqrt((Pe(1)-xi3)**2 + (Pe(2)-eta3)**2 + Pe(3)**2)
        r4=sqrt((Pe(1)-xi4)**2 + (Pe(2)-eta4)**2 + Pe(3)**2)

        e1=Pe(3)**2 + (Pe(1)-xi1)**2
        e2=Pe(3)**2 + (Pe(1)-xi2)**2
        e3=Pe(3)**2 + (Pe(1)-xi3)**2
        e4=Pe(3)**2 + (Pe(1)-xi4)**2
        ! 
        h1=(Pe(2)-eta1)*(Pe(1)-xi1)
        h2=(Pe(2)-eta2)*(Pe(1)-xi2)
        h3=(Pe(2)-eta3)*(Pe(1)-xi3)
        h4=(Pe(2)-eta4)*(Pe(1)-xi4)

        ! Velocities in element frame 
        ! --------------------------------------------------------------------------------
        ! --- Log term 
        ! --------------------------------------------------------------------------------
        ! Security - Katz Plotkin page 608 appendix D Code 11 
        if ( r1+r2-d12<=0.0_MK .or. d12 <=0.0_MK ) then
            RJ12=0._MK
        else
            RJ12=1/d12 * log((r1+r2-d12)/(r1+r2+d12))
        endif

        if ( r2+r3-d23<=0.0_MK .or. d23 <=0.0_MK ) then
            RJ23=0._MK
        else
            RJ23=1/d23 * log((r2+r3-d23)/(r2+r3+d23))
        endif

        if ( r3+r4-d34<=0.0_MK .or. d34 <=0.0_MK ) then
            RJ34=0._MK
        else
            RJ34=1/d34 * log((r3+r4-d34)/(r3+r4+d34))
        endif

        if ( r4+r1-d41<=0.0_MK .or. d41 <=0.0_MK ) then
            RJ41=0._MK
        else
            RJ41=1/d41 * log((r4+r1-d41)/(r4+r1+d41))
        endif
        ! --------------------------------------------------------------------------------
        ! --- Tan term 
        ! --------------------------------------------------------------------------------
        ! 12
        if (precision_equal(xi2,xi1)) then ! Security - Hess 1962 - page 47 - bottom
            TAN12=0._MK
        else
            m12=(eta2-eta1)/(xi2-xi1)
            if( abs(Pe(3))<eps_quadsource ) then ! case where z is too small, jumps may occur 
                TAN12=pi*aint((mysign(1.0_MK,(m12*e1-h1)) - mysign(1.0_MK,(m12*e2-h2)))/2) ! Security-Hess1962-page47-top
            else
                TAN12= atan((m12*e1-h1)/(Pe(3)*r1)) - atan((m12*e2-h2)/(Pe(3)*r2))
            endif
        endif
        ! 23
        if (precision_equal(xi3,xi2)) then ! Security - Hess 1962 - page 47 - bottom
            TAN23=0._MK
        else
            m23=(eta3-eta2)/(xi3-xi2)
            if( abs(Pe(3))<eps_quadsource ) then ! case where z is too small, jumps may occur 
                TAN23=pi*aint((mysign(1.0_MK,(m23*e2-h2)) - mysign(1.0_MK,(m23*e3-h3)))/2) ! Security-Hess1962-page47-top
            else
                TAN23= atan((m23*e2-h2)/(Pe(3)*r2)) - atan((m23*e3-h3)/(Pe(3)*r3))
            endif
        endif 
        ! 34
        if (precision_equal(xi4,xi3)) then ! Security - Hess 1962 - page 47 - bottom
            TAN34=0._MK
        else
            m34=(eta4-eta3)/(xi4-xi3)
            if( abs(Pe(3))<eps_quadsource ) then ! case where z is too small, jumps may occur 
                TAN34=pi*aint((mysign(1.0_MK,(m34*e3-h3)) - mysign(1.0_MK,(m34*e4-h4)))/2) ! Security-Hess1962-page47-top
            else
                TAN34= atan((m34*e3-h3)/(Pe(3)*r3)) - atan((m34*e4-h4)/(Pe(3)*r4))
            endif
        endif
        ! 41
        if (precision_equal(xi1,xi4)) then ! Security - Hess 1962 - page 47 - bottom
            TAN41=0._MK
        else
            m41=(eta1-eta4)/(xi1-xi4)
            if( abs(Pe(3))<eps_quadsource ) then ! case where z is too small, jumps may occur 
                TAN41=pi*aint((mysign(1.0_MK,(m41*e4-h4)) - mysign(1.0_MK,(m41*e1-h1)))/2) ! Security-Hess1962-page47-top
            else
                TAN41= atan((m41*e4-h4)/(Pe(3)*r4)) - atan((m41*e1-h1)/(Pe(3)*r1))
            endif
        endif


        ! --------------------------------------------------------------------------------
        ! --- Velocity  in Panel frame
        ! --------------------------------------------------------------------------------
        ! Vx 
        V(1)= Sigma/(fourpi)*( (eta2-eta1)*RJ12 + &
            (eta3-eta2)*RJ23 + &
            (eta4-eta3)*RJ34 + &
            (eta1-eta4)*RJ41 )
        ! vy 

        V(2)= Sigma/(fourpi)*( (xi1-xi2) *RJ12 + &
            (xi2-xi3) *RJ23 + &
            (xi3-xi4) *RJ34 + &
            (xi4-xi1) *RJ41 )

        ! Vz 
        V(3)= Sigma/(fourpi)*( ( TAN12 ) + &
            ( TAN23 ) + &
            ( TAN34 ) + &
            ( TAN41 ) )

        ! --------------------------------------------------------------------------------
        ! --- Velocity in Reference frame 
        ! --------------------------------------------------------------------------------
        UI(1:3)=matmul(tA,V(1:3)) ! summation for all panels
        !              write(*,*) 'A(1,:)',A(1,:)
        !              write(*,*) 'A(2,:)',A(2,:)
        !              write(*,*) 'A(3,:)',A(3,:)
        !              write(*,*) 'TAN12', TAN12
        !              write(*,*) 'TAN23', TAN23
        !              write(*,*) 'TAN34', TAN34
        !              write(*,*) 'TAN41', TAN41
        !              write(*,*) 'RJ', RJ12,RJ23,RJ34,RJ41
        !              write(*,*) 'eta', eta
        !              write(*,*) 'xi', xi
        !              write(*,*) 'V', V
        !              write(*,*) 'UI', UI
        !              if (ip==1) then
        !                  STOP
        !              endif
        !  

        if (bComputeGrad) then
          if (nWarning<2) then
            write(*,*)'Panels: Smooth and Grad TODO'
            nWarning=nWarning+1
          endif
            ! Just a fake test to use Viscous
            if (int(Viscous)==-1000) then
            endif
        endif
    end subroutine  fUi_QuadSourcePlaneCst_11


    !> Same as QuadDoubletCst but uses precomputed length l
    subroutine fUi_QuadDoubletCst_11_b(CP,P1,P2,P3,P4,l2,Intensity,SmoothModel,SmoothParam,bComputeGrad,Uind,Grad)
        use UISegments, only : fUi_SegmentCst_11_b
        ! Arguments declarations 
        real(MK), dimension(3), intent(in)  :: CP           !< 
        real(MK), dimension(3), intent(in)  :: P1           !< 
        real(MK), dimension(3), intent(in)  :: P2           !< 
        real(MK), dimension(3), intent(in)  :: P3           !< 
        real(MK), dimension(3), intent(in)  :: P4           !< 
        real(MK), dimension(4), intent(in)  :: l2           !< Norm2 Length of sides
        real(MK),  intent(in)               :: Intensity    !< 
        integer , intent(in)                :: SmoothModel  !< 
        real(MK) , intent(in)               :: SmoothParam  !< 
        logical, intent(in)                 :: bComputeGrad !< 
        real(MK), dimension(3), intent(out) :: Uind         !< No Side effects!!!
        real(MK), dimension(9), intent(out) :: Grad         !< No Side effects!!!
        ! Variable declarations 
        real(MK), dimension(3) :: Uindtmp            !< No Side effects!!!
        real(MK), dimension(9) :: Gradtmp            !< No Side effects!!!
        real(MK), dimension(3) :: DP1, DP2
        ! Side1 - 1-2 segment
        DP1=CP-P1; DP2=CP-P2
        call fUi_SegmentCst_11_b( DP1, DP2, l2(1), Intensity , SmoothModel , SmoothParam , bComputeGrad , Uindtmp , Gradtmp ) 
        Uind(1:3)=Uindtmp(1:3) ! Important first pass
        Grad(1:9)=Gradtmp(1:9)
        !print "(A,3F21.16)","ui :",Uindtmp
        ! Side2 - 2-3 segment
        DP1=CP-P2; DP2=CP-P3
        call fUi_SegmentCst_11_b( DP1 , DP2, l2(2), Intensity , SmoothModel , SmoothParam , bComputeGrad , Uindtmp , Gradtmp ) 
        Uind(1:3)=Uind(1:3)+Uindtmp(1:3)
        Grad(1:9)=Grad(1:9)+Gradtmp(1:9)
        !print "(A,3F21.16)","ui :",Uindtmp
        ! Side3 -3-4 segment
        DP1=CP-P3; DP2=CP-P4
        call fUi_SegmentCst_11_b( DP1 , DP2, l2(3), Intensity , SmoothModel , SmoothParam , bComputeGrad , Uindtmp , Gradtmp ) 
        Uind(1:3)=Uind(1:3)+Uindtmp(1:3)
        Grad(1:9)=Grad(1:9)+Gradtmp(1:9)
        !print "(A,3F21.16)","ui :",Uindtmp
        ! Side4 -4-1 segment
        DP1=CP-P4; DP2=CP-P1
        call fUi_SegmentCst_11_b( DP1, DP2, l2(4), Intensity , SmoothModel , SmoothParam , bComputeGrad , Uindtmp , Gradtmp ) 
        Uind(1:3)=Uind(1:3)+Uindtmp(1:3)
        Grad(1:9)=Grad(1:9)+Gradtmp(1:9)
        !print "(A,3F21.16)","ui :",Uindtmp
    end subroutine  fUi_QuadDoubletCst_11_b

    !> Same as HorseShoe_11 but uses precomputed length l
    subroutine fUi_HorseShoe_11_b(CP,P1,P2,P3,P4,l2,Intensity,SmoothModel,SmoothParam,bComputeGrad,Uind,Grad)
        use UISegments, only : fUi_SegmentCst_11_b
        ! Arguments declarations 
        real(MK), dimension(3), intent(in)  :: CP           !< 
        real(MK), dimension(3), intent(in)  :: P1           !< 
        real(MK), dimension(3), intent(in)  :: P2           !< 
        real(MK), dimension(3), intent(in)  :: P3           !< 
        real(MK), dimension(3), intent(in)  :: P4           !< 
        real(MK), dimension(4), intent(in)  :: l2           !< Norm2 Length of sides
        real(MK),  intent(in)               :: Intensity    !< 
        integer , intent(in)                :: SmoothModel  !< 
        real(MK) , intent(in)               :: SmoothParam  !< 
        logical, intent(in)                 :: bComputeGrad !< 
        real(MK), dimension(3), intent(out) :: Uind         !< No Side effects!!!
        real(MK), dimension(9), intent(out) :: Grad         !< No Side effects!!!
        ! Variable declarations 
        real(MK), dimension(3) :: Uindtmp            !< No Side effects!!!
        real(MK), dimension(9) :: Gradtmp            !< No Side effects!!!
        real(MK), dimension(3) :: DP1          !< 
        real(MK), dimension(3) :: DP2          !< 
        ! Side1 - 1-2 segment
        DP1=CP-P1 ; DP2= CP-P2;
        call fUi_SegmentCst_11_b( DP1, DP2, l2(1), Intensity , SmoothModel , SmoothParam , bComputeGrad , Uindtmp , Gradtmp ) 
        Uind(1:3)=Uindtmp(1:3) ! Important first pass
        Grad(1:9)=Gradtmp(1:9)
        !print "(A,3F21.16)","ui :",Uindtmp
        ! Side2 - 2-3 segment
        DP1=CP-P2 ; DP2= CP-P3;
        call fUi_SegmentCst_11_b( DP1 , DP2, l2(2), Intensity , SmoothModel , SmoothParam , bComputeGrad , Uindtmp , Gradtmp ) 
        Uind(1:3)=Uind(1:3)+Uindtmp(1:3)
        Grad(1:9)=Grad(1:9)+Gradtmp(1:9)
        !print "(A,3F21.16)","ui :",Uindtmp
        ! Side4 -4-1 segment
        DP1=CP-P4 ; DP2= CP-P1;
        call fUi_SegmentCst_11_b( DP1, DP2, l2(4), Intensity , SmoothModel , SmoothParam , bComputeGrad , Uindtmp , Gradtmp ) 
        Uind(1:3)=Uind(1:3)+Uindtmp(1:3)
        Grad(1:9)=Grad(1:9)+Gradtmp(1:9)
        !print "(A,3F21.16)","ui :",Uindtmp
    end subroutine  fUi_HorseShoe_11_b
    !>  
    subroutine fUi_QuadDoubletCst_11(CP,P1,P2,P3,P4,Intensity,SmoothModel,SmoothParam,bComputeGrad,Uind,Grad)
        use UISegments, only : fUi_SegmentCst_11
        ! Arguments declarations 
        real(MK), dimension(3), intent(in)  :: CP           !< 
        real(MK), dimension(3), intent(in)  :: P1           !< 
        real(MK), dimension(3), intent(in)  :: P2           !< 
        real(MK), dimension(3), intent(in)  :: P3           !< 
        real(MK), dimension(3), intent(in)  :: P4           !< 
        real(MK),  intent(in)               :: Intensity    !< 
        integer , intent(in)                :: SmoothModel  !< 
        real(MK) , intent(in)               :: SmoothParam  !< 
        logical, intent(in)                 :: bComputeGrad !< 
        real(MK), dimension(3), intent(out) :: Uind         !< No Side effects!!!
        real(MK), dimension(9), intent(out) :: Grad         !< No Side effects!!!
        ! Variable declarations 
        real(MK), dimension(3) :: Uindtmp    !< No Side effects!!!
        real(MK), dimension(9) :: Gradtmp    !< No Side effects!!!
        real(MK), dimension(3) :: DP1, DP2
        ! Side1 - 1-2 segment
        DP1=CP-P1; DP2=CP-P2
        call fUi_SegmentCst_11 ( DP1, DP2, Intensity , SmoothModel , SmoothParam , bComputeGrad , Uindtmp , Gradtmp ) 
        Uind(1:3)=Uindtmp(1:3) ! Important first pass
        Grad(1:9)=Gradtmp(1:9)
!         print "(A,3F35.30)","ui :",Uindtmp
        ! Side2 - 2-3 segment
        DP1=CP-P2; DP2=CP-P3
        call fUi_SegmentCst_11 ( DP1, DP2, Intensity , SmoothModel , SmoothParam , bComputeGrad , Uindtmp , Gradtmp ) 
        Uind(1:3)=Uind(1:3)+Uindtmp(1:3)
        Grad(1:9)=Grad(1:9)+Gradtmp(1:9)
!         print "(A,3F35.30)","ui :",Uindtmp
        ! Side3 - 3-4 segment
        DP1=CP-P3; DP2=CP-P4
        call fUi_SegmentCst_11 ( DP1, DP2, Intensity , SmoothModel , SmoothParam , bComputeGrad , Uindtmp , Gradtmp ) 
        Uind(1:3)=Uind(1:3)+Uindtmp(1:3)
        Grad(1:9)=Grad(1:9)+Gradtmp(1:9)
!         print "(A,3F35.30)","ui :",Uindtmp
        ! Side4 - 4-1 segment
        DP1=CP-P4; DP2=CP-P1
        call fUi_SegmentCst_11 ( DP1, DP2, Intensity , SmoothModel , SmoothParam , bComputeGrad , Uindtmp , Gradtmp ) 
        Uind(1:3)=Uind(1:3)+Uindtmp(1:3)
        Grad(1:9)=Grad(1:9)+Gradtmp(1:9)
!         print "(A,3F35.30)","ui :",Uindtmp
    end subroutine  fUi_QuadDoubletCst_11


    !> Like constantdoublet /vortex ring but without the 34 segment
    subroutine fUi_HorseShoe_11(CP,P1,P2,P3,P4,Intensity,SmoothModel,SmoothParam,bComputeGrad,Uind,Grad)
        use UISegments, only : fUi_SegmentCst_11
        ! Arguments declarations 
        real(MK), dimension(3), intent(in)  :: CP           !< 
        real(MK), dimension(3), intent(in)  :: P1           !< 
        real(MK), dimension(3), intent(in)  :: P2           !< 
        real(MK), dimension(3), intent(in)  :: P3           !< 
        real(MK), dimension(3), intent(in)  :: P4           !< 
        real(MK),  intent(in)               :: Intensity    !< 
        integer , intent(in)                :: SmoothModel  !< 
        real(MK) , intent(in)               :: SmoothParam  !< 
        logical, intent(in)                 :: bComputeGrad !< 
        real(MK), dimension(3), intent(out) :: Uind         !< No Side effects!!!
        real(MK), dimension(9), intent(out) :: Grad         !< No Side effects!!!
        ! Variable declarations 
        real(MK), dimension(3) :: Uindtmp    !< No Side effects!!!
        real(MK), dimension(9) :: Gradtmp    !< No Side effects!!!
        real(MK), dimension(3) :: DP1, DP2
        ! Side1 - 1-2 segment
        DP1=CP-P1; DP2=CP-P2
        call fUi_SegmentCst_11 ( DP1, DP2, Intensity , SmoothModel , SmoothParam , bComputeGrad , Uindtmp , Gradtmp ) 
        Uind(1:3)=Uindtmp(1:3) ! Important first pass
        Grad(1:9)=Gradtmp(1:9)
!         print "(A,3F35.30)","ui :",Uindtmp
        ! Side2 - 2-3 segment
        DP1=CP-P2; DP2=CP-P3
        call fUi_SegmentCst_11 ( DP1, DP2, Intensity , SmoothModel , SmoothParam , bComputeGrad , Uindtmp , Gradtmp ) 
        Uind(1:3)=Uind(1:3)+Uindtmp(1:3)
        Grad(1:9)=Grad(1:9)+Gradtmp(1:9)
!         print "(A,3F35.30)","ui :",Uindtmp
        ! Side4 - 4-1 segment
        DP1=CP-P4; DP2=CP-P1
        call fUi_SegmentCst_11 ( DP1, DP2, Intensity , SmoothModel , SmoothParam , bComputeGrad , Uindtmp , Gradtmp ) 
        Uind(1:3)=Uind(1:3)+Uindtmp(1:3)
        Grad(1:9)=Grad(1:9)+Gradtmp(1:9)
!         print "(A,3F35.30)","ui :",Uindtmp
    end subroutine  fUi_HorseShoe_11



    !> 
    subroutine fUi_QuadSourcePlaneCst(CPs,Sigmas,xi,eta,RefPoint,TransfoMat,Viscous,bComputeGrad,UI,Grad,nCPs,nPanels)
        use MathConstants,only : pi,fourpi
        use PrecisionMod, only : precision_equal
        real(MK),parameter :: eps_quadsource=1e-6_MK !!!!!!!!!!!!!!!!!! !< Used if z coordinate close to zero
        ! Arguments declarations 
        logical, intent(in)                          :: bComputeGrad !< 
        integer, intent(in)                          :: nCPs
        integer, intent(in)                          :: nPanels
        real(MK), dimension(3,nCPs), intent(in)      :: CPs            !< 3 x n
        real(MK), dimension(3,nCPs), intent(inout)   :: UI             !< 3 x n
        real(MK), dimension(9,nCPs), intent(inout)   :: Grad           !< 9 x n
        real(MK), dimension(nPanels), intent(in)     :: Sigmas         !< 
        real(MK), dimension(3,nPanels), intent(in)   :: RefPoint       !< 3 x n !coordinate of panel origin in ref coordinates
        real(MK), dimension(4,nPanels), intent(in)   :: xi             !< 4 x n
        real(MK), dimension(4,nPanels), intent(in)   :: eta            !< 4 x n
        real(MK), dimension(3,3,nPanels), intent(in) :: TransfoMat     !< 3 x 3 x n
        real(MK), dimension(:), intent(in)           :: Viscous ! TODO
        ! Variable declarations 
        real(MK), dimension(3,3) :: A  !<  
        real(MK), dimension(3,3) :: tA  !<  
        real(MK) :: d12  !<  
        real(MK) :: d23  !<  
        real(MK) :: d34  !<  
        real(MK) :: d41  !<  
        real(MK) :: m12  !<  
        real(MK) :: m23  !<  
        real(MK) :: m34  !<  
        real(MK) :: m41  !<  
        real(MK) :: e1    !< 
        real(MK) :: e2    !< 
        real(MK) :: e3    !< 
        real(MK) :: e4    !< 
        real(MK) :: h1    !< 
        real(MK) :: h2    !< 
        real(MK) :: h3    !< 
        real(MK) :: h4    !< 
        real(MK) :: r1    !< 
        real(MK) :: r2    !< 
        real(MK) :: r3    !< 
        real(MK) :: r4    !< 
        real(MK) :: RJ12  !< 
        real(MK) :: RJ23  !< 
        real(MK) :: RJ34  !< 
        real(MK) :: RJ41  !< 
        real(MK) :: TAN12 !< 
        real(MK) :: TAN23 !< 
        real(MK) :: TAN34 !< 
        real(MK) :: TAN41 !< 
        real(MK), dimension(3) :: V     !< 
        real(MK), dimension(3) :: P    !< 
        real(MK), dimension(3) :: Pe   !< 
        ! 
        integer :: ip,icp ! loop index

        ! [ ui vi wi] = fUi_QuadSourcePlaneConstantN(1:4,4:7,7:10,0:5,0:5,0:5,0:5,0:5,0:5,0) 


        ! TODO viscous is not used yet 

        V = 0.0_MK 
        do ip=1,nPanels !loop on panels 
            !param that are constant for each panel, distances and slopes - The slopes can be if divided by zero NaN => security required 
            d12=sqrt((xi(2,ip)-xi(1,ip))**2+(eta(2,ip)-eta(1,ip))**2)
            d23=sqrt((xi(3,ip)-xi(2,ip))**2+(eta(3,ip)-eta(2,ip))**2)
            d34=sqrt((xi(4,ip)-xi(3,ip))**2+(eta(4,ip)-eta(3,ip))**2)
            d41=sqrt((xi(1,ip)-xi(4,ip))**2+(eta(1,ip)-eta(4,ip))**2)

            do icp=1,nCPs ! loop on Control Points
                ! transform control points in panel coordinate system using matrix - vectorial way 
                A=TransfoMat(1:3,1:3,ip)
                tA=transpose(A)

                P(1:3)=CPs(1:3,icp)-RefPoint(1:3,ip)
                !             write(*,*) 'Px', P(1,:)
                !             write(*,*) 'Py', P(2,:)
                !             write(*,*) 'Pz', P(3,:)
                Pe = matmul(A,P) ! transfo in element coordinate system, noted x,y,z, but in fact xi eta zeta
                !             write(*,*) 'Pex', Pe(1,:)
                !             write(*,*) 'Pey', Pe(2,:)
                !             write(*,*) 'Pez', Pe(3,:)
                !             write(*,*) 'P', shape(MP)
                !             write(*,*) 'Pe', shape(MPe)
                !             write(*,*) 'r1', shape(r1)
                !             write(*,*) 'expr', shape(sqrt((MPe(1,:)-xi(1,ip))**2))

                ! scalars
                r1=sqrt((Pe(1)-xi(1,ip))**2 + (Pe(2)-eta(1,ip))**2 + Pe(3)**2)
                r2=sqrt((Pe(1)-xi(2,ip))**2 + (Pe(2)-eta(2,ip))**2 + Pe(3)**2)
                r3=sqrt((Pe(1)-xi(3,ip))**2 + (Pe(2)-eta(3,ip))**2 + Pe(3)**2)
                r4=sqrt((Pe(1)-xi(4,ip))**2 + (Pe(2)-eta(4,ip))**2 + Pe(3)**2)

                e1=Pe(3)**2 + (Pe(1)-xi(1,ip))**2
                e2=Pe(3)**2 + (Pe(1)-xi(2,ip))**2
                e3=Pe(3)**2 + (Pe(1)-xi(3,ip))**2
                e4=Pe(3)**2 + (Pe(1)-xi(4,ip))**2
                ! 
                h1=(Pe(2)-eta(1,ip))*(Pe(1)-xi(1,ip))
                h2=(Pe(2)-eta(2,ip))*(Pe(1)-xi(2,ip))
                h3=(Pe(2)-eta(3,ip))*(Pe(1)-xi(3,ip))
                h4=(Pe(2)-eta(4,ip))*(Pe(1)-xi(4,ip))

                ! Velocities in element frame 
                !             write(*,*) 'RJbefore', RJ12,RJ23,RJ34,RJ41
                ! --------------------------------------------------------------------------------
                ! --- Log term 
                ! --------------------------------------------------------------------------------
                ! Security - Katz Plotkin page 608 appendix D Code 11 
                if ( r1+r2-d12<=0.0_MK .or. d12 <=0.0_MK ) then
                    RJ12=0._MK
                else
                    RJ12=1/d12 * log((r1+r2-d12)/(r1+r2+d12))
                endif

                if ( r2+r3-d23<=0.0_MK .or. d23 <=0.0_MK ) then
                    RJ23=0._MK
                else
                    RJ23=1/d23 * log((r2+r3-d23)/(r2+r3+d23))
                endif

                if ( r3+r4-d34<=0.0_MK .or. d34 <=0.0_MK ) then
                    RJ34=0._MK
                else
                    RJ34=1/d34 * log((r3+r4-d34)/(r3+r4+d34))
                endif

                if ( r4+r1-d41<=0.0_MK .or. d41 <=0.0_MK ) then
                    RJ41=0._MK
                else
                    RJ41=1/d41 * log((r4+r1-d41)/(r4+r1+d41))
                endif
                !             write(*,*) 'RJafter', RJ12,RJ23,RJ34,RJ41
                ! 
                ! --------------------------------------------------------------------------------
                ! --- Tan term 
                ! --------------------------------------------------------------------------------
                ! 12
                if (precision_equal(xi(2,ip),xi(1,ip))) then ! Security - Hess 1962 - page 47 - bottom
                    TAN12=0._MK
                else
                    m12=(eta(2,ip)-eta(1,ip))/(xi(2,ip)-xi(1,ip))
                    if( abs(Pe(3))<eps_quadsource ) then ! case where z is too small, jumps may occur 
                        TAN12=pi*aint((mysign(1.0_MK,(m12*e1-h1)) - mysign(1.0_MK,(m12*e2-h2)))/2) ! Security-Hess1962-page47-top
                    else
                        TAN12= atan((m12*e1-h1)/(Pe(3)*r1)) - atan((m12*e2-h2)/(Pe(3)*r2))
                    endif
                endif
                ! 23
                if (precision_equal(xi(3,ip),xi(2,ip))) then ! Security - Hess 1962 - page 47 - bottom
                    TAN23=0._MK
                else
                    m23=(eta(3,ip)-eta(2,ip))/(xi(3,ip)-xi(2,ip))
                    if( abs(Pe(3))<eps_quadsource ) then ! case where z is too small, jumps may occur 
                        TAN23=pi*aint((mysign(1.0_MK,(m23*e2-h2)) - mysign(1.0_MK,(m23*e3-h3)))/2) ! Security-Hess1962-page47-top
                    else
                        TAN23= atan((m23*e2-h2)/(Pe(3)*r2)) - atan((m23*e3-h3)/(Pe(3)*r3))
                    endif
                endif 
                ! 34
                if (precision_equal(xi(4,ip),xi(3,ip))) then ! Security - Hess 1962 - page 47 - bottom
                    TAN34=0._MK
                else
                    m34=(eta(4,ip)-eta(3,ip))/(xi(4,ip)-xi(3,ip))
                    if( abs(Pe(3))<eps_quadsource ) then ! case where z is too small, jumps may occur 
                        TAN34=pi*aint((mysign(1.0_MK,(m34*e3-h3)) - mysign(1.0_MK,(m34*e4-h4)))/2) ! Security-Hess1962-page47-top
                    else
                        TAN34= atan((m34*e3-h3)/(Pe(3)*r3)) - atan((m34*e4-h4)/(Pe(3)*r4))
                    endif
                endif
                ! 41
                if (precision_equal(xi(1,ip),xi(4,ip))) then ! Security - Hess 1962 - page 47 - bottom
                    TAN41=0._MK
                else
                    m41=(eta(1,ip)-eta(4,ip))/(xi(1,ip)-xi(4,ip))
                    if( abs(Pe(3))<eps_quadsource ) then ! case where z is too small, jumps may occur 
                        TAN41=pi*aint((mysign(1.0_MK,(m41*e4-h4)) - mysign(1.0_MK,(m41*e1-h1)))/2) ! Security-Hess1962-page47-top
                    else
                        TAN41= atan((m41*e4-h4)/(Pe(3)*r4)) - atan((m41*e1-h1)/(Pe(3)*r1))
                    endif
                endif


                ! --------------------------------------------------------------------------------
                ! --- Velocity  in Panel frame
                ! --------------------------------------------------------------------------------
                ! Vx 
                V(1)= Sigmas(ip)/(fourpi)*( (eta(2,ip)-eta(1,ip))*RJ12 + &
                    (eta(3,ip)-eta(2,ip))*RJ23 + &
                    (eta(4,ip)-eta(3,ip))*RJ34 + &
                    (eta(1,ip)-eta(4,ip))*RJ41 )

                !          write(*,*)                     (eta(2,ip)-eta(1,ip))*RJ12 
                !          write(*,*)                     (eta(3,ip)-eta(2,ip))*RJ23 
                !          write(*,*)                     (eta(4,ip)-eta(3,ip))*RJ34 
                !          write(*,*)                     (eta(1,ip)-eta(4,ip))*RJ41 
                !                  write(*,*) Sigmas(ip)/(fourpi) 
                !                  write(*,*) Sigmas(ip)
                !                  write(*,*) fourpi

                ! vy 

                V(2)= Sigmas(ip)/(fourpi)*( (xi(1,ip)-xi(2,ip)) *RJ12 + &
                    (xi(2,ip)-xi(3,ip)) *RJ23 + &
                    (xi(3,ip)-xi(4,ip)) *RJ34 + &
                    (xi(4,ip)-xi(1,ip)) *RJ41 )

                ! Vz 
                V(3)= Sigmas(ip)/(fourpi)*( ( TAN12 ) + &
                    ( TAN23 ) + &
                    ( TAN34 ) + &
                    ( TAN41 ) )

                ! --------------------------------------------------------------------------------
                ! --- Velocity in Reference frame 
                ! --------------------------------------------------------------------------------
                UI(1:3,icp)=UI(1:3,icp)+matmul(tA,V(1:3)) ! summation for all panels
                !              write(*,*) 'A(1,:)',A(1,:)
                !              write(*,*) 'A(2,:)',A(2,:)
                !              write(*,*) 'A(3,:)',A(3,:)
                !              write(*,*) 'TAN12', TAN12
                !              write(*,*) 'TAN23', TAN23
                !              write(*,*) 'TAN34', TAN34
                !              write(*,*) 'TAN41', TAN41
                !              write(*,*) 'RJ', RJ12,RJ23,RJ34,RJ41
                !              write(*,*) 'eta', eta
                !              write(*,*) 'xi', xi
                !              write(*,*) 'V', V
                !              write(*,*) 'UI', UI
                !              if (ip==1) then
                !                  STOP
                !              endif
                !  

            end do ! control points
        end do !end panel loop
        if (bComputeGrad) then
            Grad=Grad+0.0_MK
            write(*,*)'Panels: Smooth and Grad TODO'
            ! Just a fake test to use Viscous
            if (int(Viscous(1))==-1000) then
            endif
        endif
        !         write(*,*) UI
        !         ui=Ui(:,transpose(1))
        !         vi=Ui(:,transpose(2))
        !         wi=Ui(:,transpose(3))
    end subroutine  fUi_QuadSourcePlaneCst 






    elemental real(MK) function mysign(ref,val)
        use PrecisionMod, only: PRECISION_EPS
        real(MK),intent(in) ::ref
        real(MK),intent(in) ::val
        if ( abs(val)>PRECISION_EPS ) then
            mysign=sign(ref,val)
        else
            mysign=1.0_MK
        endif
    endfunction

end module UIPanels


