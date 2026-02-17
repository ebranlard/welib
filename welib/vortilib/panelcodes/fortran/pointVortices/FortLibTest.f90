!>  Module intended to facilitate the creation of test files
module FortlibTest
    use PrecisionMod, only: MK
    implicit none

    character(len=255),save :: testname
    logical :: bTestStop=.true.
    logical :: bTestPrint=.true.
    logical :: bTestPlot=.false.
            

    interface test_almost_equal; module procedure &
          test_almost_equal_1, &
          test_almost_equal_2, &
          test_almost_equal_4
    end interface
    interface test_equal; module procedure &
          test_equal_i1, &
          test_equal_i0
    end interface
    interface test_vector_almost_equal; module procedure &
          test_vector_almost_equal_2, &
          test_vector_almost_equal_4
    end interface

contains

    subroutine test_success(info,bPrint_in)
        character(len=*), intent(in) :: info
        logical, intent(in), optional  ::  bPrint_in
        if(present(bPrint_in)) then
            if(bPrint_in) then
                write(*,'(A)')'[ OK ] '//trim(testname)//': '//trim(Info)
            endif
        else
            write(*,'(A)')'[ OK ] '//trim(testname)//': '//trim(Info)
        endif
    end subroutine

    subroutine test_fail(info,bPrint_in,bStop_in)
        character(len=*), intent(in) :: info
        logical, intent(in), optional  ::  bPrint_in
        logical, intent(in), optional  ::  bStop_in
        if(present(bPrint_in)) then
            if(bPrint_in) then
                write(*,'(A)')'[FAIL] '//trim(testname)//': '//trim(Info)
            endif
        else
            write(*,'(A)')'[FAIL] '//trim(testname)//': '//trim(Info)
        endif
        if(present(bStop_in)) then
            if(bStop_in) then
                STOP -1 !OTHER-COMPILER
                STOP ! COMPAQ-COMPILER
            endif
        else
            STOP -1 !OTHER-COMPILER
            STOP ! COMPAQ-COMPILER
        endif
    end subroutine


    subroutine  test_vector_almost_equal_2(Var,VecRef,VecTry,MINNORM,bStop,bPrint,bDEBUG,bPassed)
!	    use UtilsOldFortran !COMPAQ-COMPILER
        ! Arguments
        character(len=*), intent(in) :: Var
        real(MK), dimension(:,:), intent(in) :: VecRef         !< 
        real(MK), dimension(:,:), intent(in) :: VecTry         !< 
        real(MK), intent(in) :: MINNORM
        logical, intent(in) :: bStop
        logical, intent(in) :: bPrint
        logical, intent(in) :: bDEBUG
        logical, intent(out),optional :: bPassed
        ! Variables
        real(MK), dimension(:),allocatable :: AbsError      !< 
        real(MK), dimension(:),allocatable :: RelError      !<
        real(MK), dimension(:),allocatable :: NormRef      !<
        integer :: cpt
        integer :: nD, nCPs
        real(MK) :: Max_abs_error, Min_abs_error, Mean_abs_error
        real(MK) :: Max_avg_error, Min_avg_error, Mean_avg_error
        real(MK) :: Max_rel_error, Min_rel_error, Mean_rel_error, MeanNorm
        real(MK) :: FilledRatio
        character(len=20) :: sfmt
        character(len=255) :: InfoRel
        character(len=255) :: InfoAbs
        character(len=255) :: InfoAvg
        integer :: i
        character(len=*),parameter :: RFMTSHORT='F12.5'
        ! 
        nD   = size(VecRef,1)
        nCPs = size(VecRef,2)
        !
        allocate ( NormRef  ( nCPs )  ) 
        allocate ( AbsError ( nCPs )  ) 
        allocate ( RelError ( nCPs )  ) 
        cpt=0
        RelError=0.
        AbsError=0.
        !  For output
        write(sfmt,*) nD
        ! First get the mean norm
        do i=1,nCPs
            NormRef(i)  = norm2(VecRef(1:nD,i))
        enddo
        MeanNorm = sum(NormRef)/nCPs
        !
        do i=1,nCPs
            if(i<10 .and. bDEBUG) then
                print"(A,"// adjustl(sfmt)//RFMTSHORT//")",trim(Var)//' ref:',VecRef(1:nD,i)
                print"(A,"// adjustl(sfmt)//RFMTSHORT//")",trim(Var)//' num:',VecTry(1:nD,i)
                print"(A,"// adjustl(sfmt)//RFMTSHORT//")",trim(Var)//' rel:',VecTry(1:nD,i)/VecRef(1:nD,i)
            endif
            AbsError(i) = norm2(VecRef(1:nD,i)-VecTry(1:nD,i))
            if(NormRef(i)>MeanNorm/100) then
                cpt=cpt+1
                RelError(cpt)=AbsError(i)/NormRef(i)
                if(RelError(cpt)>2 .and. bDEBUG) then
                    !print*,'Huge Rel error',RelError(cpt),norm2(CPs(1:3,i))
                    ! print*,trim(Var)//' ref:',VecRef(1:nD,i)
                    ! print*,trim(Var)//' num:',VecTry(1:nD,i)
                    ! print*,trim(Var)//' rel:',VecTry(1:nD,i)/VecRef(1:nD,i)
                    ! print*,trim(Var)//' abs:',AbsError(i)
                    ! print*,trim(Var)//' Nrm:',NormRef(i)
                    !STOP
                endif
                !print*,'ErrorNorm: ',ErrorNorm(i)
                !print*,'Error : ',VecRef(1:3,i)-VecTry(1:3,i)
                !print*,'Uind bs:',VecRef(1:3,i)
                !print*,'Uind th  :',VecTry  (1:3,i)
            endif
        enddo
        FilledRatio=real(cpt)/nCPs*100 !
        if(FilledRatio<1) then
            print'(A,F8.3)','[WARN] TOO many values around 0 for: '//trim(Var)//' Filled ratio:',FilledRatio
        endif
        !if(FilledRatio<90) then
        !    print*,'A lot of values around 0 for: '//trim(Var)//' Filled ratio:',FilledRatio
        !endif
        Max_Rel_Error  = maxval(RelError(1:cpt))
        Min_Rel_Error  = minval(RelError(1:cpt))
        Mean_Rel_Error = sum(RelError(1:cpt))/cpt
        !
        Max_Abs_Error  = maxval(AbsError)
        Min_Abs_Error  = minval(AbsError)
        Mean_Abs_Error = sum(AbsError)/nCPs
        !
        Max_Avg_Error  =Max_Abs_Error  /MeanNorm*100
        Min_Avg_Error  =Min_Abs_Error  /MeanNorm*100
        Mean_Avg_Error =Mean_Abs_Error /MeanNorm*100


        if (FilledRatio<90) then
            write(InfoRel,'(A,F7.3,A,F7.3,A,F7.3,A,F5.1,A)') trim(Var)//' Rel+ Error, Mean:',Mean_Rel_error*100,&
                & '%, Max:', Max_Rel_error*100,'%, Min:',Min_Rel_error*100,'% - (',FilledRatio,'%)'
        else
            write(InfoRel,'(A,F7.3,A,F7.3,A,F7.3,A)') trim(Var)//' Rel+ Error, Mean:',Mean_Rel_error*100,&
                & '%, Max:', Max_Rel_error*100,'%, Min:',Min_Rel_error*100,'%'
        endif
        write(InfoAvg,'(A,F7.3,A,F7.3,A,F7.3,A)') trim(Var)//' Rel. Error, Mean:',Mean_avg_error,&
            & '% - Max:', Max_avg_error,'% - Min:',Min_avg_error,'%'
        write(InfoAbs,'(A,EN10.2E2,A,EN10.2E2,A,EN10.2E2,A,EN10.2E2)') trim(Var)//' Abs. Error, Mean:',Mean_Abs_error,&
            & ', Max: ', Max_abs_error,', Min: ',Min_abs_error,', Norm: ',MeanNorm
        if(bDEBUG) then
            print'(A)',trim(InfoRel)
            print'(A)',trim(InfoAvg)
            if(bDEBUG) then
                print'(A)',trim(InfoAbs)
            endif
        endif
        if(Mean_Rel_Error>MINNORM) then
            call test_fail(InfoRel)
            if(FilledRatio<90) then
                call test_fail(InfoAvg,bPrint,.false.)
                call test_fail(InfoAbs,bPrint,bStop)
            endif
        else
            call test_success(InfoRel,bPrint)
            if(FilledRatio<90) then
                call test_success(InfoAvg,bPrint)
                call test_success(InfoAbs,bPrint)
            endif
        endif
        !
        if(present(bPassed)) then
            bPassed=.true.
            if(Mean_Rel_Error>MINNORM) then
                bPassed=.false.
                print'(A,A)','[FAIL] '//trim(testname),trim(InfoRel)
                print'(A,A)','[FAIL] '//trim(testname),trim(InfoAvg)
                print'(A,A)','[FAIL] '//trim(testname),trim(InfoAbs)
            endif
        endif
    end subroutine

    subroutine  test_vector_almost_equal_4(Var,VecRef,VecTry,MINNORM,bStop,bPrint,bDEBUG)
!	    use UtilsOldFortran !COMPAQ-COMPILER
        ! Arguments
        character(len=*), intent(in) :: Var
        real(MK), dimension(:,:,:,:), intent(in) :: VecRef         !< 
        real(MK), dimension(:,:,:,:), intent(in) :: VecTry         !< 
        real(MK), intent(in) :: MINNORM
        logical, intent(in) :: bStop,bPrint,bDEBUG
        ! Variables
        real(MK), dimension(:,:),allocatable :: VecRef2    !< 
        real(MK), dimension(:,:),allocatable :: VecTry2   !<
        integer :: p, i,j,k,nD,n1,n2,n3,nCPs
        ! 
        nD = size(VecRef,1); n1 = size(VecRef,2); n2 = size(VecRef,3); n3 = size(VecRef,4);
        nCPs=n1*n2*n3
        allocate ( VecRef2 (nD,nCPs)  ) ; allocate ( VecTry2 (nD,nCPs)  ) 
        p=0
        do k=1,n3; do j=1,n2; do i=1,n1
            p=p+1
            VecRef2(1:nD,p)=VecRef(1:nD,i,j,k)
            VecTry2(1:nD,p)=VecTry(1:nD,i,j,k)
        enddo; enddo; enddo
        call  test_vector_almost_equal_2(Var,VecRef2,VecTry2,MINNORM,bStop,bPrint,bDEBUG)
    end subroutine

    subroutine  test_almost_equal_1(Var,VecRef,VecTry,MINNORM,bStop,bPrint,bPassed)
        ! Arguments
        character(len=*), intent(in) :: Var
        real(MK), dimension(:), intent(in) :: VecRef         !< 
        real(MK), dimension(:), intent(in) :: VecTry         !< 
        real(MK), intent(in) :: MINNORM
        logical, intent(in) :: bStop
        logical, intent(in) :: bPrint
        logical, intent(out),optional :: bPassed
        ! Variables
        character(len=255) :: InfoAbs
        integer :: i,cpt
        real(MK) :: delta
        real(MK) :: delta_cum
        ! 
        cpt=0
        delta_cum=0.0_MK
        do i=1,size(VecRef,1)
            delta=abs(VecRef(i)-VecTry(i))
            delta_cum=delta_cum+delta
            if(delta>MINNORM) then
                cpt=cpt+1
            endif
        enddo
        delta_cum=delta_cum/size(VecRef)

        if(cpt>0) then
            write(InfoAbs,'(A,ES8.1E2,A,ES8.1E2,A,I0)') trim(Var)//' tol: ',MINNORM,', mean: ',delta_cum,' - Failed:',cpt
            call test_fail(InfoAbs,bPrint,bStop)
        else
            write(InfoAbs,'(A,ES8.1E2,A,ES8.1E2)') trim(Var)//' tol: ',MINNORM,', mean: ',delta_cum
            call test_success(InfoAbs,bPrint)
        endif
        if(present(bPassed)) then
            bPassed=(cpt==0)
        endif
    end subroutine
    subroutine  test_almost_equal_2(Var,VecRef,VecTry,MINNORM,bStop,bPrint,bPassed)
        ! Arguments
        character(len=*), intent(in) :: Var
        real(MK), dimension(:,:), intent(in) :: VecRef         !< 
        real(MK), dimension(:,:), intent(in) :: VecTry         !< 
        real(MK), intent(in) :: MINNORM
        logical, intent(in) :: bStop
        logical, intent(in) :: bPrint
        logical, intent(out),optional :: bPassed
        ! Variables
        real(MK), dimension(:),allocatable :: VecRef2    !< 
        real(MK), dimension(:),allocatable :: VecTry2   !<
        integer :: p, i,j,n1,n2,nCPs
        ! 
        n1 = size(VecRef,1); n2 = size(VecRef,2); nCPs=n1*n2
        allocate ( VecRef2 (n1*n2)  ) ; allocate ( VecTry2 (n1*n2)  ) 
        p=0
        do j=1,n2; do i=1,n1
            p=p+1
            VecRef2(p)=VecRef(i,j)
            VecTry2(p)=VecTry(i,j)
        enddo; enddo;
        call  test_almost_equal(Var,VecRef2,VecTry2,MINNORM,bStop,bPrint,bPassed)
    end subroutine

    subroutine  test_almost_equal_4(Var,VecRef,VecTry,MINNORM,bStop,bPrint,bPassed)
!	    use UtilsOldFortran !COMPAQ-COMPILER
        ! Arguments
        character(len=*), intent(in) :: Var
        real(MK), dimension(:,:,:,:), intent(in) :: VecRef         !< 
        real(MK), dimension(:,:,:,:), intent(in) :: VecTry         !< 
        real(MK), intent(in) :: MINNORM
        logical, intent(in) :: bStop,bPrint
        logical, optional, intent(out) :: bPassed
        ! Variables
        real(MK), dimension(:),allocatable :: VecRef2    !< 
        real(MK), dimension(:),allocatable :: VecTry2   !<
        integer :: p, i,j,k,l,nD,n1,n2,n3,nCPs
        ! 
        nD = size(VecRef,1); n1 = size(VecRef,2); n2 = size(VecRef,3); n3 = size(VecRef,4);
        nCPs=n1*n2*n3
        allocate ( VecRef2 (nD*nCPs)  ) ; allocate ( VecTry2 (nD*nCPs)  ) 
        p=0
        do l=1,nD;do k=1,n3; do j=1,n2; do i=1,n1
            p=p+1
            VecRef2(p)=VecRef(l,i,j,k)
            VecTry2(p)=VecTry(l,i,j,k)
        enddo; enddo; enddo; enddo
        call  test_almost_equal(Var,VecRef2,VecTry2,MINNORM,bStop,bPrint,bPassed)
    end subroutine



    subroutine  test_equal_i0(Var,iRef,iTry)
        ! Arguments
        character(len=*), intent(in) :: Var
        integer, intent(in) :: iRef         !< 
        integer, intent(in) :: iTry         !< 
        ! Variables
        character(len=255) :: InfoAbs
        if(iRef/=iTry) then
            write(InfoAbs,'(A,I0,A,I0)') trim(Var),iRef,'/',iTry
            call test_fail(InfoAbs)
            STOP -1 !OTHER-COMPILER
            STOP ! COMPAQ-COMPILER
        else
            write(InfoAbs,'(A,A,I0)') trim(Var),' ok ',iRef
            call test_success(InfoAbs)
        endif
    end subroutine

    subroutine  test_equal_i1(Var,VecRef,VecTry,bTest,bPrintOnly,bPassed)
        ! Arguments
        character(len=*), intent(in) :: Var
        integer, dimension(:), intent(in) :: VecRef         !< 
        integer, dimension(:), intent(in) :: VecTry         !< 
        logical, intent(in) :: bTest
        logical, intent(in) :: bPrintOnly
        logical, intent(out),optional :: bPassed
        ! Variables
        character(len=255) :: InfoAbs
        integer :: i,cpt
        ! 
        cpt=0
        do i=1,size(VecRef)
            if(VecRef(i)/=VecTry(i)) then
                cpt=cpt+1
            endif
        enddo
        if(cpt>0) then
            write(InfoAbs,'(A,I0)') trim(Var)//' Elements different: ',cpt
            if(present(bPassed)) then
                bPassed=.false.
            endif
        else
            write(InfoAbs,'(A)') trim(Var)//' reproduced to identity'
            if(present(bPassed)) then
                bPassed=.true.
            endif
        endif
        if(bPrintOnly) then
            print'(A)',trim(InfoAbs)
        endif
        if(bTest) then
            if(cpt>0) then
                call test_fail(InfoAbs)
                STOP -1 !OTHER-COMPILER
	         	STOP ! COMPAQ-COMPILER
            else
                call test_success(InfoAbs)
            endif
        endif
    end subroutine






    !> This needs debugging
!     subroutine vector_comparison_vtk(Var,VecRef,VecTry,n1,n2,n3)
!         use GridTools, only: natural_order_flat_to_reverse_flat
!         use VTK
!         ! Arguments
!         character(len=*), intent(in)         :: Var
!         real(MK), dimension(:,:), intent(in) :: VecRef             !< 
!         real(MK), dimension(:,:), intent(in) :: VecTry             !< 
!         integer                              :: n1,n2,n3
!         !
!         real(MK),dimension(:,:), pointer     :: Reordered =>null()
!         real(MK),dimension(:,:), pointer     :: Delta =>null()
!         !
!         allocate(Reordered(1:size(VecRef,1),1:size(VecRef,2)))
!         allocate(Delta(1:size(VecRef,1),1:size(VecRef,2)))
!         call natural_order_flat_to_reverse_flat(VecRef,Reordered,n1,n2,n3)
!         call vtk_point_data_vector(Reordered,trim(Var)//'_ref')
!         call natural_order_flat_to_reverse_flat(VecTry,Reordered,n1,n2,n3)
!         call vtk_point_data_vector(Reordered,trim(Var)//'_tst')
!         Delta=VecRef-VecTry
!         call natural_order_flat_to_reverse_flat(Delta,Reordered,n1,n2,n3)
!         call vtk_point_data_vector(Reordered,trim(Var)//'_del')
!         if(associated(Delta))deallocate(Delta)
!         if(associated(Reordered))deallocate(Reordered)
!     end subroutine 


end module FortlibTest
