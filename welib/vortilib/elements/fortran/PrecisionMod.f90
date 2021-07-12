module PrecisionMod
    use iso_c_binding!OTHER-COMPILER
    implicit none
    integer, parameter :: MK = C_DOUBLE  !< MK stand for MyKind : The kind used throughout the program.
    ! Iso c binding pseudo replacement... compiler and platform dependent!
    ! Gfortran and intel have those values, but better to comment and use iso_c_binding above)
    !integer, parameter :: C_INT    = 4 !COMPAQ-COMPILER
    !integer, parameter :: C_FLOAT  = 4 !COMPAQ-COMPILER
    !integer, parameter :: C_DOUBLE = 8 !COMPAQ-COMPILER
    !integer, parameter :: C_CHAR   = 1 !COMPAQ-COMPILER

    integer, parameter :: SP = kind( 1.0e0 ) !< Single Precision
    integer, parameter :: DP = kind( 1.0d0 ) !< Double Precision
    !integer, parameter :: DP = selected_real_kind( 16 ) !< Double Precision
!     integer, parameter :: QP = selected_real_kind (32) !< Quad Precision
    integer, parameter :: QP = selected_real_kind(33, 4931)
    !integer, parameter :: MK = C_DOUBLE       !< MK stand for MyKind : The kind used throughout the program. MODIFIED BY MAKEFILE
    integer, parameter :: CDK = C_DOUBLE       !< C Double Kind
    integer, parameter :: CIK = C_INT          !< C Integer Kind
    integer, parameter :: CIC = C_CHAR         !< C Char Kind

    real(MK),parameter :: PRECISION_EPS =  epsilon(1.0_MK) !< Machine Precision For the given MyKind for problems of scale 1!
!     real(MK),parameter :: PRECISION_UI  = 2.2204460e-10_MK  !< 
    real(MK),parameter :: PRECISION_UI  = PRECISION_EPS/100  !< 
    ! #define FLT_EPSILON 1.19209290e-07F
    ! #define DBL_EPSILON 2.2204460492503131e-16
    ! #define LDBL_EPSILON 1.0842021724855044340075E-19L
    ! #define PRECISION_EPS DBL_EPSILON //< Machine Precision For the given MyKind for problems of scale 1!
   
    ! Pointers
    !integer, parameter :: IPTR_KIND = int_ptr_kind() !INTEL-COMPILER !For compatibility 32/64bits 
    integer, parameter :: IPTR_KIND = 4               !OTHER-COMPILER !For compatibility 32/64bits 
    !integer, parameter :: IPTR_KIND = 4              !COMPAQ-COMPILER!For compatibility 32/64bits


    integer, parameter :: I4 = 4
    
contains

    !> Checks if two reals are equal given machine precision
    logical function precision_equal(x,y)
        real(MK), intent(in) :: x
        real(MK), intent(in) :: y
        precision_equal=abs(x-y)<PRECISION_EPS
    end function

    !> Checks if two reals are different
    logical function precision_different(x,y)
        real(MK), intent(in) :: x
        real(MK), intent(in) :: y
        precision_different=.not.precision_equal(x,y)
    end function
    
    !> Checks if two reals are equal given machine precision
    logical function precision_equal_soft(x,y)
        real(MK), intent(in) :: x
        real(MK), intent(in) :: y
        precision_equal_soft=abs(x-y)<1000*PRECISION_EPS
    end function

    !> Checks if two reals are equal given machine precision
    logical function precision_equal_r4(x,y)
        real(4), intent(in) :: x
        real(4), intent(in) :: y
        precision_equal_r4=abs(x-y)<epsilon(real(1.0,4))
    end function

    !> Checks if two reals are different
    logical function precision_different_r4(x,y)
        real(4), intent(in) :: x
        real(4), intent(in) :: y
        precision_different_r4=.not.precision_equal_r4(x,y)
    end function

    !> Checks if two reals are equal given machine precision
    logical function precision_equal_r8(x,y)
        real(8), intent(in) :: x
        real(8), intent(in) :: y
        precision_equal_r8=abs(x-y)<epsilon(real(1.0,8))
    end function

    !> Checks if two reals are different
    logical function precision_different_r8(x,y)
        real(8), intent(in) :: x
        real(8), intent(in) :: y
        precision_different_r8=.not.precision_equal_r8(x,y)
    end function

    !> Checks if two reals are equal given machine precision
    logical function precision_equal_rdp(x,y)
        double precision, intent(in) :: x
        double precision, intent(in) :: y
        precision_equal_rdp=abs(x-y)<epsilon(1.0d0)
    end function

    !> Checks if two reals are different
    logical function precision_different_rdp(x,y)
        double precision, intent(in) :: x
        double precision, intent(in) :: y
        precision_different_rdp=.not. precision_equal_rdp(x,y)
    end function

    !> Print precision 
    subroutine precision_print()
        print'(A,A)'     ,'Precision : ',precision_string()
        print'(A,5I3)'   ,'Precision types: ',MK,SP,DP,C_FLOAT,C_DOUBLE
        print'(a,E80.64)','Precision eps: ',PRECISION_EPS
        print'(a,E80.64)','Precision ui : ',PRECISION_UI
        print'(A,I3)','Bytes MK: ',sizeof(1.0_MK)!OTHER-COMPILER
        print'(A,I3)','Bytes SP: ',sizeof(1.0_SP)!OTHER-COMPILER
        print'(A,I3)','Bytes DP: ',sizeof(1.0_DP)!OTHER-COMPILER
    end subroutine
    
    integer function sizeOfReal() result(num_bytes)
        !num_bytes=0             !COMPAQ-COMPILER
        num_bytes=int(sizeof(1.0_MK))!OTHER-COMPILER
!         implicit none
!         character, dimension(:), allocatable :: c
!         real(MK) :: r
!         r=0.0_MK
! !         intrinsic size
! !         intrinsic transfer
! !         sizeOfReal = size(transfer(r,c))
!         inquire(iolength=sizeOfReal)r
!          use iso_fortran_env
!         integer :: num_file_storage_units
!         inquire(iolength=num_file_storage_units) 1.0D0
!         num_bytes = num_file_storage_units*FILE_STORAGE_SIZE/8
!         write(*,*) "double has size: ", num_bytes
!         return
    end function sizeOfReal 


    !> Print precision 
    function precision_string() result(s)
        character(len=40) :: s
        if(MK==SP) then
            if(SP==C_FLOAT) then
                s="single precision - C float"
            else
                s="single precision"
            endif
        elseif(MK==DP) then
            if(SP==C_DOUBLE) then
                s="double precision - C double"
            else
                s="double precision"
            endif
        else
            if(MK==C_FLOAT) then
            s="C float"
            elseif(MK==C_DOUBLE) then
                s="C double"
            else
                s="unknown"
            endif
        endif
    end function

end module PrecisionMod
