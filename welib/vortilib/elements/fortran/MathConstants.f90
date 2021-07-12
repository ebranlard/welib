!> 
module MathConstants
    use PrecisionMod, only: MK, SP
    implicit none
    real(MK), parameter :: eps           = 2.220446049250313e-16_MK
    ! Pi related
    real(MK), parameter :: pi            = 3.141592653589793238462643383279502884197169399375105820974944_MK
    real(MK), parameter :: pi_half       = 1.570796326794896619231321691639751442098584699687552910487472_MK
    real(MK), parameter :: twopi         = 6.283185307179586476925286766559005768394338798750211641949889_MK
    real(MK), parameter :: fourpi        = 12.56637061435917295385057353311801153678867759750042328389977_MK
    real(MK), parameter :: twopi_inv     = 0.159154943091895335768883763372514362034459645740456448747667_MK
    real(MK), parameter :: fourpi_inv    = 0.079577471545947667884441881686257181017229822870228224373833_MK
    real(MK), parameter :: four_third_pi = 4.188790204786390984616857844372670512262892532500141094633259_MK
    real(MK), parameter :: rad2deg       = 57.29577951308232087679815481410517033240547246656432154916024_MK
    real(MK), parameter :: deg2rad       = 0.017453292519943295769236907684886127134428718885417254560971_MK
    ! NaN
    integer, parameter  :: NaN=-999999
    real(MK), parameter :: rNaN=-999999._MK
    real(SP), parameter :: NaN_SP=-999999._SP

    interface is_nan
        module procedure iis_nan, ris_nan
    end interface
contains
    logical function iis_nan(x)
        integer,intent(in) ::x
        iis_nan=(x==-NaN)
    end function
    
    logical function ris_nan(x)
        use PrecisionMod, only:PRECISION_EPS
        real(MK),intent(in) ::x
        ris_nan=abs(x-rNaN)<PRECISION_EPS
    end function


end module MathConstants

