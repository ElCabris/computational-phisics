! =============================================================================
! precision_kinds.f90 — Precision constants shared by all Fortran modules
!
! Raman & Fan, PRL 104, 087401 (2010).
!
! Using ISO_C_BINDING kinds ensures that Fortran complex arrays passed to/from
! C++ match the layout of std::complex<double> (which is IEEE 754 double-
! precision, i.e., 2 × 8-byte real values, C-interoperable).
! =============================================================================

module precision_kinds
    use iso_c_binding, only: c_double, c_int
    implicit none

    ! Double-precision real kind (matches C double / std::complex<double> real part)
    integer, parameter :: dp = c_double

    ! Default integer kind for sizes and indices (matches C int)
    integer, parameter :: ip = c_int

    ! Complex kind built from dp — use KIND= to avoid precision loss warning
    integer, parameter :: zp = kind( cmplx(1.0_dp, 1.0_dp, kind=dp) )

end module precision_kinds
