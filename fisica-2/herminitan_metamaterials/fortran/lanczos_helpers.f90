! =============================================================================
! lanczos_helpers.f90 — Fortran utility routines for ARPACK post-processing
!
! Raman & Fan, PRL 104, 087401 (2010).
!
! These routines complement the C++ ARPACK driver in eigensolver.cpp by
! providing Fortran-native operations that benefit from compiler optimizations
! on array syntax and BLAS-style inner products.
!
! Routines:
!   recover_eigenvalues  — converts ARPACK Ritz values 1/(ω−σ) back to ω
!   compute_ritz_residual — computes || Ĥ ŷ − ω ŷ ||₂ / |ω|
!   normalize_eigvec      — normalizes a complex vector to unit L₂ norm
! =============================================================================

module lanczos_helpers
    use precision_kinds
    use arpack_interface
    implicit none

contains

    ! ── Recover physical eigenfrequencies from ARPACK Ritz values ────────────
    ! ARPACK in shift-invert mode returns λ = 1/(ω − σ).
    ! This routine inverts: ω = σ + 1/λ.
    !
    ! d(1:nconv)    : Ritz values from zneupd (complex)
    ! sigma         : spectral shift used during znaupd
    ! omega(1:nconv): recovered eigenfrequencies (complex, Im ≈ 0 for lossless)
    subroutine recover_eigenvalues(d, sigma, omega, nconv)
        integer(ip),  intent(in)  :: nconv
        complex(zp),  intent(in)  :: d(nconv)
        complex(zp),  intent(in)  :: sigma
        complex(zp),  intent(out) :: omega(nconv)
        integer(ip) :: j

        do j = 1, nconv
            omega(j) = sigma + (1.0_dp, 0.0_dp) / d(j)
        end do
    end subroutine recover_eigenvalues

    ! ── Compute Ritz residual for one eigenpair ───────────────────────────────
    ! Returns || r ||₂ / |ω|  where  r = Ĥ·ŷ − ω·ŷ.
    !
    ! Ĥ·ŷ must be provided pre-computed (C++ calls matvec before invoking this).
    ! This avoids passing the full sparse matrix to Fortran.
    !
    ! Hy(1:n) : Ĥ·ŷ  (pre-computed matrix-vector product)
    ! y(1:n)  : eigenvector ŷ
    ! omega   : eigenvalue
    real(dp) function ritz_residual(Hy, y, omega, n)
        integer(ip),  intent(in) :: n
        complex(zp),  intent(in) :: Hy(n), y(n)
        complex(zp),  intent(in) :: omega
        complex(zp)  :: r(n)
        real(dp)     :: norm_r, abs_omega
        integer(ip)  :: k

        do k = 1, n
            r(k) = Hy(k) - omega * y(k)
        end do

        norm_r    = sqrt( real( dot_product(r, r), dp ) )
        abs_omega = abs(omega) + 1.0e-30_dp
        ritz_residual = norm_r / abs_omega
    end function ritz_residual

    ! ── Normalize complex vector to unit L₂ norm (in-place) ──────────────────
    subroutine normalize_eigvec(y, n)
        integer(ip),  intent(in)    :: n
        complex(zp),  intent(inout) :: y(n)
        real(dp) :: nrm

        nrm = sqrt( real( dot_product(y, y), dp ) )
        if (nrm > 0.0_dp) y = y / cmplx(nrm, 0.0_dp, kind=zp)
    end subroutine normalize_eigvec

    ! ── Hermitian inner product  <x|y> = Σ conj(x_k) y_k ────────────────────
    complex(zp) function hermitian_dot(x, y, n)
        integer(ip), intent(in) :: n
        complex(zp), intent(in) :: x(n), y(n)
        hermitian_dot = dot_product(x, y)   ! Fortran dot_product conjugates x
    end function hermitian_dot

end module lanczos_helpers
