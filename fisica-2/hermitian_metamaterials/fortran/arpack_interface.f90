! =============================================================================
! arpack_interface.f90 — Fortran module wrapping ARPACK znaupd/zneupd
!
! Raman & Fan, PRL 104, 087401 (2010).
!
! This module provides:
!   (1) Explicit interfaces to znaupd and zneupd so that the compiler can
!       type-check calls from lanczos_helpers.f90.
!   (2) A helper subroutine arpack_error_check() that converts ARPACK INFO
!       codes to human-readable messages.
!
! Why Fortran wrappers around ARPACK rather than calling it directly from C++?
!   ARPACK's reverse-communication interface uses assumed-shape arrays and
!   CHARACTER arguments with Fortran semantics.  Calling from C++ requires
!   careful handling of string termination, array layout, and ABI alignment.
!   The wrappers here provide a clean, type-safe Fortran interface that the
!   C++ code calls via extern "C" declarations (eigensolver.cpp).
!
! NOTE: The actual reverse-communication loop runs in C++ (eigensolver.cpp),
!   not in Fortran.  These wrappers only provide explicit interface blocks
!   and error utilities, not the full loop.
! =============================================================================

module arpack_interface
    use precision_kinds
    implicit none

    ! ── Explicit interface to ARPACK znaupd ──────────────────────────────────
    ! znaupd performs one step of the Implicitly Restarted Arnoldi method for
    ! complex non-symmetric eigenvalue problems.
    interface
        subroutine znaupd(ido, bmat, n, which, nev, tol, resid, ncv, &
                          v, ldv, iparam, ipntr, workd, workl, lworkl, &
                          rwork, info)
            use precision_kinds
            integer(ip), intent(inout) :: ido
            character,   intent(in)    :: bmat
            integer(ip), intent(in)    :: n
            character(2),intent(in)    :: which
            integer(ip), intent(in)    :: nev
            real(dp),    intent(inout) :: tol
            complex(zp), intent(inout) :: resid(n)
            integer(ip), intent(in)    :: ncv
            complex(zp), intent(inout) :: v(ldv, ncv)
            integer(ip), intent(in)    :: ldv
            integer(ip), intent(inout) :: iparam(11)
            integer(ip), intent(inout) :: ipntr(14)
            complex(zp), intent(inout) :: workd(3*n)
            complex(zp), intent(inout) :: workl(lworkl)
            integer(ip), intent(in)    :: lworkl
            real(dp),    intent(inout) :: rwork(ncv)
            integer(ip), intent(inout) :: info
        end subroutine znaupd
    end interface

    ! ── Explicit interface to ARPACK zneupd ──────────────────────────────────
    ! zneupd extracts Ritz values and Ritz vectors after znaupd converges.
    interface
        subroutine zneupd(rvec, howmny, select, d, z, ldz, sigma, workev, &
                          bmat, n, which, nev, tol, resid, ncv, v, ldv,   &
                          iparam, ipntr, workd, workl, lworkl, rwork, info)
            use precision_kinds
            integer(ip), intent(in)    :: rvec
            character,   intent(in)    :: howmny
            logical,     intent(inout) :: select(*)
            complex(zp), intent(out)   :: d(*)
            complex(zp), intent(out)   :: z(ldz, *)
            integer(ip), intent(in)    :: ldz
            complex(zp), intent(in)    :: sigma
            complex(zp), intent(inout) :: workev(3*ncv)
            character,   intent(in)    :: bmat
            integer(ip), intent(in)    :: n
            character(2),intent(in)    :: which
            integer(ip), intent(in)    :: nev
            real(dp),    intent(inout) :: tol
            complex(zp), intent(inout) :: resid(n)
            integer(ip), intent(in)    :: ncv
            complex(zp), intent(inout) :: v(ldv, ncv)
            integer(ip), intent(in)    :: ldv
            integer(ip), intent(inout) :: iparam(11)
            integer(ip), intent(inout) :: ipntr(14)
            complex(zp), intent(inout) :: workd(3*n)
            complex(zp), intent(inout) :: workl(lworkl)
            integer(ip), intent(in)    :: lworkl
            real(dp),    intent(inout) :: rwork(*)
            integer(ip), intent(inout) :: info
        end subroutine zneupd
    end interface

contains

    ! ── Error checker ─────────────────────────────────────────────────────────
    ! Converts ARPACK INFO code to a message string.
    ! Fatal codes (|info| > 1) cause a STOP; warnings print to stderr.
    subroutine arpack_error_check(info, context)
        integer(ip),      intent(in) :: info
        character(len=*), intent(in) :: context

        if (info == 0) return

        select case (info)
        case (1)
            write(*, '(A,A)') '[ARPACK WARNING] ', trim(context) // &
                ': maxiter reached before convergence.'
        case (3)
            write(*, '(A,A)') '[ARPACK WARNING] ', trim(context) // &
                ': no shifts applied; linear combination breakdown.'
        case (-1)
            write(*, '(A)') '[ARPACK FATAL] N must be positive.'
            stop 1
        case (-3)
            write(*, '(A)') '[ARPACK FATAL] NCV must satisfy NEV+2 <= NCV <= N.'
            stop 1
        case (-6)
            write(*, '(A)') '[ARPACK FATAL] BMAT must be I or G.'
            stop 1
        case (-8)
            write(*, '(A)') '[ARPACK FATAL] Error in internal LAPACK eigensolve.'
            stop 1
        case (-9)
            write(*, '(A)') '[ARPACK FATAL] Starting vector is zero.'
            stop 1
        case (-9999)
            write(*, '(A)') '[ARPACK FATAL] Unable to build Arnoldi factorization.'
            stop 1
        case default
            write(*, '(A,I0,A,A)') '[ARPACK FATAL] INFO=', info, &
                ' in ', trim(context)
            stop 1
        end select

    end subroutine arpack_error_check

end module arpack_interface
