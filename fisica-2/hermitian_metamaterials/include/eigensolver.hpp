#pragma once
// =============================================================================
// eigensolver.hpp — ARPACK-based implicitly restarted Arnoldi eigensolver
//
// Raman & Fan, PRL 104, 087401 (2010).
//
// ─────────────────────────────────────────────────────────────────────────────
// ALGORITHM: IRLM via ARPACK  (Implicitly Restarted Lanczos/Arnoldi Method)
// ─────────────────────────────────────────────────────────────────────────────
//
// We solve   ω ŷ = Ĥ ŷ   for the nbands smallest positive eigenvalues ω
// using ARPACK's reverse-communication interface.
//
// Because Ĥ is complex Hermitian (and non-Hermitian for the lossy case),
// we use:
//   znaupd  — builds the Arnoldi factorization iteratively
//   zneupd  — extracts Ritz values and Ritz vectors after convergence
//
// NOTE: although ARPACK provides dsaupd/dseupd for real-symmetric problems,
// we cannot use them here because the Bloch-periodic Yee operators produce
// complex entries in Ĥ for generic k ∈ (0, π).  The complex variants
// znaupd/zneupd are required for the full Γ→X path.
//
// ─────────────────────────────────────────────────────────────────────────────
// SHIFT-INVERT MODE (ARPACK mode = 3)
// ─────────────────────────────────────────────────────────────────────────────
//
// Standard mode (mode=1) finds eigenvalues of LARGEST magnitude, which for
// a photonic Hamiltonian are the high-frequency modes — the opposite of what
// we want.  Shift-invert mode (mode=3) transforms the problem:
//
//   Ĥ ŷ = ω ŷ   →   (Ĥ − σI)⁻¹ ŷ = λ ŷ,   λ = 1/(ω − σ)
//
// Requesting the LARGEST |λ| retrieves eigenvalues ω closest to σ.
// Choosing 0 < σ < ω_min_physical ensures:
//   • Null-space modes (ω=0) map to λ = −1/σ < 0  (excluded)
//   • Physical modes with ω > σ map to λ > 0       (retrieved)
//
// At each reverse-communication step, ARPACK requests:
//   y = (Ĥ − σI)⁻¹ x
// This is supplied by the ShiftInvertOp from matrices.hpp.
//
// ─────────────────────────────────────────────────────────────────────────────
// ARPACK REVERSE-COMMUNICATION LOOP
// ─────────────────────────────────────────────────────────────────────────────
//
// IDO codes returned by znaupd:
//   IDO =  0  — first call, initialize
//   IDO = -1  — compute y = OP·x  (OP = (Ĥ−σI)⁻¹ for shift-invert)
//   IDO =  1  — compute y = OP·x  (same as -1 in mode 3)
//   IDO =  2  — compute y = B·x   (B = I for standard eigenvalue problem)
//   IDO = 99  — converged, exit loop and call zneupd
//
// BMAT parameter:
//   'I'  — standard eigenvalue problem (B = I, which is our case since
//           we have already reduced ω A x = B x to ω ŷ = Ĥ ŷ)
//
// WHICH parameter (for znaupd):
//   'LM' — largest magnitude of λ = 1/(ω−σ), i.e., ω closest to σ
//           This is the correct choice for shift-invert mode.
//
// ─────────────────────────────────────────────────────────────────────────────
// LANCZOS VECTOR COUNT  ncv
// ─────────────────────────────────────────────────────────────────────────────
//
// ncv is the number of Lanczos basis vectors kept in the Krylov subspace.
// It must satisfy:   nev + 2  ≤  ncv  ≤  n_dof
//
// Performance rule of thumb:
//   ncv = min(max(2*nev + 1,  20),  n_dof)
//
// A larger ncv:
//   • Reduces the number of ARPACK restarts (fewer matvec calls)
//   • Increases memory: ncv × n_dof × 16 bytes per Lanczos vector
//   • Increases per-iteration cost of the dense sub-problem (ncv × ncv)
//
// For nev = nbands = 15 and ncv = 35:
//   N=20 TM: 35 × 1600 × 16 = 896 KB  (negligible vs. LU storage)
//
// ─────────────────────────────────────────────────────────────────────────────
// CONVERGENCE CRITERION
// ─────────────────────────────────────────────────────────────────────────────
//
// ARPACK declares convergence when the Ritz residual satisfies:
//   || Ĥ ŷⱼ − ωⱼ ŷⱼ ||₂  <  tol × |ωⱼ|
//
// Default tolerance: tol = 1e-10.
// The paper's band frequencies are accurate to ~1% at N=20, so tol=1e-10
// is many orders of magnitude tighter than the discretization error.
// Setting tol=0 uses LAPACK's machine epsilon (dlamch('EPS') ≈ 2.2e-16).
//
// ─────────────────────────────────────────────────────────────────────────────
// FORTRAN INTERFACE
// ─────────────────────────────────────────────────────────────────────────────
//
// ARPACK is a Fortran library.  The C++ code calls it via the Fortran wrappers
// declared in fortran/arpack_interface.f90.  All arrays passed to ARPACK must
// be in Fortran (column-major) order.  The Lanczos vector matrix V (n×ncv)
// is stored column-major in a flat std::vector<Complex>.
// =============================================================================

#ifndef HERMETIC_BANDS_EIGENSOLVER_HPP
#define HERMETIC_BANDS_EIGENSOLVER_HPP

#include "utils.hpp"
#include "sparse.hpp"
#include "matrices.hpp"
#include <vector>
#include <string>

// =============================================================================
// § 1. ARPACK configuration
// =============================================================================

struct ARPACKConfig {
    int   nev     = 15;     // number of eigenvalues to compute (nbands)
    int   ncv     = 0;      // Lanczos vectors; 0 = auto-select (max(2*nev+1,20))
    int   maxiter = 500;    // maximum ARPACK restarts
    Real  tol     = 1e-10;  // relative residual tolerance (0 = machine eps)
    bool  verbose = false;  // print ARPACK iteration progress

    // Compute ncv if not set by user.
    int effective_ncv(int n_dof) const {
        if (ncv > 0) return ncv;
        int auto_ncv = std::max(2 * nev + 1, 20);
        return std::min(auto_ncv, n_dof);
    }
};

// =============================================================================
// § 2. Eigenresult — eigenvalues and back-transformed eigenvectors
// =============================================================================

struct EigenResult {
    RVec             eigenvalues;    // ω values, sorted ascending, length = nconv
    std::vector<CVec> eigenvectors;  // physical x = S⁻¹ ŷ, length = nconv
                                     // eigenvectors[n] has length n_dof
    int              nconv;          // number of converged eigenvalues
    int              niter;          // actual ARPACK iterations used
    Real             max_residual;   // max ||Ĥ ŷ − ω ŷ||/|ω| over converged modes

    // True if all requested nev eigenvalues converged.
    bool fully_converged() const { return nconv >= static_cast<int>(eigenvalues.size()); }
};

// Complex eigenresult for the non-Hermitian (lossy) problem.
struct ComplexEigenResult {
    std::vector<Complex> eigenvalues;    // complex ω = ω_r + i ω_i
    std::vector<CVec>    eigenvectors;
    int                  nconv;
    int                  niter;
    Real                 max_residual;
};

// =============================================================================
// § 3. Primary solver interface
// =============================================================================

// Solve the lossless Hermitian problem  ω ŷ = Ĥ₀ ŷ  using ARPACK znaupd
// in shift-invert mode.  Returns the nbands lowest physical eigenvalues.
//
// Precondition: sys.H_hat is complex Hermitian and sys.sigma > 0.
// The ShiftInvertOp in sys is built and factored inside this function.
EigenResult solve_hermitian(const HermitianSystem& sys,
                            const ARPACKConfig&    cfg);

// Solve the lossy non-Hermitian problem  ω ŷ = (Ĥ₀ + δĤ) ŷ.
// Eigenvalues are complex:  ω = ω_r + i ω_i  (ω_i < 0 for physical decay).
ComplexEigenResult solve_nonhermitian(const NonHermitianSystem& sys,
                                      const ARPACKConfig&       cfg);

// =============================================================================
// § 4. ARPACK reverse-communication driver
//
// This is the low-level loop that calls znaupd repeatedly and dispatches
// the shift-invert solve at each step.  Exposed here for testing/debugging;
// normally called internally by solve_hermitian / solve_nonhermitian.
//
// Parameters match the ARPACK znaupd Fortran interface:
//   op    : shift-invert operator y = (Ĥ−σI)⁻¹ x
//   H_hat : the full Ĥ matrix (used only to compute final residuals)
//   n     : problem dimension (n_dof)
//   cfg   : solver configuration
//   sigma : spectral shift (must match the shift used to factor op)
// =============================================================================

EigenResult run_arpack_znaupd(const ShiftInvertOp& op,
                               const CSRMatrix&     H_hat,
                               Int                  n,
                               const ARPACKConfig&  cfg,
                               Real                 sigma);

// =============================================================================
// § 5. Post-processing utilities
// =============================================================================

// Compute true residuals  ||Ĥ ŷⱼ − ωⱼ ŷⱼ||₂ / |ωⱼ|  for all converged modes.
// Uses the physical eigenvectors (already back-transformed to x = S⁻¹ ŷ).
RVec compute_residuals(const CSRMatrix&  H_hat,
                       const EigenResult& result,
                       const RVec&        sqrtA);

// Sort eigenvalues and eigenvectors by ascending Re[ω].
void sort_by_eigenvalue(EigenResult& result);
void sort_by_eigenvalue(ComplexEigenResult& result);

// Filter out null-space modes: remove any eigenvalue with |ω| < omega_threshold.
// These are the V-in-air spurious modes that survive despite the shift.
// (Normally absent in shift-invert mode, but can appear near k=Γ at low σ.)
void filter_nullspace(EigenResult& result, Real omega_threshold = 1e-6);

// =============================================================================
// § 6. Orthogonality verification  (eq. 13, Raman & Fan)
//
// The lossless eigenmodes satisfy the A-weighted orthogonality condition:
//   <xₘ | A | xₙ> = δₘₙ       (eq. 13)
//
// Equivalently, since ŷ = S x and S = √A:
//   <ŷₘ | ŷₙ> = δₘₙ            (standard L² orthogonality in transformed basis)
//
// The verification computes the Gram matrix G:
//   G_{mn} = <ŷₘ | ŷₙ>         (should be δ_{mn} up to numerical noise)
//
// This is DISTINCT from check_hermitian() in sparse.hpp, which verifies
// Ĥ_{ij} = conj(Ĥ_{ji}) — a matrix property unrelated to eigenvectors.
//
// Returns max_{m≠n} |G_{mn}| and max_m ||G_{mm}| − 1|.
// Well-converged modes should give values below ~10 × tol.
// =============================================================================

struct OrthogonalityReport {
    Real max_off_diagonal  = 0.0;  // max |<ŷₘ|ŷₙ>| for m != n
    Real max_diagonal_error= 0.0;  // max ||<ŷₙ|ŷₙ>| - 1|
    int  n_modes           = 0;    // number of modes checked
    bool passed            = false;
};

OrthogonalityReport verify_modal_orthogonality(const EigenResult& result);

// =============================================================================
// § 7. ARPACK error handling
//
// znaupd and zneupd return an INFO code.  Non-zero codes indicate errors.
// This function converts the numeric code to a human-readable message and
// throws std::runtime_error for fatal codes.
//
// Fatal codes (always throw):
//   INFO = -1  : n_dof must be positive
//   INFO = -3  : ncv must satisfy nev+2 ≤ ncv ≤ n_dof
//   INFO = -6  : bmat must be 'I' or 'G'
//   INFO = -8  : error in internal LAPACK eigenvalue calculation
//   INFO = -9  : starting vector is zero — try a different random seed
//   INFO = -9999: unable to build an Arnoldi factorization (matrix likely singular)
//
// Warning codes (print but do not throw):
//   INFO =  1  : maxiter reached before convergence — increase maxiter or ncv
//   INFO =  3  : no shifts applied — linear combination broke down (rare)
// =============================================================================

void check_arpack_info(int info, const std::string& context);

// =============================================================================
// § 8. Band convergence study utility
//
// Runs the solver at progressively finer grids (N = N_start, N_start+step, …)
// and reports how eigenvalues at a fixed k-point change with N.
// Used to verify spectral convergence as claimed in the paper.
// =============================================================================

struct ConvergenceEntry {
    int  N;
    RVec eigenvalues;
    Real cpu_seconds;
};

std::vector<ConvergenceEntry>
convergence_study(const SimParams& params,
                  Real kx, Real ky,
                  int N_start, int N_end, int N_step,
                  const ARPACKConfig& cfg);

#endif // HERMETIC_BANDS_EIGENSOLVER_HPP
