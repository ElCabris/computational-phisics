#pragma once
// =============================================================================
// matrices.hpp — High-level problem assembly: lossless Hermitian system,
//                lossy non-Hermitian system, and shift-invert preparation.
//
// Raman & Fan, PRL 104, 087401 (2010).
//
// This module sits one level above operators.hpp.  operators.hpp provides
// the elemental building blocks (curl operators, A, B, Ĥ).  matrices.hpp
// packages them into the complete structures that eigensolver.hpp consumes.
//
// ─────────────────────────────────────────────────────────────────────────────
// TWO PROBLEMS TO SOLVE (paper §III)
// ─────────────────────────────────────────────────────────────────────────────
//
//  Problem 1 — Lossless Hermitian (Γ = 0):
//    ω ŷ = Ĥ₀ ŷ
//    Ĥ₀ = S⁻¹ B₀ S⁻¹  (eq. 12, real positive eigenvalues)
//    Solved with ARPACK znaupd shift-invert.
//
//  Problem 2 — Lossy non-Hermitian (Γ ≠ 0):
//    ω ŷ = (Ĥ₀ + δĤ) ŷ
//    δĤ is anti-Hermitian → eigenvalues become complex: ω = ω_r + i ω_i
//    Solved with ARPACK znaupd (no symmetry assumed).
//
//  Problem 3 — Perturbative losses (Γ ≪ ωp):
//    Uses lossless eigenmodes {ω₀ₙ, y₀ₙ} + first-order correction (eq. 15).
//    Avoids re-solving the non-Hermitian problem at small Γ.
//
// ─────────────────────────────────────────────────────────────────────────────
// NON-HERMITIAN MODIFICATION  δĤ  (eq. 15 derivation)
// ─────────────────────────────────────────────────────────────────────────────
//
// With damping Γ ≠ 0, the Drude EOM becomes (e^{-iωt} convention):
//   (−iω + Γ) V = ε∞ ωp² E
//
// Rewriting as ω A x = B x + δB x (perturbation on B):
//   The V-row of B acquires an extra −iΓ I_metal term on the diagonal:
//     B_NH = B₀  +  δB,   δB_{VV} = −iΓ I_metal
//
// Transforming to the ŷ = S x basis (eq. 12):
//   δĤ = S⁻¹ δB S⁻¹
//   For the (V,V) block:
//     δĤ_{VV} = S⁻¹_V · (−iΓ I_metal) · S⁻¹_V
//              = (−iΓ) · S⁻¹_V² · I_metal
//              = (−iΓ) · (1/(ε∞ωp²)) · I_metal   [since S_V = 1/√(ε∞ωp²)]
//
// This is the perturbation used in eq. (15) of the paper.
// It is purely imaginary (anti-Hermitian), so Im[ω₁] < 0 for Γ > 0
// (energy decay, physical dissipation).
//
// ─────────────────────────────────────────────────────────────────────────────
// SHIFT-INVERT STRATEGY AND NULL-SPACE MANAGEMENT
// ─────────────────────────────────────────────────────────────────────────────
//
// As documented in operators.hpp §6, Ĥ has an exact null space from V-in-air
// DOFs.  The shift-invert scheme:
//
//   Replace Ĥ ŷ = ω ŷ  with  (Ĥ − σI)⁻¹ ŷ = λ ŷ,  λ = 1/(ω − σ)
//
// • Null-space modes (ω = 0) map to λ = 1/(0−σ) = −1/σ < 0.
// • Physical modes near σ map to λ ≈ 1/small > 0 (large positive).
// • ARPACK requests the LARGEST λ → finds physical modes first.
// • Requesting nbands modes never encounters the null space, provided
//   σ > 0 and σ is chosen below the lowest physical band.
//
// Automatic σ selection:
//   For a Drude metal with fill fraction f and ωp = 1,
//   the lowest physical band at Γ starts near ω ≈ 0 for TE and TM.
//   We use σ = 0.01 ωp by default, well above the null space (ω = 0)
//   but below any physical band for the typical resolution.
//   The user can override via SimParams::sigma_shift.
//
// CAUTION: σ must satisfy:
//   0 < σ < ω_min_physical
//   where ω_min_physical is the smallest non-zero photonic band frequency.
//   If σ is too large, physical modes below σ are missed.
//   If σ = 0, the shift-invert factorization encounters the singular Ĥ.
//
// ─────────────────────────────────────────────────────────────────────────────
// EIGENVECTOR BACK-TRANSFORMATION
// ─────────────────────────────────────────────────────────────────────────────
//
// ARPACK returns eigenvectors ŷ in the transformed basis (ŷ = S x).
// Physical field amplitudes are recovered as:
//   x = S⁻¹ ŷ
//
// Components:
//   TM:  H_block = ŷ_H / √μ₀
//        E_block = ŷ_E / √(ε∞(r))
//        V_block = ŷ_V × √(ε∞(r) ωp²(r))   (zero in air automatically)
//
// The inner product in the ŷ basis is Euclidean: <ŷₘ|ŷₙ> = δₘₙ.
// The inner product in the x basis is A-weighted: <xₘ|A|xₙ> = δₘₙ (eq. 13).
// =============================================================================

#ifndef HERMETIC_BANDS_MATRICES_HPP
#define HERMETIC_BANDS_MATRICES_HPP

#include "utils.hpp"
#include "sparse.hpp"
#include "yee_grid.hpp"
#include "operators.hpp"

// =============================================================================
// § 1. HermitianSystem — packaged lossless problem at one k-point
// =============================================================================

struct HermitianSystem {
    CSRMatrix H_hat;        // Ĥ₀ = S⁻¹ B₀ S⁻¹  (complex Hermitian, eq. 12)
    RVec      sqrtA;        // S  diagonal (length = total_dof)
    RVec      sqrtA_inv;    // S⁻¹ diagonal (length = total_dof)
    Int       n_dof;        // total DOF = 4N² (TM) or 5N² (TE)
    Int       N;            // grid resolution
    Real      kx, ky;       // Bloch wave-vector (normalized, in [0,π])
    Real      sigma;        // spectral shift for shift-invert ARPACK

    // Number of air-V null-space modes embedded in Ĥ.
    // These are NOT physical eigenmodes.  They appear at ω = 0 and are
    // excluded from the requested spectrum by shift-invert (see §6 operators.hpp).
    Int       n_nullspace;

    // Convenience: total degrees of freedom per field block (= N²).
    Int n2() const noexcept { return N * N; }
};

// =============================================================================
// § 2. NonHermitianSystem — lossy problem with Γ ≠ 0
// =============================================================================

struct NonHermitianSystem {
    CSRMatrix H_hat_nh;     // Ĥ₀ + δĤ  (complex, NOT Hermitian for Γ≠0)
    RVec      sqrtA;
    RVec      sqrtA_inv;
    Int       n_dof;
    Int       N;
    Real      kx, ky;
    Real      gamma;        // damping rate Γ (in units of ωp)
    Real      sigma;        // spectral shift for shift-invert ARPACK
};

// =============================================================================
// § 3. PerturbationOperator — δĤ for eq. (15)
//
// Stores only the (V,V) diagonal block of δĤ (the only non-zero part).
// The full matrix δĤ is sparse: nnz = n_metal (number of metallic nodes).
//
// δĤ_{VV}[k] = −iΓ / (ε∞(r_k) ωp²(r_k))   for k in metal
//             = 0                              for k in air
// =============================================================================

struct PerturbationOperator {
    CVec delta_H_diag_V;   // diagonal of (V,V) block of δĤ (length = N²)
                            // zero entries for air nodes
    Int  n_metal;           // number of metallic nodes (non-zero entries)
    Real gamma;             // Γ used to build this operator
};

// =============================================================================
// § 4. Assembly functions
// =============================================================================

// Assemble the lossless Hermitian system at wave-vector (kx, ky).
// sigma_shift: shift for shift-invert; if ≤ 0, uses default 0.01 * omega_p.
HermitianSystem assemble_hermitian(const YeeGrid& g,
                                   Real kx, Real ky,
                                   Real sigma_shift = 0.0);

// Assemble the lossy non-Hermitian system at (kx, ky) with damping Γ = gamma.
NonHermitianSystem assemble_nonhermitian(const YeeGrid& g,
                                         Real kx, Real ky,
                                         Real gamma);

// Build the perturbation operator δĤ (used in eq. 15 of the paper).
// Takes the assembled lossless system to reuse S⁻¹ diagonal.
PerturbationOperator build_perturbation_op(const HermitianSystem& sys,
                                           const YeeGrid& g,
                                           Real gamma);

// =============================================================================
// § 5. Eigenvector back-transformation
//
// ARPACK returns ŷ in the S-transformed basis.  Recover physical amplitudes:
//   x = S⁻¹ ŷ    (in-place, overwrites ŷ)
//
// After calling this, columns of eigvecs contain physical field amplitudes.
// =============================================================================

// Back-transform a single eigenvector ŷ → x = S⁻¹ ŷ  (in-place).
void backtransform_eigvec(CVec& y_hat, const RVec& sqrtA_inv);

// Back-transform all eigenvectors stored column-major in a flat vector.
// eigvecs: flat array of n_dof × n_vecs complex values (column-major).
void backtransform_all(std::vector<CVec>& eigvecs, const RVec& sqrtA_inv);

// =============================================================================
// § 6. Field-block extraction helpers
//
// After back-transformation, extract individual field components from x.
// These return const CVec views (as copies) into the relevant block.
//
// TM block layout:  x = [Hx | Hy | Ez | Vz],  each block has N² entries.
// TE block layout:  x = [Hz | Ex | Ey | Vx | Vy], each block has N² entries.
// =============================================================================

// Extract one block from a physical eigenvector x (zero-copy view as span).
// blk: block index (use TM::BLK_* or TE::BLK_* constants).
CVec extract_block(const CVec& x, int blk, Int N2);

// Compute |field|² on the N×N grid from a block (for field output).
// Returns a real N²-vector of squared magnitudes.
RVec field_intensity(const CVec& block);

// =============================================================================
// § 7. Shift-invert operator  (Ĥ − σI)⁻¹
//
// During ARPACK reverse communication, every Arnoldi step requires the action
// of (Ĥ − σI)⁻¹ on a vector.  This is the computationally dominant cost.
//
// ─────────────────────────────────────────────────────────────────────────────
// MEMORY SCALING
// ─────────────────────────────────────────────────────────────────────────────
//
// Dense LU storage for the n_dof × n_dof complex matrix:
//   TM (n_dof = 4N²):  mem = (4N²)² × 16 bytes = 256 N⁴  bytes
//   TE (n_dof = 5N²):  mem = (5N²)² × 16 bytes = 400 N⁴  bytes
//
//  N=16  TM ≈  17 MB   TE ≈  26 MB    ← fast prototyping
//  N=20  TM ≈  41 MB   TE ≈  64 MB    ← paper validation default
//  N=30  TM ≈ 207 MB   TE ≈ 324 MB    ← comfortable on 4 GB
//  N=40  TM ≈ 655 MB   TE ≈ 1.0 GB    ← feasible on 8 GB
//  N=50  TM ≈ 1.6 GB   TE ≈ 2.5 GB    ← borderline on 16 GB
//  N=64  TM ≈ 4.3 GB   TE ≈ 6.7 GB    ← switch to sparse LU
//
// ─────────────────────────────────────────────────────────────────────────────
// COMPUTATIONAL SCALING
// ─────────────────────────────────────────────────────────────────────────────
//
// zgetrf (factorization): O(n_dof³) flops — done ONCE before the ARPACK loop.
//   N=20 TM: (1600)³/3 ≈ 1.4 × 10⁹ flops  ≈  0.5 s  (single-core, LAPACK)
//   N=30 TM: (3600)³/3 ≈ 1.6 × 10¹⁰ flops ≈  5–15 s
//
// zgetrs (solve): O(n_dof²) flops — called once per Arnoldi iteration.
//   Typical ARPACK: ~200 iterations for nbands=15, ncv=45 Lanczos vectors.
//   N=20 TM: 200 × (1600)² ≈ 5.1 × 10⁸ flops ≈  0.2 s  total solve time
//   N=30 TM: 200 × (3600)² ≈ 2.6 × 10⁹ flops ≈  1–3 s  total solve time
//
// Total wall time at default N=20 is dominated by factorization (< 1 s).
// At N=30, factorization and iteration are comparable (< 20 s total).
//
// ─────────────────────────────────────────────────────────────────────────────
// PRACTICAL RESOLUTION LIMITS
// ─────────────────────────────────────────────────────────────────────────────
//
//  N ≤ 30  — Dense shift-invert is optimal: no preconditioning, machine-
//             precision accuracy, deterministic convergence.  Covers the full
//             convergence study of the paper (N=16,20,25,30).
//
//  30 < N ≤ 50  — Dense LU still feasible on 16 GB RAM but factorization
//             time grows rapidly (O(N⁶) total).  Practical but slow (minutes).
//
//  N > 50  — Dense LU is infeasible.  Requires replacing ShiftInvertOp with
//             a sparse direct solver (UMFPACK, PARDISO) or a preconditioned
//             iterative solve (GMRES + ILU).  This is outside the scope of
//             the paper's validation but can be added as a future extension.
//
// ─────────────────────────────────────────────────────────────────────────────
// WHY DENSE SHIFT-INVERT IS CORRECT FOR THE PAPER'S VALIDATION REGIME
// ─────────────────────────────────────────────────────────────────────────────
//
// Raman & Fan (2010) validate their method at N ≈ 20 grid points per period.
// At this resolution:
//
// (a) Spectral accuracy: the photonic bands converge to within ~1% of their
//     continuum values at N=20 for smooth Drude materials.  There is no
//     benefit from N > 40 until the geometry is near-singular (sharp corners).
//
// (b) Conditioning: Ĥ − σI with σ = 0.01 is well-conditioned for all
//     k-points on Γ→X.  LAPACK's dense LU achieves backward error < 10⁻¹³,
//     far below the band discretization error (~10⁻²).  An iterative inner
//     solver would require a very tight tolerance to avoid contaminating the
//     Arnoldi vectors — which would itself be expensive.
//
// (c) Reproducibility: dense LU is deterministic.  An iterative inner solve
//     with loose tolerance produces non-deterministic Arnoldi recurrences,
//     making it harder to reproduce Figure 1 exactly across platforms.
//
// (d) Implementation complexity: dense LU requires only zgetrf + zgetrs from
//     LAPACK (already a dependency).  A sparse solver would add a mandatory
//     third dependency (UMFPACK / SuperLU), violating the project constraint
//     of BLAS + LAPACK + ARPACK only.
//
// Conclusion: for N ≤ 30 (the paper's validation regime), dense shift-invert
// with LAPACK is numerically superior, simpler, and dependency-free.
// =============================================================================

struct ShiftInvertOp {
    // Dense LU factorization of (Ĥ − σI) for small systems.
    // For large systems this will be replaced by a sparse iterative solve.
    // mutable: LAPACK zgetrs_ takes non-const pointers even for read-only solve.
    mutable std::vector<Complex> LU_dense;   // column-major dense matrix (n_dof × n_dof)
    mutable std::vector<int>     ipiv;       // LAPACK pivot indices (int, not Int: LAPACK ABI)
    Int                          n_dof = 0;
    bool                         factored = false;

    // Factor (Ĥ − σI) using LAPACK zgetrf (dense LU).
    // Called once before the ARPACK iteration loop.
    void factor(const CSRMatrix& H_hat, Real sigma);

    // Apply (Ĥ − σI)⁻¹ x → y using LAPACK zgetrs.
    // Called at every ARPACK reverse-communication step.
    void apply(const CVec& x, CVec& y) const;
};

// Build and factor the shift-invert operator for the given system.
ShiftInvertOp build_shift_invert(const HermitianSystem& sys);
ShiftInvertOp build_shift_invert_nh(const NonHermitianSystem& sys);

// =============================================================================
// § 8. System diagnostics
// =============================================================================

// Print a concise summary of the assembled system to stdout.
void print_system_info(const HermitianSystem& sys, int verbosity);

// Estimate the spectral radius of Ĥ (Gershgorin bound).
// Useful for choosing a good initial shift σ.
Real gershgorin_radius(const CSRMatrix& H_hat);

// Count the exact number of structural zeros on the diagonal of Ĥ
// (= n_nullspace, the air-V modes).
Int count_nullspace_dof(const CSRMatrix& H_hat);

#endif // HERMETIC_BANDS_MATRICES_HPP
