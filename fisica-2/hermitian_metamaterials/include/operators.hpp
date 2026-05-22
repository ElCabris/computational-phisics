#pragma once
// =============================================================================
// operators.hpp — Discrete curl operators, mass matrix A, operator B,
//                 and Hermitian Hamiltonian Ĥ
//
// Raman & Fan, PRL 104, 087401 (2010).
//
// ─────────────────────────────────────────────────────────────────────────────
// PHYSICAL EQUATIONS (lossless Drude, e^{-iωt} convention)
// ─────────────────────────────────────────────────────────────────────────────
// Maxwell + Drude auxiliary-field formulation (ω₀ = 0, Γ = 0):
//
//   iω μ₀ H  =  ∇×E                           (Faraday)
//   -iω ε∞ E  =  ∇×H − V                       (Ampère − polarization)
//   -iω V    =  ε∞ ωp² E                        (Drude EOM for V = ∂P/∂t)
//
// Rewritten as  ω A x = B x  (eq. 6, Raman & Fan):
//
//   ω (μ₀ H)           =  −i ∇×E               → A_H = μ₀
//   ω (ε∞ E)           =  +i ∇×H − i V         → A_E = ε∞
//   ω ((ε∞ωp²)⁻¹ V)   =  +i E                 → A_V = 1/(ε∞ωp²)
//
// NOTE: units ε₀ = 1 throughout (lengths in a, frequencies in ωp).
//
// ─────────────────────────────────────────────────────────────────────────────
// BLOCK STRUCTURE OF THE STATE VECTOR
// ─────────────────────────────────────────────────────────────────────────────
//
//  TM mode — Drude (ω₀ = 0):
//  ┌──────────┬──────────────────────┬────────────────┐
//  │  Block   │  Field component     │  DOF range     │
//  ├──────────┼──────────────────────┼────────────────┤
//  │  0       │  Hx  (N² entries)    │  [0,    N²)    │
//  │  1       │  Hy  (N² entries)    │  [N²,  2N²)   │
//  │  2       │  Ez  (N² entries)    │  [2N², 3N²)   │
//  │  3       │  Vz  (N² entries)    │  [3N², 4N²)   │
//  └──────────┴──────────────────────┴────────────────┘
//  Total DOF (TM-Drude) = 4N²
//
//  TE mode — Drude (ω₀ = 0):
//  ┌──────────┬──────────────────────┬────────────────┐
//  │  Block   │  Field component     │  DOF range     │
//  ├──────────┼──────────────────────┼────────────────┤
//  │  0       │  Hz  (N² entries)    │  [0,    N²)    │
//  │  1       │  Ex  (N² entries)    │  [N²,  2N²)   │
//  │  2       │  Ey  (N² entries)    │  [2N², 3N²)   │
//  │  3       │  Vx  (N² entries)    │  [3N², 4N²)   │
//  │  4       │  Vy  (N² entries)    │  [4N², 5N²)   │
//  └──────────┴──────────────────────┴────────────────┘
//  Total DOF (TE-Drude) = 5N²
//
// ─────────────────────────────────────────────────────────────────────────────
// BLOCK STRUCTURE OF B (TM-Drude, Hermitian)
// ─────────────────────────────────────────────────────────────────────────────
//
//  x = (Hx, Hy, Ez, Vz)^T
//
//        Hx    Hy       Ez          Vz
//  Hx [   0     0    −i ∂y_f        0  ]
//  Hy [   0     0    +i ∂x_f        0  ]
//  Ez [−i∂y_b† +i∂x_b†   0         −iI ]
//  Vz [   0     0       +iI          0  ]
//
//  where:
//    ∂y_f : forward finite difference in y (maps Ez → Hx-curl contribution)
//    ∂x_f : forward finite difference in x (maps Ez → Hy-curl contribution)
//    ∂y_b†: backward diff in y (adjoint of ∂y_f) → maps Hx → Ez-curl contribution
//    ∂x_b†: backward diff in x (adjoint of ∂x_f) → maps Hy → Ez-curl contribution
//    I    : identity on metallic nodes (0 in air, where Vz is absent)
//
//  The (Ez,Vz) coupling:
//    Ez  row: ... − iI Vz  (Ampère: −iV term)
//    Vz  row: ... + iI Ez  (Drude: +iε∞ωp² E, absorbed into A_V normalization)
//
//  B is Hermitian iff ∂y_b† = (∂y_f)†  and  ∂x_b† = (∂x_f)†.
//  This is guaranteed by the Yee staggering with forward (curl-E) and
//  backward (curl-H) difference conventions.  See yee_grid.hpp.
//
// ─────────────────────────────────────────────────────────────────────────────
// BLOCK STRUCTURE OF B (TE-Drude, Hermitian)
// ─────────────────────────────────────────────────────────────────────────────
//
//  x = (Hz, Ex, Ey, Vx, Vy)^T
//
//        Hz        Ex       Ey       Vx    Vy
//  Hz [   0      −i∂y_b  +i∂x_b      0     0  ]
//  Ex [ −i∂y_f†    0       0         −iI    0  ]
//  Ey [ +i∂x_f†    0       0          0    −iI ]
//  Vx [   0       +iI      0          0     0  ]
//  Vy [   0        0      +iI         0     0  ]
//
//  Notes:
//    TE-Ampère: ω(ε∞ Ex) = i ∂Hz/∂y  →  B_{Ex,Hz} = −i ∂y_b (backward diff)
//              ω(ε∞ Ey) = −i ∂Hz/∂x  →  B_{Ey,Hz} = +i ∂x_b
//    TE-Faraday: ω(μ₀ Hz) = −i(∂Ey/∂x − ∂Ex/∂y)
//              →  B_{Hz,Ex} = −i(−∂y_f†) = +i∂y_f... actually see impl.
//
// ─────────────────────────────────────────────────────────────────────────────
// HERMITIAN HAMILTONIAN  Ĥ  (eq. 12, Raman & Fan)
// ─────────────────────────────────────────────────────────────────────────────
//
//  Define  y = S x  where  S = √A  (block-diagonal, real positive diagonal).
//  Then the eigenvalue problem ω A x = B x becomes:
//
//      ω y = Ĥ y,    with  Ĥ = S⁻¹ B S⁻¹   (eq. 12)
//
//  Ĥ is Hermitian because B is Hermitian and S is real diagonal positive.
//
//  Diagonal entries of S = √A:
//    Block Hx, Hy : √μ₀  per entry
//    Block Ez, Ex, Ey : √ε∞  per entry (position-dependent: ε∞=1 in air)
//    Block Vz, Vx, Vy : 1/√(ε∞ωp²)  per entry (nonzero only in metal)
//
//  Diagonal entries of S⁻¹ = (√A)⁻¹:
//    Block H : 1/√μ₀
//    Block E : 1/√ε∞(r)
//    Block V : √(ε∞(r)ωp²(r))     (= 0 in air since ωp=0 there)
//
//  Ĥ is built by: for each non-zero B_{ij},
//      Ĥ_{ij} = B_{ij} · (S⁻¹)_i · (S⁻¹)_j
//  i.e., left-multiply each row i by (S⁻¹)_i and right-multiply each
//  column j by (S⁻¹)_j.
//
// =============================================================================

#ifndef HERMETIC_BANDS_OPERATORS_HPP
#define HERMETIC_BANDS_OPERATORS_HPP

#include "utils.hpp"
#include "sparse.hpp"
#include "yee_grid.hpp"

// =============================================================================
// § 1. Block-index constants
// =============================================================================

// TM block indices into the state vector x
namespace TM {
    inline constexpr int BLK_HX = 0;
    inline constexpr int BLK_HY = 1;
    inline constexpr int BLK_EZ = 2;
    inline constexpr int BLK_VZ = 3;
    inline constexpr int N_BLOCKS_DRUDE = 4;   // Hx, Hy, Ez, Vz

    inline Int total_dof(Int N2) noexcept { return N_BLOCKS_DRUDE * N2; }
}

// TE block indices into the state vector x
namespace TE {
    inline constexpr int BLK_HZ = 0;
    inline constexpr int BLK_EX = 1;
    inline constexpr int BLK_EY = 2;
    inline constexpr int BLK_VX = 3;
    inline constexpr int BLK_VY = 4;
    inline constexpr int N_BLOCKS_DRUDE = 5;   // Hz, Ex, Ey, Vx, Vy

    inline Int total_dof(Int N2) noexcept { return N_BLOCKS_DRUDE * N2; }
}

// =============================================================================
// § 2. Elemental 1D finite-difference operators (N×N sparse matrices)
//
// Df_x: forward difference in x with Bloch phase e^{+ikx} at wrap.
//   (Df_x)_{flat(i,j), flat(i+1,j)} = +1/Δ  (× e^{+ikx} if i+1 wraps)
//   (Df_x)_{flat(i,j), flat(i,  j)} = −1/Δ
//
// Df_y: forward difference in y with Bloch phase e^{+iky} at wrap.
//   (Df_y)_{flat(i,j), flat(i,j+1)} = +1/Δ  (× e^{+iky} if j+1 wraps)
//   (Df_y)_{flat(i,j), flat(i,  j)} = −1/Δ
//
// Db_x = −(Df_x)†: backward difference in x (adjoint of forward, k → −k).
//   (Db_x)_{flat(i,j), flat(i-1,j)} = −1/Δ  (× e^{-ikx} if i-1 wraps)
//   (Db_x)_{flat(i,j), flat(i,  j)} = +1/Δ
//
// Db_y = −(Df_y)†: backward difference in y.
//
// Adjoint relation (critical for Hermiticity of B):
//   Db_x = −Df_x†,   Db_y = −Df_y†
// This ensures the curl-H operator is the negative adjoint of curl-E.
// =============================================================================

CSRMatrix build_Df_x(const YeeGrid& g, Real kx);
CSRMatrix build_Df_y(const YeeGrid& g, Real ky);

// These are derived from the above via conj_transpose + sign flip.
// Provided as convenience wrappers; internally call build_Df_{x,y}.
CSRMatrix build_Db_x(const YeeGrid& g, Real kx);
CSRMatrix build_Db_y(const YeeGrid& g, Real ky);

// =============================================================================
// § 3. TM curl operators
//
// curl-E for TM: maps Ez (N²) → (Hx, Hy) (2N²)
//   Dimension: (2N²) × (N²)
//
//   Row block 0 (Hx): rows 0 … N²−1
//     C_E_TM[i, j] = (Df_y)_{i,j}       ← ∂Ez/∂y at Hx position
//   Row block 1 (Hy): rows N² … 2N²−1
//     C_E_TM[i+N², j] = −(Df_x)_{i,j}   ← −∂Ez/∂x at Hy position
//                                            (minus from (∇×E)_y = −∂Ez/∂x)
//
// curl-H for TM: maps (Hx, Hy) (2N²) → Ez (N²)
//   Dimension: (N²) × (2N²)
//   = C_E_TM†  (adjoint, with k→−k in Bloch phases)
//
//   Row i (Ez):
//     Col block 0 (Hx): C_H_TM[i, j]     = (Db_y)†_{i,j} = −(Df_y†)_{i,j}
//     Col block 1 (Hy): C_H_TM[i, j+N²]  = −(Db_x)†_{i,j} = +(Df_x†)_{i,j}
//   Combined: (∂Hy/∂x)_backward − (∂Hx/∂y)_backward  (backward diffs)
// =============================================================================

// Build C_E for TM: (2N²) × (N²)
CSRMatrix build_CE_TM(const YeeGrid& g, Real kx, Real ky);

// Build C_H for TM: (N²) × (2N²)
// Implemented as −C_E^† (guarantees exact adjoint relationship).
CSRMatrix build_CH_TM(const YeeGrid& g, Real kx, Real ky);

// =============================================================================
// § 4. TE curl operators
//
// curl-E for TE: maps (Ex, Ey) (2N²) → Hz (N²)
//   Dimension: (N²) × (2N²)
//
//   Row i (Hz):  ∂Ey/∂x − ∂Ex/∂y  at Hz cell-center position
//     Col block 0 (Ex): −(Df_y)_{i,j}   ← −∂Ex/∂y (forward diff in y)
//     Col block 1 (Ey): +(Df_x)_{i,j}   ← +∂Ey/∂x (forward diff in x)
//
// curl-H for TE: maps Hz (N²) → (Ex, Ey) (2N²)
//   Dimension: (2N²) × (N²)
//   = C_E_TE†
//
//   Row block 0 (Ex): (Db_y)_{i,j} = +∂Hz/∂y (backward in y)
//   Row block 1 (Ey): (Db_x)_{i,j} = −∂Hz/∂x (backward in x, with sign)
// =============================================================================

// Build C_E for TE: (N²) × (2N²)
CSRMatrix build_CE_TE(const YeeGrid& g, Real kx, Real ky);

// Build C_H for TE: (2N²) × (N²)
CSRMatrix build_CH_TE(const YeeGrid& g, Real kx, Real ky);

// =============================================================================
// § 5. Metal-region identity operators
//
// In the Drude model, the V field exists ONLY inside the metal (ωp ≠ 0).
// Outside the metal, the V degrees of freedom are absent but retained in
// the vector for uniform block structure; their rows/columns in B and Ĥ
// are simply zero.
//
// I_metal(N², N²): diagonal matrix with 1 at nodes inside the metal, 0 outside.
// I_air(N², N²)  : diagonal matrix with 1 at nodes in air, 0 in metal.
//
// I_metal is used to build the (Ez, Vz) and (Ex/Ey, Vx/Vy) coupling blocks.
// =============================================================================

CSRMatrix build_I_metal(const YeeGrid& g);
CSRMatrix build_I_air(const YeeGrid& g);

// =============================================================================
// § 6. Diagonal mass matrix A  (eq. 6, Raman & Fan)
//
// A is block-diagonal with spatially varying entries in the E and V blocks
// (because ε∞ and ωp are position-dependent: different in metal vs. air).
//
// TM-Drude diagonal entries (4N² total):
//   A[BLK_HX * N² + k] = μ₀             ∀k
//   A[BLK_HY * N² + k] = μ₀             ∀k
//   A[BLK_EZ * N² + k] = ε∞(r_k)        (1 in air, ε∞ in metal)
//   A[BLK_VZ * N² + k] = 1/(ε∞(r_k) ωp²(r_k))   (0 in air → regularized)
//
// TE-Drude diagonal entries (5N² total):
//   A[BLK_HZ * N² + k] = μ₀             ∀k
//   A[BLK_EX * N² + k] = ε∞_x(r_k)     (harmonic avg at Ex position)
//   A[BLK_EY * N² + k] = ε∞_y(r_k)     (harmonic avg at Ey position)
//   A[BLK_VX * N² + k] = 1/(ε∞_x ωp²)  (metal only)
//   A[BLK_VY * N² + k] = 1/(ε∞_y ωp²)  (metal only)
//
// ─────────────────────────────────────────────────────────────────────────────
// NULL-SPACE TREATMENT FOR V DEGREES OF FREEDOM IN AIR
// ─────────────────────────────────────────────────────────────────────────────
//
// In the Drude model, ωp(r) = 0 in air, so A_V = 1/(ε∞ωp²) diverges there.
// Naively setting those entries to ∞ (or leaving them undefined) would corrupt
// the entire S = √A transformation.  The following explains the correct
// treatment and why it preserves the physical spectrum.
//
// STEP 1 — Why V-in-air DOFs are exactly zero in B:
//   The B operator couples V to E via:
//     B_{V←E} = +i I_metal    (Drude EOM: −iωV = ε∞ωp² E, only in metal)
//     B_{E←V} = −i I_metal    (Ampère: −iωε∞E = ∇×H − V,  V≠0 only in metal)
//   In air, I_metal = 0, so both the (V,E) and (E,V) block entries are
//   identically zero for air nodes.  The V-in-air rows AND columns of B
//   are structurally zero — not approximately, but exactly by construction.
//
// STEP 2 — Consequence for Ĥ:
//   Ĥ = S⁻¹ B S⁻¹.  For any air-V index k:
//     row k of Ĥ = S⁻¹[k] × (row k of B) × S⁻¹[·] = S⁻¹[k] × 0 = 0
//   regardless of what value we assign to S⁻¹[k].
//   Identically, column k of Ĥ = 0.
//   Therefore: Ĥ has an exact null space of dimension n_air
//   (number of V DOFs in air), with corresponding eigenvectors that are
//   purely in the air-V subspace and eigenvalue ω = 0.
//
// STEP 3 — Why these are NOT physical modes:
//   Any null-space eigenvector y₀ has all components zero except the
//   air-V block.  Back-transforming to x = S⁻¹ y₀, the physical fields
//   (H, E) are zero and only V-in-air is nonzero.  But V-in-air has no
//   physical meaning (ωp = 0 means no oscillating dipole density there).
//   These are purely numerical artifacts of retaining unused DOFs.
//
// STEP 4 — Preventing contamination of ARPACK's output:
//   ARPACK in standard mode (no shift) would find the ω = 0 null-space
//   modes FIRST when searching for the smallest eigenvalues, burying the
//   physical photonic bands.
//
//   The remedy is SHIFT-INVERT mode with shift σ > 0 (see matrices.hpp):
//     Solve   (Ĥ − σI) z = r   at each Arnoldi step
//     Request the LARGEST eigenvalues of (Ĥ − σI)⁻¹
//   Under this transform, a null-space mode at ω = 0 maps to:
//     λ_null = 1/(0 − σ) = −1/σ  < 0
//   Physical modes near σ map to λ ≈ +1/small, i.e., LARGE POSITIVE values.
//   ARPACK retrieves the largest λ first → physical modes first.
//   The n_air null-space modes are the smallest (most negative) eigenvalues
//   of (Ĥ−σI)⁻¹ and are NEVER requested by the solver.
//
// STEP 5 — The regularization choice:
//   We set A_V[air] = 1 (any finite positive constant) for two reasons:
//   (a) It keeps S = √A well-defined everywhere (no division by zero).
//   (b) It has NO effect on Ĥ or on the physical spectrum, as shown in
//       Step 2: zero rows/columns in B remain zero in Ĥ regardless of S.
//   A_V[air] = 1 is chosen for simplicity; any value > 0 gives the same Ĥ.
// =============================================================================

// Returns the full A diagonal as a real vector (length = total_dof).
RVec build_A_diagonal_TM(const YeeGrid& g);
RVec build_A_diagonal_TE(const YeeGrid& g);

// Build A as a CSR diagonal matrix.
CSRMatrix build_A_TM(const YeeGrid& g);
CSRMatrix build_A_TE(const YeeGrid& g);

// =============================================================================
// § 7. Square-root transformation vectors  S = √A,  S⁻¹ = (√A)⁻¹
//
// Since A is diagonal, √A and its inverse are just element-wise sqrt.
// Stored as real vectors (length = total_dof) for efficiency:
// applying S or S⁻¹ to a vector is an O(N) element-wise multiply.
//
// S⁻¹[k] = 1/√(A[k]),   with the regularization that if A[k]=0, S⁻¹[k]=0.
// =============================================================================

RVec build_sqrtA_TM(const YeeGrid& g);
RVec build_sqrtA_TE(const YeeGrid& g);

RVec build_sqrtA_inv_TM(const YeeGrid& g);
RVec build_sqrtA_inv_TE(const YeeGrid& g);

// Apply S⁻¹ in-place:  x[k] *= s_inv[k]
void apply_sqrtA_inv(CVec& x, const RVec& s_inv);

// Apply S in-place:  x[k] *= s[k]
void apply_sqrtA(CVec& x, const RVec& s);

// =============================================================================
// § 8. Full operator B  (Hermitian, eq. 10-11, Raman & Fan)
//
// Assembly strategy:
//   1. Build elemental curl operators (§ 3–4).
//   2. Embed them into the full (total_dof × total_dof) block matrix B
//      at the correct block offsets (§ 1).
//   3. Add the ±iI coupling blocks between E and V.
//   4. Verify Hermiticity via check_hermitian (optional, --verify-orthogonality).
//
// The B matrix is COMPLEX even for real k (because of the ±i prefactors).
// For k=0 and k=π (Γ and X points), the Bloch phases are real (±1), and B
// becomes a real antisymmetric + purely imaginary symmetric combination —
// still complex Hermitian.
// =============================================================================

CSRMatrix build_B_TM(const YeeGrid& g, Real kx, Real ky);
CSRMatrix build_B_TE(const YeeGrid& g, Real kx, Real ky);

// =============================================================================
// § 9. Hermitian Hamiltonian  Ĥ = S⁻¹ B S⁻¹  (eq. 12, Raman & Fan)
//
// Implementation: for each stored non-zero B[i,j], set
//   Ĥ[i,j] = s_inv[i] * B[i,j] * s_inv[j]
// Since B is already assembled in CSR, this is an O(nnz) operation.
//
// Ĥ has the same sparsity pattern as B.
// Ĥ is Hermitian (by construction from Hermitian B and real positive S⁻¹).
// =============================================================================

CSRMatrix build_H_hat_TM(const YeeGrid& g, Real kx, Real ky);
CSRMatrix build_H_hat_TE(const YeeGrid& g, Real kx, Real ky);

// Dispatch by polarization mode.
CSRMatrix build_H_hat(const YeeGrid& g, Real kx, Real ky);

// =============================================================================
// § 10. Adjoint verification
//
// Verifies that C_H = C_E† holds numerically for the given (kx, ky).
// Returns the max absolute element-wise error.  Should be < 1e-12 for
// any well-assembled curl pair.
// =============================================================================

Real verify_curl_adjoint_TM(const YeeGrid& g, Real kx, Real ky);
Real verify_curl_adjoint_TE(const YeeGrid& g, Real kx, Real ky);

#endif // HERMETIC_BANDS_OPERATORS_HPP
