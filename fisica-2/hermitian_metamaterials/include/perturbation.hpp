#pragma once
// =============================================================================
// perturbation.hpp — Perturbative loss correction (eq. 15, Raman & Fan)
//                    and modal orthogonality verification (eq. 13)
//
// Raman & Fan, PRL 104, 087401 (2010).
//
// ─────────────────────────────────────────────────────────────────────────────
// PHYSICAL CONTEXT
// ─────────────────────────────────────────────────────────────────────────────
//
// For small damping Γ ≪ ωp, the complex eigenfrequency can be expanded as:
//
//   ω = ω₀ + ω₁ + O(Γ²)
//
// where ω₀ is the lossless (real) eigenfrequency and ω₁ is the first-order
// complex correction due to Γ.
//
// ─────────────────────────────────────────────────────────────────────────────
// EQUATION (15) — PERTURBATIVE IMAGINARY FREQUENCY SHIFT
// ─────────────────────────────────────────────────────────────────────────────
//
// From first-order perturbation theory in the ŷ = S x basis:
//
//   ω₁ = <ŷ₀ | δĤ | ŷ₀>                                     (first order)
//
// where δĤ is the anti-Hermitian perturbation from losses (see matrices.hpp).
// For the Drude model with rate Γ, δĤ has only a (V,V) block:
//
//   δĤ_{VV}[k] = −iΓ / (ε∞(r_k) ωp²(r_k))    in metal
//              =  0                              in air
//
// Substituting and expanding in terms of physical fields (x = S⁻¹ ŷ₀):
//
//   ω₁ = −(iΓ/2) × [ ∫_metal |V₀(r)|² / (ε∞ ωp²) dr ]
//                  / [ ∫ W₀(r) dr ]                        (eq. 15)
//
// SIGN CONVENTION (e^{-iωt}):
//   Im[ω₁] < 0 ↔ exponential decay e^{Im[ω₁]t} ↔ physical damping.
//   The (−iΓ/2) factor ensures Im[ω₁] = −(Γ/2)·ratio < 0 for Γ > 0.
//
// where the modal energy density is:
//
//   W₀(r) = ½ [ ε∞(r)|E₀(r)|² + μ₀|H₀(r)|² + |V₀(r)|²/(ε∞(r) ωp²(r)) ]
//                                                           (eq. 14)
//
// IMPORTANT:
//   The factor ½ in W₀ is the time-averaged electromagnetic energy density.
//   The integral over metal only in the numerator reflects that V₀ ≠ 0 only
//   where ωp ≠ 0.
//   The factor (iΓ/2) versus (iΓ) difference between δĤ and eq. (15)
//   comes from the normalization convention: the eigenvectors ŷ₀ satisfy
//   <ŷ₀|ŷ₀> = 1 (unit norm in ARPACK), while eq. (15) uses the energy
//   normalization ∫ W₀ dr = ½ (for a mode with unit peak amplitude).
//
// ─────────────────────────────────────────────────────────────────────────────
// DISCRETE IMPLEMENTATION
// ─────────────────────────────────────────────────────────────────────────────
//
// On the N×N Yee grid with cell volume Δ² = (1/N)²:
//
//   ∫ f(r) dr  ≈  Δ² × Σ_k f[k]        (quadrature by cell area)
//
// The factor Δ² cancels between numerator and denominator in eq. (15),
// so in practice:
//
//   ω₁ = (iΓ/2) × Σ_k∈metal |V₀[k]|²/(ε∞[k] ωp²[k])
//                / Σ_k W₀[k]
//
// where W₀[k] = ½(ε∞[k]|E₀[k]|² + μ₀|H₀[k]|² + |V₀[k]|²/(ε∞[k]ωp²[k]))
//
// For TM: E₀ = Ez, H₀² = |Hx|² + |Hy|²
// For TE: E₀² = |Ex|² + |Ey|², H₀ = Hz
//
// ─────────────────────────────────────────────────────────────────────────────
// MODAL ORTHOGONALITY VERIFICATION  (eq. 13)
// ─────────────────────────────────────────────────────────────────────────────
//
// The lossless eigenmodes satisfy the A-weighted inner product (eq. 13):
//
//   <xₘ | A | xₙ>  =  δₘₙ
//
// where x = S⁻¹ ŷ and A = S² (diagonal).  This is equivalent to:
//
//   <ŷₘ | ŷₙ>  =  δₘₙ    (standard Euclidean in the ŷ basis)
//
// The function verify_modal_orthogonality() computes the full Gram matrix
// G_{mn} = <ŷₘ|ŷₙ> and returns a report on its deviation from the identity.
//
// This is a DISTINCT check from check_hermitian() (sparse.hpp), which tests
// the matrix operator Ĥ, not the computed eigenvectors.
// =============================================================================

#ifndef HERMETIC_BANDS_PERTURBATION_HPP
#define HERMETIC_BANDS_PERTURBATION_HPP

#include "utils.hpp"
#include "yee_grid.hpp"
#include "eigensolver.hpp"

// =============================================================================
// § 1. Per-mode energy decomposition  (eq. 14)
// =============================================================================

struct ModeEnergy {
    Real E_field;    // ½ ∫ ε∞(r)|E₀(r)|² dr   (electric energy)
    Real H_field;    // ½ ∫ μ₀|H₀(r)|² dr        (magnetic energy)
    Real V_field;    // ½ ∫ |V₀(r)|²/(ε∞ωp²) dr  (oscillator energy, metal only)
    Real total;      // E_field + H_field + V_field  = ∫ W₀ dr

    // Fraction of energy stored in the metallic oscillators.
    Real metal_fraction() const {
        return (total > 0.0) ? V_field / total : 0.0;
    }
};

// Compute modal energy from a physical eigenvector x = S⁻¹ ŷ.
// n2 = N² (DOF per block), mode = TM or TE.
ModeEnergy compute_mode_energy(const CVec&    x,
                                const YeeGrid& g);

// =============================================================================
// § 2. Perturbative loss rate  (eq. 15)
// =============================================================================

struct PerturbResult {
    Real    omega0;         // lossless eigenfrequency ω₀  (real)
    Complex omega1;         // first-order correction ω₁ = Im[ω₁] (imaginary)
    Real    Q_factor;       // quality factor Q = −ω₀ / (2 Im[ω₁])
    ModeEnergy energy;      // energy decomposition for this mode
};

// Compute the perturbative loss correction for one lossless eigenmode.
//
// Parameters:
//   omega0  : lossless eigenfrequency
//   x       : physical eigenvector (back-transformed, length = total_dof)
//   g       : Yee grid (provides ε∞, ωp masks)
//   gamma   : damping rate Γ
//
// Returns ω₁ such that ω ≈ ω₀ + ω₁ with Im[ω₁] < 0 for physical decay.
PerturbResult compute_perturbative_loss(Real           omega0,
                                         const CVec&    x,
                                         const YeeGrid& g,
                                         Real           gamma);

// Compute perturbative corrections for all converged modes.
std::vector<PerturbResult>
compute_all_perturbative(const EigenResult& lossless,
                          const YeeGrid&     g,
                          Real               gamma);

// =============================================================================
// § 3. Exact vs. perturbative comparison  (Figure 2 of the paper)
//
// Compares exact Im[ω] from the non-Hermitian solve with the perturbative
// estimate from eq. (15).
//
// For Γ = 0.01 ωp (as in the paper), the two should agree to within ~1%
// for modes whose energy is primarily in the metallic oscillator (high
// metal_fraction), and may differ more for modes with low metal_fraction
// (weakly lossy modes where higher-order corrections matter).
// =============================================================================

struct LossComparison {
    Real omega_real;          // Re[ω] (common to both methods up to O(Γ²))
    Real im_exact;            // Im[ω] from full non-Hermitian solve
    Real im_perturbative;     // Im[ω] from eq. (15)
    Real relative_error;      // |im_exact − im_perturb| / |im_exact|
};

std::vector<LossComparison>
compare_loss_methods(const EigenResult&        lossless,
                     const ComplexEigenResult&  exact_lossy,
                     const std::vector<PerturbResult>& perturbative);

// =============================================================================
// § 4. Modal orthogonality verification  (eq. 13)
//
// Computes G_{mn} = <ŷₘ|ŷₙ> for all pairs of converged eigenvectors.
// Note: eigenvectors here are ŷ (ARPACK output, unit-norm in L² sense),
//       NOT the physical x = S⁻¹ ŷ.  Back-transforming before calling this
//       function would test a different (physically meaningful) inner product.
//
// Two overloads:
//   (a) From raw ARPACK eigenvectors ŷ  → tests <ŷ|ŷ> = I
//   (b) From physical eigenvectors x   → tests <x|A|x> = I  (eq. 13)
//
// Both should give the same result numerically since ŷ = S x and S is unitary
// up to scaling, but testing both separately isolates where errors arise.
// =============================================================================

// OrthogonalityResult is defined in utils.hpp (included above via eigensolver.hpp).

// Test <ŷₘ|ŷₙ> = δₘₙ  (standard inner product on ARPACK eigenvectors).
OrthogonalityResult
verify_transformed_orthogonality(const EigenResult& result);

// Test <xₘ|A|xₙ> = δₘₙ  (eq. 13, A-weighted inner product on physical fields).
OrthogonalityResult
verify_physical_orthogonality(const EigenResult& result,
                               const RVec&        A_diagonal);

void print_orthogonality_report(const OrthogonalityResult& r, int verbosity);

// =============================================================================
// § 5. Numerator and denominator integrals for eq. (15)
//
// Exposed separately for debugging and for reproducing the paper's plots.
// =============================================================================

// ∫_metal |V₀(r)|² / (ε∞(r) ωp²(r)) dr   — numerator of eq. (15) (÷ iΓ/2)
Real integral_V_metal(const CVec& x, const YeeGrid& g);

// ∫ W₀(r) dr = E_field + H_field + V_field  — denominator of eq. (15)
Real integral_W0(const CVec& x, const YeeGrid& g);

#endif // HERMETIC_BANDS_PERTURBATION_HPP
