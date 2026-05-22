// =============================================================================
// perturbation.cpp — Perturbative loss correction (eq. 15) and
//                    modal orthogonality verification (eq. 13)
//
// Raman & Fan, PRL 104, 087401 (2010).
// =============================================================================

#include "perturbation.hpp"
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <iostream>
#include <iomanip>

// =============================================================================
// § 1. Numerator integral: ∫_metal |V₀|² / (ε∞ ωp²) dr
//
// Discrete: Δ² × Σ_{k∈metal} |V[k]|² / (ε∞[k] ωp²[k])
// The Δ² factor cancels with the denominator, so we return the bare sum.
// =============================================================================

Real integral_V_metal(const CVec& x, const YeeGrid& g) {
    const Int N2  = g.n_dof;
    Real sum = 0.0;

    if (g.params.mode == Polarization::TM) {
        // V block = BLK_VZ
        for (Int k = 0; k < N2; ++k) {
            if (!g.metal_mask[static_cast<std::size_t>(k)]) continue;
            const Real eps_i = g.eps_inf_mask[static_cast<std::size_t>(k)];
            const Real wp_i  = g.omega_p_mask[static_cast<std::size_t>(k)];
            if (wp_i <= 0.0) continue;
            const std::size_t idx =
                static_cast<std::size_t>(block_offset(TM::BLK_VZ, N2) + k);
            sum += std::norm(x[idx]) / (eps_i * wp_i * wp_i);
        }
    } else {
        // TE: two V blocks (Vx and Vy)
        for (Int k = 0; k < N2; ++k) {
            if (!g.metal_mask[static_cast<std::size_t>(k)]) continue;
            const Real eps_x = g.eps_x[static_cast<std::size_t>(k)];
            const Real eps_y = g.eps_y[static_cast<std::size_t>(k)];
            const Real wp    = g.omega_p_mask[static_cast<std::size_t>(k)];
            if (wp <= 0.0) continue;
            const std::size_t vx_idx =
                static_cast<std::size_t>(block_offset(TE::BLK_VX, N2) + k);
            const std::size_t vy_idx =
                static_cast<std::size_t>(block_offset(TE::BLK_VY, N2) + k);
            sum += std::norm(x[vx_idx]) / (eps_x * wp * wp)
                 + std::norm(x[vy_idx]) / (eps_y * wp * wp);
        }
    }
    return sum;
}

// =============================================================================
// § 2. Denominator integral: ∫ W₀ dr  (eq. 14)
//
// W₀ = ½ (ε∞|E|² + μ₀|H|² + |V|²/(ε∞ωp²))
// Returns the bare sum (Δ² cancels with numerator).
// =============================================================================

Real integral_W0(const CVec& x, const YeeGrid& g) {
    const Int N2 = g.n_dof;
    Real W = 0.0;

    if (g.params.mode == Polarization::TM) {
        for (Int k = 0; k < N2; ++k) {
            const std::size_t ks = static_cast<std::size_t>(k);
            const Real eps_k = g.eps_z[ks];
            const Real mu0   = PhysConst::mu0;

            // Electric energy: ½ ε∞ |Ez|²
            W += 0.5 * eps_k
               * std::norm(x[static_cast<std::size_t>(block_offset(TM::BLK_EZ, N2) + k)]);

            // Magnetic energy: ½ μ₀ (|Hx|² + |Hy|²)
            W += 0.5 * mu0
               * std::norm(x[static_cast<std::size_t>(block_offset(TM::BLK_HX, N2) + k)]);
            W += 0.5 * mu0
               * std::norm(x[static_cast<std::size_t>(block_offset(TM::BLK_HY, N2) + k)]);

            // Oscillator energy: ½ |Vz|²/(ε∞ωp²)  (metal only)
            if (g.metal_mask[ks]) {
                const Real wp = g.omega_p_mask[ks];
                if (wp > 0.0)
                    W += 0.5
                       * std::norm(x[static_cast<std::size_t>(block_offset(TM::BLK_VZ, N2) + k)])
                       / (eps_k * wp * wp);
            }
        }
    } else {
        for (Int k = 0; k < N2; ++k) {
            const std::size_t ks = static_cast<std::size_t>(k);
            const Real eps_x = g.eps_x[ks];
            const Real eps_y = g.eps_y[ks];
            const Real mu0   = PhysConst::mu0;

            W += 0.5 * mu0
               * std::norm(x[static_cast<std::size_t>(block_offset(TE::BLK_HZ, N2) + k)]);
            W += 0.5 * eps_x
               * std::norm(x[static_cast<std::size_t>(block_offset(TE::BLK_EX, N2) + k)]);
            W += 0.5 * eps_y
               * std::norm(x[static_cast<std::size_t>(block_offset(TE::BLK_EY, N2) + k)]);

            if (g.metal_mask[ks]) {
                const Real wp = g.omega_p_mask[ks];
                if (wp > 0.0) {
                    W += 0.5
                       * std::norm(x[static_cast<std::size_t>(block_offset(TE::BLK_VX, N2) + k)])
                       / (eps_x * wp * wp);
                    W += 0.5
                       * std::norm(x[static_cast<std::size_t>(block_offset(TE::BLK_VY, N2) + k)])
                       / (eps_y * wp * wp);
                }
            }
        }
    }
    return W;
}

// =============================================================================
// § 3. Modal energy decomposition
// =============================================================================

ModeEnergy compute_mode_energy(const CVec& x, const YeeGrid& g) {
    ModeEnergy e{};
    const Int N2 = g.n_dof;

    if (g.params.mode == Polarization::TM) {
        for (Int k = 0; k < N2; ++k) {
            const std::size_t ks = static_cast<std::size_t>(k);
            const Real eps_k = g.eps_z[ks];
            const Real mu0   = PhysConst::mu0;

            e.E_field += 0.5 * eps_k
                * std::norm(x[static_cast<std::size_t>(block_offset(TM::BLK_EZ, N2) + k)]);
            e.H_field += 0.5 * mu0
                * (std::norm(x[static_cast<std::size_t>(block_offset(TM::BLK_HX, N2) + k)])
                +  std::norm(x[static_cast<std::size_t>(block_offset(TM::BLK_HY, N2) + k)]));

            if (g.metal_mask[ks]) {
                const Real wp = g.omega_p_mask[ks];
                if (wp > 0.0)
                    e.V_field += 0.5
                        * std::norm(x[static_cast<std::size_t>(block_offset(TM::BLK_VZ, N2) + k)])
                        / (eps_k * wp * wp);
            }
        }
    } else {
        for (Int k = 0; k < N2; ++k) {
            const std::size_t ks = static_cast<std::size_t>(k);
            e.H_field += 0.5 * PhysConst::mu0
                * std::norm(x[static_cast<std::size_t>(block_offset(TE::BLK_HZ, N2) + k)]);
            e.E_field += 0.5 * g.eps_x[ks]
                * std::norm(x[static_cast<std::size_t>(block_offset(TE::BLK_EX, N2) + k)]);
            e.E_field += 0.5 * g.eps_y[ks]
                * std::norm(x[static_cast<std::size_t>(block_offset(TE::BLK_EY, N2) + k)]);

            if (g.metal_mask[ks]) {
                const Real wp = g.omega_p_mask[ks];
                if (wp > 0.0) {
                    e.V_field += 0.5
                        * std::norm(x[static_cast<std::size_t>(block_offset(TE::BLK_VX, N2) + k)])
                        / (g.eps_x[ks] * wp * wp);
                    e.V_field += 0.5
                        * std::norm(x[static_cast<std::size_t>(block_offset(TE::BLK_VY, N2) + k)])
                        / (g.eps_y[ks] * wp * wp);
                }
            }
        }
    }
    e.total = e.E_field + e.H_field + e.V_field;
    return e;
}

// =============================================================================
// § 4. Perturbative loss correction  (eq. 15)
// =============================================================================

PerturbResult compute_perturbative_loss(Real        omega0,
                                         const CVec& x,
                                         const YeeGrid& g,
                                         Real        gamma) {
    const Real num = integral_V_metal(x, g);
    const Real den = integral_W0(x, g);

    check(den > 0.0, "compute_perturbative_loss: zero energy density.");

    // ω₁ = (iΓ/2) × (∫|V|²/(ε∞ωp²)) / (∫W₀)
    const Complex omega1 = Complex{0.0, gamma / 2.0} * (num / den);

    PerturbResult r;
    r.omega0  = omega0;
    r.omega1  = omega1;
    r.Q_factor = (omega1.imag() < 0.0)
                ? -omega0 / (2.0 * omega1.imag())
                : std::numeric_limits<Real>::infinity();
    r.energy  = compute_mode_energy(x, g);
    return r;
}

std::vector<PerturbResult>
compute_all_perturbative(const EigenResult& lossless,
                          const YeeGrid&     g,
                          Real               gamma) {
    std::vector<PerturbResult> results;
    results.reserve(static_cast<std::size_t>(lossless.nconv));
    for (int j = 0; j < lossless.nconv; ++j) {
        results.push_back(
            compute_perturbative_loss(
                lossless.eigenvalues[static_cast<std::size_t>(j)],
                lossless.eigenvectors[static_cast<std::size_t>(j)],
                g, gamma));
    }
    return results;
}

// =============================================================================
// § 5. Exact vs. perturbative comparison
// =============================================================================

std::vector<LossComparison>
compare_loss_methods(const EigenResult&                lossless,
                     const ComplexEigenResult&          exact_lossy,
                     const std::vector<PerturbResult>& perturbative) {
    const int n = std::min({lossless.nconv,
                             exact_lossy.nconv,
                             static_cast<int>(perturbative.size())});
    std::vector<LossComparison> cmp;
    cmp.reserve(static_cast<std::size_t>(n));

    for (int j = 0; j < n; ++j) {
        LossComparison c;
        c.omega_real       = lossless.eigenvalues[static_cast<std::size_t>(j)];
        c.im_exact         = exact_lossy.eigenvalues[static_cast<std::size_t>(j)].imag();
        c.im_perturbative  = perturbative[static_cast<std::size_t>(j)].omega1.imag();
        const Real denom   = std::abs(c.im_exact) + 1e-30;
        c.relative_error   = std::abs(c.im_exact - c.im_perturbative) / denom;
        cmp.push_back(c);
    }
    return cmp;
}

// =============================================================================
// § 6. Modal orthogonality verification (A-weighted, eq. 13)
// =============================================================================

OrthogonalityResult
verify_transformed_orthogonality(const EigenResult& result) {
    OrthogonalityResult r{};
    r.n_modes = result.nconv;

    for (int m = 0; m < result.nconv; ++m) {
        for (int p = 0; p < result.nconv; ++p) {
            const Complex G = dot_herm(result.eigenvectors[static_cast<std::size_t>(m)],
                                       result.eigenvectors[static_cast<std::size_t>(p)]);
            if (m == p) {
                r.max_diagonal_error = std::max(r.max_diagonal_error,
                    std::abs(G.real() - 1.0) + std::abs(G.imag()));
            } else {
                r.max_off_diagonal = std::max(r.max_off_diagonal, std::abs(G));
            }
        }
    }
    const Real thr = 1e-6;
    r.passed = (r.max_off_diagonal < thr) && (r.max_diagonal_error < thr);
    return r;
}

OrthogonalityResult
verify_physical_orthogonality(const EigenResult& result,
                               const RVec&        A_diag) {
    // Compute <xₘ|A|xₙ> = Σ_k A[k] conj(xₘ[k]) xₙ[k]
    OrthogonalityResult r{};
    r.n_modes = result.nconv;

    for (int m = 0; m < result.nconv; ++m) {
        for (int p = 0; p < result.nconv; ++p) {
            Complex G{0.0, 0.0};
            const CVec& xm = result.eigenvectors[static_cast<std::size_t>(m)];
            const CVec& xp = result.eigenvectors[static_cast<std::size_t>(p)];
            for (std::size_t k = 0; k < xm.size(); ++k)
                G += A_diag[k] * std::conj(xm[k]) * xp[k];

            if (m == p) {
                r.max_diagonal_error = std::max(r.max_diagonal_error,
                    std::abs(G.real() - 1.0) + std::abs(G.imag()));
            } else {
                r.max_off_diagonal = std::max(r.max_off_diagonal, std::abs(G));
            }
        }
    }
    const Real thr = 1e-6;
    r.passed = (r.max_off_diagonal < thr) && (r.max_diagonal_error < thr);
    return r;
}

void print_orthogonality_report(const OrthogonalityResult& r, int verbosity) {
    if (verbosity < 1) return;
    std::cout << std::scientific << std::setprecision(3)
        << "  Orthogonality (eq. 13):  "
        << "max|G_mn| (m≠n) = " << r.max_off_diagonal
        << "  |G_nn−1| = "       << r.max_diagonal_error
        << "  [" << (r.passed ? "PASS" : "FAIL") << "]\n";
}
