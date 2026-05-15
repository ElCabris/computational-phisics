#pragma once
// =============================================================================
// output.hpp — File output: band data, loss comparison, field profiles
//
// Raman & Fan, PRL 104, 087401 (2010).
//
// Output formats are designed for direct use with gnuplot or matplotlib.
// All files are plain ASCII; header lines begin with '#'.
// =============================================================================

#ifndef HERMETIC_BANDS_OUTPUT_HPP
#define HERMETIC_BANDS_OUTPUT_HPP

#include "utils.hpp"
#include "yee_grid.hpp"
#include "eigensolver.hpp"
#include "perturbation.hpp"
#include <string>
#include <vector>

// =============================================================================
// § 1. Band structure output  →  bands.dat
//
// Format (one line per k-point):
//   k_index   k_norm   ω₁   ω₂   …   ωₙ
//
// k_norm = kx / (π/a) ∈ [0, 1]  along Γ→X
// Frequencies in units of ωp.
// Columns separated by spaces; header explains each column.
// =============================================================================

struct BandPoint {
    int  k_index;           // integer index along the k-path
    Real k_norm;            // kx·a/π ∈ [0,1]
    RVec frequencies;       // converged eigenfrequencies (ascending)
};

void write_bands(const std::string&           filename,
                 const std::vector<BandPoint>& bands,
                 const SimParams&              params);

// =============================================================================
// § 2. Loss comparison output  →  loss_comparison.dat
//
// Format:
//   band_index   omega_real   im_exact   im_perturbative   relative_error
//
// omega_real    : Re[ω] (lossless band frequency, in units of ωp)
// im_exact      : Im[ω] from full non-Hermitian solve (negative, decay)
// im_perturbative : Im[ω] from eq. (15) (negative)
// relative_error  : |im_exact − im_perturb| / |im_exact|
// =============================================================================

void write_loss_comparison(const std::string&                filename,
                            const std::vector<LossComparison>& data,
                            const SimParams&                   params);

// =============================================================================
// § 3. Spatial field profile  →  field_n_k.dat
//
// Written when --field-mode n k is specified.
//
// TM format (one line per Yee node):
//   x   y   |Ez|²   |Hx|²   |Hy|²   |Vz|²
//
// TE format:
//   x   y   |Hz|²   |Ex|²   |Ey|²   |Vx|²   |Vy|²
//
// x, y : physical coordinates in units of a ∈ [0, 1)
// |F|² : field intensity (square magnitude) at the node's Yee position
//
// Nodes ordered row-major: inner loop over ix, outer loop over iy.
// Two blank lines between iy blocks (gnuplot splot format).
// =============================================================================

struct FieldProfile {
    int     mode_index;     // band index n (0-based)
    int     k_index;        // k-point index
    Real    kx_norm;        // kx·a/π
    Real    omega;          // eigenfrequency ω/ωp
    CVec    eigenvector;    // physical x = S⁻¹ ŷ (back-transformed)
};

void write_field_profile(const std::string&  filename,
                          const FieldProfile& fp,
                          const YeeGrid&      g);

// =============================================================================
// § 4. Energy density profile  →  energy_n_k.dat  (optional diagnostic)
//
// Format:
//   x   y   W_E   W_H   W_V   W_total
//
// W_E(r) = ½ ε∞(r)|E₀(r)|²
// W_H(r) = ½ μ₀|H₀(r)|²
// W_V(r) = ½ |V₀(r)|²/(ε∞(r) ωp²(r))
// =============================================================================

void write_energy_profile(const std::string&  filename,
                           const FieldProfile& fp,
                           const YeeGrid&      g);

// =============================================================================
// § 5. Convergence study output  →  convergence.dat
//
// Format:
//   N   band_1   band_2   …   band_nbands   cpu_s
// =============================================================================

void write_convergence(const std::string& filename,
                        const std::vector<ConvergenceEntry>& data,
                        int nbands);

// =============================================================================
// § 6. Console progress reporting
// =============================================================================

// Print one-line summary after solving a k-point.
void report_kpoint(int k_index, Real kx_norm, const EigenResult& result,
                   int verbosity);

// Print final summary after the full k-path is solved.
void report_summary(const std::vector<BandPoint>& bands,
                    double elapsed_seconds,
                    int verbosity);

#endif // HERMETIC_BANDS_OUTPUT_HPP
