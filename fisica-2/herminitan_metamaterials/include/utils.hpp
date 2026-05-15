#pragma once
// =============================================================================
// utils.hpp — Fundamental types, physical constants, and utility functions
//
// Raman & Fan, PRL 104, 087401 (2010).
//
// Convention throughout the code:
//   - Lengths are in units of the lattice constant  a.
//   - Frequencies are in units of ωp (plasma frequency).
//   - The dispersion relation is ε(ω) = ε∞ (1 - ωp²/ω²)  [Drude, lossless].
//   - All matrices are built in SI-like normalized form; see matrices.hpp.
// =============================================================================

#ifndef HERMETIC_BANDS_UTILS_HPP
#define HERMETIC_BANDS_UTILS_HPP

#include <complex>
#include <cstddef>
#include <cstdint>
#include <cmath>
#include <string>
#include <vector>
#include <array>
#include <stdexcept>
#include <cassert>
#include <iostream>
#include <iomanip>
#include <chrono>

// =============================================================================
// § 1. Scalar type aliases
// =============================================================================

// All floating-point work uses double precision.
using Real    = double;
using Complex = std::complex<double>;

// Integer type for matrix dimension indices (signed, 64-bit).
// Using int32 would overflow for N > 46340 on a square matrix.
using Int     = std::int64_t;
using UInt    = std::uint64_t;

// Short alias for the imaginary unit i = sqrt(-1).
inline constexpr Complex I_UNIT {0.0, 1.0};

// =============================================================================
// § 2. Physical and mathematical constants (SI, normalized)
// =============================================================================

namespace PhysConst {

    // Vacuum permeability  μ₀  [H/m]
    inline constexpr Real mu0     = 1.2566370614359173e-6;

    // Vacuum permittivity  ε₀  [F/m]
    inline constexpr Real eps0    = 8.854187817620389e-12;

    // Speed of light in vacuum  c = 1/√(μ₀ε₀)  [m/s]
    inline constexpr Real c_light = 2.99792458e8;

    // Mathematical π
    inline constexpr Real PI      = 3.141592653589793238462643383279502884L;

    // 2π
    inline constexpr Real TWO_PI  = 2.0 * PI;

} // namespace PhysConst

// =============================================================================
// § 3. Simulation parameters — all parameters the user can set via CLI
// =============================================================================

// Polarization mode: TM has scalar Ez; TE has scalar Hz.
enum class Polarization { TE, TM };

// Shape of the metallic inclusion within the unit cell.
enum class InclusionShape { SQUARE, CIRCLE };

struct SimParams {
    // --- Physical model ---
    Polarization    mode          = Polarization::TE;
    InclusionShape  shape         = InclusionShape::SQUARE;

    // Drude model parameters (normalized to ωp = 1).
    Real            eps_inf       = 1.0;    // ε∞ high-frequency permittivity
    Real            omega_p       = 1.0;    // plasma frequency [normalized]
    Real            gamma         = 0.0;    // damping rate Γ [normalized to ωp]
    Real            omega_0       = 0.0;    // resonance frequency (0 = Drude)

    // Fill fraction: side length of metallic square / lattice constant a.
    // For SQUARE: s = fill_fraction * a  (paper uses s/a = 0.25 as default).
    // For CIRCLE: r = fill_fraction * a.
    Real            fill_fraction = 0.25;

    // --- Numerical parameters ---
    int             resolution    = 20;     // Yee grid points per unit cell side
    int             nk            = 20;     // k-points along Γ→X
    int             nbands        = 15;     // number of eigenfrequencies to compute

    // --- Output control ---
    std::string     output_file   = "bands.dat";
    bool            perturb       = false;  // compute perturbative losses
    bool            verify_ortho  = false;  // verify eq. (13) orthogonality

    // For spatial field export: mode index and k-point index.
    int             field_mode    = -1;     // -1 = don't export
    int             field_kindex  = -1;

    // Verbosity level: 0 = silent, 1 = normal, 2 = verbose/debug.
    int             verbosity     = 1;
};

// =============================================================================
// § 4. Dense complex vector — thin wrapper to avoid std::vector<Complex>
//      boilerplate in inner loops.
// =============================================================================

// We keep this as a plain typedef; operators are free functions below.
using CVec = std::vector<Complex>;
using RVec = std::vector<Real>;

// Allocate a zero-initialized complex vector of length n.
inline CVec make_cvec(Int n) { return CVec(static_cast<std::size_t>(n), 0.0); }
inline RVec make_rvec(Int n) { return RVec(static_cast<std::size_t>(n), 0.0); }

// ---------------------------------------------------------------------------
// § 4.1  BLAS-style dot products (conjugate and non-conjugate)
// ---------------------------------------------------------------------------

// Hermitian inner product  <x|y> = Σ conj(x_i) y_i
inline Complex dot_herm(const CVec& x, const CVec& y) {
    assert(x.size() == y.size());
    Complex acc{0.0, 0.0};
    for (std::size_t i = 0; i < x.size(); ++i)
        acc += std::conj(x[i]) * y[i];
    return acc;
}

// Euclidean norm  ||x||₂
inline Real norm2(const CVec& x) {
    return std::sqrt(std::real(dot_herm(x, x)));
}

// ---------------------------------------------------------------------------
// § 4.2  Scale and axpy   (y += alpha * x)
// ---------------------------------------------------------------------------
inline void axpy(Complex alpha, const CVec& x, CVec& y) {
    assert(x.size() == y.size());
    for (std::size_t i = 0; i < x.size(); ++i)
        y[i] += alpha * x[i];
}

// =============================================================================
// § 5. Yee grid index utilities
//
// The 2D Yee cell has resolution N×N points.  We use row-major (C) ordering:
//   flat_index(ix, iy) = iy * N + ix,   0 ≤ ix,iy < N.
//
// Periodicity is handled by modular wrap:
//   wrap(i, N) = ((i % N) + N) % N
// =============================================================================

// Flat index into an N×N array, row-major.
inline Int flat(Int ix, Int iy, Int N) noexcept {
    return iy * N + ix;
}

// Periodic wrap (handles negative indices correctly).
inline Int wrap(Int i, Int N) noexcept {
    return ((i % N) + N) % N;
}

// Flat index with periodic wrap.
inline Int flat_wrap(Int ix, Int iy, Int N) noexcept {
    return flat(wrap(ix, N), wrap(iy, N), N);
}

// =============================================================================
// § 6. Bloch-periodic phase factor
//
// When taking a finite-difference step that crosses the unit cell boundary,
// the Bloch condition supplies a phase factor:
//
//   f(r + a ê_μ) = e^{i k_μ a} f(r)
//
// In normalized units (a = 1), k_x and k_y are in [0, π/a].
// =============================================================================

// Returns the Bloch phase e^{i k_mu * a} for a step in direction mu ∈ {x,y}.
// kx, ky: wave-vector components in units of 1/a (so k ∈ [0, π]).
inline Complex bloch_phase_x(Real kx) noexcept {
    return std::exp(I_UNIT * kx);
}

inline Complex bloch_phase_y(Real ky) noexcept {
    return std::exp(I_UNIT * ky);
}

// =============================================================================
// § 7. Drude permittivity
//
// ε(ω) = ε∞ (1 − ωp² / (ω² + iΓω))
//
// For Γ = 0 (lossless):  ε(ω) = ε∞ (1 − ωp²/ω²)
// ω, ωp, Γ all in the same normalized units.
// =============================================================================

inline Complex drude_epsilon(Real omega, const SimParams& p) {
    Complex denom = Complex(omega * omega - p.omega_0 * p.omega_0,
                            p.gamma * omega);
    return Complex(p.eps_inf, 0.0)
           * (1.0 - p.omega_p * p.omega_p / denom);
}

// =============================================================================
// § 8. Material mask on the Yee grid
//
// Returns true if the grid point (ix, iy) lies inside the metallic inclusion.
// Coordinates: cell is [0, 1) × [0, 1) normalized to the lattice constant a.
//   SQUARE: |x - 0.5| < s/2 AND |y - 0.5| < s/2   with s = fill_fraction
//   CIRCLE: (x-0.5)² + (y-0.5)² < r²               with r = fill_fraction
// =============================================================================

inline bool in_metal(Int ix, Int iy, Int N, const SimParams& p) noexcept {
    const Real dx = 1.0 / static_cast<Real>(N);
    // Cell-centered coordinates in [0,1)
    const Real cx = (static_cast<Real>(ix) + 0.5) * dx;
    const Real cy = (static_cast<Real>(iy) + 0.5) * dx;
    const Real fx = cx - 0.5;
    const Real fy = cy - 0.5;

    switch (p.shape) {
        case InclusionShape::SQUARE:
            return (std::abs(fx) < 0.5 * p.fill_fraction) &&
                   (std::abs(fy) < 0.5 * p.fill_fraction);
        case InclusionShape::CIRCLE:
            return (fx * fx + fy * fy) < (p.fill_fraction * p.fill_fraction);
    }
    return false; // unreachable
}

// =============================================================================
// § 9. k-path generation   Γ → X   in the first Brillouin zone
//
// For a 2D square lattice with lattice constant a = 1:
//   Γ = (0,   0  )
//   X = (π/a, 0  ) = (π, 0) in normalized units
//
// Returns a vector of nk evenly-spaced k_x values in [0, π].
// k_y = 0 throughout (the paper's Figure 1 uses this path).
// =============================================================================

inline RVec kpath_GammaX(int nk) {
    RVec kx(static_cast<std::size_t>(nk));
    for (int i = 0; i < nk; ++i)
        kx[static_cast<std::size_t>(i)] =
            PhysConst::PI * static_cast<Real>(i) / static_cast<Real>(nk - 1);
    return kx;
}

// =============================================================================
// § 10. Wall-clock timer (useful for benchmarking ARPACK convergence)
// =============================================================================

struct Timer {
    using Clock = std::chrono::high_resolution_clock;
    Clock::time_point t0;

    void start() { t0 = Clock::now(); }

    // Returns elapsed time in seconds.
    double elapsed() const {
        auto dt = Clock::now() - t0;
        return std::chrono::duration<double>(dt).count();
    }
};

// =============================================================================
// § 11. Logging helpers
// =============================================================================

// Log at given verbosity level.
//   level 0: always printed (errors)
//   level 1: normal progress messages
//   level 2: debug/verbose
inline void log(int verbosity_threshold, int verbosity,
                const std::string& msg) {
    if (verbosity >= verbosity_threshold)
        std::cout << msg << "\n";
}

// Print a section banner for readable terminal output.
inline void banner(const std::string& title, int verbosity) {
    if (verbosity >= 1) {
        const int fill = static_cast<int>(60 - title.size() - 4);
        std::cout << "\n-- " << title << " "
                  << std::string(static_cast<std::size_t>(fill > 0 ? fill : 0), '-')
                  << "\n";
    }
}

// =============================================================================
// § 12. Error handling — throw on numerical or logical errors
// =============================================================================

// Throw a std::runtime_error with a formatted message.
[[noreturn]] inline void fatal(const std::string& msg) {
    throw std::runtime_error("[hermetic-bands] FATAL: " + msg);
}

// Assert with message (active in all build types, unlike <cassert>).
inline void check(bool cond, const std::string& msg) {
    if (!cond) fatal(msg);
}

// =============================================================================
// § 12b. CLI argument parsing and parameter validation (implemented in utils.cpp)
// =============================================================================

SimParams parse_args(int argc, char** argv);
void print_usage(const char* prog);
void validate_params(const SimParams& p);
void print_params(const SimParams& p, int verbosity);

// =============================================================================
// § 13. Small utility: pretty-print a vector of eigenfrequencies
// =============================================================================

inline void print_eigenvalues(const RVec& evals, int verbosity) {
    if (verbosity < 1) return;
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "  Eigenfrequencies [ω/ωp]:";
    for (auto v : evals)
        std::cout << "  " << v;
    std::cout << "\n";
}

#endif // HERMETIC_BANDS_UTILS_HPP
