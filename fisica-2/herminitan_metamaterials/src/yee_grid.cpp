// =============================================================================
// yee_grid.cpp — YeeGrid construction and material mask assembly
//
// Raman & Fan, PRL 104, 087401 (2010).
//
// Field positions and curl orientations follow the conventions in yee_grid.hpp.
// =============================================================================

#include "yee_grid.hpp"
#include <iostream>
#include <iomanip>
#include <cmath>

// =============================================================================
// YeeGrid constructor
// =============================================================================

YeeGrid::YeeGrid(const SimParams& p)
    : N(p.resolution),
      delta(1.0 / static_cast<Real>(p.resolution)),
      n_dof(static_cast<Int>(p.resolution) * static_cast<Int>(p.resolution)),
      params(p)
{
    check(N >= 4, "YeeGrid: resolution must be >= 4");

    const std::size_t sz = static_cast<std::size_t>(n_dof);
    metal_mask  .assign(sz, false);
    eps_z       .assign(sz, 1.0);
    eps_x       .assign(sz, 1.0);
    eps_y       .assign(sz, 1.0);
    eps_inf_mask.assign(sz, 1.0);
    omega_p_mask.assign(sz, 0.0);

    build_material_masks(*this);

    if (p.mode == Polarization::TM)
        build_eps_tm(*this);
    else
        build_eps_te(*this);
}

// =============================================================================
// Material masks — primary-grid nodes (i·Δ, j·Δ)
// =============================================================================

void build_material_masks(YeeGrid& g) {
    const int N = g.N;
    for (int j = 0; j < N; ++j) {
        for (int i = 0; i < N; ++i) {
            const std::size_t k = static_cast<std::size_t>(flat(i, j, N));
            const bool metal = in_metal(i, j, N, g.params);
            g.metal_mask  [k] = metal;
            g.eps_inf_mask[k] = metal ? g.params.eps_inf : 1.0;
            g.omega_p_mask[k] = metal ? g.params.omega_p : 0.0;
        }
    }
}

// =============================================================================
// TM permittivity: eps_z at primary-grid nodes (i·Δ, j·Δ)
//
// In the Hermitian formulation, the Drude contribution to ε(ω) is carried
// by the auxiliary V field.  Therefore eps_z stores only ε∞(r):
//   eps_z[k] = ε∞(r_k)  =  eps_inf in metal, 1.0 in air.
// The dynamic Drude correction −ε∞ωp²/ω² enters via the A_V diagonal block.
// =============================================================================

void build_eps_tm(YeeGrid& g) {
    const int N = g.N;
    for (int j = 0; j < N; ++j)
        for (int i = 0; i < N; ++i) {
            const std::size_t k = static_cast<std::size_t>(flat(i, j, N));
            g.eps_z[k] = g.eps_inf_mask[k];
        }
}

// =============================================================================
// TE permittivity: harmonic averaging at edge positions
//
// Ex[i,j] sits at ((i+½)·Δ, j·Δ) — between cells (i,j) and (i,j−1).
// Ey[i,j] sits at (i·Δ, (j+½)·Δ) — between cells (i,j) and (i−1,j).
//
// Harmonic average of ε∞ (= arithmetic average of 1/ε∞) is correct for
// the tangential field component at a dielectric interface.
// See Zhao & Cangellaris, IEEE TMTT 44 (1996).
// =============================================================================

Real eps_inf_Ex(Int i, Int j, const YeeGrid& g) noexcept {
    // Cells adjacent to the Ex edge at ((i+½)·Δ, j·Δ):
    //   upper: (i, j)    lower: (i, j−1)
    const Int N = static_cast<Int>(g.N);
    const Real e_upper = g.eps_inf_mask[static_cast<std::size_t>(flat_wrap(i,     j,     N))];
    const Real e_lower = g.eps_inf_mask[static_cast<std::size_t>(flat_wrap(i, j - 1, N))];
    // Harmonic average: 1/ε̄ = ½(1/ε_upper + 1/ε_lower)
    return 2.0 / (1.0 / e_upper + 1.0 / e_lower);
}

Real eps_inf_Ey(Int i, Int j, const YeeGrid& g) noexcept {
    // Cells adjacent to the Ey edge at (i·Δ, (j+½)·Δ):
    //   right: (i, j)    left: (i−1, j)
    const Int N = static_cast<Int>(g.N);
    const Real e_right = g.eps_inf_mask[static_cast<std::size_t>(flat_wrap(i,     j, N))];
    const Real e_left  = g.eps_inf_mask[static_cast<std::size_t>(flat_wrap(i - 1, j, N))];
    return 2.0 / (1.0 / e_right + 1.0 / e_left);
}

void build_eps_te(YeeGrid& g) {
    const Int N = static_cast<Int>(g.N);
    for (Int j = 0; j < N; ++j) {
        for (Int i = 0; i < N; ++i) {
            const std::size_t k = static_cast<std::size_t>(flat(i, j, N));
            g.eps_x[k] = eps_inf_Ex(i, j, g);
            g.eps_y[k] = eps_inf_Ey(i, j, g);
        }
    }
}

// =============================================================================
// Grid summary
// =============================================================================

void print_grid_summary(const YeeGrid& g) {
    if (g.params.verbosity < 1) return;

    // Count metallic nodes.
    int n_metal = 0;
    for (auto m : g.metal_mask) if (m) ++n_metal;
    const double fill_actual =
        static_cast<double>(n_metal) / static_cast<double>(g.n_dof);

    std::cout << std::fixed << std::setprecision(4)
        << "┌─ Yee grid ───────────────────────────────────────────────┐\n"
        << "│  Resolution      : " << g.N << " × " << g.N
        <<   " = " << g.n_dof << " nodes\n"
        << "│  Cell size Δ     : " << g.delta << " (units of a)\n"
        << "│  Metallic nodes  : " << n_metal
        <<   " / " << g.n_dof
        <<   " (fill = " << fill_actual << ")\n"
        << "│  Requested fill  : " << g.params.fill_fraction << "\n"
        << "│  Total DOF       : "
        << ((g.params.mode == Polarization::TM) ? 4 : 5) * g.n_dof
        << "\n└──────────────────────────────────────────────────────────┘\n";
}
