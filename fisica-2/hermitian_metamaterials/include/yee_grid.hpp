#pragma once
// =============================================================================
// yee_grid.hpp — 2D Yee staggered-grid geometry and material mask
//
// Raman & Fan, PRL 104, 087401 (2010).
//
// ─────────────────────────────────────────────────────────────────────────────
// TIME CONVENTION
// ─────────────────────────────────────────────────────────────────────────────
// All time-harmonic fields are written as  f(r,t) = f(r) e^{-iωt}.
// Consequently:
//   ∂/∂t  →  -iω
//   ∇×E = -μ₀ ∂H/∂t  →  ∇×E = +iω μ₀ H    (Faraday)
//   ∇×H =  ε  ∂E/∂t  →  ∇×H = -iω ε  E    (Ampère)
//
// The eigenvalue problem ωAx = Bx (Raman & Fan, eq. 6) produces real
// positive ω for propagating lossless modes, consistent with this convention.
//
// ─────────────────────────────────────────────────────────────────────────────
// YEE GRID GEOMETRY
// ─────────────────────────────────────────────────────────────────────────────
// The 2D unit cell [0,a) × [0,a) is divided into N×N uniform cells.
// Cell side: Δ = a/N.  All lengths are normalized to a, so Δ = 1/N.
// Integer indices:  i ∈ {0,…,N-1} (x-direction),  j ∈ {0,…,N-1} (y-direction).
//
// FIELD POSITIONS — TM mode  (active: Ez, Hx, Hy)
// ┌─────────────────────────────────────────────────────────────────────────┐
// │  Field   │  Position                  │  Physical meaning              │
// │──────────┼────────────────────────────┼────────────────────────────────│
// │  Ez[i,j] │  (  i·Δ,    j·Δ  )        │  E-field on primary lattice    │
// │  Hx[i,j] │  (  i·Δ,   (j+½)·Δ )      │  H-field on y-face (↕ edge)   │
// │  Hy[i,j] │  ( (i+½)·Δ,  j·Δ  )       │  H-field on x-face (↔ edge)   │
// └─────────────────────────────────────────────────────────────────────────┘
//
// FIELD POSITIONS — TE mode  (active: Hz, Ex, Ey)
// ┌─────────────────────────────────────────────────────────────────────────┐
// │  Field   │  Position                  │  Physical meaning              │
// │──────────┼────────────────────────────┼────────────────────────────────│
// │  Hz[i,j] │  ( (i+½)·Δ, (j+½)·Δ )     │  H-field at cell center        │
// │  Ex[i,j] │  ( (i+½)·Δ,   j·Δ   )     │  E-field on bottom x-edge      │
// │  Ey[i,j] │  (   i·Δ,  (j+½)·Δ  )     │  E-field on left y-edge        │
// └─────────────────────────────────────────────────────────────────────────┘
//
// ─────────────────────────────────────────────────────────────────────────────
// CURL OPERATOR ORIENTATION
// ─────────────────────────────────────────────────────────────────────────────
// Finite-difference approximations of spatial derivatives use FORWARD
// differences for curl-E (Faraday) and BACKWARD differences for curl-H
// (Ampère).  This is the standard Yee choice that keeps curl-E and curl-H
// as adjoint transposes of each other.
//
//  TM — curl-E  →  H  (Faraday, forward differences):
//
//    [curl_E]_x at Hx[i,j] = (i·Δ, (j+½)·Δ):
//        ∂Ez/∂y  ≈  [Ez[i, j+1] − Ez[i, j]] / Δ
//        wrap: if j+1 == N  →  Ez[i, 0] × e^{+i ky}   (Bloch phase)
//
//    [curl_E]_y at Hy[i,j] = ((i+½)·Δ, j·Δ):
//        ∂Ez/∂x  ≈  [Ez[i+1, j] − Ez[i, j]] / Δ
//        wrap: if i+1 == N  →  Ez[0, j] × e^{+i kx}   (Bloch phase)
//
//  TM — curl-H  →  E  (Ampère, backward differences):
//
//    [curl_H]_z at Ez[i,j] = (i·Δ, j·Δ):
//        ∂Hy/∂x − ∂Hx/∂y
//            ≈  [Hy[i, j] − Hy[i-1, j]] / Δ
//             − [Hx[i, j] − Hx[i, j-1]] / Δ
//        wrap: if i-1 == -1  →  Hy[N-1, j] × e^{-i kx}  (Bloch phase)
//              if j-1 == -1  →  Hx[i, N-1] × e^{-i ky}  (Bloch phase)
//
//  TE — curl-E  →  H  (Faraday, backward differences):
//
//    [curl_E]_z at Hz[i,j] = ((i+½)·Δ, (j+½)·Δ):
//        ∂Ey/∂x − ∂Ex/∂y
//            ≈  [Ey[i+1, j] − Ey[i, j]] / Δ
//             − [Ex[i, j+1] − Ex[i, j]] / Δ
//        wrap: if i+1 == N  →  Ey[0,  j] × e^{+i kx}
//              if j+1 == N  →  Ex[i,  0] × e^{+i ky}
//
//  TE — curl-H  →  E  (Ampère, backward differences):
//
//    [curl_H]_x at Ex[i,j] = ((i+½)·Δ, j·Δ):
//        ∂Hz/∂y  ≈  [Hz[i, j] − Hz[i, j-1]] / Δ
//        wrap: if j-1 == -1  →  Hz[i, N-1] × e^{-i ky}
//
//    [curl_H]_y at Ey[i,j] = (i·Δ, (j+½)·Δ):
//        ∂Hz/∂x  ≈  [Hz[i, j] − Hz[i-1, j]] / Δ
//        wrap: if i-1 == -1  →  Hz[N-1, j] × e^{-i kx}
//
// Summary of Bloch phases:
//   Forward step i → i+1 crossing N: multiply by  e^{+i kx}
//   Forward step j → j+1 crossing N: multiply by  e^{+i ky}
//   Backward step i → i-1 crossing 0: multiply by e^{-i kx}
//   Backward step j → j-1 crossing 0: multiply by e^{-i ky}
//
// ─────────────────────────────────────────────────────────────────────────────
// MATERIAL AVERAGING
// ─────────────────────────────────────────────────────────────────────────────
// For TM, the E-field Ez lives on primary lattice nodes (i·Δ, j·Δ).
// The material ε is sampled at those same points.
//
// For TE, Ex lives at ((i+½)·Δ, j·Δ) and Ey at (i·Δ, (j+½)·Δ).
// These are on edges, not cell centers.  The material property there is
// approximated by harmonic averaging of the two adjacent cells:
//
//   1/ε_Ex[i,j] = ½ [ 1/ε_cell[i,j] + 1/ε_cell[i,j-1] ]   (avg over y)
//   1/ε_Ey[i,j] = ½ [ 1/ε_cell[i,j] + 1/ε_cell[i-1,j] ]   (avg over x)
//
// Harmonic averaging is used (not arithmetic) because ε enters via 1/ε in
// the Ampère curl-H equation, and harmonic averaging of ε corresponds to
// arithmetic averaging of 1/ε — the physically correct interpolation for
// interface-normal components (Zhao & Cangellaris 1996).
// =============================================================================

#ifndef HERMETIC_BANDS_YEE_GRID_HPP
#define HERMETIC_BANDS_YEE_GRID_HPP

#include "utils.hpp"
#include <vector>
#include <string>

// =============================================================================
// § 1. YeeGrid — core grid descriptor
//
// A YeeGrid holds all geometry and material information for one unit cell.
// It is constructed once per simulation and passed const-by-reference to
// every operator-building function.
// =============================================================================

struct YeeGrid {
    // Grid resolution: N×N cells, each of side Δ = 1/N (in units of a).
    int  N;
    Real delta;          // Δ = 1/N

    // Total degrees of freedom per field component.
    Int  n_dof;          // = N * N

    // Simulation parameters (mode, fill fraction, Drude parameters, …).
    SimParams params;

    // ── Material masks ────────────────────────────────────────────────────
    // metal_mask[flat(i,j,N)] = true  if the primary-grid point (i·Δ, j·Δ)
    // lies inside the metallic inclusion.  Used for TM Ez and TE Hz/Ex/Ey.
    std::vector<bool> metal_mask;

    // ── Permittivity arrays ───────────────────────────────────────────────
    // For TM: eps_z[flat(i,j,N)] = ε at Ez node (i·Δ, j·Δ).
    //         In air: ε = 1.  In metal: ε = ε_Drude(ω) (frequency-dependent).
    //         The Drude part is handled via the V auxiliary field; at ω=0
    //         (initial construction) eps_z stores only ε∞ in the metal.
    std::vector<Real> eps_z;    // TM: at (i·Δ, j·Δ)

    // For TE: harmonic-averaged permittivities at Ex and Ey positions.
    //   eps_x[flat(i,j,N)] = harmonic avg of eps_cell at (i+½,j) edge.
    //   eps_y[flat(i,j,N)] = harmonic avg of eps_cell at (i,j+½) edge.
    std::vector<Real> eps_x;    // TE: at ((i+½)·Δ, j·Δ)
    std::vector<Real> eps_y;    // TE: at (i·Δ, (j+½)·Δ)

    // eps_inf_mask: ε∞ value at each primary-grid node.
    //   = 1.0 in air, = params.eps_inf in metal.
    std::vector<Real> eps_inf_mask;  // at (i·Δ, j·Δ)

    // ── Auxiliary-field mask ─────────────────────────────────────────────
    // omega_p_mask[k] = ωp at primary-grid node k.
    //   = 0.0 in air, = params.omega_p in metal.
    // Used to build the A^{-1} diagonal entries.
    std::vector<Real> omega_p_mask;

    // ── Constructor ───────────────────────────────────────────────────────
    explicit YeeGrid(const SimParams& p);

    // ── Query helpers ─────────────────────────────────────────────────────

    // Is primary-grid node (i,j) inside the metal?
    bool is_metal(Int i, Int j) const noexcept {
        return metal_mask[static_cast<std::size_t>(flat_wrap(i, j, N))];
    }

    // ε∞ at primary-grid node (i,j).
    Real eps_inf_at(Int i, Int j) const noexcept {
        return eps_inf_mask[static_cast<std::size_t>(flat_wrap(i, j, N))];
    }

    // ωp at primary-grid node (i,j).  Zero in air.
    Real omegap_at(Int i, Int j) const noexcept {
        return omega_p_mask[static_cast<std::size_t>(flat_wrap(i, j, N))];
    }
};

// =============================================================================
// § 2. Construction helpers (implemented in yee_grid.cpp)
// =============================================================================

// Build the metal_mask, eps_inf_mask, omega_p_mask arrays.
void build_material_masks(YeeGrid& g);

// Build eps_z for TM:
//   eps_z[k] = eps_inf_mask[k]   (static ε∞; Drude correction enters via V).
void build_eps_tm(YeeGrid& g);

// Build eps_x, eps_y for TE using harmonic averaging.
void build_eps_te(YeeGrid& g);

// =============================================================================
// § 3. Field-position coordinate queries
//
// These return the physical (x,y) coordinate of a given field component,
// useful for output and validation.
// =============================================================================

// TM field positions (in units of a = 1).
inline Real Ez_x(Int i, const YeeGrid& g) noexcept {
    return static_cast<Real>(i) * g.delta;
}
inline Real Ez_y(Int j, const YeeGrid& g) noexcept {
    return static_cast<Real>(j) * g.delta;
}

inline Real Hx_x(Int i, const YeeGrid& g) noexcept {
    return static_cast<Real>(i) * g.delta;
}
inline Real Hx_y(Int j, const YeeGrid& g) noexcept {
    return (static_cast<Real>(j) + 0.5) * g.delta;
}

inline Real Hy_x(Int i, const YeeGrid& g) noexcept {
    return (static_cast<Real>(i) + 0.5) * g.delta;
}
inline Real Hy_y(Int j, const YeeGrid& g) noexcept {
    return static_cast<Real>(j) * g.delta;
}

// TE field positions.
inline Real Hz_x(Int i, const YeeGrid& g) noexcept {
    return (static_cast<Real>(i) + 0.5) * g.delta;
}
inline Real Hz_y(Int j, const YeeGrid& g) noexcept {
    return (static_cast<Real>(j) + 0.5) * g.delta;
}

inline Real Ex_x(Int i, const YeeGrid& g) noexcept {
    return (static_cast<Real>(i) + 0.5) * g.delta;
}
inline Real Ex_y(Int j, const YeeGrid& g) noexcept {
    return static_cast<Real>(j) * g.delta;
}

inline Real Ey_x(Int i, const YeeGrid& g) noexcept {
    return static_cast<Real>(i) * g.delta;
}
inline Real Ey_y(Int j, const YeeGrid& g) noexcept {
    return (static_cast<Real>(j) + 0.5) * g.delta;
}

// =============================================================================
// § 4. Harmonic permittivity averaging for TE edges
//
// See header comment on material averaging above.
// These functions compute the harmonic average of ε∞ at the two cells
// adjacent to an edge position.  They are used in build_eps_te().
//
// IMPORTANT: harmonic average of ε is used (NOT arithmetic), because the
// E-field component normal to a material interface must satisfy continuity
// of D = εE, which translates to harmonic averaging of ε along the
// interpolation path.  Using arithmetic averaging introduces O(Δ) errors
// at interfaces that do not vanish with mesh refinement.
// =============================================================================

// ε∞ at Ex position ((i+½)·Δ, j·Δ): harmonic avg over cells (i,j) and (i,j-1).
Real eps_inf_Ex(Int i, Int j, const YeeGrid& g) noexcept;

// ε∞ at Ey position (i·Δ, (j+½)·Δ): harmonic avg over cells (i,j) and (i-1,j).
Real eps_inf_Ey(Int i, Int j, const YeeGrid& g) noexcept;

// =============================================================================
// § 5. Grid summary printout
// =============================================================================

void print_grid_summary(const YeeGrid& g);

#endif // HERMETIC_BANDS_YEE_GRID_HPP
