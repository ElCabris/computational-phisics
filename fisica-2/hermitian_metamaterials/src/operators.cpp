// =============================================================================
// operators.cpp — Discrete curl operators, A matrix, B operator, and Ĥ
//
// Raman & Fan, PRL 104, 087401 (2010).
//
// Yee grid conventions (all documented in yee_grid.hpp):
//   Time convention  : e^{-iωt}
//   curl-E           : forward differences (→ produces H from E via Faraday)
//   curl-H           : backward differences (→ produces E from H via Ampère)
//   Bloch BC         : forward wrap multiplies by e^{+ik}, backward by e^{-ik}
// =============================================================================

#include "operators.hpp"
#include <cmath>
#include <stdexcept>

// =============================================================================
// § 1. Elemental finite-difference operators
//
// Df_x:  (Df_x)_{row,col} where row = flat(i,j), col = flat(i',j')
//   col = flat(i+1,j,N) with value +1/Δ (Bloch phase if i+1 wraps)
//   col = flat(i,  j,N) with value -1/Δ
//
// All N² × N² operators on the primary grid.
// =============================================================================

CSRMatrix build_Df_x(const YeeGrid& g, Real kx) {
    const Int N  = static_cast<Int>(g.N);
    const Real dxinv = 1.0 / g.delta;
    const Complex phase_p = bloch_phase_x(kx);   // e^{+ikx}

    COOBuffer buf;
    coo_reserve(buf, N * N, 2);

    for (Int j = 0; j < N; ++j) {
        for (Int i = 0; i < N; ++i) {
            const Int row = flat(i, j, N);
            // Diagonal: −1/Δ
            coo_push(buf, row, row, Complex{-dxinv, 0.0});
            // Off-diagonal: +1/Δ at (i+1, j) with Bloch phase if wrap
            if (i + 1 < N) {
                coo_push(buf, row, flat(i + 1, j, N), Complex{dxinv, 0.0});
            } else {
                // i+1 wraps to 0: multiply by Bloch phase e^{+ikx}
                coo_push(buf, row, flat(0, j, N), phase_p * dxinv);
            }
        }
    }
    return coo_to_csr(buf, N * N, N * N);
}

CSRMatrix build_Df_y(const YeeGrid& g, Real ky) {
    const Int N = static_cast<Int>(g.N);
    const Real dyinv = 1.0 / g.delta;
    const Complex phase_p = bloch_phase_y(ky);   // e^{+iky}

    COOBuffer buf;
    coo_reserve(buf, N * N, 2);

    for (Int j = 0; j < N; ++j) {
        for (Int i = 0; i < N; ++i) {
            const Int row = flat(i, j, N);
            coo_push(buf, row, row, Complex{-dyinv, 0.0});
            if (j + 1 < N) {
                coo_push(buf, row, flat(i, j + 1, N), Complex{dyinv, 0.0});
            } else {
                // j+1 wraps to 0: multiply by e^{+iky}
                coo_push(buf, row, flat(i, 0, N), phase_p * dyinv);
            }
        }
    }
    return coo_to_csr(buf, N * N, N * N);
}

// Backward differences: Db_x = −Df_x†  (adjoint with conjugated Bloch phase)
//   (Db_x)_{row,col}:
//     col = flat(i,  j,N) with value +1/Δ
//     col = flat(i-1,j,N) with value −1/Δ (× e^{-ikx} if i-1 wraps)

CSRMatrix build_Db_x(const YeeGrid& g, Real kx) {
    const Int N = static_cast<Int>(g.N);
    const Real dxinv = 1.0 / g.delta;
    const Complex phase_m = std::conj(bloch_phase_x(kx));  // e^{-ikx}

    COOBuffer buf;
    coo_reserve(buf, N * N, 2);

    for (Int j = 0; j < N; ++j) {
        for (Int i = 0; i < N; ++i) {
            const Int row = flat(i, j, N);
            coo_push(buf, row, row, Complex{dxinv, 0.0});
            if (i - 1 >= 0) {
                coo_push(buf, row, flat(i - 1, j, N), Complex{-dxinv, 0.0});
            } else {
                // i-1 wraps to N-1: multiply by e^{-ikx}
                coo_push(buf, row, flat(N - 1, j, N), phase_m * (-dxinv));
            }
        }
    }
    return coo_to_csr(buf, N * N, N * N);
}

CSRMatrix build_Db_y(const YeeGrid& g, Real ky) {
    const Int N = static_cast<Int>(g.N);
    const Real dyinv = 1.0 / g.delta;
    const Complex phase_m = std::conj(bloch_phase_y(ky));  // e^{-iky}

    COOBuffer buf;
    coo_reserve(buf, N * N, 2);

    for (Int j = 0; j < N; ++j) {
        for (Int i = 0; i < N; ++i) {
            const Int row = flat(i, j, N);
            coo_push(buf, row, row, Complex{dyinv, 0.0});
            if (j - 1 >= 0) {
                coo_push(buf, row, flat(i, j - 1, N), Complex{-dyinv, 0.0});
            } else {
                // j-1 wraps to N-1: multiply by e^{-iky}
                coo_push(buf, row, flat(i, N - 1, N), phase_m * (-dyinv));
            }
        }
    }
    return coo_to_csr(buf, N * N, N * N);
}

// =============================================================================
// § 2. TM curl operators
//
// C_E_TM : (2N²) × (N²)
//   Upper N² rows (Hx block): +Df_y  (∂Ez/∂y  at Hx position)
//   Lower N² rows (Hy block): −Df_x  (−∂Ez/∂x at Hy position)
//     Sign from (∇×E)_y = −∂Ez/∂x  with e^{-iωt} convention.
//
// C_H_TM : (N²) × (2N²)  = C_E_TM†
// =============================================================================

CSRMatrix build_CE_TM(const YeeGrid& g, Real kx, Real ky) {
    const Int N2 = g.n_dof;
    const CSRMatrix Dfy = build_Df_y(g, ky);
    const CSRMatrix Dfx = build_Df_x(g, kx);

    // Assemble (2N²) × (N²) by embedding Dfy and (−Dfx) at row offsets 0 and N².
    COOBuffer buf;
    coo_reserve(buf, 2 * N2, 2);

    // Upper block (Hx): rows [0, N²)  ←  +Df_y
    for (Int r = 0; r < N2; ++r) {
        const Int s = Dfy.row_ptr[static_cast<std::size_t>(r)];
        const Int e = Dfy.row_ptr[static_cast<std::size_t>(r + 1)];
        for (Int k = s; k < e; ++k)
            coo_push(buf, r, Dfy.col_ind[static_cast<std::size_t>(k)],
                     Dfy.val[static_cast<std::size_t>(k)]);
    }
    // Lower block (Hy): rows [N², 2N²)  ←  −Df_x
    for (Int r = 0; r < N2; ++r) {
        const Int s = Dfx.row_ptr[static_cast<std::size_t>(r)];
        const Int e = Dfx.row_ptr[static_cast<std::size_t>(r + 1)];
        for (Int k = s; k < e; ++k)
            coo_push(buf, N2 + r, Dfx.col_ind[static_cast<std::size_t>(k)],
                     -Dfx.val[static_cast<std::size_t>(k)]);   // minus sign
    }
    return coo_to_csr(buf, 2 * N2, N2);
}

// C_H_TM = C_E_TM†  — built from the adjoint, guarantees exact adjoint relation.
CSRMatrix build_CH_TM(const YeeGrid& g, Real kx, Real ky) {
    return conj_transpose(build_CE_TM(g, kx, ky));
}

// =============================================================================
// § 3. TE curl operators
//
// C_E_TE : (N²) × (2N²)
//   Row i (Hz): ∂Ey/∂x − ∂Ex/∂y  at Hz cell-center
//   Col block 0 (Ex):  −Df_y  (−∂Ex/∂y, forward diff in y)
//   Col block 1 (Ey):  +Df_x  (+∂Ey/∂x, forward diff in x)
//
// C_H_TE : (2N²) × (N²) = C_E_TE†
// =============================================================================

CSRMatrix build_CE_TE(const YeeGrid& g, Real kx, Real ky) {
    const Int N2 = g.n_dof;
    const CSRMatrix Dfx = build_Df_x(g, kx);
    const CSRMatrix Dfy = build_Df_y(g, ky);

    // Assemble (N²) × (2N²) with column offset N² for the Ey block.
    COOBuffer buf;
    coo_reserve(buf, N2, 4);

    // Col block 0 (Ex): −Df_y at columns [0, N²)
    for (Int r = 0; r < N2; ++r) {
        const Int s = Dfy.row_ptr[static_cast<std::size_t>(r)];
        const Int e = Dfy.row_ptr[static_cast<std::size_t>(r + 1)];
        for (Int k = s; k < e; ++k)
            coo_push(buf, r, Dfy.col_ind[static_cast<std::size_t>(k)],
                     -Dfy.val[static_cast<std::size_t>(k)]);   // minus sign
    }
    // Col block 1 (Ey): +Df_x at columns [N², 2N²)
    for (Int r = 0; r < N2; ++r) {
        const Int s = Dfx.row_ptr[static_cast<std::size_t>(r)];
        const Int e = Dfx.row_ptr[static_cast<std::size_t>(r + 1)];
        for (Int k = s; k < e; ++k)
            coo_push(buf, r, N2 + Dfx.col_ind[static_cast<std::size_t>(k)],
                     Dfx.val[static_cast<std::size_t>(k)]);
    }
    return coo_to_csr(buf, N2, 2 * N2);
}

CSRMatrix build_CH_TE(const YeeGrid& g, Real kx, Real ky) {
    return conj_transpose(build_CE_TE(g, kx, ky));
}

// =============================================================================
// § 4. Metal-region identity operators
// =============================================================================

CSRMatrix build_I_metal(const YeeGrid& g) {
    const Int N2 = g.n_dof;
    COOBuffer buf;
    buf.reserve(static_cast<std::size_t>(N2));
    for (Int k = 0; k < N2; ++k)
        if (g.metal_mask[static_cast<std::size_t>(k)])
            coo_push(buf, k, k, Complex{1.0, 0.0});
    return coo_to_csr(buf, N2, N2);
}

CSRMatrix build_I_air(const YeeGrid& g) {
    const Int N2 = g.n_dof;
    COOBuffer buf;
    buf.reserve(static_cast<std::size_t>(N2));
    for (Int k = 0; k < N2; ++k)
        if (!g.metal_mask[static_cast<std::size_t>(k)])
            coo_push(buf, k, k, Complex{1.0, 0.0});
    return coo_to_csr(buf, N2, N2);
}

// =============================================================================
// § 5. Mass matrix A diagonal
//
// TM-Drude state vector: x = (Hx, Hy, Ez, Vz)
//   Block 0 (Hx): μ₀
//   Block 1 (Hy): μ₀
//   Block 2 (Ez): ε∞(r)          [position-dependent]
//   Block 3 (Vz): 1/(ε∞(r)ωp²)   [metal only; regularized to 1 in air]
//
// TE-Drude state vector: x = (Hz, Ex, Ey, Vx, Vy)
//   Block 0 (Hz): μ₀
//   Block 1 (Ex): ε∞_x(r)        [harmonic avg at Ex edge]
//   Block 2 (Ey): ε∞_y(r)        [harmonic avg at Ey edge]
//   Block 3 (Vx): 1/(ε∞_x ωp²)   [regularized to 1 in air]
//   Block 4 (Vy): 1/(ε∞_y ωp²)   [regularized to 1 in air]
// =============================================================================

RVec build_A_diagonal_TM(const YeeGrid& g) {
    const Int N2 = g.n_dof;
    RVec A(static_cast<std::size_t>(TM::N_BLOCKS_DRUDE * N2));

    for (Int k = 0; k < N2; ++k) {
        const std::size_t ks = static_cast<std::size_t>(k);
        // Blocks Hx, Hy: μ₀ = 1 in normalized (lattice) units
        A[static_cast<std::size_t>(TM::BLK_HX * N2 + k)] = 1.0;
        A[static_cast<std::size_t>(TM::BLK_HY * N2 + k)] = 1.0;
        // Block Ez: ε∞(r)
        A[static_cast<std::size_t>(TM::BLK_EZ * N2 + k)] = g.eps_z[ks];
        // Block Vz: 1/(ε∞ ωp²), regularized to 1 in air
        const Real eps_i = g.eps_inf_mask[ks];
        const Real wp_i  = g.omega_p_mask[ks];
        if (wp_i > 0.0)
            A[static_cast<std::size_t>(TM::BLK_VZ * N2 + k)] = 1.0 / (eps_i * wp_i * wp_i);
        else
            A[static_cast<std::size_t>(TM::BLK_VZ * N2 + k)] = 1.0;   // regularization
    }
    return A;
}

RVec build_A_diagonal_TE(const YeeGrid& g) {
    const Int N2 = g.n_dof;
    RVec A(static_cast<std::size_t>(TE::N_BLOCKS_DRUDE * N2));

    for (Int k = 0; k < N2; ++k) {
        const std::size_t ks = static_cast<std::size_t>(k);
        A[static_cast<std::size_t>(TE::BLK_HZ * N2 + k)] = 1.0;  // μ₀ = 1 in normalized units
        A[static_cast<std::size_t>(TE::BLK_EX * N2 + k)] = g.eps_x[ks];
        A[static_cast<std::size_t>(TE::BLK_EY * N2 + k)] = g.eps_y[ks];

        // Vx block: edge at Ex position ((i+½)Δ, jΔ)
        // ωp at that edge: use metal_mask at primary node (i,j) as proxy.
        // (More accurate: average ωp from two adjacent cells as done for ε∞.)
        const Real eps_x_k = g.eps_x[ks];
        const Real wp_k    = g.omega_p_mask[ks];   // proxy: primary node ωp
        if (wp_k > 0.0) {
            A[static_cast<std::size_t>(TE::BLK_VX * N2 + k)] = 1.0 / (eps_x_k * wp_k * wp_k);
            A[static_cast<std::size_t>(TE::BLK_VY * N2 + k)] = 1.0 / (g.eps_y[ks] * wp_k * wp_k);
        } else {
            A[static_cast<std::size_t>(TE::BLK_VX * N2 + k)] = 1.0;
            A[static_cast<std::size_t>(TE::BLK_VY * N2 + k)] = 1.0;
        }
    }
    return A;
}

CSRMatrix build_A_TM(const YeeGrid& g) {
    RVec d = build_A_diagonal_TM(g);
    CVec cd(d.size());
    for (std::size_t i = 0; i < d.size(); ++i) cd[i] = Complex{d[i], 0.0};
    return diag_csr(cd);
}

CSRMatrix build_A_TE(const YeeGrid& g) {
    RVec d = build_A_diagonal_TE(g);
    CVec cd(d.size());
    for (std::size_t i = 0; i < d.size(); ++i) cd[i] = Complex{d[i], 0.0};
    return diag_csr(cd);
}

// =============================================================================
// § 6. √A and (√A)⁻¹ vectors
// =============================================================================

RVec build_sqrtA_TM(const YeeGrid& g) {
    RVec A = build_A_diagonal_TM(g);
    for (auto& v : A) v = std::sqrt(v);
    return A;
}

RVec build_sqrtA_TE(const YeeGrid& g) {
    RVec A = build_A_diagonal_TE(g);
    for (auto& v : A) v = std::sqrt(v);
    return A;
}

RVec build_sqrtA_inv_TM(const YeeGrid& g) {
    RVec A = build_A_diagonal_TM(g);
    for (auto& v : A) v = (v > 0.0) ? 1.0 / std::sqrt(v) : 0.0;
    return A;
}

RVec build_sqrtA_inv_TE(const YeeGrid& g) {
    RVec A = build_A_diagonal_TE(g);
    for (auto& v : A) v = (v > 0.0) ? 1.0 / std::sqrt(v) : 0.0;
    return A;
}

void apply_sqrtA_inv(CVec& x, const RVec& s_inv) {
    for (std::size_t i = 0; i < x.size(); ++i)
        x[i] *= s_inv[i];
}

void apply_sqrtA(CVec& x, const RVec& s) {
    for (std::size_t i = 0; i < x.size(); ++i)
        x[i] *= s[i];
}

// =============================================================================
// § 7. Full operator B (Hermitian, TM-Drude)
//
// State vector: x = (Hx, Hy, Ez, Vz)^T  with block size N²
//
// B block structure:
//
//        Hx         Hy        Ez         Vz
//  Hx [   0          0      −i Df_y      0  ]
//  Hy [   0          0      +i Df_x*     0  ]   (* Df_x is preceded by -1 in CE_TM)
//  Ez [ −i Db_y†   +i Db_x†   0         −iI ]
//  Vz [   0          0       +iI          0  ]
//
// The (Ez,Hx) and (Ez,Hy) coupling entries match (Hx,Ez)† and (Hy,Ez)†
// by construction since C_H_TM = C_E_TM†.
//
// The coupling (Ez,Vz) = −iI_metal and (Vz,Ez) = +iI_metal
// are restricted to metallic nodes (V = 0 in air).
// =============================================================================

CSRMatrix build_B_TM(const YeeGrid& g, Real kx, Real ky) {
    const Int N2  = g.n_dof;
    const Int Ntot = TM::N_BLOCKS_DRUDE * N2;

    const CSRMatrix CE = build_CE_TM(g, kx, ky);   // (2N²) × (N²)
    const CSRMatrix CH = build_CH_TM(g, kx, ky);   // (N²)  × (2N²) = CE†
    const CSRMatrix Im = build_I_metal(g);          // (N²)  × (N²)

    COOBuffer buf;
    coo_reserve(buf, Ntot, 5);

    // ── (H, Ez) block: −i CE  ────────────────────────────────────────────
    // CE has 2N² rows; its upper N² rows correspond to Hx, lower N² to Hy.
    for (Int r = 0; r < 2 * N2; ++r) {
        const Int row_B = block_offset(TM::BLK_HX, N2) + r;   // Hx then Hy rows
        const Int col_off = block_offset(TM::BLK_EZ, N2);
        const Int s = CE.row_ptr[static_cast<std::size_t>(r)];
        const Int e = CE.row_ptr[static_cast<std::size_t>(r + 1)];
        for (Int k = s; k < e; ++k)
            coo_push(buf, row_B,
                     col_off + CE.col_ind[static_cast<std::size_t>(k)],
                     Complex{0, -1} * CE.val[static_cast<std::size_t>(k)]);
    }

    // ── (Ez, H) block: +i CH  ────────────────────────────────────────────
    // CH has N² rows; its 2N² columns span Hx and Hy blocks.
    for (Int r = 0; r < N2; ++r) {
        const Int row_B  = block_offset(TM::BLK_EZ, N2) + r;
        const Int s = CH.row_ptr[static_cast<std::size_t>(r)];
        const Int e = CH.row_ptr[static_cast<std::size_t>(r + 1)];
        for (Int k = s; k < e; ++k) {
            // col_ind into the (2N²) CH col space → add Hx block offset 0
            const Int c_local = CH.col_ind[static_cast<std::size_t>(k)];
            // c_local < N² → Hx block; c_local >= N² → Hy block
            const Int col_B = block_offset(TM::BLK_HX, N2) + c_local;
            coo_push(buf, row_B, col_B,
                     Complex{0, 1} * CH.val[static_cast<std::size_t>(k)]);
        }
    }

    // ── (Ez, Vz) block: −iI_metal  ───────────────────────────────────────
    for (Int k = 0; k < N2; ++k) {
        if (g.metal_mask[static_cast<std::size_t>(k)]) {
            const Int rB = block_offset(TM::BLK_EZ, N2) + k;
            const Int cB = block_offset(TM::BLK_VZ, N2) + k;
            coo_push(buf, rB, cB, Complex{0, -1});
        }
    }

    // ── (Vz, Ez) block: +iI_metal  ───────────────────────────────────────
    for (Int k = 0; k < N2; ++k) {
        if (g.metal_mask[static_cast<std::size_t>(k)]) {
            const Int rB = block_offset(TM::BLK_VZ, N2) + k;
            const Int cB = block_offset(TM::BLK_EZ, N2) + k;
            coo_push(buf, rB, cB, Complex{0, 1});
        }
    }

    return coo_to_csr(buf, Ntot, Ntot);
}

// =============================================================================
// § 8. Full operator B (Hermitian, TE-Drude)
//
// State vector: x = (Hz, Ex, Ey, Vx, Vy)^T
//
// B block structure:
//
//        Hz        Ex       Ey       Vx    Vy
//  Hz [   0      −i·CE_x  −i·CE_y    0     0  ]
//  Ex [ −i·CH_x    0        0       −iI    0  ]
//  Ey [ −i·CH_y    0        0        0    −iI ]
//  Vx [   0       +iI       0        0     0  ]
//  Vy [   0        0       +iI       0     0  ]
//
// where CE has two row blocks (Hz row is split into Ex and Ey contributions)
// and CH has two row blocks (Ex and Ey rows).
// =============================================================================

CSRMatrix build_B_TE(const YeeGrid& g, Real kx, Real ky) {
    const Int N2   = g.n_dof;
    const Int Ntot = TE::N_BLOCKS_DRUDE * N2;

    // CE_TE: (N²) × (2N²) — maps (Ex, Ey) to Hz row
    const CSRMatrix CE = build_CE_TE(g, kx, ky);
    // CH_TE: (2N²) × (N²) = CE† — maps Hz to (Ex, Ey) rows
    const CSRMatrix CH = build_CH_TE(g, kx, ky);

    COOBuffer buf;
    coo_reserve(buf, Ntot, 5);

    // ── (Hz, E) block: −i CE  ─────────────────────────────────────────────
    // CE: N² rows × 2N² cols; row r → Hz block row r; cols 0..N²−1 → Ex, N²..2N²-1 → Ey
    for (Int r = 0; r < N2; ++r) {
        const Int rB = block_offset(TE::BLK_HZ, N2) + r;
        const Int s = CE.row_ptr[static_cast<std::size_t>(r)];
        const Int e = CE.row_ptr[static_cast<std::size_t>(r + 1)];
        for (Int k = s; k < e; ++k) {
            const Int c_local = CE.col_ind[static_cast<std::size_t>(k)];
            // c_local < N² → Ex block; c_local >= N² → Ey block
            const Int col_blk_off = (c_local < N2)
                ? block_offset(TE::BLK_EX, N2) + c_local
                : block_offset(TE::BLK_EY, N2) + (c_local - N2);
            coo_push(buf, rB, col_blk_off,
                     Complex{0, -1} * CE.val[static_cast<std::size_t>(k)]);
        }
    }

    // ── (E, Hz) block: −i CH  ─────────────────────────────────────────────
    // CH: 2N² rows × N² cols; rows 0..N²-1 → Ex block, N²..2N²-1 → Ey block
    for (Int r = 0; r < 2 * N2; ++r) {
        const Int blk  = (r < N2) ? TE::BLK_EX : TE::BLK_EY;
        const Int r_in_blk = (r < N2) ? r : r - N2;
        const Int rB = block_offset(blk, N2) + r_in_blk;
        const Int s = CH.row_ptr[static_cast<std::size_t>(r)];
        const Int e = CH.row_ptr[static_cast<std::size_t>(r + 1)];
        for (Int k = s; k < e; ++k) {
            const Int cB = block_offset(TE::BLK_HZ, N2)
                         + CH.col_ind[static_cast<std::size_t>(k)];
            coo_push(buf, rB, cB,
                     Complex{0, -1} * CH.val[static_cast<std::size_t>(k)]);
        }
    }

    // ── (Ex, Vx) and (Ey, Vy) coupling: −iI_metal  ───────────────────────
    for (Int k = 0; k < N2; ++k) {
        if (g.metal_mask[static_cast<std::size_t>(k)]) {
            coo_push(buf, block_offset(TE::BLK_EX, N2) + k,
                          block_offset(TE::BLK_VX, N2) + k, Complex{0, -1});
            coo_push(buf, block_offset(TE::BLK_EY, N2) + k,
                          block_offset(TE::BLK_VY, N2) + k, Complex{0, -1});
        }
    }

    // ── (Vx, Ex) and (Vy, Ey) coupling: +iI_metal  ───────────────────────
    for (Int k = 0; k < N2; ++k) {
        if (g.metal_mask[static_cast<std::size_t>(k)]) {
            coo_push(buf, block_offset(TE::BLK_VX, N2) + k,
                          block_offset(TE::BLK_EX, N2) + k, Complex{0, 1});
            coo_push(buf, block_offset(TE::BLK_VY, N2) + k,
                          block_offset(TE::BLK_EY, N2) + k, Complex{0, 1});
        }
    }

    return coo_to_csr(buf, Ntot, Ntot);
}

// =============================================================================
// § 9. Hermitian Hamiltonian  Ĥ = S⁻¹ B S⁻¹  (eq. 12, Raman & Fan)
//
// For each stored non-zero B[i,j]:   Ĥ[i,j] = s_inv[i] * B[i,j] * s_inv[j]
// This is O(nnz) and preserves the sparsity pattern of B.
// =============================================================================

static CSRMatrix apply_similarity(const CSRMatrix& B, const RVec& s_inv) {
    CSRMatrix H = B;   // copy sparsity pattern
    for (Int r = 0; r < B.nrows; ++r) {
        const Real si = s_inv[static_cast<std::size_t>(r)];
        const Int start = B.row_ptr[static_cast<std::size_t>(r)];
        const Int end   = B.row_ptr[static_cast<std::size_t>(r + 1)];
        for (Int k = start; k < end; ++k) {
            const Int   c  = B.col_ind[static_cast<std::size_t>(k)];
            const Real sj  = s_inv[static_cast<std::size_t>(c)];
            H.val[static_cast<std::size_t>(k)] = si * B.val[static_cast<std::size_t>(k)] * sj;
        }
    }
    return H;
}

CSRMatrix build_H_hat_TM(const YeeGrid& g, Real kx, Real ky) {
    const CSRMatrix B      = build_B_TM(g, kx, ky);
    const RVec      s_inv  = build_sqrtA_inv_TM(g);
    return apply_similarity(B, s_inv);
}

CSRMatrix build_H_hat_TE(const YeeGrid& g, Real kx, Real ky) {
    const CSRMatrix B      = build_B_TE(g, kx, ky);
    const RVec      s_inv  = build_sqrtA_inv_TE(g);
    return apply_similarity(B, s_inv);
}

CSRMatrix build_H_hat(const YeeGrid& g, Real kx, Real ky) {
    if (g.params.mode == Polarization::TM)
        return build_H_hat_TM(g, kx, ky);
    else
        return build_H_hat_TE(g, kx, ky);
}

// =============================================================================
// § 10. Adjoint verification
// =============================================================================

Real verify_curl_adjoint_TM(const YeeGrid& g, Real kx, Real ky) {
    const CSRMatrix CE  = build_CE_TM(g, kx, ky);
    const CSRMatrix CH  = build_CH_TM(g, kx, ky);
    const CSRMatrix CEd = conj_transpose(CE);
    return relative_diff(CH, CEd);
}

Real verify_curl_adjoint_TE(const YeeGrid& g, Real kx, Real ky) {
    const CSRMatrix CE  = build_CE_TE(g, kx, ky);
    const CSRMatrix CH  = build_CH_TE(g, kx, ky);
    const CSRMatrix CEd = conj_transpose(CE);
    return relative_diff(CH, CEd);
}
