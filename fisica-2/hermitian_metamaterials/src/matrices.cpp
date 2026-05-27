// =============================================================================
// matrices.cpp — High-level problem assembly
//
// Raman & Fan, PRL 104, 087401 (2010).
// =============================================================================

#include "matrices.hpp"
#include <cmath>
#include <stdexcept>
#include <iostream>
#include <algorithm>

// LAPACK declarations (Fortran calling convention, underscored symbols).
extern "C" {
    // LU factorization of a general complex matrix.
    void zgetrf_(int* m, int* n, Complex* A, int* lda,
                 int* ipiv, int* info);
    // Solve using factored LU (single RHS or multiple).
    void zgetrs_(char* trans, int* n, int* nrhs, Complex* A, int* lda,
                 int* ipiv, Complex* B, int* ldb, int* info);
}

// =============================================================================
// § 1. Null-space DOF count
// =============================================================================

Int count_nullspace_dof(const CSRMatrix& H_hat) {
    Int count = 0;
    for (Int i = 0; i < H_hat.nrows; ++i) {
        // A DOF is in the null space if its entire row is zero.
        const Int s = H_hat.row_ptr[static_cast<std::size_t>(i)];
        const Int e = H_hat.row_ptr[static_cast<std::size_t>(i + 1)];
        bool row_zero = true;
        for (Int k = s; k < e; ++k) {
            if (std::abs(H_hat.val[static_cast<std::size_t>(k)]) > 0.0) {
                row_zero = false;
                break;
            }
        }
        if (row_zero) ++count;
    }
    return count;
}

// =============================================================================
// § 2. Gershgorin spectral radius estimate
// =============================================================================

Real gershgorin_radius(const CSRMatrix& H_hat) {
    Real max_bound = 0.0;
    for (Int i = 0; i < H_hat.nrows; ++i) {
        const Int s = H_hat.row_ptr[static_cast<std::size_t>(i)];
        const Int e = H_hat.row_ptr[static_cast<std::size_t>(i + 1)];
        Real diag = 0.0, offsum = 0.0;
        for (Int k = s; k < e; ++k) {
            const Int   c = H_hat.col_ind[static_cast<std::size_t>(k)];
            const Real  v = std::abs(H_hat.val[static_cast<std::size_t>(k)]);
            if (c == i) diag = v; else offsum += v;
        }
        max_bound = std::max(max_bound, diag + offsum);
    }
    return max_bound;
}

// =============================================================================
// § 3. Automatic sigma selection
// =============================================================================

static Real choose_sigma(const CSRMatrix& H_hat, Real omega_p, Real user_sigma) {
    if (user_sigma > 0.0) return user_sigma;
    // Default: 1% of ωp.  Well above ω=0 null space and below physical bands.
    (void)H_hat;
    return 0.01 * omega_p;
}

// =============================================================================
// § 4. Assemble lossless Hermitian system
// =============================================================================

HermitianSystem assemble_hermitian(const YeeGrid& g,
                                    Real kx, Real ky,
                                    Real sigma_shift) {
    HermitianSystem sys;
    sys.N    = g.N;
    sys.kx   = kx;
    sys.ky   = ky;

    if (g.params.mode == Polarization::TM) {
        sys.H_hat    = build_H_hat_TM(g, kx, ky);
        sys.sqrtA    = build_sqrtA_TM(g);
        sys.sqrtA_inv = build_sqrtA_inv_TM(g);
        sys.n_dof    = TM::total_dof(g.n_dof);
    } else {
        sys.H_hat    = build_H_hat_TE(g, kx, ky);
        sys.sqrtA    = build_sqrtA_TE(g);
        sys.sqrtA_inv = build_sqrtA_inv_TE(g);
        sys.n_dof    = TE::total_dof(g.n_dof);
    }

    sys.n_nullspace = count_nullspace_dof(sys.H_hat);
    sys.sigma       = choose_sigma(sys.H_hat, g.params.omega_p, sigma_shift);

    // ── Hermiticity diagnostic ────────────────────────────────────────────────
    // Verifies ||Ĥ - Ĥ†||_max / ||Ĥ||_F ≈ 0.
    // A correct assembly gives values below O(N · ε_machine) ≈ N × 2.2e-16.
    // Active only at verbosity >= 2 to avoid overhead in production runs.
    // NOTE: this is a MATRIX property check (eq. 12 of Raman & Fan), distinct
    //       from modal orthogonality (eq. 13) which requires eigenvectors.
    if (g.params.verbosity >= 2) {
        const Real herm_err = check_hermitian(sys.H_hat);
        const Real frob     = frobenius_norm(sys.H_hat);
        const Real rel_err  = (frob > 0.0) ? herm_err / frob : herm_err;
        std::cout << "  [hermiticity] ||Ĥ - Ĥ†||_max = " << herm_err
                  << "  (rel/Frob = " << rel_err << ")\n";
        if (herm_err > 1e-10 * frob)
            std::cerr << "  [WARNING] Hermiticity error exceeds 1e-10 × ||Ĥ||_F. "
                         "Check Bloch phase signs in operators.\n";
    }

    return sys;
}

// =============================================================================
// § 5. Assemble lossy non-Hermitian system
//
// B_NH = B₀ + δB,   δB_{VV} = −iΓ I_metal
// Ĥ_NH = S⁻¹ B_NH S⁻¹ = Ĥ₀ + δĤ
// δĤ_{VV}[k] = s_inv[V,k]² × (−iΓ) = −iΓ/(ε∞(k) ωp²(k))   in metal
// =============================================================================

NonHermitianSystem assemble_nonhermitian(const YeeGrid& g,
                                          Real kx, Real ky,
                                          Real gamma) {
    // Start from the lossless system.
    HermitianSystem h0 = assemble_hermitian(g, kx, ky);

    NonHermitianSystem sys;
    sys.N        = g.N;
    sys.kx       = kx;
    sys.ky       = ky;
    sys.gamma    = gamma;
    sys.n_dof    = h0.n_dof;
    sys.sqrtA    = h0.sqrtA;
    sys.sqrtA_inv = h0.sqrtA_inv;
    sys.sigma    = h0.sigma;

    // Copy Ĥ₀ and add δĤ_{VV} = −iΓ/(ε∞ ωp²) on metallic V nodes.
    sys.H_hat_nh = h0.H_hat;   // copy

    const Int N2 = g.n_dof;

    // ── add_vv_diag: add δĤ_{VV} = −iΓ·s_inv² to ALL metallic V nodes ────────
    //
    // FIX (was: bug that only added the perturbation to the FIRST metallic node):
    //
    // In the lossless system, B[VV] = 0 for all metallic nodes, therefore
    // Ĥ₀[VV] = S⁻¹·0·S⁻¹ = 0 — no diagonal entry exists in the CSR structure.
    // The old code searched for the diagonal, and when not found (!found), it
    // rebuilt the COO with ONE new entry and then `break`-ed out of the loop,
    // silently skipping all remaining metallic nodes.
    //
    // Correct approach: collect δ for ALL metallic nodes first, then do a
    // single COO rebuild that inserts all of them at once.  coo_to_csr() sums
    // duplicate (row,col) pairs, so this is safe even if a diagonal already
    // exists (future-proof).
    auto add_vv_diag = [&](int blk) {
        // Pass 1: collect (global_row, delta) for every metallic V node.
        COOBuffer new_diag;
        for (Int k = 0; k < N2; ++k) {
            if (!g.metal_mask[static_cast<std::size_t>(k)]) continue;
            const Int global_k = block_offset(blk, N2) + k;
            const Real si = h0.sqrtA_inv[static_cast<std::size_t>(global_k)];
            // δĤ_{kk} = s_inv[k]² × (−iΓ) = −iΓ/(ε∞(r_k) ωp²(r_k))
            const Complex delta = Complex{0.0, -gamma} * (si * si);
            new_diag.push_back({global_k, global_k, delta});
        }
        if (new_diag.empty()) return;   // no metallic nodes in this block

        // Pass 2: single COO rebuild — copy all existing non-zeros, then
        //         append all V-V delta entries.  One coo_to_csr call handles
        //         both insertion (structural zeros) and accumulation (if entry
        //         already existed, which cannot happen in the lossless H_hat₀
        //         but is handled correctly by the duplicate-sum rule).
        COOBuffer buf;
        coo_reserve(buf, sys.H_hat_nh.nrows, 5);
        for (Int r = 0; r < sys.H_hat_nh.nrows; ++r) {
            const Int s2 = sys.H_hat_nh.row_ptr[static_cast<std::size_t>(r)];
            const Int e2 = sys.H_hat_nh.row_ptr[static_cast<std::size_t>(r + 1)];
            for (Int ki = s2; ki < e2; ++ki)
                buf.push_back({r,
                    sys.H_hat_nh.col_ind[static_cast<std::size_t>(ki)],
                    sys.H_hat_nh.val[static_cast<std::size_t>(ki)]});
        }
        for (const auto& e : new_diag)
            buf.push_back(e);

        sys.H_hat_nh = coo_to_csr(buf, sys.H_hat_nh.nrows, sys.H_hat_nh.ncols);
    };

    if (g.params.mode == Polarization::TM) {
        add_vv_diag(TM::BLK_VZ);
    } else {
        add_vv_diag(TE::BLK_VX);
        add_vv_diag(TE::BLK_VY);
    }

    return sys;
}

// =============================================================================
// § 6. Perturbation operator  δĤ  (diagonal of V,V block only)
// =============================================================================

PerturbationOperator build_perturbation_op(const HermitianSystem& sys,
                                            const YeeGrid& g,
                                            Real gamma) {
    PerturbationOperator op;
    op.gamma = gamma;

    const Int N2 = g.n_dof;
    op.delta_H_diag_V.assign(static_cast<std::size_t>(N2), Complex{0.0, 0.0});
    op.n_metal = 0;

    const int V_blk = (g.params.mode == Polarization::TM)
                    ? TM::BLK_VZ : TE::BLK_VX;

    for (Int k = 0; k < N2; ++k) {
        if (!g.metal_mask[static_cast<std::size_t>(k)]) continue;
        const Int gk = block_offset(V_blk, N2) + k;
        const Real si = sys.sqrtA_inv[static_cast<std::size_t>(gk)];
        // δĤ_{kk} = s_inv[k]² × (−iΓ) = −iΓ/(ε∞ ωp²)
        op.delta_H_diag_V[static_cast<std::size_t>(k)] = Complex{0.0, -gamma} * (si * si);
        ++op.n_metal;
    }
    return op;
}

// =============================================================================
// § 7. Eigenvector back-transformation
// =============================================================================

void backtransform_eigvec(CVec& y_hat, const RVec& sqrtA_inv) {
    for (std::size_t i = 0; i < y_hat.size(); ++i)
        y_hat[i] *= sqrtA_inv[i];
}

void backtransform_all(std::vector<CVec>& eigvecs, const RVec& sqrtA_inv) {
    for (auto& v : eigvecs)
        backtransform_eigvec(v, sqrtA_inv);
}

// =============================================================================
// § 8. Field-block extraction
// =============================================================================

CVec extract_block(const CVec& x, int blk, Int N2) {
    const std::size_t start = static_cast<std::size_t>(block_offset(blk, N2));
    const std::size_t sz    = static_cast<std::size_t>(N2);
    return CVec(x.begin() + start, x.begin() + start + sz);
}

RVec field_intensity(const CVec& block) {
    RVec out(block.size());
    for (std::size_t i = 0; i < block.size(); ++i)
        out[i] = std::norm(block[i]);   // |z|² = z.real²+z.imag²
    return out;
}

// =============================================================================
// § 9. ShiftInvertOp — dense LU factorization of (Ĥ − σI)
// =============================================================================

void ShiftInvertOp::factor(const CSRMatrix& H_hat, Real sigma) {
    n_dof = H_hat.nrows;
    const std::size_t n = static_cast<std::size_t>(n_dof);

    // Convert sparse Ĥ to dense column-major layout for LAPACK.
    LU_dense.assign(n * n, Complex{0.0, 0.0});
    for (Int r = 0; r < H_hat.nrows; ++r) {
        const Int s = H_hat.row_ptr[static_cast<std::size_t>(r)];
        const Int e = H_hat.row_ptr[static_cast<std::size_t>(r + 1)];
        for (Int k = s; k < e; ++k) {
            const Int c = H_hat.col_ind[static_cast<std::size_t>(k)];
            // Column-major: index = col * n + row
            LU_dense[static_cast<std::size_t>(c) * n + static_cast<std::size_t>(r)]
                = H_hat.val[static_cast<std::size_t>(k)];
        }
    }

    // Subtract σI.
    for (std::size_t i = 0; i < n; ++i)
        LU_dense[i * n + i] -= Complex{sigma, 0.0};

    // Factor with LAPACK zgetrf.
    ipiv.resize(n);          // ipiv is std::vector<int> — matches LAPACK int*
    int n_int = static_cast<int>(n);
    int info  = 0;
    zgetrf_(&n_int, &n_int, LU_dense.data(), &n_int,
            ipiv.data(), &info);

    if (info != 0)
        throw std::runtime_error(
            "[ShiftInvertOp::factor] zgetrf failed with INFO=" +
            std::to_string(info) + ". Matrix may be singular (check sigma > 0).");

    factored = true;
}

void ShiftInvertOp::apply(const CVec& x, CVec& y) const {
    if (!factored)
        throw std::runtime_error("[ShiftInvertOp::apply] Not yet factored.");

    const std::size_t n = static_cast<std::size_t>(n_dof);
    y = x;   // LAPACK solves in-place

    char trans = 'N';
    int  n_int = static_cast<int>(n);
    int  nrhs  = 1;
    int  info  = 0;

    // zgetrs solves (Ĥ − σI) y = x using the stored LU factorization.
    // LU_dense and ipiv are mutable members; no const_cast needed.
    zgetrs_(&trans, &n_int, &nrhs,
            LU_dense.data(), &n_int,
            ipiv.data(),
            y.data(), &n_int, &info);

    if (info != 0)
        throw std::runtime_error(
            "[ShiftInvertOp::apply] zgetrs failed with INFO=" +
            std::to_string(info));
}

ShiftInvertOp build_shift_invert(const HermitianSystem& sys) {
    ShiftInvertOp op;
    op.factor(sys.H_hat, sys.sigma);
    return op;
}

ShiftInvertOp build_shift_invert_nh(const NonHermitianSystem& sys) {
    ShiftInvertOp op;
    op.factor(sys.H_hat_nh, sys.sigma);
    return op;
}

// =============================================================================
// § 10. System diagnostics
// =============================================================================

void print_system_info(const HermitianSystem& sys, int verbosity) {
    if (verbosity < 1) return;
    const Real r = gershgorin_radius(sys.H_hat);
    std::cout
        << "  n_dof      = " << sys.n_dof << "\n"
        << "  n_nullspace= " << sys.n_nullspace
        <<   " (air-V DOFs, excluded by shift-invert)\n"
        << "  sigma      = " << sys.sigma << " ωp\n"
        << "  Gershgorin ≈ " << r << " ωp  (rough spectral radius)\n"
        << "  LU memory  ≈ "
        << (static_cast<double>(sys.n_dof) * sys.n_dof * 16) / 1e6
        << " MB\n";
}
