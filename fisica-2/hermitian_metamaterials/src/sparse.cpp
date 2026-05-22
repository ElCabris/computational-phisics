// =============================================================================
// sparse.cpp — CSR matrix engine implementation
//
// Raman & Fan, PRL 104, 087401 (2010).
// =============================================================================

#include "sparse.hpp"
#include <algorithm>
#include <numeric>
#include <cmath>
#include <stdexcept>
#include <iostream>
#include <iomanip>

// =============================================================================
// § 1. COO → CSR conversion
// =============================================================================

CSRMatrix coo_to_csr(COOBuffer& entries, Int nrows, Int ncols) {
    if (entries.empty()) {
        CSRMatrix M(nrows, ncols);
        return M;
    }

    // Sort by (row, col) lexicographically — O(nnz log nnz).
    std::sort(entries.begin(), entries.end(),
              [](const COOEntry& a, const COOEntry& b) {
                  return (a.row != b.row) ? (a.row < b.row) : (a.col < b.col);
              });

    // Merge duplicates (same row,col) by summing their values.
    COOBuffer merged;
    merged.reserve(entries.size());
    for (auto& e : entries) {
        if (!merged.empty()
            && merged.back().row == e.row
            && merged.back().col == e.col) {
            merged.back().val += e.val;
        } else {
            merged.push_back(e);
        }
    }

    // Remove explicitly-zero entries that result from cancellation.
    // This keeps the sparsity pattern tight.
    merged.erase(
        std::remove_if(merged.begin(), merged.end(),
                       [](const COOEntry& e) {
                           return std::abs(e.val) == 0.0;
                       }),
        merged.end());

    const Int nnz = static_cast<Int>(merged.size());
    CSRMatrix M(nrows, ncols);
    M.nnz = nnz;
    M.val.resize(static_cast<std::size_t>(nnz));
    M.col_ind.resize(static_cast<std::size_t>(nnz));
    M.row_ptr.assign(static_cast<std::size_t>(nrows + 1), 0);

    // Count non-zeros per row.
    for (auto& e : merged)
        M.row_ptr[static_cast<std::size_t>(e.row + 1)]++;

    // Exclusive prefix sum → row_ptr.
    for (Int r = 0; r < nrows; ++r)
        M.row_ptr[static_cast<std::size_t>(r + 1)] +=
            M.row_ptr[static_cast<std::size_t>(r)];

    // Fill val and col_ind.
    for (Int k = 0; k < nnz; ++k) {
        M.val    [static_cast<std::size_t>(k)] = merged[static_cast<std::size_t>(k)].val;
        M.col_ind[static_cast<std::size_t>(k)] = merged[static_cast<std::size_t>(k)].col;
    }

    return M;
}

// =============================================================================
// § 2. Matrix-vector product:  y = alpha * A * x + beta * y
// =============================================================================

void matvec(const CSRMatrix& A,
            Complex           alpha,
            const CVec&       x,
            Complex           beta,
            CVec&             y) {
    const std::size_t n = static_cast<std::size_t>(A.nrows);
    if (y.size() != n)
        throw std::invalid_argument("matvec: y size mismatch");
    if (x.size() != static_cast<std::size_t>(A.ncols))
        throw std::invalid_argument("matvec: x size mismatch");

    for (std::size_t i = 0; i < n; ++i) {
        Complex acc{0.0, 0.0};
        const Int start = A.row_ptr[i];
        const Int end   = A.row_ptr[i + 1];
        for (Int k = start; k < end; ++k)
            acc += A.val[static_cast<std::size_t>(k)]
                 * x[static_cast<std::size_t>(A.col_ind[static_cast<std::size_t>(k)])];
        y[i] = alpha * acc + beta * y[i];
    }
}

// =============================================================================
// § 3. Conjugate transpose  C = A†
// =============================================================================

CSRMatrix conj_transpose(const CSRMatrix& A) {
    // Build COO with conjugated values and swapped row/col.
    COOBuffer buf;
    buf.reserve(static_cast<std::size_t>(A.nnz));
    for (Int r = 0; r < A.nrows; ++r) {
        const Int start = A.row_ptr[static_cast<std::size_t>(r)];
        const Int end   = A.row_ptr[static_cast<std::size_t>(r + 1)];
        for (Int k = start; k < end; ++k) {
            const Int   c  = A.col_ind[static_cast<std::size_t>(k)];
            const Complex v = std::conj(A.val[static_cast<std::size_t>(k)]);
            buf.push_back({c, r, v});    // swapped: new row = old col
        }
    }
    return coo_to_csr(buf, A.ncols, A.nrows);
}

// =============================================================================
// § 4. Matrix arithmetic
// =============================================================================

CSRMatrix csr_add(const CSRMatrix& A, const CSRMatrix& B) {
    if (A.nrows != B.nrows || A.ncols != B.ncols)
        throw std::invalid_argument("csr_add: dimension mismatch");

    COOBuffer buf;
    buf.reserve(static_cast<std::size_t>(A.nnz + B.nnz));

    auto dump = [&](const CSRMatrix& M) {
        for (Int r = 0; r < M.nrows; ++r) {
            const Int s = M.row_ptr[static_cast<std::size_t>(r)];
            const Int e = M.row_ptr[static_cast<std::size_t>(r + 1)];
            for (Int k = s; k < e; ++k)
                buf.push_back({r,
                               M.col_ind[static_cast<std::size_t>(k)],
                               M.val[static_cast<std::size_t>(k)]});
        }
    };
    dump(A);
    dump(B);
    return coo_to_csr(buf, A.nrows, A.ncols);
}

void csr_scale(CSRMatrix& A, Complex s) {
    for (auto& v : A.val) v *= s;
}

// =============================================================================
// § 5. Diagonal and identity matrices
// =============================================================================

CSRMatrix diag_csr(const CVec& d) {
    const Int n = static_cast<Int>(d.size());
    COOBuffer buf;
    buf.reserve(static_cast<std::size_t>(n));
    for (Int i = 0; i < n; ++i)
        if (d[static_cast<std::size_t>(i)] != Complex{0.0, 0.0})
            buf.push_back({i, i, d[static_cast<std::size_t>(i)]});
    return coo_to_csr(buf, n, n);
}

CSRMatrix identity_csr(Int n, Complex s) {
    COOBuffer buf;
    buf.reserve(static_cast<std::size_t>(n));
    for (Int i = 0; i < n; ++i)
        buf.push_back({i, i, s});
    return coo_to_csr(buf, n, n);
}

// =============================================================================
// § 6. Verification utilities
// =============================================================================

// Check Hermiticity: max |Ĥ_{ij} − conj(Ĥ_{ji})|
// This tests only the operator structure (eq. 12), NOT modal orthogonality.
Real check_hermitian(const CSRMatrix& H_hat) {
    Real max_err = 0.0;
    for (Int r = 0; r < H_hat.nrows; ++r) {
        const Int s = H_hat.row_ptr[static_cast<std::size_t>(r)];
        const Int e = H_hat.row_ptr[static_cast<std::size_t>(r + 1)];
        for (Int k = s; k < e; ++k) {
            const Int c = H_hat.col_ind[static_cast<std::size_t>(k)];
            const Complex H_ij = H_hat.val[static_cast<std::size_t>(k)];
            const Complex H_ji = H_hat.get(c, r);      // O(row_nnz) lookup
            const Real err = std::abs(H_ij - std::conj(H_ji));
            if (err > max_err) max_err = err;
        }
    }
    return max_err;
}

Real relative_diff(const CSRMatrix& A, const CSRMatrix& B) {
    const Real norm_B = frobenius_norm(B);
    if (norm_B == 0.0) return frobenius_norm(A);
    const CSRMatrix diff = csr_add(A, csr_scaled(B, Complex{-1.0, 0.0}));
    return frobenius_norm(diff) / norm_B;
}

Real frobenius_norm(const CSRMatrix& A) {
    Real sum = 0.0;
    for (auto& v : A.val) sum += std::norm(v);   // norm = |v|²
    return std::sqrt(sum);
}

void print_dense(const CSRMatrix& A, int max_dim) {
    const Int nr = std::min(A.nrows, static_cast<Int>(max_dim));
    const Int nc = std::min(A.ncols, static_cast<Int>(max_dim));
    std::cout << std::fixed << std::setprecision(4);
    for (Int i = 0; i < nr; ++i) {
        for (Int j = 0; j < nc; ++j) {
            Complex v = A.get(i, j);
            std::cout << "(" << std::setw(7) << v.real()
                      << "," << std::setw(7) << v.imag() << ") ";
        }
        std::cout << "\n";
    }
}
