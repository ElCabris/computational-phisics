#pragma once
// =============================================================================
// sparse.hpp — Manually implemented CSR (Compressed Sparse Row) matrix engine
//
// Raman & Fan, PRL 104, 087401 (2010).
//
// ALL matrices in this project (curl operators, A, B, Ĥ) are stored in CSR
// format.  No external sparse-matrix library is used.
//
// Why CSR (not COO or CSC)?
//   - Row-major traversal in matvec matches CPU cache lines.
//   - ARPACK's reverse-communication interface requires y = A*x in place;
//     CSR performs this in one cache-friendly pass over val[].
//   - Column-major CSC would be equally fast but requires transposing the
//     operator ordering, which is less natural for curl-based assembly.
//
// Memory layout (0-based indexing):
//
//   row_ptr[i]   : index of the first non-zero of row i in val[] and col_ind[].
//   row_ptr[i+1] : one-past-last non-zero of row i.
//   col_ind[k]   : column index of non-zero k.
//   val[k]       : value of non-zero k.
//
// Example (3×3 matrix):
//   [ 1+2i   0    3  ]       val     = [1+2i, 3, 4, 5+i, 6]
//   [ 0      4    0  ]       col_ind = [0,    2, 1, 0,   2]
//   [ 5+i    0    6  ]       row_ptr = [0, 2, 3, 5]
//
// Assembly workflow:
//   1. Collect (row, col, value) triplets in a COO buffer.
//   2. Call coo_to_csr() to sort, deduplicate (sum duplicates), and compress.
//   3. Use the resulting CSRMatrix in matvec() calls.
//
// Hermitian-symmetry note:
//   The Hermitian operator Ĥ = (√A)^{-1} B (√A)^{-1}  (eq. 12, Raman & Fan)
//   is stored as a FULL complex matrix (both upper and lower triangles).
//   We do NOT exploit symmetry in storage because znaupd requires the full
//   operator and because assembling only the upper triangle while maintaining
//   correct Bloch phases would be error-prone.
// =============================================================================

#ifndef HERMETIC_BANDS_SPARSE_HPP
#define HERMETIC_BANDS_SPARSE_HPP

#include "utils.hpp"

#include <vector>
#include <algorithm>
#include <numeric>
#include <stdexcept>
#include <cmath>

// =============================================================================
// § 1. COO (Coordinate) triplet — intermediate assembly format
// =============================================================================

struct COOEntry {
    Int     row;
    Int     col;
    Complex val;
};

// A COO buffer is simply a list of triplets.
// Duplicate (row,col) pairs are SUMMED during coo_to_csr().
using COOBuffer = std::vector<COOEntry>;

// =============================================================================
// § 2. CSR matrix
// =============================================================================

struct CSRMatrix {
    Int               nrows{0};
    Int               ncols{0};
    Int               nnz{0};       // number of stored non-zeros
    std::vector<Complex> val;       // non-zero values,   length nnz
    std::vector<Int>  col_ind;      // column indices,    length nnz
    std::vector<Int>  row_ptr;      // row pointers,      length nrows+1

    // Default-constructed matrix is empty.
    CSRMatrix() = default;

    CSRMatrix(Int nr, Int nc)
        : nrows(nr), ncols(nc), nnz(0),
          row_ptr(static_cast<std::size_t>(nr + 1), 0) {}

    bool empty() const { return nnz == 0; }

    // Return number of non-zeros in row i.
    Int row_nnz(Int i) const {
        return row_ptr[static_cast<std::size_t>(i + 1)]
             - row_ptr[static_cast<std::size_t>(i)];
    }

    // Read element (i,j) — O(row_nnz) search.
    // Returns 0 if the element is not stored.
    Complex get(Int i, Int j) const;

    // Memory usage in bytes.
    std::size_t memory_bytes() const {
        return val.size()     * sizeof(Complex)
             + col_ind.size() * sizeof(Int)
             + row_ptr.size() * sizeof(Int);
    }
};

// =============================================================================
// § 3. CSR construction from COO
//
// Algorithm (two-pass, O(nnz log nnz)):
//   Pass 1 — sort triplets lexicographically by (row, col).
//   Pass 2 — merge duplicates by summing values; build row_ptr, col_ind, val.
//
// Complexity: O(nnz log nnz) sort + O(nnz) fill.
// =============================================================================

// Build a CSRMatrix from an unsorted COO buffer.
// Duplicate (row,col) pairs are summed (additive assembly, standard FEM/FDFD).
CSRMatrix coo_to_csr(COOBuffer& entries, Int nrows, Int ncols);

// Reserve estimated space in a COO buffer to avoid repeated reallocations.
// For Yee operators, each row has at most 5 non-zeros (self + 4 neighbors).
inline void coo_reserve(COOBuffer& buf, Int nrows, Int max_per_row = 5) {
    buf.reserve(static_cast<std::size_t>(nrows * max_per_row));
}

// Convenience: append one entry.
inline void coo_push(COOBuffer& buf, Int row, Int col, Complex val) {
    buf.push_back({row, col, val});
}

// =============================================================================
// § 4. Matrix-vector product   y = alpha * A * x + beta * y
//
// This is the primary kernel called by ARPACK's reverse-communication loop
// at every Lanczos iteration.  It must be as fast as possible.
//
// Inner loop accesses val[] and col_ind[] sequentially (cache-friendly).
// =============================================================================

void matvec(const CSRMatrix& A,
            Complex           alpha,
            const CVec&       x,
            Complex           beta,
            CVec&             y);

// Simplified form: y = A * x  (overwrites y).
inline void matvec(const CSRMatrix& A, const CVec& x, CVec& y) {
    matvec(A, Complex{1.0, 0.0}, x, Complex{0.0, 0.0}, y);
}

// =============================================================================
// § 5. Conjugate transpose (Hermitian adjoint)   C = A†
//
// Used to verify Hermitian symmetry: Ĥ = Ĥ† (eq. 12, Raman & Fan).
// Also used to construct the B operator from its sub-blocks.
//
// Algorithm: transpose COO and conjugate values, then re-sort.
// =============================================================================

CSRMatrix conj_transpose(const CSRMatrix& A);

// =============================================================================
// § 6. Matrix arithmetic
// =============================================================================

// C = A + B  (must have the same sparsity pattern or will be merged).
CSRMatrix csr_add(const CSRMatrix& A, const CSRMatrix& B);

// A *= s  (in-place scale of all non-zeros).
void csr_scale(CSRMatrix& A, Complex s);

// C = s * A  (out-of-place).
inline CSRMatrix csr_scaled(const CSRMatrix& A, Complex s) {
    CSRMatrix C = A;
    csr_scale(C, s);
    return C;
}

// =============================================================================
// § 7. Block matrix assembly
//
// The Hermitian operator Ĥ is built from sub-blocks:
//
//   For TM (Drude, ω₀=0, lossless):
//
//       Ĥ = (√A)^{-1} B (√A)^{-1}
//
//   The state vector is x = (Hz, Ez, Vz)^T  [or (Hz, Ex, Ey, Vx, Vy) for TE].
//   Each field block occupies N² consecutive rows/columns.
//
// block_offset(block_id, N) returns the row/column offset of a given field
// block.  Block ordering is defined in matrices.hpp and referenced here for
// index calculations.
// =============================================================================

// Row/column offset of block b when each block has n_dof entries (n_dof = N²).
inline Int block_offset(Int b, Int n_dof) noexcept { return b * n_dof; }

// =============================================================================
// § 8. Diagonal matrix (stored as CSR for uniform interface)
//
// A diagonal matrix D with entries d[i] is trivially sparse: exactly N
// non-zeros, one per row.  Storing it in CSR is slightly wasteful (row_ptr
// costs N+1 integers) but maintains a uniform interface for all operators.
// =============================================================================

// Build a CSR diagonal matrix from a vector of diagonal entries.
CSRMatrix diag_csr(const CVec& d);

// Build a CSR scalar multiple of the identity: D = s * I_n.
CSRMatrix identity_csr(Int n, Complex s = {1.0, 0.0});

// =============================================================================
// § 9. Matrix verification utilities
//
// IMPORTANT DISTINCTION:
//
//   check_hermitian()  — property of the MATRIX Ĥ alone.
//     Verifies Ĥ_{ij} = conj(Ĥ_{ji}) for all stored pairs.
//     This is a structural check on the operator, independent of any
//     eigenvectors.  Relevant to eq. (12) of Raman & Fan.
//
//   Modal orthogonality — property of the EIGENVECTORS under the A-metric.
//     Condition:  <x_m | A | x_n> = δ_{mn}  (eq. 13, Raman & Fan).
//     This is NOT a matrix property; it requires computed eigenvectors and
//     the weight matrix A.  Implemented in perturbation.hpp::verify_modal_orthogonality().
//
// Do NOT use check_hermitian() to test modal orthogonality or vice versa.
// =============================================================================

// Returns max_{i,j} |Ĥ_{ij} − conj(Ĥ_{ji})|  over all stored non-zeros.
// Purely a matrix-symmetry check (eq. 12).  A well-assembled Ĥ should give
// a value below O(N · ε_machine) where ε_machine ≈ 2.2e-16 for double.
Real check_hermitian(const CSRMatrix& H_hat);

// Returns max |A_{ij} − B_{ij}| / ||B||_F  (relative Frobenius deviation).
Real relative_diff(const CSRMatrix& A, const CSRMatrix& B);

// Print a dense representation of A to stdout (only for small N, debugging).
void print_dense(const CSRMatrix& A, int max_dim = 16);

// Compute the Frobenius norm: ||A||_F = sqrt(sum |a_ij|^2).
Real frobenius_norm(const CSRMatrix& A);

// =============================================================================
// § 10. CSRMatrix::get() — inline definition
// =============================================================================

inline Complex CSRMatrix::get(Int i, Int j) const {
    const Int start = row_ptr[static_cast<std::size_t>(i)];
    const Int end   = row_ptr[static_cast<std::size_t>(i + 1)];
    for (Int k = start; k < end; ++k) {
        if (col_ind[static_cast<std::size_t>(k)] == j)
            return val[static_cast<std::size_t>(k)];
    }
    return Complex{0.0, 0.0};
}

#endif // HERMETIC_BANDS_SPARSE_HPP
