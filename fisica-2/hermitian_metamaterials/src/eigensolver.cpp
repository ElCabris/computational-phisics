// =============================================================================
// eigensolver.cpp — ARPACK reverse-communication Arnoldi eigensolver
//
// Raman & Fan, PRL 104, 087401 (2010).
//
// ARPACK Fortran routines used:
//   znaupd_  — Arnoldi iteration (implicitly restarted)
//   zneupd_  — eigenvalue/eigenvector extraction
// Both declared via arpack_interface.f90.
// =============================================================================

#include "eigensolver.hpp"
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <numeric>
#include <random>

// =============================================================================
// ARPACK Fortran declarations
// Named with trailing underscore (GNU Fortran ABI).
// =============================================================================
extern "C" {
    void znaupd_(int*     ido,
                 char*    bmat,
                 int*     n,
                 char*    which,
                 int*     nev,
                 double*  tol,
                 Complex* resid,
                 int*     ncv,
                 Complex* v,
                 int*     ldv,
                 int*     iparam,
                 int*     ipntr,
                 Complex* workd,
                 Complex* workl,
                 int*     lworkl,
                 double*  rwork,
                 int*     info);

    void zneupd_(int*     rvec,
                 char*    howmny,
                 int*     select,
                 Complex* d,
                 Complex* z,
                 int*     ldz,
                 Complex* sigma,
                 Complex* workev,
                 char*    bmat,
                 int*     n,
                 char*    which,
                 int*     nev,
                 double*  tol,
                 Complex* resid,
                 int*     ncv,
                 Complex* v,
                 int*     ldv,
                 int*     iparam,
                 int*     ipntr,
                 Complex* workd,
                 Complex* workl,
                 int*     lworkl,
                 double*  rwork,
                 int*     info);
}

// =============================================================================
// § 1. ARPACK error handling
// =============================================================================

void check_arpack_info(int info, const std::string& ctx) {
    if (info == 0) return;

    // Warning codes: print but do not throw.
    if (info == 1) {
        std::cerr << "[ARPACK/" << ctx << "] WARNING INFO=1: "
                     "Maximum iterations reached before convergence. "
                     "Increase --max-iter or ncv.\n";
        return;
    }
    if (info == 3) {
        std::cerr << "[ARPACK/" << ctx << "] WARNING INFO=3: "
                     "No shifts applied; linear combination breakdown.\n";
        return;
    }

    // Fatal codes — each branch throws, so no fallthrough possible.
    std::string errmsg = std::string("[ARPACK/") + ctx + "] FATAL INFO="
                       + std::to_string(info) + ": ";
    switch (info) {
        case -1:    errmsg += "N (problem size) must be positive."; break;
        case -3:    errmsg += "NCV must satisfy NEV+2 <= NCV <= N."; break;
        case -6:    errmsg += "BMAT must be 'I' or 'G'."; break;
        case -8:    errmsg += "Error in internal LAPACK eigensolve."; break;
        case -9:    errmsg += "Starting vector is zero; try a different random seed."; break;
        case -9999: errmsg += "Unable to build Arnoldi factorization (matrix singular?)."; break;
        default:    errmsg += "Unspecified ARPACK error."; break;
    }
    throw std::runtime_error(errmsg);
}

// =============================================================================
// § 2. Core ARPACK reverse-communication loop
// =============================================================================

EigenResult run_arpack_znaupd(const ShiftInvertOp& op,
                               const CSRMatrix&     H_hat,
                               Int                  n,
                               const ARPACKConfig&  cfg,
                               Real                 sigma) {
    // Non-const: Fortran ABI requires non-const int* for all scalar arguments.
    int N      = static_cast<int>(n);
    int nev    = cfg.nev;
    int ncv    = cfg.effective_ncv(N);

    if (ncv > N)
        throw std::invalid_argument("ncv > n_dof: reduce nbands or increase resolution.");

    // ── ARPACK work arrays ─────────────────────────────────────────────────
    std::vector<Complex> resid(static_cast<std::size_t>(N));
    std::vector<Complex> v(static_cast<std::size_t>(N * ncv));
    int lworkl = 3 * ncv * ncv + 5 * ncv;
    std::vector<Complex> workl(static_cast<std::size_t>(lworkl));
    std::vector<Complex> workd(static_cast<std::size_t>(3 * N));
    std::vector<double>  rwork(static_cast<std::size_t>(ncv));

    // Seed starting residual with a fixed random vector for reproducibility.
    {
        std::mt19937_64 rng(42);
        std::uniform_real_distribution<double> dist(-1.0, 1.0);
        for (auto& r : resid) r = Complex{dist(rng), dist(rng)};
    }

    // ── ARPACK parameters ──────────────────────────────────────────────────
    int iparam[11] = {};
    iparam[0] = 1;              // ishift = 1: exact shifts
    iparam[2] = cfg.maxiter;    // maxiter
    iparam[6] = 3;              // mode 3: shift-invert  (OP = (Ĥ−σI)⁻¹)

    int ipntr[14] = {};

    char bmat[]  = "I";         // standard eigenvalue problem
    char which[] = "LM";        // largest magnitude of 1/(ω−σ) → closest to σ
    double tol   = cfg.tol;

    int ido  = 0;
    int info = 1;               // info=1: use supplied starting vector (resid)

    // ── Reverse-communication loop ─────────────────────────────────────────
    int n_matvec = 0;
    while (true) {
        znaupd_(&ido, bmat, &N, which, &nev, &tol,
                resid.data(), &ncv,
                v.data(), &N,
                iparam, ipntr,
                workd.data(), workl.data(), &lworkl,
                rwork.data(), &info);

        if (info < 0) check_arpack_info(info, "znaupd");

        if (ido == -1 || ido == 1) {
            // Compute y = (Ĥ − σI)⁻¹ x
            // x is in workd[ipntr[0]-1 .. ipntr[0]+N-2]  (1-based Fortran)
            // y goes to workd[ipntr[1]-1 .. ipntr[1]+N-2]
            const std::size_t x_off = static_cast<std::size_t>(ipntr[0] - 1);
            const std::size_t y_off = static_cast<std::size_t>(ipntr[1] - 1);

            CVec x_in(workd.begin() + x_off,
                      workd.begin() + x_off + N);
            CVec y_out(static_cast<std::size_t>(N));

            op.apply(x_in, y_out);

            std::copy(y_out.begin(), y_out.end(),
                      workd.begin() + y_off);
            ++n_matvec;
        }
        else if (ido == 2) {
            // B = I: y = x (copy workd[ipntr[0]] → workd[ipntr[1]])
            const std::size_t x_off = static_cast<std::size_t>(ipntr[0] - 1);
            const std::size_t y_off = static_cast<std::size_t>(ipntr[1] - 1);
            std::copy(workd.begin() + x_off,
                      workd.begin() + x_off + N,
                      workd.begin() + y_off);
        }
        else if (ido == 99) {
            break;   // converged
        }
        else {
            throw std::runtime_error("[ARPACK] Unexpected IDO=" + std::to_string(ido));
        }
    }

    // Print convergence info.
    if (cfg.verbose)
        std::cout << "  ARPACK: " << iparam[8] << " iterations, "
                  << n_matvec << " matvec calls, "
                  << iparam[4] << " converged Ritz values.\n";

    check_arpack_info(info, "znaupd (final)");

    // ── Extract eigenvalues and eigenvectors with zneupd ──────────────────
    const int nconv = iparam[4];

    std::vector<Complex> d(static_cast<std::size_t>(nev + 1));
    std::vector<Complex> z(static_cast<std::size_t>(N * nev));
    std::vector<Complex> workev(static_cast<std::size_t>(3 * ncv));
    std::vector<int>     select(static_cast<std::size_t>(ncv), 0);

    int  rvec     = 1;          // compute eigenvectors
    char howmny[] = "A";        // all nev Ritz vectors
    Complex sigma_c{sigma, 0.0};

    zneupd_(&rvec, howmny, select.data(),
            d.data(), z.data(), &N,
            &sigma_c, workev.data(),
            bmat, &N, which, &nev, &tol,
            resid.data(), &ncv,
            v.data(), &N,
            iparam, ipntr,
            workd.data(), workl.data(), &lworkl,
            rwork.data(), &info);

    check_arpack_info(info, "zneupd");

    // ── Build EigenResult ─────────────────────────────────────────────────
    EigenResult result;
    result.nconv = nconv;
    result.niter = iparam[2];

    // d[] contains the eigenvalues of (Ĥ−σI)⁻¹: λ = 1/(ω−σ)
    // Recover ω = σ + 1/λ
    for (int j = 0; j < nconv; ++j) {
        const Complex lambda = d[static_cast<std::size_t>(j)];
        const Complex omega  = Complex{sigma, 0.0} + 1.0 / lambda;
        // For the lossless case, Im[ω] should be < 1e-8; take real part.
        result.eigenvalues.push_back(omega.real());
    }

    // Store eigenvectors (columns of z, column-major, length N each).
    result.eigenvectors.resize(static_cast<std::size_t>(nconv));
    for (int j = 0; j < nconv; ++j) {
        const std::size_t off = static_cast<std::size_t>(j * N);
        result.eigenvectors[static_cast<std::size_t>(j)] =
            CVec(z.begin() + off, z.begin() + off + N);
    }

    sort_by_eigenvalue(result);
    filter_nullspace(result);

    // Compute max residual.
    const RVec resids = compute_residuals(H_hat, result, {});
    result.max_residual = resids.empty() ? 0.0
        : *std::max_element(resids.begin(), resids.end());

    return result;
}

// =============================================================================
// § 3. Public solver interfaces
// =============================================================================

EigenResult solve_hermitian(const HermitianSystem& sys,
                             const ARPACKConfig&    cfg) {
    if (cfg.verbose) {
        std::cout << "  Building shift-invert LU (sigma=" << sys.sigma << ")...\n";
    }
    Timer t;
    t.start();
    const ShiftInvertOp op = build_shift_invert(sys);
    if (cfg.verbose)
        std::cout << "  LU factored in " << t.elapsed() << " s\n";

    EigenResult result = run_arpack_znaupd(op, sys.H_hat, sys.n_dof, cfg, sys.sigma);

    // Back-transform eigenvectors: x = S⁻¹ ŷ
    backtransform_all(result.eigenvectors, sys.sqrtA_inv);

    return result;
}

ComplexEigenResult solve_nonhermitian(const NonHermitianSystem& sys,
                                       const ARPACKConfig&       cfg) {
    // For the non-Hermitian problem, run the same ARPACK loop but without
    // assuming real eigenvalues.
    ShiftInvertOp op = build_shift_invert_nh(sys);

    int N   = static_cast<int>(sys.n_dof);
    int nev = cfg.nev;
    int ncv = cfg.effective_ncv(N);

    std::vector<Complex> resid(static_cast<std::size_t>(N));
    std::vector<Complex> v(static_cast<std::size_t>(N * ncv));
    int lworkl = 3 * ncv * ncv + 5 * ncv;
    std::vector<Complex> workl(static_cast<std::size_t>(lworkl));
    std::vector<Complex> workd(static_cast<std::size_t>(3 * N));
    std::vector<double>  rwork(static_cast<std::size_t>(ncv));

    {
        std::mt19937_64 rng(42);
        std::uniform_real_distribution<double> dist(-1.0, 1.0);
        for (auto& r : resid) r = Complex{dist(rng), dist(rng)};
    }

    int iparam[11] = {};
    iparam[0] = 1;
    iparam[2] = cfg.maxiter;
    iparam[6] = 3;
    int ipntr[14] = {};

    char bmat[]  = "I";
    char which[] = "LM";
    double tol   = cfg.tol;
    int ido = 0, info = 1;

    while (true) {
        znaupd_(&ido, bmat, &N, which, &nev, &tol,
                resid.data(), &ncv, v.data(), &N,
                iparam, ipntr,
                workd.data(), workl.data(), &lworkl,
                rwork.data(), &info);

        if (info < 0) check_arpack_info(info, "znaupd_nh");

        if (ido == -1 || ido == 1) {
            const std::size_t x_off = static_cast<std::size_t>(ipntr[0] - 1);
            const std::size_t y_off = static_cast<std::size_t>(ipntr[1] - 1);
            CVec x_in(workd.begin() + x_off, workd.begin() + x_off + N);
            CVec y_out(static_cast<std::size_t>(N));
            op.apply(x_in, y_out);
            std::copy(y_out.begin(), y_out.end(), workd.begin() + y_off);
        }
        else if (ido == 2) {
            const std::size_t x_off = static_cast<std::size_t>(ipntr[0] - 1);
            const std::size_t y_off = static_cast<std::size_t>(ipntr[1] - 1);
            std::copy(workd.begin() + x_off, workd.begin() + x_off + N,
                      workd.begin() + y_off);
        }
        else if (ido == 99) break;
    }

    check_arpack_info(info, "znaupd_nh (final)");

    const int nconv = iparam[4];
    std::vector<Complex> d(static_cast<std::size_t>(nev + 1));
    std::vector<Complex> z(static_cast<std::size_t>(N * nev));
    std::vector<Complex> workev(static_cast<std::size_t>(3 * ncv));
    std::vector<int>     select(static_cast<std::size_t>(ncv), 0);

    int  rvec = 1;
    char howmny[] = "A";
    Complex sigma_c{sys.sigma, 0.0};

    zneupd_(&rvec, howmny, select.data(),
            d.data(), z.data(), &N,
            &sigma_c, workev.data(),
            bmat, &N, which, &nev, &tol,
            resid.data(), &ncv, v.data(), &N,
            iparam, ipntr,
            workd.data(), workl.data(), &lworkl,
            rwork.data(), &info);

    check_arpack_info(info, "zneupd_nh");

    ComplexEigenResult result;
    result.nconv = nconv;
    result.niter = iparam[2];

    for (int j = 0; j < nconv; ++j) {
        const Complex lambda = d[static_cast<std::size_t>(j)];
        result.eigenvalues.push_back(Complex{sys.sigma, 0.0} + 1.0 / lambda);
    }
    result.eigenvectors.resize(static_cast<std::size_t>(nconv));
    for (int j = 0; j < nconv; ++j) {
        const std::size_t off = static_cast<std::size_t>(j * N);
        result.eigenvectors[static_cast<std::size_t>(j)] =
            CVec(z.begin() + off, z.begin() + off + N);
        backtransform_eigvec(result.eigenvectors[static_cast<std::size_t>(j)],
                             sys.sqrtA_inv);
    }

    sort_by_eigenvalue(result);
    result.max_residual = 0.0;   // TODO: compute complex residuals
    return result;
}

// =============================================================================
// § 4. Post-processing
// =============================================================================

RVec compute_residuals(const CSRMatrix&   H_hat,
                        const EigenResult& result,
                        const RVec&        sqrtA) {
    // Residual in the ŷ basis: ||Ĥ ŷ − ω ŷ||₂ / |ω|
    // sqrtA is currently unused here (residuals are in the ŷ basis since
    // eigenvectors are back-transformed; re-transform if needed).
    RVec res(static_cast<std::size_t>(result.nconv), 0.0);
    (void)sqrtA;

    for (int j = 0; j < result.nconv; ++j) {
        const Real omega = result.eigenvalues[static_cast<std::size_t>(j)];
        const CVec& y    = result.eigenvectors[static_cast<std::size_t>(j)];
        CVec Hy(y.size(), Complex{0.0, 0.0});
        matvec(H_hat, y, Hy);

        Real norm_r = 0.0;
        for (std::size_t k = 0; k < y.size(); ++k) {
            const Complex r = Hy[k] - omega * y[k];
            norm_r += std::norm(r);
        }
        res[static_cast<std::size_t>(j)] =
            std::sqrt(norm_r) / (std::abs(omega) + 1e-30);
    }
    return res;
}

void sort_by_eigenvalue(EigenResult& result) {
    const int n = result.nconv;
    std::vector<int> idx(static_cast<std::size_t>(n));
    std::iota(idx.begin(), idx.end(), 0);
    std::sort(idx.begin(), idx.end(), [&](int a, int b) {
        return result.eigenvalues[static_cast<std::size_t>(a)]
             < result.eigenvalues[static_cast<std::size_t>(b)];
    });
    RVec   sorted_eval(static_cast<std::size_t>(n));
    std::vector<CVec> sorted_evec(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) {
        sorted_eval[static_cast<std::size_t>(i)] =
            result.eigenvalues[static_cast<std::size_t>(idx[static_cast<std::size_t>(i)])];
        sorted_evec[static_cast<std::size_t>(i)] =
            result.eigenvectors[static_cast<std::size_t>(idx[static_cast<std::size_t>(i)])];
    }
    result.eigenvalues  = std::move(sorted_eval);
    result.eigenvectors = std::move(sorted_evec);
}

void sort_by_eigenvalue(ComplexEigenResult& result) {
    const int n = result.nconv;
    std::vector<int> idx(static_cast<std::size_t>(n));
    std::iota(idx.begin(), idx.end(), 0);
    std::sort(idx.begin(), idx.end(), [&](int a, int b) {
        return result.eigenvalues[static_cast<std::size_t>(a)].real()
             < result.eigenvalues[static_cast<std::size_t>(b)].real();
    });
    std::vector<Complex> sorted_eval(static_cast<std::size_t>(n));
    std::vector<CVec>    sorted_evec(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) {
        sorted_eval[static_cast<std::size_t>(i)] =
            result.eigenvalues[static_cast<std::size_t>(idx[static_cast<std::size_t>(i)])];
        sorted_evec[static_cast<std::size_t>(i)] =
            result.eigenvectors[static_cast<std::size_t>(idx[static_cast<std::size_t>(i)])];
    }
    result.eigenvalues  = std::move(sorted_eval);
    result.eigenvectors = std::move(sorted_evec);
}

void filter_nullspace(EigenResult& result, Real omega_threshold) {
    EigenResult filtered;
    for (int j = 0; j < result.nconv; ++j) {
        if (std::abs(result.eigenvalues[static_cast<std::size_t>(j)]) > omega_threshold) {
            filtered.eigenvalues.push_back(
                result.eigenvalues[static_cast<std::size_t>(j)]);
            filtered.eigenvectors.push_back(
                std::move(result.eigenvectors[static_cast<std::size_t>(j)]));
        }
    }
    filtered.nconv        = static_cast<int>(filtered.eigenvalues.size());
    filtered.niter        = result.niter;
    filtered.max_residual = result.max_residual;
    result = std::move(filtered);
}

// =============================================================================
// § 5. Modal orthogonality verification  (eq. 13)
// =============================================================================

OrthogonalityReport verify_modal_orthogonality(const EigenResult& result) {
    const int n = result.nconv;
    OrthogonalityReport r{};
    r.n_modes = n;

    for (int m = 0; m < n; ++m) {
        for (int p = 0; p < n; ++p) {
            const Complex g_mp = dot_herm(result.eigenvectors[static_cast<std::size_t>(m)],
                                          result.eigenvectors[static_cast<std::size_t>(p)]);
            if (m == p) {
                const Real dev = std::abs(g_mp.real() - 1.0) + std::abs(g_mp.imag());
                r.max_diagonal_error = std::max(r.max_diagonal_error, dev);
            } else {
                r.max_off_diagonal = std::max(r.max_off_diagonal, std::abs(g_mp));
            }
        }
    }
    const Real threshold = 1e-6;
    r.passed = (r.max_off_diagonal < threshold) && (r.max_diagonal_error < threshold);
    return r;
}

// =============================================================================
// § 6. Convergence study
// =============================================================================

std::vector<ConvergenceEntry>
convergence_study(const SimParams& params,
                  Real kx, Real ky,
                  int N_start, int N_end, int N_step,
                  const ARPACKConfig& cfg) {
    std::vector<ConvergenceEntry> entries;
    for (int Ni = N_start; Ni <= N_end; Ni += N_step) {
        SimParams p = params;
        p.resolution = Ni;
        YeeGrid g(p);

        Timer t;
        t.start();
        HermitianSystem sys = assemble_hermitian(g, kx, ky);
        EigenResult     res = solve_hermitian(sys, cfg);
        const double    cpu = t.elapsed();

        ConvergenceEntry e;
        e.N           = Ni;
        e.eigenvalues = res.eigenvalues;
        e.cpu_seconds = cpu;
        entries.push_back(std::move(e));
    }
    return entries;
}
