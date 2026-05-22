// =============================================================================
// output.cpp — File output routines
//
// Raman & Fan, PRL 104, 087401 (2010).
// =============================================================================

#include "output.hpp"
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <iostream>
#include <cmath>

// =============================================================================
// § 1. Band structure
// =============================================================================

void write_bands(const std::string&            filename,
                 const std::vector<BandPoint>& bands,
                 const SimParams&              params) {
    std::ofstream f(filename);
    if (!f)
        throw std::runtime_error("Cannot open output file: " + filename);

    f << "# hermetic-bands: photonic band structure\n"
      << "# Raman & Fan, PRL 104, 087401 (2010)\n"
      << "# Mode       : " << (params.mode == Polarization::TE ? "TE" : "TM") << "\n"
      << "# Resolution : " << params.resolution << "\n"
      << "# Fill frac  : " << params.fill_fraction << "\n"
      << "# eps_inf    : " << params.eps_inf << "\n"
      << "# gamma/omegap: " << params.gamma << "\n"
      << "# Columns    : k_index  k*a/pi  omega_1/omegap  omega_2  ...\n";

    f << std::scientific << std::setprecision(8);
    for (const auto& bp : bands) {
        f << bp.k_index << "  " << bp.k_norm;
        for (auto w : bp.frequencies)
            f << "  " << w;
        f << "\n";
    }
}

// =============================================================================
// § 2. Loss comparison
// =============================================================================

void write_loss_comparison(const std::string&                 filename,
                            const std::vector<LossComparison>& data,
                            const SimParams&                   params) {
    std::ofstream f(filename);
    if (!f)
        throw std::runtime_error("Cannot open: " + filename);

    f << "# Loss comparison: exact (non-Hermitian) vs. perturbative (eq. 15)\n"
      << "# gamma/omegap = " << params.gamma << "\n"
      << "# Columns: band_idx  omega_real  im_omega_exact  im_omega_perturb  rel_err\n";

    f << std::scientific << std::setprecision(8);
    for (std::size_t j = 0; j < data.size(); ++j) {
        const auto& d = data[j];
        f << j
          << "  " << d.omega_real
          << "  " << d.im_exact
          << "  " << d.im_perturbative
          << "  " << d.relative_error
          << "\n";
    }
}

// =============================================================================
// § 3. Spatial field profile
// =============================================================================

void write_field_profile(const std::string&  filename,
                          const FieldProfile& fp,
                          const YeeGrid&      g) {
    std::ofstream f(filename);
    if (!f)
        throw std::runtime_error("Cannot open: " + filename);

    const Int N   = static_cast<Int>(g.N);
    const Int N2  = g.n_dof;
    const bool tm = (g.params.mode == Polarization::TM);

    f << "# Field profile: mode " << fp.mode_index
      << "  k-point " << fp.k_index
      << "  kx*a/pi = " << fp.kx_norm
      << "  omega/omegap = " << fp.omega << "\n";

    if (tm)
        f << "# x  y  |Ez|^2  |Hx|^2  |Hy|^2  |Vz|^2\n";
    else
        f << "# x  y  |Hz|^2  |Ex|^2  |Ey|^2  |Vx|^2  |Vy|^2\n";

    f << std::scientific << std::setprecision(6);

    for (Int j = 0; j < N; ++j) {
        for (Int i = 0; i < N; ++i) {
            const std::size_t k = static_cast<std::size_t>(flat(i, j, N));

            // Physical x,y coordinates of the primary-grid node.
            const Real x = Ez_x(i, g);
            const Real y = Ez_y(j, g);

            f << x << "  " << y;

            if (tm) {
                f << "  " << std::norm(fp.eigenvector[static_cast<std::size_t>(
                                       block_offset(TM::BLK_EZ, N2) + k)])
                  << "  " << std::norm(fp.eigenvector[static_cast<std::size_t>(
                                       block_offset(TM::BLK_HX, N2) + k)])
                  << "  " << std::norm(fp.eigenvector[static_cast<std::size_t>(
                                       block_offset(TM::BLK_HY, N2) + k)])
                  << "  " << std::norm(fp.eigenvector[static_cast<std::size_t>(
                                       block_offset(TM::BLK_VZ, N2) + k)]);
            } else {
                f << "  " << std::norm(fp.eigenvector[static_cast<std::size_t>(
                                       block_offset(TE::BLK_HZ, N2) + k)])
                  << "  " << std::norm(fp.eigenvector[static_cast<std::size_t>(
                                       block_offset(TE::BLK_EX, N2) + k)])
                  << "  " << std::norm(fp.eigenvector[static_cast<std::size_t>(
                                       block_offset(TE::BLK_EY, N2) + k)])
                  << "  " << std::norm(fp.eigenvector[static_cast<std::size_t>(
                                       block_offset(TE::BLK_VX, N2) + k)])
                  << "  " << std::norm(fp.eigenvector[static_cast<std::size_t>(
                                       block_offset(TE::BLK_VY, N2) + k)]);
            }
            f << "\n";
        }
        f << "\n";   // blank line between rows for gnuplot splot
    }
}

// =============================================================================
// § 4. Energy density profile
// =============================================================================

void write_energy_profile(const std::string&  filename,
                           const FieldProfile& fp,
                           const YeeGrid&      g) {
    std::ofstream f(filename);
    if (!f)
        throw std::runtime_error("Cannot open: " + filename);

    f << "# Energy density profile: mode " << fp.mode_index << "\n"
      << "# x  y  W_E  W_H  W_V  W_total\n";

    const Int N  = static_cast<Int>(g.N);
    const Int N2 = g.n_dof;
    f << std::scientific << std::setprecision(6);

    for (Int j = 0; j < N; ++j) {
        for (Int i = 0; i < N; ++i) {
            const std::size_t k = static_cast<std::size_t>(flat(i, j, N));
            const Real x = Ez_x(i, g);
            const Real y = Ez_y(j, g);

            Real WE = 0.0, WH = 0.0, WV = 0.0;

            if (g.params.mode == Polarization::TM) {
                WE = 0.5 * g.eps_z[k]
                   * std::norm(fp.eigenvector[static_cast<std::size_t>(
                                block_offset(TM::BLK_EZ, N2) + k)]);
                WH = 0.5 * PhysConst::mu0
                   * (std::norm(fp.eigenvector[static_cast<std::size_t>(
                                block_offset(TM::BLK_HX, N2) + k)])
                   +  std::norm(fp.eigenvector[static_cast<std::size_t>(
                                block_offset(TM::BLK_HY, N2) + k)]));
                if (g.metal_mask[k]) {
                    const Real wp = g.omega_p_mask[k];
                    if (wp > 0.0)
                        WV = 0.5 * std::norm(fp.eigenvector[static_cast<std::size_t>(
                                              block_offset(TM::BLK_VZ, N2) + k)])
                                 / (g.eps_z[k] * wp * wp);
                }
            } else {
                WH = 0.5 * PhysConst::mu0
                   * std::norm(fp.eigenvector[static_cast<std::size_t>(
                                block_offset(TE::BLK_HZ, N2) + k)]);
                WE = 0.5 * g.eps_x[k]
                   * std::norm(fp.eigenvector[static_cast<std::size_t>(
                                block_offset(TE::BLK_EX, N2) + k)])
                   + 0.5 * g.eps_y[k]
                   * std::norm(fp.eigenvector[static_cast<std::size_t>(
                                block_offset(TE::BLK_EY, N2) + k)]);
                if (g.metal_mask[k]) {
                    const Real wp = g.omega_p_mask[k];
                    if (wp > 0.0) {
                        WV = 0.5 * std::norm(fp.eigenvector[static_cast<std::size_t>(
                                              block_offset(TE::BLK_VX, N2) + k)])
                                 / (g.eps_x[k] * wp * wp)
                           + 0.5 * std::norm(fp.eigenvector[static_cast<std::size_t>(
                                              block_offset(TE::BLK_VY, N2) + k)])
                                 / (g.eps_y[k] * wp * wp);
                    }
                }
            }

            f << x << "  " << y << "  " << WE << "  " << WH
              << "  " << WV << "  " << (WE + WH + WV) << "\n";
        }
        f << "\n";
    }
}

// =============================================================================
// § 5. Convergence study
// =============================================================================

void write_convergence(const std::string& filename,
                        const std::vector<ConvergenceEntry>& data,
                        int nbands) {
    std::ofstream f(filename);
    if (!f)
        throw std::runtime_error("Cannot open: " + filename);

    f << "# Convergence study: eigenvalues vs. grid resolution\n"
      << "# N  band_1  ...  band_" << nbands << "  cpu_s\n";
    f << std::scientific << std::setprecision(8);

    for (const auto& e : data) {
        f << e.N;
        for (std::size_t j = 0; j < static_cast<std::size_t>(nbands)
                                 && j < e.eigenvalues.size(); ++j)
            f << "  " << e.eigenvalues[j];
        f << "  " << e.cpu_seconds << "\n";
    }
}

// =============================================================================
// § 6. Console progress
// =============================================================================

void report_kpoint(int k_index, Real kx_norm, const EigenResult& result,
                   int verbosity) {
    if (verbosity < 1) return;
    std::cout << std::fixed << std::setprecision(4)
              << "  k[" << std::setw(3) << k_index << "]"
              << "  k·a/π=" << kx_norm
              << "  nconv=" << result.nconv
              << "  ω₁=" << (result.eigenvalues.empty() ? 0.0 : result.eigenvalues[0])
              << "  res=" << std::scientific << std::setprecision(2)
              << result.max_residual << "\n";
}

void report_summary(const std::vector<BandPoint>& bands,
                    double elapsed_seconds,
                    int verbosity) {
    if (verbosity < 1) return;
    std::cout << "\n  Done. " << bands.size() << " k-points, "
              << (bands.empty() ? 0 : (int)bands[0].frequencies.size())
              << " bands, "
              << std::fixed << std::setprecision(2)
              << elapsed_seconds << " s total.\n";
}
