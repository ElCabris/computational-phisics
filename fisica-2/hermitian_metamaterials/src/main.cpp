// =============================================================================
// main.cpp — hermetic-bands CLI entry point
//
// Raman & Fan, PRL 104, 087401 (2010).
//
// Execution flow:
//   1. Parse CLI arguments → SimParams
//   2. Build YeeGrid (geometry + material masks)
//   3. For each k-point on Γ→X:
//      a. Assemble lossless Hermitian system Ĥ₀
//      b. Solve with ARPACK (shift-invert)
//      c. Back-transform eigenvectors
//      d. Optionally verify orthogonality (eq. 13)
//   4. Write bands.dat
//   5. If --loss Γ ≠ 0:
//      a. Run non-Hermitian solve
//      b. Compute perturbative corrections (eq. 15) if --perturb
//      c. Write loss_comparison.dat
//   6. If --field-mode n k: write spatial field profile
// =============================================================================

#include "utils.hpp"
#include "yee_grid.hpp"
#include "operators.hpp"
#include "matrices.hpp"
#include "eigensolver.hpp"
#include "perturbation.hpp"
#include "output.hpp"

#include <iostream>
#include <stdexcept>
#include <vector>

int main(int argc, char** argv) {
    try {
        // ── 1. Parse and validate parameters ─────────────────────────────
        SimParams params = parse_args(argc, argv);
        validate_params(params);
        print_params(params, params.verbosity);

        // ── 2. Build grid ─────────────────────────────────────────────────
        banner("Building Yee grid", params.verbosity);
        YeeGrid g(params);
        print_grid_summary(g);

        // ── 3. ARPACK solver configuration ───────────────────────────────
        ARPACKConfig arpack_cfg;
        arpack_cfg.nev     = params.nbands;
        arpack_cfg.tol     = 1e-10;
        arpack_cfg.maxiter = 500;
        arpack_cfg.verbose = (params.verbosity >= 2);

        // ── 4. k-path Γ→X ─────────────────────────────────────────────────
        const RVec kx_path = kpath_GammaX(params.nk);

        banner("Computing photonic bands", params.verbosity);

        std::vector<BandPoint> all_bands;
        all_bands.reserve(static_cast<std::size_t>(params.nk));

        // Storage for the field profile k-point (if requested)
        EigenResult    field_result;
        HermitianSystem field_sys;
        bool field_collected = false;

        Timer total_timer;
        total_timer.start();

        for (int ki = 0; ki < params.nk; ++ki) {
            const Real kx     = kx_path[static_cast<std::size_t>(ki)];
            const Real ky     = 0.0;
            const Real kx_norm = kx / PhysConst::PI;

            // Assemble lossless Hermitian system at this k-point
            HermitianSystem sys = assemble_hermitian(g, kx, ky);

            if (params.verbosity >= 2)
                print_system_info(sys, params.verbosity);

            // Solve
            EigenResult res = solve_hermitian(sys, arpack_cfg);
            report_kpoint(ki, kx_norm, res, params.verbosity);

            // Optional orthogonality check
            if (params.verify_ortho) {
                const RVec A_diag = (params.mode == Polarization::TM)
                                  ? build_A_diagonal_TM(g)
                                  : build_A_diagonal_TE(g);
                OrthogonalityReport ortho =
                    verify_physical_orthogonality(res, A_diag);
                print_orthogonality_report(ortho, params.verbosity);
            }

            // Collect field profile data
            if (ki == params.field_kindex && params.field_mode >= 0) {
                field_result   = res;
                field_sys      = sys;
                field_collected = true;
            }

            // Store band data
            BandPoint bp;
            bp.k_index    = ki;
            bp.k_norm     = kx_norm;
            bp.frequencies = res.eigenvalues;
            // Pad to nbands if fewer converged
            while ((int)bp.frequencies.size() < params.nbands)
                bp.frequencies.push_back(0.0);
            bp.frequencies.resize(static_cast<std::size_t>(params.nbands));

            all_bands.push_back(std::move(bp));
        }

        report_summary(all_bands, total_timer.elapsed(), params.verbosity);

        // ── 5. Write band structure ───────────────────────────────────────
        banner("Writing output", params.verbosity);
        write_bands(params.output_file, all_bands, params);
        if (params.verbosity >= 1)
            std::cout << "  Bands → " << params.output_file << "\n";

        // ── 6. Loss calculation ────────────────────────────────────────────
        if (params.gamma > 0.0) {
            banner("Loss calculation (Γ = " + std::to_string(params.gamma) + " ωp)",
                   params.verbosity);

            // Use the middle k-point for loss comparison (reproduces Fig. 2).
            const int ki_loss = params.nk / 2;
            const Real kx_loss = kx_path[static_cast<std::size_t>(ki_loss)];

            // Lossless reference at the chosen k-point
            HermitianSystem h0  = assemble_hermitian(g, kx_loss, 0.0);
            EigenResult lossless = solve_hermitian(h0, arpack_cfg);

            // Exact non-Hermitian solve
            NonHermitianSystem nh = assemble_nonhermitian(g, kx_loss, 0.0,
                                                           params.gamma);
            ComplexEigenResult exact_lossy = solve_nonhermitian(nh, arpack_cfg);

            // Perturbative correction
            std::vector<PerturbResult> perturb_results;
            if (params.perturb) {
                perturb_results = compute_all_perturbative(lossless, g, params.gamma);
            }

            // Comparison
            if (params.perturb) {
                const std::vector<LossComparison> cmp =
                    compare_loss_methods(lossless, exact_lossy, perturb_results);
                const std::string loss_file = "loss_comparison.dat";
                write_loss_comparison(loss_file, cmp, params);
                if (params.verbosity >= 1)
                    std::cout << "  Loss comparison → " << loss_file << "\n";
            }
        }

        // ── 7. Field profile export ────────────────────────────────────────
        if (params.field_mode >= 0 && field_collected) {
            const int n = params.field_mode;
            if (n < field_result.nconv) {
                FieldProfile fp;
                fp.mode_index  = n;
                fp.k_index     = params.field_kindex;
                fp.kx_norm     = kx_path[static_cast<std::size_t>(params.field_kindex)]
                                 / PhysConst::PI;
                fp.omega       = field_result.eigenvalues[static_cast<std::size_t>(n)];
                fp.eigenvector = field_result.eigenvectors[static_cast<std::size_t>(n)];

                const std::string field_file =
                    "field_" + std::to_string(n) + "_k"
                    + std::to_string(params.field_kindex) + ".dat";
                write_field_profile(field_file, fp, g);
                write_energy_profile("energy_" + std::to_string(n) + "_k"
                                     + std::to_string(params.field_kindex) + ".dat",
                                     fp, g);
                if (params.verbosity >= 1)
                    std::cout << "  Field profile → " << field_file << "\n";
            } else {
                std::cerr << "[warning] --field-mode " << n
                          << " exceeds number of converged modes ("
                          << field_result.nconv << ").\n";
            }
        }

        return 0;

    } catch (const std::exception& e) {
        std::cerr << "\n[hermetic-bands] Error: " << e.what() << "\n";
        return 1;
    }
}
