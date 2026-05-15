// =============================================================================
// utils.cpp — Implementation of utility functions declared in utils.hpp
//
// Raman & Fan, PRL 104, 087401 (2010).
// =============================================================================

#include "utils.hpp"
#include <iostream>
#include <iomanip>
#include <sstream>
#include <stdexcept>

// =============================================================================
// CLI argument parsing
//
// Parses argv into a SimParams struct.  Unknown flags trigger a usage message
// and exit.  Invalid numeric values throw std::invalid_argument.
// =============================================================================

SimParams parse_args(int argc, char** argv) {
    SimParams p;

    auto next_arg = [&](int& i, const char* flag) -> std::string {
        if (i + 1 >= argc) {
            throw std::invalid_argument(
                std::string("Flag ") + flag + " requires an argument.");
        }
        return argv[++i];
    };

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];

        if (arg == "--help" || arg == "-h") {
            print_usage(argv[0]);
            std::exit(0);
        }
        else if (arg == "--mode") {
            std::string m = next_arg(i, "--mode");
            if (m == "TE" || m == "te") p.mode = Polarization::TE;
            else if (m == "TM" || m == "tm") p.mode = Polarization::TM;
            else throw std::invalid_argument("--mode must be TE or TM");
        }
        else if (arg == "--resolution") {
            p.resolution = std::stoi(next_arg(i, "--resolution"));
            if (p.resolution < 4)
                throw std::invalid_argument("--resolution must be >= 4");
        }
        else if (arg == "--fill-fraction") {
            p.fill_fraction = std::stod(next_arg(i, "--fill-fraction"));
            if (p.fill_fraction <= 0.0 || p.fill_fraction >= 1.0)
                throw std::invalid_argument("--fill-fraction must be in (0,1)");
        }
        else if (arg == "--nk") {
            p.nk = std::stoi(next_arg(i, "--nk"));
            if (p.nk < 1)
                throw std::invalid_argument("--nk must be >= 1");
        }
        else if (arg == "--nbands") {
            p.nbands = std::stoi(next_arg(i, "--nbands"));
            if (p.nbands < 1)
                throw std::invalid_argument("--nbands must be >= 1");
        }
        else if (arg == "--loss") {
            p.gamma = std::stod(next_arg(i, "--loss"));
            if (p.gamma < 0.0)
                throw std::invalid_argument("--loss (gamma) must be >= 0");
        }
        else if (arg == "--perturb") {
            p.perturb = true;
        }
        else if (arg == "--output") {
            p.output_file = next_arg(i, "--output");
        }
        else if (arg == "--field-mode") {
            p.field_mode   = std::stoi(next_arg(i, "--field-mode"));
            p.field_kindex = std::stoi(next_arg(i, "--field-mode (k)"));
        }
        else if (arg == "--verify-orthogonality") {
            p.verify_ortho = true;
        }
        else if (arg == "--verbosity") {
            p.verbosity = std::stoi(next_arg(i, "--verbosity"));
        }
        else if (arg == "--shape") {
            std::string s = next_arg(i, "--shape");
            if (s == "square" || s == "SQUARE")
                p.shape = InclusionShape::SQUARE;
            else if (s == "circle" || s == "CIRCLE")
                p.shape = InclusionShape::CIRCLE;
            else throw std::invalid_argument("--shape must be square or circle");
        }
        else if (arg == "--eps-inf") {
            p.eps_inf = std::stod(next_arg(i, "--eps-inf"));
            if (p.eps_inf <= 0.0)
                throw std::invalid_argument("--eps-inf must be > 0");
        }
        else {
            throw std::invalid_argument("Unknown flag: " + arg);
        }
    }
    return p;
}

// =============================================================================
// Usage message
// =============================================================================

void print_usage(const char* prog) {
    std::cout <<
        "Usage: " << prog << " [OPTIONS]\n\n"
        "Compute photonic band structure of a 2D plasmonic crystal\n"
        "using the Hermitian formulation of Raman & Fan, PRL 104, 087401 (2010).\n\n"
        "Options:\n"
        "  --mode [TE|TM]          Polarization (default: TE)\n"
        "  --resolution N          Yee grid points per unit cell (default: 20)\n"
        "  --fill-fraction f       Metallic bar size / lattice constant (default: 0.25)\n"
        "  --shape [square|circle] Inclusion geometry (default: square)\n"
        "  --eps-inf VALUE         High-frequency permittivity ε∞ (default: 1.0)\n"
        "  --nk N                  k-points along Gamma->X (default: 20)\n"
        "  --nbands N              Number of bands to compute (default: 15)\n"
        "  --loss GAMMA            Damping rate in units of ωp (default: 0.0)\n"
        "  --perturb               Compute perturbative loss correction (eq. 15)\n"
        "  --output FILE           Output file (default: bands.dat)\n"
        "  --field-mode n k        Export spatial field profile for mode n at k-point k\n"
        "  --verify-orthogonality  Verify modal orthogonality (eq. 13)\n"
        "  --verbosity [0|1|2]     Output verbosity (default: 1)\n"
        "  --help                  Show this message\n\n"
        "Examples:\n"
        "  " << prog << " --mode TE --resolution 20 --nbands 15\n"
        "  " << prog << " --mode TM --loss 0.01 --perturb --output bands_lossy.dat\n"
        "  " << prog << " --field-mode 3 10 --output field.dat\n";
}

// =============================================================================
// Parameter validation and summary
// =============================================================================

void validate_params(const SimParams& p) {
    if (p.fill_fraction <= 0.0 || p.fill_fraction >= 1.0)
        fatal("fill_fraction must be in (0,1)");
    if (p.resolution < 4)
        fatal("resolution must be >= 4");
    if (p.nbands < 1)
        fatal("nbands must be >= 1");
    if (p.nk < 1)
        fatal("nk must be >= 1");
    if (p.gamma < 0.0)
        fatal("gamma must be >= 0");
    if (p.eps_inf <= 0.0)
        fatal("eps_inf must be > 0");

    // For the Drude model (omega_0 = 0), the V field is always present.
    // For a Lorentz oscillator (omega_0 > 0), the P field would also be
    // needed — but that extension is outside the scope of this implementation.
    if (p.omega_0 != 0.0) {
        std::cerr << "[warning] omega_0 != 0: Lorentz case not yet implemented; "
                     "using Drude (omega_0 = 0).\n";
    }

    // Warn if nbands is close to the Krylov subspace limit for the given N.
    const int n2 = p.resolution * p.resolution;
    const int ndof = (p.mode == Polarization::TM) ? 4 * n2 : 5 * n2;
    const int ncv_min = std::max(2 * p.nbands + 1, 20);
    if (ncv_min >= ndof) {
        std::cerr << "[warning] nbands=" << p.nbands
                  << " is too large for resolution=" << p.resolution
                  << " (n_dof=" << ndof << "). Reduce nbands or increase resolution.\n";
    }
}

void print_params(const SimParams& p, int verbosity) {
    if (verbosity < 1) return;
    auto yn = [](bool b) { return b ? "yes" : "no"; };
    std::cout << std::fixed << std::setprecision(4);
    std::cout
        << "\n┌─ Simulation parameters ─────────────────────────────────┐\n"
        << "│  Mode              : " << (p.mode == Polarization::TE ? "TE" : "TM") << "\n"
        << "│  Resolution N      : " << p.resolution << "\n"
        << "│  Fill fraction s/a : " << p.fill_fraction << "\n"
        << "│  Inclusion shape   : " << (p.shape == InclusionShape::SQUARE ? "square" : "circle") << "\n"
        << "│  ε∞                : " << p.eps_inf << "\n"
        << "│  ωp (normalized)   : " << p.omega_p << "\n"
        << "│  Γ  (normalized)   : " << p.gamma << "\n"
        << "│  k-points (Γ→X)    : " << p.nk << "\n"
        << "│  Bands             : " << p.nbands << "\n"
        << "│  Perturbative      : " << yn(p.perturb) << "\n"
        << "│  Verify ortho      : " << yn(p.verify_ortho) << "\n"
        << "│  Output file       : " << p.output_file << "\n"
        << "└──────────────────────────────────────────────────────────┘\n\n";
}
