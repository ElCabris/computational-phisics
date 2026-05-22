# hermetic-bands

Photonic band structure of 2D plasmonic metamaterials via Hermitian eigenvalue
reformulation.

Implements the method of **Raman & Fan, Phys. Rev. Lett. 104, 087401 (2010)**:
the dispersive Maxwell eigenvalue problem with a Drude metallic inclusion is
cast as a standard Hermitian eigenproblem and solved with ARPACK shift-invert
(IRLM).

---

## Physics background

### Material model

A 2D square-lattice photonic crystal with lattice constant *a*.  The unit cell
contains a metallic inclusion (square or circle) embedded in air.  The metal
follows the lossless Drude model:

    ε(ω) = ε∞ (1 − ωp² / ω²)

All frequencies are normalized to the plasma frequency ωp.

### Hermitian reformulation (eq. 12, Raman & Fan)

The time-harmonic Maxwell + Drude system

    ω A x = B x

where A is a positive-definite material matrix, is converted to a standard
Hermitian problem via the congruence transformation ŷ = √A x:

    ω ŷ = Ĥ ŷ,    Ĥ = (√A)⁻¹ B (√A)⁻¹

Because √A is diagonal and positive, this is always numerically stable.

### State vector

| Mode | DOFs | Size  |
|------|------|-------|
| TM   | Hx, Hy, Ez, Vz | 4 N² |
| TE   | Hz, Ex, Ey, Vx, Vy | 5 N² |

V = ∂P/∂t is the Drude auxiliary polarization velocity.

### Discretization

2D Yee staggered grid, N×N cells, Bloch-periodic boundary conditions.
Forward differences (curl-E) acquire phase e^{+ik} at the upper boundary;
backward differences (curl-H) acquire e^{−ik}.  The curl operators satisfy
C_H = C_E†, guaranteeing Hermiticity of Ĥ for real σ (the lossless case).

### Eigensolver

ARPACK shift-invert mode (mode = 3) via `znaupd`/`zneupd`:

    (Ĥ − σI)⁻¹ ŷ = λ ŷ,   λ = 1/(ω − σ)

Choosing 0 < σ < ω_min_physical maps null-space modes (ω = 0) to λ = −1/σ < 0,
which ARPACK never retrieves when requesting largest |λ|.

The LU factorization of (Ĥ − σI) is computed once via LAPACK `zgetrf`/`zgetrs`.
Memory: 16 N⁴ bytes for the dense factor (TM: 256 N⁴ bytes total; N=20 → ~41 MB).

### Perturbative losses (eq. 15)

For small damping Γ, the first-order frequency shift is:

    ω₁ = (iΓ/2) × ∫_metal |V₀|² / (ε∞ ωp²) dV  /  ∫ W₀ dV

where W₀ is the mode energy density of the lossless mode.

---

## Dependencies

| Library | Purpose | Arch Linux | Ubuntu/Debian |
|---------|---------|-----------|---------------|
| ARPACK  | IRLM eigensolver | `sudo pacman -S arpack` | `sudo apt install libarpack2-dev` |
| LAPACK  | Dense LU (`zgetrf`/`zgetrs`) | `sudo pacman -S lapack` | `sudo apt install liblapack-dev` |
| BLAS    | Level-1/2 kernels | bundled with LAPACK | `sudo apt install libblas-dev` |
| gfortran | Fortran wrappers | `sudo pacman -S gcc-fortran` | `sudo apt install gfortran` |
| CMake ≥ 3.18 | Build system | `sudo pacman -S cmake` | `sudo apt install cmake` |

---

## Building

### CMake (recommended)

```bash
mkdir build && cd build
cmake ..
make -j$(nproc)
```

The executable is `build/hermetic-bands`.

### Makefile

```bash
make -j$(nproc)
```

Debug build with AddressSanitizer and UBSan:

```bash
make debug -j$(nproc)
```

---

## Usage

```
hermetic-bands [OPTIONS]
```

### Options

| Flag | Default | Description |
|------|---------|-------------|
| `--mode TE\|TM` | `TE` | Polarization |
| `--resolution N` | `20` | Yee grid points per unit cell side |
| `--fill-fraction f` | `0.25` | Metal inclusion side/lattice constant |
| `--shape square\|circle` | `square` | Inclusion geometry |
| `--eps-inf eps` | `1.0` | ε∞ of metal |
| `--nk nk` | `20` | k-points along Γ→X |
| `--nbands nb` | `15` | Number of bands to compute |
| `--gamma Γ` | `0.0` | Drude damping rate (sets non-Hermitian problem) |
| `--perturb` | off | Compute perturbative loss correction |
| `--sigma σ` | auto | Spectral shift (default: 0.01 ωp) |
| `--output FILE` | `bands.dat` | Output file for band data |
| `--field-mode n` | off | Export field profile for band n |
| `--field-kindex i` | off | k-point index for field export |
| `--verify-ortho` | off | Verify modal orthogonality (eq. 13) |
| `--verbosity 0\|1\|2` | `1` | Output verbosity |

---

## Reproducing paper figures

### Figure 1 — TM bands, metallic squares, s/a = 0.25

```bash
hermetic-bands --mode TM --resolution 20 --fill-fraction 0.25 \
               --shape square --nk 40 --nbands 10 \
               --output bands_TM_s025.dat
```

Expected: two photonic bands below and above the plasma frequency, with a gap
centered near ω ≈ 0.5 ωp.

### Figure 2 — TE bands

```bash
hermetic-bands --mode TE --resolution 20 --fill-fraction 0.25 \
               --shape square --nk 40 --nbands 10 \
               --output bands_TE_s025.dat
```

### Convergence check

```bash
for N in 10 15 20 25 30; do
    hermetic-bands --mode TM --resolution $N --fill-fraction 0.25 \
                   --nk 1 --nbands 5 \
                   --output bands_N${N}.dat
done
```

Eigenfrequencies at the X point (kx = π) should converge at second order in 1/N.

### Perturbative loss estimate (Γ = 0.01 ωp)

```bash
hermetic-bands --mode TM --resolution 20 --fill-fraction 0.25 \
               --gamma 0.01 --perturb \
               --output bands_lossy.dat
```

Outputs both ω_r (real part) and the perturbative Im[ω₁] for each band.

### Field profile export

```bash
hermetic-bands --mode TM --resolution 20 --fill-fraction 0.25 \
               --nk 20 --nbands 6 \
               --field-mode 3 --field-kindex 10 \
               --output bands.dat
```

Writes `field_mode3_k10.dat` with columns: ix, iy, |Ez|², |Hx|², |Hy|², |Vz|².

---

## Output format

### `bands.dat`

```
# hermetic-bands output
# mode=TM N=20 fill=0.25 shape=square
# kx/pi   band_1   band_2   ...   band_nb   [omega/omega_p]
0.000000   ...
```

Suitable for direct plotting with gnuplot:

```gnuplot
plot for [i=2:11] 'bands.dat' using 1:i with lines title ''
```

### `loss_comparison.dat` (with `--perturb`)

```
# kx/pi  band  omega_r  gamma_pert  gamma_direct
```

### `field_mode<n>_k<i>.dat`

Grid data in gnuplot `splot` format (blank line between iy blocks).

---

## Mathematical notes

### Why `znaupd` instead of `dsaupd`

The Bloch-periodic Yee operators produce complex entries in Ĥ for generic
k ∈ (0, π) because the Bloch phase e^{ik} is complex.  The real-symmetric
ARPACK variants (`dsaupd`/`dseupd`) require a real matrix and cannot be used
here.  `znaupd`/`zneupd` (complex non-symmetric Arnoldi) are used throughout,
even though Ĥ is Hermitian.

### Null-space modes

V-in-air DOFs are decoupled from all other DOFs because the Drude coupling
factor I_metal = 0 for air cells.  Their rows and columns in Ĥ are exactly
zero, producing a null space of dimension equal to the number of air V-DOFs.
The spectral shift σ > 0 maps these modes to λ = −1/σ < 0; ARPACK's request
for largest |λ| among the physical (positive) eigenvalues never retrieves them.
No projection or subspace removal is required.

### Dense vs sparse shift-invert

For the paper's validation regime (N ≤ 30, TM n_dof ≤ 3600), the dense LAPACK
LU factorization is:
- O(n_dof³) build cost: computed once, amortized over all ARPACK matvec calls
- Exact (no iterative refinement, no preconditioning tuning)
- Condition-independent: no CG breakdown near the spectral shift

Memory scaling:

| N  | TM n_dof | LU memory |
|----|----------|-----------|
| 20 | 1600     | 41 MB     |
| 30 | 3600     | 207 MB    |
| 40 | 6400     | 655 MB    |

For N > 40 (beyond the paper's validation target), a sparse direct solver
(PARDISO, MUMPS) or iterative preconditioned Krylov method would be preferable.

---

## Reference

A. Raman and S. Fan,
"Photonic Band Structure of Dispersive Metamaterials Formulated as a
Hermitian Eigenvalue Problem,"
*Phys. Rev. Lett.* **104**, 087401 (2010).
DOI: [10.1103/PhysRevLett.104.087401](https://doi.org/10.1103/PhysRevLett.104.087401)
