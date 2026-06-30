<!--
  Authorship:
    Davide Batic  (Mathematics Department, Khalifa University of Science and
                   Technology, Abu Dhabi, UAE)
    Denys Dutykh  (Mathematics Department, Khalifa University of Science and
                   Technology, Abu Dhabi, UAE)
-->

# Global and dynamical viability of the string T-duality regular black hole

Official code repository accompanying the manuscript

> **Global and dynamical viability of the string T-duality regular black hole**
> Davide Batic and Denys Dutykh, Mathematics Department, Khalifa University of
> Science and Technology, Abu Dhabi, UAE.

The full manuscript is included as
[`pdf/T_Dual_BHs_Geodesic_Completeness.pdf`](pdf/T_Dual_BHs_Geodesic_Completeness.pdf).

## Overview

The paper analyses the global and dynamical viability of the
string-T-duality–Bardeen regular black hole, whose static, spherically
symmetric line element is

```
ds² = -f(r) dt² + dr²/f(r) + r² dΩ²,   f(r) = 1 - 2 m(r)/r,
m(r) = M r³ / (r² + ℓ₀²)^(3/2),        f(r) = 1 - 2 M r² / (r² + ℓ₀²)^(3/2),
```

where `ℓ₀` is the zero-point length induced by the T-duality–modified
propagator. The analysis separates three logically distinct questions: local
curvature regularity, global geodesic completeness, and dynamical stability of
the inner Cauchy horizon. The codes in this repository produce the manuscript
figures and the numerical convergence and stability checks used in that
analysis:

1. a set of Python (Matplotlib) scripts that render the analytic figures of the
   manuscript: the effective potentials for causal geodesics, the three local
   radial diagnostics of the core, the light-ring/horizon thresholds, the
   horizon phase diagram with the surface gravities, the zero-point-length form
   factor with the classical-versus-filtered mass inflation, and the exact
   versus endpoint proper time along the Painlevé–Gullstrand congruence;
2. a Julia reproduction of the proper-time comparison
   (geodesic-completeness diagnostic), and
3. a Chebyshev spectral quasinormal-mode (QNM) solver for the scalar and
   electromagnetic perturbation sectors, including a three-resolution
   convergence check and a positive-imaginary-part (growing-mode) test of
   Cauchy-horizon stability.

## Repository layout

```
tduality-bh/
├── README.md
├── LICENSE                                          GNU LGPL v2.1
├── python/
│   ├── ku_style.py                                  shared Matplotlib style/helpers
│   ├── fig_geodesic_potentials.py                   effective potentials (timelike/null)
│   ├── fig_three_scale_diagnostic.py                three local radial diagnostics
│   ├── fig_lightring_structure.py                   light-ring/horizon thresholds
│   ├── fig_horizon_phase_kappa.py                   horizon phase diagram + surface gravities
│   ├── fig_tduality_filter.py                       form factor + filtered mass inflation
│   └── fig_proper_time_comparison.py                exact vs endpoint proper time
├── figures/                                         rendered figure PDFs
├── julia/
│   ├── Project.toml                                 pinned Julia environment
│   ├── tau_comparison.jl                            geodesic figure
│   └── TD_spectral_solver_with_unstable_check.jl    QNM spectral solver
├── maple/
│   ├── 01-CASE-s0_s1-SUPPLEMENTAL MATERIAL.mw        symbolic derivation
│   └── Matrix_Assembler_all_s_TD-Final Version.mw    Chebyshev matrix assembler
└── pdf/
    └── T_Dual_BHs_Geodesic_Completeness.pdf          the manuscript
```

## Requirements

**Python** 3 with `numpy`, `scipy` and `matplotlib`, used by the figure scripts
in `python/`. No virtual environment is pinned; install the three packages with

```bash
pip install numpy scipy matplotlib
```

**Julia** (tested with 1.12; the environment declares `julia = "1.10"` as the
lower bound). The non-standard-library dependencies are pinned in
`julia/Project.toml`:

| Package        | Role                                                        |
| -------------- | ----------------------------------------------------------- |
| `Plots`        | figure and scatter-plot rendering                           |
| `PlotlyJS`     | interactive HTML backend for the QNM scatter plots          |
| `GenericSchur` | generalized Schur / `eigvals` for `BigFloat` matrices       |

The scripts also use the standard libraries `DelimitedFiles`, `Printf`,
`LinearAlgebra` and `Base.Threads`, which ship with Julia.

Instantiate the environment once:

```bash
julia --project=julia -e 'using Pkg; Pkg.instantiate()'
```

**Maple** 2024.0, used by the two worksheets in `maple/` (packages
`LinearAlgebra`, `FileTools`, `ExportMatrix`, `PDEtools`). The assembler runs at
200 decimal digits, so `kernelopts(maxdigits)` must be at least 200.

## Codes

### 1. Analytic figures (Python)

The scripts in `python/` render the analytic figures of the manuscript. They
share `ku_style.py`, which fixes a consistent publication style (serif text,
Computer Modern math fonts, a fixed colour palette) and exposes a `savefig`
helper that writes a vector PDF. Each script is standalone and run from inside
the `python/` directory:

```bash
cd python
python fig_geodesic_potentials.py
python fig_three_scale_diagnostic.py
python fig_lightring_structure.py
python fig_horizon_phase_kappa.py
python fig_tduality_filter.py
python fig_proper_time_comparison.py
```

| Script                          | Figure                                                                                              |
| ------------------------------- | --------------------------------------------------------------------------------------------------- |
| `fig_geodesic_potentials.py`    | timelike and null effective potentials; the centrifugal barrier blocks non-radial geodesics at `r=0`|
| `fig_three_scale_diagnostic.py` | Newtonian, Ricci-focusing and Raychaudhuri radial diagnostics, with their three distinct scales     |
| `fig_lightring_structure.py`    | horizon and circular-null thresholds in the `(x, μ)` plane; the horizonless ultracompact band       |
| `fig_horizon_phase_kappa.py`    | horizon phase diagram and the inner/outer surface gravities `ℓ₀ κ_∓`                                 |
| `fig_tduality_filter.py`        | zero-point-length form factor and the classical-versus-filtered Cauchy-horizon mass inflation       |
| `fig_proper_time_comparison.py` | exact proper time versus the endpoint prescription along the Painlevé–Gullstrand congruence          |

The rendered PDFs are collected in `figures/`. Several scripts embed `assert`
checks against the analytic special radii and asymptotics, so a successful run
also verifies those closed-form results.

### 2. Geodesic-completeness figure (Julia)

[`julia/tau_comparison.jl`](julia/tau_comparison.jl) is the Julia counterpart of
`fig_proper_time_comparison.py`. Along the selected Painlevé–Gullstrand radial infalling congruence
it compares the exact proper time `τ_exact(r)`, obtained by integrating
`dr/dτ = -v(r)` in closed form, against the endpoint expression `τ_JN(r)`
inferred from replacing the Raychaudhuri integral by `τ·F(r)`. The two curves
share the same congruence but differ in their small-`r` asymptotics
(logarithmic versus power law), which is the point of the figure. The parameters
are fixed in the file (`M = 5`, `ℓ₀ = 1`, `r_i = 4 ℓ₀`).

Run it:

```bash
julia --project=julia julia/tau_comparison.jl
```

This writes `tau_comparison.pdf` in the current working directory. It depends
only on `Plots`.

### 3. Quasinormal-mode spectral pipeline

The QNM analysis is a three-stage pipeline: symbolic derivation in Maple,
matrix assembly in Maple, and numerical eigenvalue solution in Julia.

**Stage A, symbolic derivation.**
[`maple/01-CASE-s0_s1-SUPPLEMENTAL MATERIAL.mw`](maple/01-CASE-s0_s1-SUPPLEMENTAL%20MATERIAL.mw)
is the supplemental material. For the `s = 0` (scalar) and `s = 1`
(electromagnetic) sectors it derives the QNM master equation, performs the
change of variables to the compact coordinate, and produces the `3×3` pencil
coefficient functions `L00, …, L22` together with the boundary limits tabulated
in the manuscript. It is purely symbolic and writes nothing to disk. The
worksheet asks for a Maple restart at two marked points so that symbol
assignments stay clean; follow those in-worksheet instructions when executing
it top to bottom.

**Stage B, matrix assembly.**
[`maple/Matrix_Assembler_all_s_TD-Final Version.mw`](maple/Matrix_Assembler_all_s_TD-Final%20Version.mw)
discretises the master equation with a Chebyshev spectral collocation method and
exports the matrices of the quadratic pencil

```
M0 + i Ω M1 + Ω² M2 = 0
```

at the three resolutions `N = 190, 195, 200` (`nmin/nstep/nmax = 190/5/200`).
Set the sector (`s`, `L`, `q`) and the output directory near the top of the
worksheet, then run it. It exports MATLAB-ASCII matrix files

```
M0-<N>-<filetag>.mat,  M1-<N>-<filetag>.mat,  M2-<N>-<filetag>.mat
```

into a `data/` folder, where `filetag = s<S>_L<L>_q<Q>` (for example
`s0_L0_q0`). The `N` range here must match `N_VALS` in the Julia solver.

**Stage C, numerical solution.**
[`julia/TD_spectral_solver_with_unstable_check.jl`](julia/TD_spectral_solver_with_unstable_check.jl)
reads the `M0/M1/M2` files, linearises the quadratic pencil into a generalized
eigenvalue problem, and solves it in `BigFloat` arithmetic at 220 decimal digits
(`GenericSchur` provides the `BigFloat` eigensolver). It then:

- filters the raw spectra to physical QNM candidates (discarding the static
  `Ω = 0` mode and keeping damped modes under the `exp(-i ω t)` convention);
- matches modes across the three resolutions `N = 190, 195, 200` and keeps only
  those that agree to the tolerance `1e-4` (the convergence check);
- separates oscillatory from purely imaginary converged modes;
- scans the unfiltered spectra for roots with `Im Ω > 0`, which under the same
  convention would be growing modes, as a direct test of Cauchy-horizon
  stability (the report states `No unstable mode detected.` when none survive).

Select the sector by editing the constants at the top of the file
(`SPIN`, `ELL`, `DEFORMATION_TAG`); these must match the Maple `filetag`. Point
the solver at the matrix files and run it:

```bash
export TD_DATA_DIR=/path/to/data          # folder containing the M0/M1/M2 files
export TD_OUT_DIR=/path/to/output         # optional; default is $TD_DATA_DIR/Julia/<filetag>
julia --project=julia -t auto julia/TD_spectral_solver_with_unstable_check.jl
```

If `TD_DATA_DIR` is unset, the solver looks for a `data/` folder next to the
repository's `julia/` directory. The `-t auto` flag enables multithreading.

#### Output files

All output is written into `TD_OUT_DIR` (default `$TD_DATA_DIR/Julia/<filetag>`),
with `<filetag> = s<S>_L<L>_q<Q>` and `N ∈ {190, 195, 200}`:

| File                                         | Contents                                            |
| -------------------------------------------- | --------------------------------------------------- |
| `spectrum_finite-N<N>-<filetag>.dat`         | raw finite generalized-eigenvalue spectrum          |
| `spectrum_candidates-N<N>-<filetag>.dat`     | spectrum after the physical-mode filters            |
| `stable_modes-<filetag>.txt`                 | converged damped QNMs (summary table with errors)   |
| `stable_modes-<filetag>_details.txt`         | readable per-mode details                           |
| `stable_modes-<filetag>_full_precision.txt`  | full-precision converged modes                      |
| `stable_modes_oscillatory-<filetag>.dat`     | converged oscillatory modes                         |
| `stable_modes_pure_imag-<filetag>.dat`       | converged purely imaginary modes                    |
| `unstable_candidates-N<N>-<filetag>.txt`     | per-resolution roots with `Im Ω > 0`                |
| `unstable_modes-<filetag>.txt` (+ variants)  | converged growing-mode report (or the "none" notice)|
| `qnm_plot*-<filetag>.png` / `.html`          | scatter plots of the spectra (Re Ω, Im Ω)           |

## Reproducing the results

```bash
# 1. Analytic figures (Python).
pip install numpy scipy matplotlib
cd python && for f in fig_*.py; do python "$f"; done && cd ..

# 2. Set up the Julia environment (once).
julia --project=julia -e 'using Pkg; Pkg.instantiate()'

# 3. Geodesic-completeness figure.
julia --project=julia julia/tau_comparison.jl       # -> tau_comparison.pdf

# 4. QNM spectral pipeline.
#    a) In Maple 2024, run the assembler worksheet to export M0/M1/M2 .mat
#       files into a data/ folder (set s, L, q and the output directory first).
#    b) Solve and run the convergence and stability checks:
export TD_DATA_DIR=/path/to/data
julia --project=julia -t auto julia/TD_spectral_solver_with_unstable_check.jl
```

## Citation

If you use these codes, please cite the manuscript:

```bibtex
@unpublished{BaticDutykh2026Tduality,
  author = {Batic, Davide and Dutykh, Denys},
  title  = {Global and dynamical viability of the string {T}-duality regular black hole},
  note   = {Preprint},
  year   = {2026}
}
```

## License

The code in this repository is released under the GNU Lesser General Public
License v2.1. See [`LICENSE`](LICENSE) for the full text.

## Authors

- Davide Batic, Mathematics Department, Khalifa University of Science and
  Technology, Abu Dhabi, UAE. `davide.batic@ku.ac.ae`
- Denys Dutykh, Mathematics Department, Khalifa University of Science and
  Technology, Abu Dhabi, UAE. `denys.dutykh@ku.ac.ae`
