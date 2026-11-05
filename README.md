# SPECTRE

**Parallel Pseudo-Spectral Navier–Stokes Solver for Turbulent Channel Flow**

SPECTRE is a two-dimensional, parallel, pseudo-spectral solver for incompressible turbulent channel flow. The streamwise direction is treated spectrally via FFTW, while the wall-normal direction uses second-order finite differences. Time integration follows a fractional-step (projection) method with an implicit treatment of viscous diffusion and explicit treatment of convection. Turbulence is modelled using Spectral Vanishing Viscosity (SVV) combined with a Chollet–Lesieur subgrid-scale viscosity.

---

## Table of Contents

1. [Overview](#overview)
2. [Requirements](#requirements)
3. [Directory Layout](#directory-layout)
4. [Compilation](#compilation)
5. [Execution](#execution)
6. [Configuration File](#configuration-file)
7. [Output](#output)
8. [Numerical Methods](#numerical-methods)
9. [Code Architecture](#code-architecture)
10. [Module Reference](#module-reference)
11. [Parallelization Strategy](#parallelization-strategy)
12. [Extending the Solver](#extending-the-solver)
13. [References](#references)

---

## Overview

SPECTRE solves the incompressible Navier–Stokes equations in a 2D channel geometry:

- **x-direction:** periodic, Fourier spectral discretization
- **y-direction:** wall-bounded, second-order finite differences
- **Turbulence model:** SVV (7 selectable kernels) + Chollet–Lesieur SGS viscosity
- **Time integration:** fractional-step projection method with adaptive CFL-based time stepping
- **Parallelization:** MPI with 1D domain decomposition in the y-direction

The solver is designed for direct numerical simulation (DNS) and large-eddy simulation (LES) of canonical turbulent channel flow at moderate Reynolds numbers.

---

## Requirements

| Component | Requirement |
|-----------|-------------|
| Fortran compiler | `gfortran`, `ifort`, or equivalent with Fortran 2003 support |
| MPI | MPI-2 or later (e.g., OpenMPI, MPICH, Intel MPI) |
| FFTW | FFTW 3.x (`libfftw3`) with Fortran interface `fftw3.f03` |
| Build system | `make` |
| Environment | POSIX shell (tested on Linux and NCI Gadi) |

The `fftw3.f03` include file must be discoverable by the compiler. On most systems it is installed alongside the FFTW development package.

---

## Directory Layout

```
spectre_bin/
├── main.f90              Program entry point; reads namelist, broadcasts, launches solver
├── solver.f90            Time integration loop and fractional-step orchestration
├── navier.f90            Convection, viscous, body force, pressure, correction steps
├── poisson.f90           Thomas TDMA solver and matrix assemblers
├── field.f90             Derivative operators, divergence, Rhie–Chow
├── turbulence.f90        SVV kernels, Chollet–Lesieur viscosity, perturbation seeding
├── fft.f90               FFTW plans, Fourier transforms, wavenumbers, dealiasing
├── comm.f90              MPI init, halo exchange, gather/scatter, mode transposes
├── domain.f90            Domain decomposition, neighbour ranks, interior indices
├── diagnostics.f90       Energy, enstrophy, divergence, CFL diagnostics
├── vti_writer.f90        VTK ImageData (.vti) output writer
├── spectre_config.nml    Runtime configuration (namelist)
├── Makefile              Build rules
└── output/               Directory for .vti output files
```

---

## Compilation

### New build

```bash
cd spectre_bin
make
```

### Clean rebuild

```bash
make clean
make
```

The build produces the executable `spectre`.

---

## Execution

### First run (after compiling)

```bash
nano spectre_config.nml       # edit simulation parameters
mpirun -np <num_procs> ./spectre
```

### Re-run without recompilation

Simulation parameters may be changed in `spectre_config.nml` at any time. No recompilation is required.

```bash
nano spectre_config.nml
make init                     # resets the output directory
mpirun -np <num_procs> ./spectre
```

### Running on NCI Gadi

```bash
make
nano spectre_config.nml
qsub myNCIJob.sh
```

The number of MPI ranks must not exceed `Ny_g`. For the mode transposes to be well balanced, `num_ranks_mpi` should divide `Nx` evenly (or nearly evenly; `Decompose_Domain` handles remainder distribution).

---

## Configuration File

All runtime parameters are supplied via the Fortran namelist file `spectre_config.nml`, which is read by rank 0 and broadcast to all other ranks.

```nml
&spectre_parameters
    Nx = 256,
    Ny_g = 257,
    Re = 2000.0,
    epsilon_svv = 3.0e-3,
    Cs = 0.25,
    cfl_target = 0.1,
    n_steps = 80000,
    dpdx = -1,
    eps_pert = 2,
    k_mode_seed = 2,
    kc_ratio = 0.1,
    svv_model = 4
/
```

### Parameter Reference

| Parameter | Type | Description |
|-----------|------|-------------|
| `Nx` | integer | Number of grid points in the periodic x-direction |
| `Ny_g` | integer | Number of grid points in the wall-normal y-direction (global) |
| `Re` | real | Reynolds number, `Re = U*h/ν` |
| `epsilon_svv` | real | SVV amplitude, controls the strength of spectral vanishing viscosity |
| `Cs` | real | Chollet–Lesieur constant for the SGS viscosity |
| `cfl_target` | real | Target CFL number for adaptive time stepping |
| `n_steps` | integer | Total number of time steps |
| `dpdx` | real | Constant driving pressure gradient (only used if mass-flux control is disabled) |
| `eps_pert` | real | Amplitude of the initial velocity perturbation (turbulence seeding) |
| `k_mode_seed` | integer | Fourier mode used for the initial seeding perturbation |
| `kc_ratio` | real | High-pass cutoff ratio `kc/k_Nyq` for the turbulence models |
| `svv_model` | integer | SVV kernel selector (see below) |

### SVV Kernel Selection

| `svv_model` | Kernel |
|:-----------:|--------|
| `1` | Exponential |
| `2` | Kolmogorov Sharp |
| `3` | Kolmogorov Smooth |
| `4` | Plateau |
| `5` | Power (Exponent 2) |
| `6` | Power (Exponent 3) |
| `7` | Power (Exponent 4) |
| _default_ | Exponential |

---

## Output

### Console output

At startup, rank 0 prints the full parameter summary:

```
============================================================================================
     SPECTRE - Parallel, Pseudo-Spectral Navier-Stokes Solver for Turbulent Channel Flow
=============================================================================================
 Grid Size Nx:                  256
 Grid Size Ny_g:                257
 Reynolds Number Re:           2000.00
 ...
=============================================================================================
```

During the run, the simulation time is printed every 500 steps. At termination, the total wall-clock time is reported.

### Field output

Every 500 steps (and at the final step), the solver writes the velocity and pressure fields as a VTK ImageData file:

```
output/flow_field_<step>.vti
```

The file contains scalar arrays `u`, `v`, `p` and a vector array `velocity`. It can be opened directly in ParaView.

### Enabling diagnostics

Per-step solver diagnostics (kinetic energy, enstrophy, maximum velocity, divergence, CFL, effective viscosity range) are available via `Solver_Diagnostics` in `diagnostics.f90`. They are currently commented out inside `Fractional_Step_Solver`; uncomment the block to enable them.

---

## Numerical Methods

### Governing equations

Incompressible Navier–Stokes in rotational form with an effective viscosity:

```
∂u/∂t + (u·∇)u = -∇p + ∇·[(ν + ν_sgs + ν_svv) ∇u] + f
∇·u = 0
```

where:

- `ν` = molecular viscosity (`1/Re`)
- `ν_sgs` = Chollet–Lesieur SGS viscosity
- `ν_svv` = spectral vanishing viscosity from the selected kernel
- `f` = body force from the driving pressure gradient

### Fractional-step scheme

Each time step executes the following sequence (see `navier.f90`):

1. **Explicit convection** — non-linear term evaluated in physical space with 2/3-rule dealiasing (`Explicit_Convection_Step`)
2. **Body force** — pressure-gradient forcing added to the x-momentum (`Body_Force_Step`)
3. **Implicit viscous diffusion** — tridiagonal solve per Fourier mode (`Implicit_Viscous_Step`)
4. **Pressure Poisson** — tridiagonal solve per Fourier mode using a precomputed matrix (`Pressure_Poisson_Step`)
5. **Velocity correction** — pressure gradient removed from the intermediate velocity (`Velocity_Correction_Step`)
6. **Adaptive time step** — recalculated from the CFL constraint (`Adaptive_Timestep`)

### Discretization

| Term | Method |
|------|--------|
| ∂/∂x | Spectral (FFTW complex-to-complex along x) |
| ∂/∂y | Second-order centered finite differences |
| Convection | Explicit (one-step) |
| Diffusion | Implicit (tridiagonal per mode) |
| Pressure | Poisson solve per mode with precomputed matrix |

### Turbulence models

**Chollet–Lesieur SGS viscosity** (`Chollet_Lesieur_Viscosity`):

```
ν_sgs(k) = Cs * sqrt(E(kc)) / kc^{3/2} * f(k/kc)
f(x)     = sqrt(1 + 34.5 * exp(-3.03/x))
```

The energy spectrum `E(k)` is computed via an MPI all-reduce across all ranks.

**SVV kernels** (`turbulence.f90`): each kernel returns a spectral viscosity `ν_svv(k)` that activates only for `|k| > kc`. Selectable via `svv_model` in the namelist.

### Dealiasing

The convection term is dealiased using the 2/3 rule (`Dealias` in `fft.f90`). Stronger variants (`Dealias_Strong`, `Dealias_Stronger`) are available for experimentation.

---

## Code Architecture

The solver is organized into modular Fortran 90 units, each with a single responsibility:

```
main.f90
   └── solver.Fractional_Step_Solver
          ├── turbulence.Seed_Perturbation
          ├── poisson.Pressure_Matrix_Precompute
          └── per-step:
                ├── turbulence.Chollet_Lesieur_Viscosity
                ├── turbulence.*_SVV (per mode)
                ├── navier.Explicit_Convection_Step
                ├── navier.Body_Force_Step
                ├── navier.Implicit_Viscous_Step
                │     ├── comm.Transpose_To_Modes
                │     ├── poisson.Viscous_Poisson_Assembler
                │     ├── poisson.Thomas_TDMA_Solver
                │     └── comm.Transpose_From_Modes
                ├── navier.Pressure_Poisson_Step
                ├── navier.Velocity_Correction_Step
                ├── fft.Dealias
                ├── vti_writer.Write_VTI
                └── solver.Adaptive_Timestep
```

### Data flow

Velocity and pressure fields are stored in **Fourier space** throughout the time loop:

```
u(Ny+2, Nx)  complex  — x-Fourier, y-physical
v(Ny+2, Nx)  complex  — x-Fourier, y-physical
p(Ny+2, Nx)  complex  — x-Fourier, y-physical
```

The two ghost rows (`j=1` and `j=Ny+2`) hold halo values exchanged between ranks. Physical-space evaluation is performed only where non-linear products are required (convection, diagnostics, output).

---

## Module Reference

### `main.f90`
Program entry. Initializes MPI, reads `spectre_config.nml` on rank 0, broadcasts all parameters, prints the banner, calls `Fractional_Step_Solver`, reports timing, and finalizes MPI.

### `solver.f90`
- `Fractional_Step_Solver` — the main time loop and orchestration of all physical steps.
- `Adaptive_Timestep` — recomputes `dt` from the CFL constraint based on the current global maximum velocity.
- `Mass_Flux_Body_Force` — optional mass-flux controller using momentum balance and a soft relaxation.
- `Filter_Y` — optional 3-point smoothing filter in the y-direction (currently unused).

### `navier.f90`
Individual sub-steps of the fractional-step method:
- `Explicit_Convection_Step`
- `Implicit_Viscous_Step`
- `Body_Force_Step`
- `Pressure_Poisson_Step`
- `Velocity_Correction_Step`

### `poisson.f90`
- `Thomas_TDMA_Solver` — complex tridiagonal solver.
- `Viscous_Poisson_Assembler` — assembles the viscous tridiagonal per mode (inside the time loop, since `ν_eff` changes each step).
- `Pressure_Poisson_Assembler` — assembles the pressure Poisson tridiagonal per mode.
- `Pressure_Matrix_Precompute` — assembles and stores the pressure matrix once, before the time loop.

### `field.f90`
Generic spectral and finite-difference operators:
- `Dx_Spectral`, `Dy_Forward`, `Dy_Backward`, `Dy_Center`
- `Divergence`
- `Rhie_Chow_Face_Velocity` (available but currently unused)

### `fft.f90`
FFTW wrappers:
- `Forward_Plan`, `Backward_Plan` — batch 1D complex FFT plans.
- `Fourier_Transform` — execute a plan on an array.
- `Wavelengths` — fills the `kx` array (positive modes then negative).
- `Dealias`, `Dealias_Strong`, `Dealias_Stronger` — 2/3, 1/3, 1/6 dealiasing rules.

### `comm.f90`
MPI collectives:
- `Init_MPI`
- `Exchange_Halos` — exchange ghost rows between neighbours using vector datatypes and non-blocking `Irecv`/`Isend`.
- `Gather_Column_To_0` / `Scatter_Column_From_0` — column-wise gather/scatter (used by the VTI writer).
- `Transpose_To_Modes` / `Transpose_From_Modes` — one `MPI_Alltoallv` per field per step so that each rank owns all rows of its own subset of Fourier modes.

### `domain.f90`
- `Decompose_Domain` — splits `Ny_g` rows (or `Nx` modes) across ranks, returning local size and local-to-global map.
- `Get_Neighbour_Ranks` — north/south neighbours, with `MPI_PROC_NULL` at boundaries.
- `Get_Interior_Indices` — returns the first and last non-ghost, non-boundary row on this rank.

### `turbulence.f90`
- Seven SVV kernels: `Exponential_SVV`, `Kolmogorov_Sharp_SVV`, `Kolmogorov_Smooth_SVV`, `Plateau_SVV`, `Power_2_SVV`, `Power_3_SVV`, `Power_4_SVV`.
- `Chollet_Lesieur_Viscosity` — SGS viscosity based on the global energy spectrum.
- `Seed_Perturbation` — initializes the velocity field with a single-mode perturbation: `ψ = ε cos(kx) sin²(πy/Ly)`.

### `diagnostics.f90`
`Solver_Diagnostics` — computes kinetic energy, enstrophy, maximum velocity components, maximum divergence, CFL number, and effective viscosity range via MPI all-reduce. Prints a formatted block on rank 0.

### `vti_writer.f90`
`Write_VTI` — gathers the real-space velocity and pressure to rank 0 and writes an ASCII VTK `.vti` file.

---

## Parallelization Strategy

SPECTRE uses a **1D row-wise decomposition** in the y-direction:

- Each rank owns `Ny` interior rows plus two ghost rows (`j=1` and `j=Ny+2`).
- Halo exchange is done with a vector datatype of stride `Ny+2` and count `Nx`.
- The Thomas solver requires all `Ny_g` rows of a Fourier mode, which a single rank does not have. Rather than gather every column to rank 0, SPECTRE **transposes** the data once per solve:
  - `Transpose_To_Modes`: after the transpose, rank `r` holds all `Ny_g` rows of the modes `m_start..m_end` it owns (the same range produced by `Decompose_Domain(Nx, ...)`).
  - `Transpose_From_Modes`: sends the solved modes back to the row-wise layout.
- Each rank therefore solves its own share of the tridiagonal systems locally, and no rank becomes a serialization bottleneck.

The two transposes per step are the dominant communication cost. Their payload is a single `MPI_Alltoallv` of size `Ny × Nx` complex values per field per direction.

### Guidance on rank counts

- For the mode decomposition to be balanced, `num_ranks_mpi` should evenly divide (or nearly divide) both `Ny_g` and `Nx`.
- Typical production values: 8–64 ranks for `Nx = 256`, `Ny_g = 257`.
- More ranks than `Nx` provides no benefit for the mode transpose (a rank will own zero modes).

---

## Extending the Solver

### Adding a new SVV kernel

1. Add a subroutine to `turbulence.f90` following the signature of the existing kernels.
2. Add a new `case` in the `select case (svv_model)` block in `Fractional_Step_Solver` (`solver.f90`).
3. Add a new `case` in the `select case (svv_model)` block in `main.f90` (banner printing).
4. Document the new `svv_model` value in the configuration table above.

### Adding a new diagnostic

Extend `Solver_Diagnostics` in `diagnostics.f90`. All per-rank values must be reduced with `MPI_Allreduce` before printing. Enable the diagnostic block in `Fractional_Step_Solver` (currently commented out).

### Changing the time-stepping scheme

The explicit and implicit steps are isolated in `navier.f90`. To swap the time integrator (e.g., Adams–Bashforth 2, RK3), modify `Explicit_Convection_Step` and its call site in `Fractional_Step_Solver`; the implicit viscous and pressure steps are unaffected.

### Changing the dealiasing rule

`fft.f90` provides `Dealias` (2/3), `Dealias_Strong` (1/3) and `Dealias_Stronger` (1/6). Swap the call sites in `Explicit_Convection_Step` and at the end of `Fractional_Step_Solver`.

### Running without SVV

Set `epsilon_svv = 0` in the namelist. The SVV contribution `ν_svv` will be identically zero and the solver reduces to LES with Chollet–Lesieur only. The `svv_model` value is still read but has no effect.

### Running pure DNS

Set `epsilon_svv = 0` and `Cs = 0`. Note that without either model the solver may become numerically unstable at high Reynolds numbers; the 2/3-rule dealiasing provides only limited stabilization.

---

## References

The numerical methods implemented in SPECTRE are based on the following publications. Users are encouraged to cite the relevant works when publishing results obtained with this solver.

### Fractional-Step / Projection Method

1. **Chorin, A. J.** (1968). Numerical solution of the Navier–Stokes equations. *Mathematics of Computation*, **22**(104), 745–762.  
   The original fractional-step projection method for incompressible flow.

2. **Temam, R.** (1969). Sur l'approximation de la solution des équations de Navier–Stokes par la méthode des pas fractionnaires (II). *Archive for Rational Mechanics and Analysis*, **33**(5), 377–385.  
   The independent formulation and analysis of the projection method.

### Chollet–Lesieur Subgrid-Scale Model

3. **Chollet, J.-P., & Lesieur, M.** (1981). Parameterization of small scales of three-dimensional isotropic turbulence utilizing spectral closures. *Journal of the Atmospheric Sciences*, **38**(12), 2747–2757.  
   The spectral eddy-viscosity model used in `turbulence.f90` for `ν_sgs`.

### Spectral Vanishing Viscosity (SVV)

4. **Tadmor, E.** (1989). Convergence of spectral methods for nonlinear conservation laws. *SIAM Journal on Numerical Analysis*, **26**(1), 30–44.  
   The original introduction of the spectral vanishing viscosity concept.

5. **Karamanos, G.-S., & Karniadakis, G. E.** (2000). A spectral vanishing viscosity method for large-eddy simulations. *Journal of Computational Physics*, **163**(1), 22–50.  
   The formulation of SVV for LES of the Navier–Stokes equations, including the high-pass cutoff and kernel construction used in the `svv_model` selection.

### Pseudo-Spectral Channel-Flow Algorithm

6. **Kim, J., Moin, P., & Moser, R.** (1987). Turbulence statistics in fully developed channel flow at low Reynolds number. *Journal of Fluid Mechanics*, **177**, 133–166.  
   The classical pseudo-spectral channel-flow algorithm (Fourier in the periodic directions, Chebyshev/finite differences in the wall-normal direction) that serves as the methodological predecessor for wall-bounded spectral solvers such as SPECTRE.

### FFTW

7. **Frigo, M., & Johnson, S. G.** (2005). The design and implementation of FFTW3. *Proceedings of the IEEE*, **93**(2), 216–231.  
   The FFT library used for the streamwise spectral transforms.

---

## License and Attribution

If you use SPECTRE in published work, please cite the method references listed in [References](#references), in particular:

- the fractional-step scheme (Chorin 1968; Temam 1969),
- the Chollet–Lesieur SGS model (Chollet & Lesieur 1981),
- the SVV model and the specific kernel selected via `svv_model` (Tadmor 1989; Karamanos & Karniadakis 2000).

Contact the maintainers for the preferred citation format.