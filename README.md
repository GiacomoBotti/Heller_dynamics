# FRIEZA

**FRIEZA** is a Fortran program for computing Gaussian wavepacket autocorrelation functions and power spectra from *ab initio* molecular dynamics trajectories. The code provides two complementary formulations:

- a **normal-mode implementation** supporting Frozen Gaussian (FGWP),  Single-Hessian Gaussian (SH-GWP), and Heller Thawed Gaussian (TGWP)  propagation;
- a **Cartesian implementation** of the Frozen Gaussian method, intended  in particular for large systems where it is useful to select active atoms and freeze spectator atoms directly in Cartesian space.

``` text
   ███████╗██████╗ ██╗███████╗███████╗ █████╗
   ██╔════╝██╔══██╗██║██╔════╝╚══███╔╝██╔══██╗
   █████╗  ██████╔╝██║█████╗    ███╔╝ ███████║
   ██╔══╝  ██╔══██╗██║██╔══╝   ███╔╝  ██╔══██║
   ██║     ██║  ██║██║███████╗███████╗██║  ██║
   ╚═╝     ╚═╝  ╚═╝╚═╝╚══════╝╚══════╝╚═╝  ╚═╝

        Gaussian Wavepacket Spectroscopy
               from ab initio data
```

## Implementations

### Normal-mode implementation

The normal-mode code works in the vibrational subspace of the molecule, with

``` text
nvib = 3 * nat - rotrasl
```

and normally uses `rotrasl = 6` for nonlinear molecules. The equilibrium Cartesian Hessian is mass-weighted, diagonalized, and transformed into a fixed equilibrium normal-mode basis.

The normal-mode implementation supports three propagation schemes:
- **Frozen Gaussian (FGWP)**: the Gaussian width is kept fixed at its equilibrium harmonic value.
- **Single-Hessian Gaussian (SH-GWP)**: the Gaussian width is propagated using the equilibrium Hessian throughout the trajectory.
- **Thawed Gaussian (TGWP)**: the Gaussian width is propagated using the instantaneous Hessian supplied for each trajectory snapshot.

The vibrational mask selects normal modes to follow the trajectory. Modes with `mask(i) = 0` remain at their initial coordinate and momentum.

### Cartesian implementation

The Cartesian code works directly with

``` text
ncart = 3 * nat
```

and does not transform trajectory quantities into normal modes. Cartesian coordinates and momenta are mass-weighted in the same way as in the normal-mode implementation.

The Cartesian implementation currently provides **Frozen Gaussian (FGWP) only**. The `calculation` variable and `rotrasl` are retained in the input interface for compatibility, but are ignored by the Cartesian executable. The trajectory-Hessian file is likewise retained for compatibility and is not used in the Cartesian FGWP calculation.

The initial Cartesian Gaussian width is obtained from the matrix square root of the full mass-weighted Cartesian Hessian,

``` text
A0 = sqrt(H)
```

where the square root is evaluated from the eigen-decomposition 
``` text
H = C diag(lambda) C^T
A0 = C diag(sqrt(lambda)) C^T
```

Small negative Hessian eigenvalues can be treated as numerical noise through an explicit threshold. The six translational/rotational zero modes of a nonlinear molecule are regularized with a small positive width so that the full Cartesian Gaussian matrices remain numerically  invertible.

The Cartesian mask is specified per atom rather than per normal mode. If an atom is active (`mask(i) = 1`), its three Cartesian components follow the trajectory. If an atom is inactive (`mask(i) = 0`), its Cartesian center and momentum remain at their initial values. The initial momentum is masked as well, so the correlation function is normalized at `t = 0` for a frozen
atom subset.

This makes the Cartesian formulation particularly useful when many spectator atoms are present and only a selected group of atoms should be allowed to follow the trajectory.

## Correlation function and numerical stabilization

The Gaussian overlap/correlation function is evaluated by `correlation_mod.f90` using complex matrix inversion with LAPACK.

The determinant factors entering the Gaussian overlap are evaluated in logarithmic form rather than by explicitly multiplying all determinant factors. In particular, the code uses a complex logarithmic determinant for the matrix

``` text
W = At + transpose(A0)
```

and combines the logarithmic prefactor with the exponential part before evaluating the final complex exponential. This avoids overflow/underflow problems in high-dimensional Cartesian calculations.

The initial Gaussian determinant is also stored as `logdetA0` rather than as a direct determinant. This is important for large systems because the Cartesian dimension scales as `3 * nat`.

## Large-system Cartesian phase shift

The Cartesian implementation includes an optional fixed phase convention in the Fourier transform to keep high-frequency spectra inside the usable FFT range. The current implementation applies a shift of

``` text
phase_shift = 0.75 * Eref
```

where `Eref` is the coherent Cartesian ZPE-like reference constructed from the initial width matrix.

With `FFTW_BACKWARD`, the correlation function is multiplied by

``` text
C_shift(t) = C(t) exp(+i phase_shift t)
```

so a spectral feature at energy `E` is moved numerically to

``` text
E_FFT = E - phase_shift
```

The frequency grid is then labeled with

``` text
omega = omega_FFT + phase_shift
```

Therefore the final spectrum is only **relabelled**; no interpolation or post-processing of the spectral intensities is required.

## Files

The source is divided into the following components:

| File | Purpose |
| --- | --- |
| `main_nm.f90` | Normal-mode main program: input parsing, mass weighting, normal-mode transformation, FGWP/SH-GWP/TGWP propagation, correlation function, and FFT. |
| `main_cartesian.f90` | Cartesian Frozen-Gaussian main program: Cartesian mass weighting, full Cartesian Hessian square root, atom masking, correlation function, phase shifting, and FFT. |
| `correlation_mod.f90` | Gaussian overlap/correlation-function evaluation, complex matrix inversion, determinant evaluation, and logarithmic determinant evaluation. |
| `constants.f90` | Atomic-unit conversion constants and numerical constants. |
| `diagonalizer.f` | Symmetric-matrix diagonalization used for the equilibrium Hessian. |
| `Makefile` | Compilation and linking rules for both executables. |

The Makefile produces two executables:

``` text
frieza_nm.x
frieza_cart.x
```

## Requirements

A Fortran compiler and the following numerical libraries are required:

- `gfortran`
- FFTW3, including the Fortran 2003 interface (`fftw3.f03`) and library   `libfftw3`
- BLAS
- LAPACK

The supplied Makefile links with:

``` text
-lfftw3 -lm -llapack -lblas
```

and expects the FFTW include file to be available through `/usr/include`.

## Compilation

The Makefile provides separate targets for the two implementations:

``` bash
make compile_nm
make compile_cart
```

which produce:

``` text
frieza_nm.x
frieza_cart.x
```

To remove object files, modules, and executables:

``` bash
make .clean
```

## Running the program

Both executables read a file named exactly:

``` text
input
```

and expect a `banner.txt` file because the banner is printed at startup.

Run the normal-mode implementation with:

``` bash
./frieza_nm.x
```

Run the Cartesian implementation with:

``` bash
./frieza_cart.x
```

The same namelist structure is used by both implementations:

1. `input_files`
2. `options`
3. `trajectory`

The interpretation of `mask` differs between the two executables and is
described below.

## Input file

A representative input is:

``` fortran
&input_files
    geom_eq="h2o_opt.xyz"
    traj="parsed_log_traj.xyz"
    velocity="velocity.xyz"
    hess="final_py.out"
    output="correlation.dat"
/
&options
    calculation = 0
    eta = 0.d0
    rotrasl = 6
/
&trajectory
    steps=2500
    dt=8.2682749151502d0
    padding=170000
    mask(1) = 0
    mask(2) = 1
    mask(3) = 1
/
```

The defaults defined by the current programs are:

``` text
geom_eq    = "equilibrium_geometry.xyz"
hess_eq    = "Hessian_flat.out"
traj       = "parsed_log_traj.xyz"
hess       = "Hessian_traj.out"
velocity   = "velocity.xyz"
output     = "correlation.dat"
fourierout = "fourier.dat"
powerout   = "power.dat"
calculation = 0
eta         = 0.d0
rotrasl     = 6
steps       = 2500
dt          = 8.2682749151502d0
padding     = 0
mask        = 1
```

The final values of `steps`, `dt`, `padding`, and `mask` are read from the
`trajectory` namelist after these defaults are assigned.

### Namelist: `input_files`

| Variable | Meaning |
| --- | --- |
| `geom_eq` | Equilibrium geometry file. |
| `hess_eq` | Equilibrium Hessian file. |
| `traj` | *Ab initio* Cartesian trajectory. |
| `hess` | Trajectory Hessian file. Used by the normal-mode TGWP calculation; retained but unused by the Cartesian FGWP implementation. |
| `velocity` | Initial velocity file. |
| `output` | Correlation-function output file. |
| `fourierout` | Complex Fourier-transform output file. |
| `powerout` | Power-spectrum output file. |

### Namelist: `options`

| Variable | Normal-mode implementation | Cartesian implementation |
| --- | --- | --- |
| `calculation` | `0`: FGWP, `1`: SH-GWP, `2` or greater: TGWP. | Ignored; Cartesian code is FGWP only. |
| `eta` | Exponential damping coefficient. | Exponential damping coefficient. |
| `rotrasl` | Number of rotational/translational modes removed from the Cartesian Hessian; normally `6` for nonlinear molecules. | Retained for input compatibility; not used in the Cartesian dimensionality, which is always `3 * nat`. |

## Input files

### 1. Equilibrium geometry: `geom_eq`

The equilibrium geometry is read as an XYZ-like file:

``` text
N
comment line
Atom1   x1   y1   z1
Atom2   x2   y2   z2
...
```

The first line gives the number of atoms. The second line is skipped. Cartesian coordinates are interpreted as **Angstrom** and converted internally to atomic units.

Atomic masses are assigned internally from the element symbol. The current supported species are:

``` text
H  D  C  N  O  S  P  I
```

`D` is treated as deuterium and internally relabelled as `H1` after the mass is assigned.

### 2. Equilibrium Hessian: `hess_eq`

The equilibrium Hessian is expected to contain two header lines followed by the lower triangle of the Hessian, with one matrix element per line:

``` text
header line 1
header line 2
H(1,1)
H(1,2)
H(2,2)
H(1,3)
H(2,3)
H(3,3)
...
```

The full symmetric Hessian is reconstructed and mass-weighted according to

``` text
H_mw(i,j) = H(i,j) / sqrt(m_i m_j)
```

In the normal-mode implementation, the mass-weighted Hessian is diagonalized to obtain the equilibrium normal modes and harmonic frequencies.

In the Cartesian implementation, the mass-weighted Cartesian Hessian is diagonalized only as a numerical route to its matrix square root. The full `3 * nat` Cartesian width is then reconstructed as `A0 = sqrt(H_mw)`.

Because a full Cartesian Hessian contains rigid-body zero modes, small negative eigenvalues are tested against `hessian_neg_tol`, and zero/small modes are assigned the finite `width_floor_factor` width required to keep the Gaussian matrices invertible.

### 3. Trajectory: `traj`

Each trajectory snapshot is expected to contain two header records followed by one line per atom. Each atom line contains the symbol, Cartesian coordinates, and Cartesian velocities:

``` text
header line 1
header line 2
Atom1   x1   y1   z1   vx1   vy1   vz1
Atom2   x2   y2   z2   vx2   vy2   vz2
...
```

Coordinates are interpreted as **Angstrom** and converted internally to atomic units. Velocities are expected to already be in atomic units; no velocity conversion is applied.

The normal-mode implementation transforms the mass-weighted Cartesian coordinates and velocities into the fixed equilibrium normal-mode basis. The Cartesian implementation keeps them directly in the mass-weighted Cartesian representation.

### 4. Initial velocities: `velocity`

The initial velocity file has the same two-line header convention, followed by one line per atom containing the atomic symbol and three velocity components:

``` text
header line 1
header line 2
Atom1   vx1   vy1   vz1
Atom2   vx2   vy2   vz2
...
```

In the normal-mode implementation, these velocities are projected onto the equilibrium normal modes to construct the initial Gaussian momentum.

In the Cartesian implementation, they are mass-weighted directly. The atom mask is then applied to the initial momentum so that frozen atoms have zero initial momentum.

### 5. Trajectory Hessians: `hess`

Trajectory Hessians are used only by the normal-mode TGWP calculation. For each trajectory step, two header records are followed by the lower triangle of an `ncart x ncart` Cartesian Hessian:

``` text
header line 1
header line 2
H(1,1)
H(1,2)
H(2,2)
...
```

Each Hessian is mass-weighted and transformed into the fixed equilibrium normal-mode basis before being used to propagate the TGWP width.

The Cartesian executable does not read trajectory Hessians during its FGWP calculation.

## Calculation modes: normal-mode implementation

The normal-mode `calculation` variable selects:

| `calculation` | Method | Width propagation |
| ---: | --- | --- |
| `0` | Frozen Gaussian | Fixed equilibrium width `A0`. |
| `1` | Single-Hessian Gaussian | Width propagated using the equilibrium Hessian. |
| `2` or greater | Thawed Gaussian | Width propagated using the instantaneous trajectory Hessian. |

The initial width is diagonal in the equilibrium normal-mode basis:

``` text
A0(i,i) = sqrt(lambda_i)
```

where `lambda_i` are the positive vibrational eigenvalues of the mass-weighted
equilibrium Hessian.

## Atom masking: Cartesian implementation

The Cartesian executable allocates `mask` with `nat` entries. Each atom mask entry is expanded to its three Cartesian components internally:

``` text
mask(i) = 1  ->  x_i, y_i, z_i active
mask(i) = 0  ->  x_i, y_i, z_i frozen
```

For an active atom, the instantaneous Cartesian trajectory is used. For a frozen atom, the coordinate is held at the initial Cartesian center and the momentum is held at its initial masked value, which is zero for the initial state.

For example, a three-atom system with

``` fortran
mask(1) = 0
mask(2) = 1
mask(3) = 1
```

is internally represented as

``` text
0 0 0  1 1 1  1 1 1
```

The Cartesian mask is therefore an atom-space mask rather than a vibrational-mode mask.

## Mode masking: normal-mode implementation

In the normal-mode executable, `mask` has `nvib` entries. For

``` text
mask(i) = 1
```

the trajectory coordinate and momentum of normal mode `i` are used. For

``` text
mask(i) = 0
```

the corresponding coordinate and momentum remain at their initial values.

For a nonlinear H2O molecule with `rotrasl = 6`,

``` text
nvib = 3 * nat - rotrasl = 9 - 6 = 3
```

so three mask entries are required.

## Damping

The parameter

``` text
eta
```

introduces an exponential damping factor in the correlation function:

``` text
C(t) -> C(t) exp(-eta t)
```

with `t` in atomic time units. Setting `eta = 0` disables damping.

## Time and zero padding

The trajectory contains `steps` propagated snapshots with time increment
`dt`:

``` text
t_k = k * dt
```

The FFT length is

``` text
fftsteps = steps + padding
```

`padding` adds zero-valued samples to the FFT array. The FFT frequency grid contains positive and negative frequencies; `power.dat` writes the non-negative branch.

## Output files

### Correlation-function output

The autocorrelation function is written to the filename specified by `output`.

The file contains:

``` text
# Time          Re[C(t)]            Im[C(t)]            |C(t)|**2
```

The columns are:

| Column | Meaning |
| ---: | --- |
| 1 | Time. |
| 2 | Real part of the correlation function. |
| 3 | Imaginary part of the correlation function. |
| 4 | Squared modulus of the correlation function. |

### `fourier.dat`

The complex Fourier transform is written as:

``` text
# Freq          Ang Freq            Re[FFT]             Im[FFT]
```

The first column is the ordinary FFT frequency in inverse atomic time units; the second is the angular frequency used for output. In the Cartesian implementation the angular-frequency grid includes the fixed phase shift described above.

### `power.dat`

The power spectrum is written as:

``` text
# Freq          Ang Freq [au]       Ang Freq [cm**-1]   |FFT|**2
```

The four columns are:

| Column | Meaning                                          |
| -----: | ------------------------------------------------ |
|      1 | Ordinary frequency in inverse atomic time units. |
|      2 | Angular frequency in atomic units.               |
|      3 | Angular frequency converted to cm^-1.            |
|      4 | Power spectrum, `FFT^2`                          |

Only the non-negative FFT branch is written to `power.dat`.

## Method overview

### Normal-mode workflow

1. Read the equilibrium geometry and assign atomic masses.
2. Read and mass-weight the equilibrium Cartesian Hessian.
3. Diagonalize the equilibrium Hessian to obtain the equilibrium normal modes.
4. Transform the initial Cartesian coordinates and velocities into the fixed normal-mode basis.
5. Construct the initial Gaussian width from the harmonic frequencies.
6. For each trajectory snapshot, transform the Cartesian coordinates and velocities into the same fixed normal-mode basis.
7. Propagate the Gaussian width according to the selected FGWP, SH-GWP, or TGWP scheme.
8. Evaluate the Gaussian overlap correlation function.
9. Apply the optional exponential damping factor.
10. Fourier transform the correlation function with FFTW3.
11. Write the complex Fourier transform and the positive-frequency power spectrum.

### Cartesian workflow

1. Read the equilibrium geometry and assign atomic masses.
2. Read and mass-weight the full Cartesian equilibrium Hessian.
3. Diagonalize the Hessian and construct the full Cartesian coherent width `A0 = sqrt(H_mw)`, with numerical regularization of near-zero modes.
4. Read the initial velocities and construct the mass-weighted Cartesian initial momentum, applying the atom mask.
5. For each trajectory snapshot, mass-weight the Cartesian coordinates and velocities and apply the atom mask.
6. Keep the Gaussian width fixed at `A0`.
7. Evaluate the Cartesian Gaussian overlap correlation function using the logarithmic determinant formulation.
8. Apply exponential damping and the Cartesian phase shift.
9. Fourier transform the shifted correlation function with FFTW3.
10. Relabel the angular-frequency grid by adding the phase shift so that the reported spectrum remains on the physical frequency scale.

## Units

The internal dynamics use atomic units. The external coordinate convention is:

- geometry and trajectory coordinates: Angstrom on input;
- velocities: atomic units as supplied;
- Hessians: atomic-unit Hessian elements, followed by mass weighting;
- frequencies: converted to cm^-1 for printed spectral data.

The normal-mode executable prints the harmonic frequencies and harmonic zero-point energy. The Cartesian executable prints the coherent Cartesian ZPE-like reference used for its phase convention.

## Example workflow

A typical calculation directory contains:

``` text
.
├── banner.txt
├── input
├── equilibrium_geometry.xyz
├── Hessian_flat.out
├── parsed_log_traj.xyz
├── Hessian_traj.out
├── velocity.xyz
├── constants.f90
├── diagonalizer.f
├── correlation_mod.f90
├── main_nm.f90
├── main_cartesian.f90
└── Makefile
```

Compile the desired implementation:

``` bash
make compile_nm
# or
make compile_cart
```

Then run:

``` bash
./frieza_nm.x
```

or

``` bash
./frieza_cart.x
```

The principal output files are:

``` text
correlation.dat
fourier.dat
power.dat
```

Their names can be changed through the corresponding `input_files` namelist variables.

## Scope

FRIEZA does not perform the electronic-structure calculation itself. It operates on pre-existing *ab initio* trajectories, equilibrium geometries, equilibrium Hessians, and, for the normal-mode TGWP calculation, a sequence of trajectory Hessians.

The normal-mode implementation is the complete Gaussian-wavepacket implementation. The Cartesian implementation is intentionally restricted to Frozen Gaussian propagation and is designed primarily for Cartesian atom-based selection of active and spectator atoms in larger systems.
