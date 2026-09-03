# FRIEZA

**FRIEZA** is a Fortran program for computing Gaussian wavepacket
autocorrelation functions and power spectra from *ab initio* molecular
dynamics trajectories. It supports Frozen Gaussian (FGWP),
Single-Hessian Gaussian (SH-GWP), and Heller Thawed Gaussian (TGWP)
propagation.

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

The program supports three Gaussian propagation schemes:

-   **Frozen Gaussian (FG)**: the Gaussian width is kept fixed at its
    equilibrium harmonic value.
-   **Single-Hessian Gaussian (SH)**: the Gaussian width is propagated
    using the equilibrium Hessian throughout the trajectory.
-   **Thawed Gaussian (TG)**: the Gaussian width is propagated using a
    Hessian supplied for every trajectory snapshot.

The correlation function is Fourier transformed with FFTW3. The final
`power.dat` file contains the positive-frequency power spectrum,
`|FFT[C(t)]|^2`, with frequency also reported in cm\^-1.

## Files

The source is divided into four Fortran components:

  ---------------------------------------------------------------------
  File                               Purpose
  ---------------------------------- ----------------------------------
  `main.f90`  Main program: input parsing, normal-mode transformation, Gaussian propagation, correlation function, and FFT.

  `correlation_mod.f90`  Gaussian overlap/correlation-function evaluation and complex matrix inversion using LAPACK.

  `constants.f90`  Atomic-unit conversion constants and numerical constants.

  `diagonalizer.f` Symmetric-matrix diagonalization used to obtain the equilibrium
normal modes.

  `Makefile`  Compilation and linking rules.
  ---------------------------------------------------------------------

The executable produced by the Makefile is `frieza.x`.

## Requirements

A Fortran compiler and the following numerical libraries are required:

-   `gfortran`
-   FFTW3, including the Fortran 2003 interface (`fftw3.f03`) and library `libfftw3`
-   BLAS
-   LAPACK

The supplied Makefile links with:

``` text
-lfftw3 -lm -llapack -lblas
```

and expects the FFTW include file to be available through `/usr/include`.

## Compilation

The intended build command is:

``` bash
make compile
```

The executable will be written as:

``` text
frieza.x
```

To remove object files, modules, and executables, the Makefile currently defines the target `.clean`:

``` bash
make .clean
```

## Running the program

The program reads a file named exactly:

``` text
input
```

and is run as:

``` bash
./frieza.x
```

A `banner.txt` file is also expected by the program because it is printed at startup.

The `input` file contains three Fortran namelists:

1.  `input_files`
2.  `options`
3.  `trajectory`

A complete example input used for an H2O calculation is:

``` fortran
&input_files
    geom_eq="h2o_opt.xyz"
    traj="parsed_log_traj.xyz"
    velocity="velocity.xyz"
    hess="final_py.out"
    output="input.dat"
/
&options
    calculation = 0
    eta = 0.d0 !3.4d-5
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

Not every namelist variable needs to be specified explicitly. Before reading the namelists, the program defines defaults for several file names and options. In particular, the example above relies on the default equilibrium-Hessian filename `Hessian_flat.out`, the default Fourier-output filename `fourier.dat`, the default power-spectrum filename `power.dat`, and the default `rotrasl = 6`.

For this three-atom nonlinear H2O example,

``` text
nvib = 3 * nat - rotrasl = 9 - 6 = 3
```

so three mask entries are supplied. The first vibrational normal mode is frozen at its initial coordinate and momentum (`mask(1) = 0`), while modes 2 and 3 follow the ab initio trajectory (`mask(2) = mask(3) = 1`).

The example uses `steps = 2500` trajectory points with a time step of `dt = 8.2682749151502` atomic units of time and adds `170000` zero-valued points to the FFT array. Thus,

``` text
fftsteps = steps + padding = 172500
```

The calculation is a Frozen Gaussian calculation (`calculation = 0`) with no exponential damping (`eta = 0`). The commented value `3.4d-5` shows an alternative damping parameter that can be activated by replacing `0.d0`.

## Namelist values in the supplied example

  ------------------------------------------------------------------------------
  Namelist        Variable        Value                   Meaning
  --------------- --------------- ----------------------- ----------------------
  `input_files`   `geom_eq`       `h2o_opt.xyz`           Equilibrium geometry

  `input_files`   `traj`          `parsed_log_traj.xyz`   Ab initio trajectory

  `input_files`   `velocity`      `velocity.xyz`          Initial velocities

  `input_files`   `hess`          `final_py.out`          Trajectory Hessian  file

  `input_files`   `output`        `input.dat`             Correlation-function output

  `options`       `calculation`   `0`                     Frozen Gaussian

  `options`       `eta`           `0.d0`                  No damping

  `trajectory`    `steps`         `2500`                  Number of trajectory  steps

  `trajectory`    `dt`            `8.2682749151502d0`     Time step in atomic units

  `trajectory`    `padding`       `170000`                Number of additional  FFT points

  `trajectory`    `mask`          `0, 1, 1`               Freeze mode 1;  propagate modes 2 and 3
  ------------------------------------------------------------------------------

Important values omitted from the example input are taken from the
defaults defined in `main.f90`, including:

``` text
hess_eq    = "Hessian_flat.out"
fourierout = "fourier.dat"
powerout   = "power.dat"
rotrasl    = 6
```

Because `calculation = 0`, the trajectory Hessian file is not used for
Gaussian-width propagation in this run even though `hess="final_py.out"`
is specified. It becomes relevant for the Thawed Gaussian calculation.

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

Atomic masses are assigned internally from the element symbol. The currently supported species are:

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

The program reconstructs the full symmetric matrix from this lower-triangular sequence and mass-weights it according to

``` text
H_mw(i,j) = H(i,j) / sqrt(m_i m_j)
```

The resulting matrix is diagonalized to obtain the equilibrium normal modes and harmonic frequencies.

### 3. Trajectory: `traj`

Each trajectory snapshot is expected to contain two header records followed by one line per atom. Each atom line contains the symbol, Cartesian coordinates, and Cartesian velocities:

``` text
header line 1
header line 2
Atom1   x1   y1   z1   vx1   vy1   vz1
Atom2   x2   y2   z2   vx2   vy2   vz2
...
```

Coordinates are interpreted as **Angstrom** and converted internally to atomic units.

The velocities are used as supplied by the code and are therefore expected to already be in the velocity units used by the dynamics, i.e. atomic units. No velocity conversion is currently applied.

### 4. Initial velocities: `velocity`

The initial velocity file has the same two-line header convention, followed by one line per atom containing the atomic symbol and three velocity components:

``` text
header line 1
header line 2
Atom1   vx1   vy1   vz1
Atom2   vx2   vy2   vz2
...
```

The velocities are projected onto the equilibrium normal modes to construct the initial Gaussian momentum.

### 5. Trajectory Hessians: `hess`

This file is only read in **Thawed Gaussian** mode (`calculation = 2` or any value greater than 1).

For each trajectory step, the program expects two header records followed by the lower triangle of an `ncart x ncart` Cartesian Hessian in the same flattened format as the equilibrium Hessian:

``` text
header line 1
header line 2
H(1,1)
H(1,2)
H(2,2)
...
```

Each Hessian is mass-weighted and transformed into the equilibrium normal-mode basis before being used to propagate the Gaussian width matrix.

## Calculation modes

The `calculation` variable selects the Gaussian approximation:

     `calculation` Method                    Hessian used for width propagation
  ---------------- ------------------------- ------------------------------------
               `0` Frozen Gaussian           Fixed equilibrium width `A0`
               `1` Single Hessian Gaussian   Equilibrium Hessian
    `2` or greater Thawed Gaussian           Instantaneous Hessian from `hess`

The initial Gaussian width is diagonal in the equilibrium normal-mode basis and is constructed from the harmonic frequencies:

``` text
A0(i,i) = sqrt(lambda_i)
```

where `lambda_i` are the eigenvalues of the mass-weighted equilibrium Hessian. The harmonic frequencies are `sqrt(lambda_i)` for positive eigenvalues.

## Mode masking

`mask` controls which vibrational normal modes follow the trajectory.

For a mode with

``` text
mask(i) = 1
```

the instantaneous trajectory coordinate and momentum are used. For

``` text
mask(i) = 0
```

the corresponding coordinate and momentum are kept at their initial values.

The mask therefore provides a simple way to restrict the dynamics to a selected subset of normal modes while leaving the remaining modes at their initial values.

## Damping

The parameter

``` text
eta
```

introduces an exponential damping factor in the stored correlation function:

``` text
C(t) -> C(t) exp(-eta t)
```

with `t` in atomic time units. Setting `eta = 0` disables damping.

## Time and zero padding

The trajectory contains `steps` propagated snapshots, with time increment `dt`:

``` text
t_k = k * dt
```

The FFT length is

``` text
fftsteps = steps + padding
```

`padding` is intended to provide zero padding in the Fourier transform. The FFT frequency grid is constructed for both positive and negative frequencies; `power.dat` writes only the non-negative branch.

## Output files

### Correlation-function output

The autocorrelation function is written to the filename specified by `output`. In the supplied example this is:

``` text
input.dat
```

The program default is `correlation.dat`.

The file contains:

``` text
# Time          Re[C(t)]            Im[C(t)]            |C(t)|**2
```

The columns are:

  Column   Meaning
  -------- ---------------------------------------------
  1        Time
  2        Real part of the correlation function
  3        Imaginary part of the correlation function
  4        Squared modulus of the correlation function

### `fourier.dat`

The complex Fourier transform is written as:

``` text
# Freq          Ang Freq            Re[FFT]             Im[FFT]
```

where frequency is reported in inverse atomic time units and angular frequency is reported in atomic units.

### `power.dat`

The final power spectrum is written as:

``` text
# Freq          Ang Freq [au]       Ang Freq [cm**-1]   |FFT|**2
```

The four columns are:

  Column   Meaning
  -------- -------------------------------------------------
  1        Ordinary frequency in inverse atomic time units
  2        Angular frequency in atomic units
  3        Angular frequency converted to cm$^{-1}$
  4        Power spectrum, `|FFT|^2`

Only frequencies from zero through the Nyquist frequency are written to `power.dat`.

## Method overview

The calculation proceeds as follows:

1.  Read the equilibrium geometry and assign atomic masses.
2.  Read and mass-weight the equilibrium Cartesian Hessian.
3.  Diagonalize the equilibrium Hessian to obtain harmonic normal modes.
4.  Transform the initial geometry and velocity into mass-weighted  normal-mode coordinates.
5.  Construct the initial Gaussian width matrix from the harmonic frequencies.
6.  For each trajectory snapshot, transform the Cartesian coordinates and velocities into the equilibrium normal-mode basis.
7.  Propagate the Gaussian width matrix according to the selected Frozen, Single-Hessian, or Thawed scheme.
8.  Evaluate the Gaussian overlap correlation function.
9.  Apply the optional exponential damping factor.
10. Fourier transform the correlation function with FFTW3.
11. Write the complex Fourier transform and the positive-frequency power spectrum.

## Important implementation notes

The code uses the equilibrium normal-mode transformation throughout. Instantaneous Hessians in the thawed calculation are transformed into this fixed normal-mode basis; the normal modes themselves are not recomputed at every trajectory point.
## Units

The internal dynamics use atomic units. The main external coordinate convention is:

-   geometry and trajectory coordinates: Angstrom on input;
-   velocities: atomic units as supplied;
-   Hessians: values are used as atomic-unit Hessian elements and then mass-weighted;
-   frequencies: converted to cm$^{-1}$ for printed spectral data.

The program prints the harmonic frequencies and harmonic zero-point energy before processing the trajectory, which can be used as a basic sanity check on the equilibrium Hessian and normal-mode transformation.

## Example workflow

A typical calculation directory can contain:

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
├── main.f90
└── Makefile
```

Then:

``` bash
make compile
./frieza.x
```

For the supplied example input, the principal results are:

``` text
input.dat
fourier.dat
power.dat
```

Here `input.dat` is the correlation-function output because the namelist explicitly sets `output="input.dat"`. The other two names are inherited from the program defaults because `fourierout` and `powerout` are not specified in the input file.

## Scope

This code does not perform the electronic-structure calculation itself. It operates on a pre-existing ab initio trajectory, an equilibrium geometry, an equilibrium Hessian, and, for the thawed calculation, a sequence of trajectory Hessians.
