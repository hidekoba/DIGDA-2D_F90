# DIGDA-2D

**DIGDA-2D** is a global sediment diagenesis model based on the *Diagenesis and Geochemistry of Deep-sea sediments for climate–carbon cycle assessment* (DIGDA) framework.

The model is designed to investigate early diagenetic processes in marine sediments and their roles in long-term ocean biogeochemical and carbon cycling.

DIGDA-2D represents independent, vertically resolved sediment columns distributed across a two-dimensional horizontal ocean grid. This configuration enables computationally efficient simulations of sediment composition, biogeochemical reactions, and burial fluxes at the global scale.

The model includes four solid components (clay, particulate organic carbon, calcium carbonate, and biogenic silica) and five dissolved pore-water components (phosphate, dissolved inorganic carbon, alkalinity, oxygen, and silicate).

The represented processes include:

- Deposition of particulate material at the seafloor
- Vertical transport of solid particles through burial and bioturbation
- Organic matter remineralization
- Dissolution of calcium carbonate and biogenic silica
- Molecular diffusion of dissolved constituents in pore water
- Exchange of dissolved constituents across the sediment–water interface

DIGDA-2D is designed for application to global ocean biogeochemical simulations, including offline experiments and coupling with ocean general circulation models (OGCMs).

## Directory Structure

The repository is organized as follows:

```text
DIGDA-2D/
├── Bound/       # Boundary forcing data
├── init/        # Initial condition generation
├── RUN/         # Simulation and post-processing scripts
├── OUT/         # Simulation output and restart files
├── analysis/    # Data analysis and visualization
├── bin/         # Compiled executables
├── lib/         # Libraries and utilities
└── src/         # Fortran source code
```

### `Bound/` — Boundary Conditions

Contains boundary forcing data used to prescribe conditions at the sediment–water interface.

DIGDA-2D requires two categories of boundary conditions.

**1. Particulate fluxes reaching the seafloor**

- Clay (lithogenic particles)
- Particulate organic carbon (POC)
- Calcium carbonate (CaCO₃)
- Biogenic silica (SiO₂)

These fluxes determine the supply of solid material to the sediment surface.

**2. Bottom-water biogeochemical properties**

- Phosphate (PO₄)
- Dissolved inorganic carbon (DIC)
- Alkalinity
- Dissolved oxygen (O₂)
- Dissolved silicate

Bottom-water concentrations prescribe the upper boundary conditions for dissolved constituents in sediment pore water.

Boundary conditions can be derived from observational datasets or ocean biogeochemical model outputs, depending on the experimental configuration.

In the observation-based experiments described in the associated manuscript, boundary forcing is constructed using satellite-derived biological production, dust deposition, and observed ocean biogeochemical properties.

**Data availability:** Large boundary-forcing NetCDF datasets are not included in the GitHub repository. Users need to prepare the required forcing files separately and configure their local paths. The data sources and any distribution links will be documented separately.

### `init/` — Initial Conditions

Contains scripts for generating the initial condition files required to initialize DIGDA-2D simulations.

These scripts prepare the initial distributions of solid and dissolved sediment components.

For sequential simulations, restart files generated during a preceding simulation segment can be used as initial conditions for the subsequent segment.

### `RUN/` — Model Execution

Contains scripts for configuring and executing model simulations, including single-period and sequential simulations.

The directory also includes post-processing scripts for combining NetCDF outputs and calculating sediment composition.

Example scripts are provided as `*_dummy.sh` files. Users should modify these scripts according to their computing environment.

### `OUT/` — Model Output

Stores simulation results, including sediment biogeochemical variables, sediment composition, and restart files.

Output files can be processed and analyzed using the scripts provided in `RUN/` and `analysis/`.

### `analysis/` — Data Analysis

Contains scripts for processing, analyzing, and visualizing model outputs.

### `bin/` — Executables

Contains compiled executables used to run DIGDA-2D simulations.

### `lib/` — Libraries

Contains libraries and utility routines required for model compilation and execution.

### `src/` — Source Code

Contains the Fortran source code implementing the model's governing equations, numerical calculations, sediment transport, and biogeochemical processes.

## Requirements

DIGDA-2D is implemented in Fortran and uses MPI for parallel computation.

The example compilation environment uses:

- ****Fortran compiler:**** Intel MPI Fortran (`mpiifort`)
- ****MPI:**** Intel MPI
- ****NetCDF-C:**** Version 4.9.2
- ****NetCDF-Fortran:**** Version 4.6.0
- ****Build tools:**** Make and Bash
- ****Additional library:**** Math library (`libm`)
- ****Post-processing:**** NCO (NetCDF Operators)

The NetCDF versions listed above correspond to the example configuration files. Other compatible versions may also work. Python and NCO are required for the initialization workflow described below; CDO is needed only for the alternative initialization script.

The example execution scripts are designed for HPC environments using the PJM job scheduler. Other environments may require modifications to the compilation and execution settings.

## Compilation

### 1. Configure the Compilation Environment

The `src/` directory contains the following example configuration files:

- `Mkinclude.dummy`: Defines the model directory, target executable, and system configuration.
- `Sysdep.dummy.par`: Defines the Fortran compiler, compilation flags, and library paths.
- `compile_dummy.sh`: Provides a compilation script with selectable grid resolutions.

Copy these files to create configurations appropriate for your computing environment.

For example:

```bash
cd src
cp Mkinclude.dummy Mkinclude
cp Sysdep.dummy.par Sysdep.local.par
cp compile_dummy.sh compile.sh
```

Modify `Mkinclude` to specify the model directory and computing environment:

```makefile
DIGDADIR = /path/to/DIGDA-2D_F90
SYSTEM = local
PARALLEL = .par
```

Modify `Sysdep.local.par` to specify the compiler and the installation paths of NetCDF-C and NetCDF-Fortran.

The Makefile automatically reads `Mkinclude` and `Sysdep.$(SYSTEM)$(PARALLEL)`.

### 2. Select the Grid Resolution

Three grid resolutions are available in the example compilation script:

| Grid resolution | Example MPI processes |
|---|---|
| 128 × 64 | 8 |
| 256 × 192 | 8 |
| 2160 × 1080 | 8 |

Edit `compile.sh` and activate the desired resolution:

```bash
resolution="128x64"
\#resolution="256x192"
\#resolution="2160x1080"
```

Only one resolution should be active at a time.

The script automatically updates the corresponding grid and MPI settings in `digda_grid.f90` and `digda_mpi.f90`.

The MPI configuration should be consistent with the number of processes used during model execution. The values above represent the example configuration and can be modified for other computing environments.

### 3. Compile the Model

Run the compilation script from the `src/` directory:

```bash
bash compile.sh
```

The script performs `make clean` followed by `make`.

Compilation messages are saved to `make.log`.

The resulting executable is placed in the repository's `bin/` directory, with a filename indicating the grid resolution and MPI configuration.

For example, a 2160 × 1080 grid compiled for eight MPI processes produces:

```text
bin/digda_2160x1080.8.par
```

****Note:**** The compilation script directly modifies `digda_grid.f90` and `digda_mpi.f90`. Users should check these changes when switching between grid resolutions, particularly when managing the source code with Git.

## Model Execution

DIGDA-2D is designed to run on high-performance computing (HPC) systems using MPI.

The `RUN/` directory contains example scripts for executing model simulations and post-processing model outputs.

### 1. Execution Scripts

The following example scripts are provided:

| Script | Description |
|---|---|
| `run_digda_1p4d_dummy.sh` | Runs a single simulation period |
| `run_digda_chain_1p4d_dummy.sh` | Performs sequential simulations using restart files |
| `pp_digda_combine_SEDTRR_dummy.sh` | Combines sediment NetCDF output files |
| `pp_digda_calc_weight_percentage_dummy.sh` | Calculates sediment component weight percentages |

These scripts serve as templates and should be adapted to the user's computing environment.

### 2. Configure the Computing Environment

The example execution scripts are designed for HPC systems using the PJM job scheduler.

Users should configure the following settings:

- Model installation directory and executable path
- Boundary forcing directory
- Initial condition and restart directories
- Simulation output directory
- Job scheduler directives and resource allocation
- MPI execution command and number of processes
- Environment modules and library paths

The scripts generate a model parameter file named `PARAMET_DIGDA`, which specifies simulation settings and biogeochemical parameters using Fortran namelist input. Optional carbon-isotope calculations are controlled separately by the `calc13` and `calc14` flags in the `nmsdiso` namelist.

### 3. Running a Single Simulation

Copy the example script:

```bash
cd RUN
cp run_digda_1p4d_dummy.sh run_digda_1p4d.sh
```

Modify the script to specify the computing environment, simulation period, input directories, and model parameters.

For a properly configured PJM environment, execute:

```bash
bash run_digda_1p4d.sh
```

The script generates `PARAMET_DIGDA`, submits the job to the scheduler, monitors its status, and invokes post-processing after job completion.

### 4. Running Sequential Simulations

For long simulations, DIGDA-2D supports sequential execution using restart files.

Copy the example script:

```bash
cp run_digda_chain_1p4d_dummy.sh run_digda_chain_1p4d.sh
```

The simulation period can be configured using:

```bash
SDYRSTR="0"
SDYREND="150000"
SDYRSTP="30000"
SDYRHST="10000"
```

Here:

- `SDYRSTR`: Starting model year
- `SDYREND`: Final model year
- `SDYRSTP`: Length of each simulation segment
- `SDYRHST`: History output interval

In this example, the simulation is divided into five consecutive segments of 30,000 years each.

The restart files generated at the end of each segment are used as initial conditions for the next segment.

Sequential execution is particularly useful for long-term sediment simulations requiring extended spin-up periods.

### 5. Post-processing

The execution scripts invoke post-processing utilities after model integration.

****Combining sediment outputs****

`pp_digda_combine_SEDTRR_dummy.sh` uses NCO (`ncks`) to combine sediment output variables into a NetCDF file named `SEDTRR.nc`.

****Calculating sediment composition****

`pp_digda_calc_weight_percentage_dummy.sh` calculates the weight percentages of four sediment components:

| Variable | Sediment component |
|---|---|
| `PCLA` | Clay |
| `PPOM` | Particulate organic matter |
| `PCAR` | Calcium carbonate |
| `POPA` | Biogenic silica (opal) |

The calculated variables are saved as separate NetCDF files for subsequent analysis and visualization.

****Note:**** The execution and post-processing scripts are system-dependent examples. They may require modification for other HPC systems or local computers.

## Input and Output Data

### Boundary Forcing Data

Boundary forcing files are stored in `Bound/`.

The forcing consists of particulate fluxes reaching the seafloor and bottom-water biogeochemical properties.

For observation-based experiments, the associated manuscript describes the use of satellite-derived biological production, observational dust deposition datasets, and climatological ocean biogeochemical fields.

For OGCM-driven experiments, forcing fields can be obtained from ocean model output.

The required input file names, variable names, and file formats depend on the experimental configuration. Large NetCDF forcing files are intentionally excluded from version control; their provision and redistribution must follow the original data providers’ terms.

### Initial Conditions

The `init/` directory provides Python scripts and coordinate tables to create a NetCDF restart template, together with example shell scripts for assigning initial pore-water concentrations. Initialization must use the same horizontal grid and 10-layer sediment vertical grid as the compiled model.

**Supported grid configurations**

| Horizontal grid | Sediment layers | Python template script |
|---|---:|---|
| 128 × 64 | 10 | `mkinit_template_128064.py` |
| 256 × 192 | 10 | `mkinit_template_256192.py` |
| 2160 × 1080 | 10 | `mkinit_template_21601080.py` |

The template scripts use the corresponding `lon_values_*`, `lat_values_*`, and `dep_values_*` coordinate files to construct a ten-layer NetCDF restart template. Both `standard` and `ciso` modes include the same ten standard fields. The `ciso` mode additionally includes six carbon-isotope variables.

**Standard initial conditions (both `standard` and `ciso` modes)**

| NetCDF variable | Sediment component | Template initial value | Unit |
|---|---|---:|---|
| `Clay` | Clay | 2,600,000 | g m⁻³ |
| `OrgM` | Particulate organic matter | 0 | mmol m⁻³ |
| `CaCO3` | Calcium carbonate | 0 | mmol m⁻³ |
| `SiO2` | Biogenic silica | 0 | mmol m⁻³ |
| `Agetracer` | Age tracer | 0 | mmol m⁻³ (template attribute) |
| `Phosphate` | Dissolved phosphate | 0 | mmol m⁻³ |
| `DIC` | Dissolved inorganic carbon | 0 | mmol m⁻³ |
| `Alkalinity` | Alkalinity | 0 | mmol m⁻³ |
| `Silicate` | Dissolved silicate | 0 | mmol m⁻³ |
| `Oxygen` | Dissolved oxygen | 0 | mmol m⁻³ |

These are **template values**, not necessarily the final values used to start the simulation. The example `mkinit_template_dummy.sh` workflow replaces the five dissolved fields (`Phosphate`, `DIC`, `Alkalinity`, `Silicate`, and `Oxygen`) with the corresponding bottom-water forcing fields. The solid-component initial values remain as specified above.

A separate 128 × 64 isotope-inclusive example supplied during development initializes these five dissolved fields to 2.5, 2300, 2300, 150, and 150 mmol m⁻³, respectively. Those values belong to that example and are **not** the fixed defaults of the published dummy workflow, which uses bottom-water forcing.

**Standard and carbon-isotope (`ciso`) configurations**

The optional `ciso` configuration adds six carbon-isotope restart variables, in addition to the standard sediment and pore-water fields:

| Carbon-isotope variable | Example initial value | Unit |
|---|---:|---|
| `C13-OrgM` | 0 | mmol m⁻³ |
| `C13-CaCO3` | 0 | mmol m⁻³ |
| `C13-DIC` | 25 | mmol m⁻³ |
| `C14-OrgM` | 0 | mmol m⁻³ |
| `C14-CaCO3` | 0 | mmol m⁻³ |
| `C14-DIC` | 2.6 × 10⁻⁹ | mmol m⁻³ |

These isotope values are example initial values from the supplied initialization template and are not prescribed universally for every experiment. Enabling isotope variables in the restart file does **not** by itself enable isotope calculations in the model: the run-time `calc13` and `calc14` flags must also be configured in `PARAMET_DIGDA`.

**Generating the initial file**

The public example script `mkinit_template_dummy.sh` is designed to select the horizontal resolution and either `standard` or `ciso` initialization. When using that script, run it from the repository's `init/` directory (after adapting its paths to your environment), for example:

```bash
cd init
bash mkinit_template_dummy.sh 128x64 standard
bash mkinit_template_dummy.sh 128x64 ciso
# Alternative resolutions: 256x192 or 2160x1080
```

The initialization step requires bottom-water forcing files with the following variables:

| File | NetCDF variable | Initialized field |
|---|---|---|
| `bnit.nc` | `NIT` | `Phosphate` |
| `bdic.nc` | `DIC` | `DIC` |
| `balk.nc` | `ALK` | `Alkalinity` |
| `bsil.nc` | `SIL` | `Silicate` |
| `bdo2.nc` | `DO2` | `Oxygen` |

The supplied scripts and input datasets are computing-environment dependent. Confirm that file locations, dimensions, horizontal coordinates, and variable units match before using generated restart files. The original `mkinit_template.sh` contains the 256 × 192 and 2160 × 1080 initialization procedures; 128 × 64 needs the corresponding forcing data and a matching initialization procedure. The `mkinit_modified.sh` script is an alternative approach based on an existing model output and also uses CDO.

**Note:** The `ciso` workflow and its restart-field names should be validated against the Fortran restart reader before scientific production runs. Large generated NetCDF restart files should be kept outside Git version control.

### Model Outputs

Model results are stored in `OUT/`.

Outputs include sediment biogeochemical variables, sediment composition, and restart files.

The post-processing scripts can generate additional NetCDF files containing sediment component weight percentages.

Detailed descriptions of the input and output variables, file formats, and naming conventions will be provided separately.

## References and Citation

DIGDA-2D is based on the model described in the following manuscript:

Kobayashi, H., and Oka, A., *Global simulation of early diagenesis in marine sediments using the Diagenesis and Geochemistry of Deep-sea sediments for climate–carbon cycle assessment (DIGDA) model*. Manuscript in preparation.

The manuscript describes the model formulation, numerical configuration, evaluation against observational constraints, and applications to global sediment biogeochemical simulations.

Users of DIGDA-2D are encouraged to acknowledge the model and cite the associated scientific publication once it becomes available.

## License

License information will be provided separately.
