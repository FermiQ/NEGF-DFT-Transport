# TransOMat

## Description

TransOMat is a Fortran-based software package for performing electronic transport calculations, particularly for nanoscale systems. It appears to be designed to interface with output from Density Functional Theory (DFT) codes, potentially such as the Conquest code, to calculate properties like transmission and current density using methods like the Non-Equilibrium Green's Function (NEGF) formalism. The use of PETSc and SLEPc suggests it handles large-scale computations and eigenvalue problems common in such simulations.

## Dependencies

TransOMat relies on several external libraries and tools. You will need to have these installed and configured in your environment for TransOMat to compile and run correctly.

*   **Fortran Compiler:** A Fortran compiler compatible with the MPI wrapper (e.g., Intel Fortran `ifort` as used with `mpiifort`).
*   **MPI Library:** A Message Passing Interface library (e.g., Intel MPI, OpenMPI). The code is built using an MPI Fortran compiler wrapper (e.g., `mpiifort`).
*   **OpenMP:** For threaded parallelism. Compiler flags like `-qopenmp` suggest its use.
*   **PETSc (Portable, Extensible Toolkit for Scientific computation):** Used for solving large-scale linear algebra problems. Path to PETSc (e.g., `$PETSC_DIR`) needs to be configured in `Makefile.arch`.
*   **SLEPc (Scalable Library for Eigenvalue Problem Computations):** Used for solving large-scale eigenvalue problems. It's typically built on top of PETSc.
*   **SuperLU_DIST:** A library for solving large, sparse systems of linear equations using direct methods. Often used as a solver within PETSc.
*   **PEXSI (Parallel EXascale Solver Infrastructure):** Another library for solving Kohn-Sham DFT equations or other large sparse matrix problems. Path to PEXSI (e.g., `$HOME/prog/pexsi/pexsi_v2.0.0/build_make`) needs to be configured in `Makefile.arch`.
*   **Intel MKL (Math Kernel Library):** Used for optimized math routines. Linked via PEXSI and potentially other parts.

**Note on Versions:**
The example script `example/al-hf-hfo2-al/run_sakura_dft_ecc.sh` loads specific modules like `gcc/gcc-9.2.1`, `intel/intel-latest`, and `intelmpi/intelmpi-latest`. While these exact versions might not be strictly mandatory, they indicate a known working environment. Compatibility issues can arise from significantly different library or compiler versions. Ensure that the chosen Fortran compiler, MPI library, and all scientific libraries (PETSc, PEXSI, MKL) are compiled and linked consistently.

## Compilation

TransOMat is compiled using `make`. The build process relies on two Makefiles:

*   `Makefile`: The main Makefile that defines the build rules and dependencies.
*   `Makefile.arch`: An architecture-specific Makefile that is included by the main `Makefile`. This file contains crucial settings like compiler choices, paths to libraries (PETSc, PEXSI), and compiler/linker flags.

**Steps to Compile:**

1.  **Configure `Makefile.arch`:**
    *   Before compiling, you **must** review and customize `Makefile.arch` for your specific system environment.
    *   **Compiler:** Set `FC` (Fortran Compiler, e.g., `mpiifort`) and `LD` (Linker, e.g., `mpiifort`).
    *   **Library Paths:** Update `PETSC_ROOT` and `PEXSI_ROOT` to point to the correct locations of your PETSc and PEXSI installations.
    *   **Flags:** Adjust optimization flags (`OPT`), debugging flags (`DBG`), and other `FFLAGS` or `LDFLAGS` as needed for your system and compiler. The provided `Makefile.arch` often contains settings for a specific cluster or environment (e.g., using Intel compilers and MKL).
    *   Ensure that the include paths (`PETSC_INCLUDE`, `PEXSI_INCLUDE`) and library paths (`PETSC_LIB`, `PEXSI_LIB`) correctly correspond to your library installations.

2.  **Build the Executable:**
    *   Once `Makefile.arch` is configured, navigate to the root directory of the TransOMat source code.
    *   Run the `make` command:
        ```bash
        make
        ```
    *   This will compile the source files located in the `src/` directory.
    *   Object files (`.o`) will be placed in the directory specified by `OBJ_DIR` (default `./obj`).
    *   Module files (`.mod`) will be placed in the directory specified by `MOD_DIR` (default `./mod`).
    *   If the compilation is successful, the executable `transomat` will be created in the directory specified by `BIN_DIR` (default `./bin`).

3.  **Clean Build Files:**
    *   To remove object files, module files, and the executable, you can use:
        ```bash
        make clean
        ```
    *   To remove only temporary files (as defined by `TMPFILES` in `Makefile`), use:
        ```bash
        make neat
        ```

**Notes:**
*   The `Makefile` system is designed to automatically handle dependencies between source files.
*   If you encounter compilation errors, they are often related to incorrect paths in `Makefile.arch`, missing dependencies, or compiler incompatibilities.

## Input Files

TransOMat requires several input files to define the system and control the calculation parameters. The primary control file is typically named `transp.ini` (though this can be changed), which sets most of the operational parameters for a run.

### Main Control File (transp.ini)

This file contains key-value pairs that dictate the behavior of the TransOMat calculation. Parameters are generally specified in the format `$keyname=value`. Comments can typically be added using '#' at the beginning of the line or after the value. The `src/control.f90` module is responsible for parsing these parameters.

Below is a list of recognized parameters. Note that not all parameters may be required for every type of calculation; some have default values or are only used in specific modes. Refer to the example files and the descriptions for guidance.

*   `$ecc_dir=`: Directory for ECC (central region) data.
*   `$elec_left_dir=`: Directory for left electrode data.
*   `$elec_right_dir=`: Directory for right electrode data.
*   `$trans_file=`: Output filename for transmission data.
*   `$eta_cc=`: Broadening parameter (eta) for charge-charge interaction.
*   `$eps_geo=`: Geometric epsilon for convergence.
*   `$eps_real=`: Real part epsilon for convergence.
*   `$eta_elec=`: Broadening parameter (eta) for electron-electron interaction.
*   `$conv_dez=`: Convergence criterion for energy integration (dez).
*   `$k_mat_mode=`: Mode for k-point matrix generation.
*   `$e_start=`: Starting energy for calculations.
*   `$e_end=`: Ending energy for calculations.
*   `$n_steps=`: Number of energy steps for calculations.
*   `$i_start=`: Starting energy step index.
*   `$i_end=`: Ending energy step index.
*   `$nkx_dez=`: Number of k-points in the x-direction for energy integration.
*   `$nky_dez=`: Number of k-points in the y-direction for energy integration.
*   `$bias_l=`: Bias voltage on the left electrode.
*   `$bias_r=`: Bias voltage on the right electrode.
*   `$integrator=`: Type of integrator to use.
*   `$nint_order=`: Order of numerical integration.
*   `$maxsub=`: Maximum number of subdivisions for integrator.
*   `$double_contour=`: Logical flag to use double contour integration.
*   `$epsfermi=`: Broadening for Fermi function.
*   `$delta_imag=`: Small imaginary part for Green's functions (eta).
*   `$eps_int=`: Convergence criterion for integration.
*   `$eps_int_contour=`: Convergence criterion for contour integration.
*   `$elow=`: Lower energy bound for integration.
*   `$temperature=`: Temperature for Fermi distribution.
*   `$oneshot=`: Logical flag for one-shot calculation (no SCF).
*   `$current_density=`: Logical flag to calculate current density.
*   `$dftsigma=`: Logical flag to use DFT sigma for self-energies.
*   `$dftsigma_occ=`: Occupation broadening for DFT sigma.
*   `$dftsigma_virt=`: Virtual state broadening for DFT sigma.
*   `$dftsigma_imu=`: First mu index for DFT sigma matrix elements.
*   `$dftsigma_jmu=`: Second mu index for DFT sigma matrix elements.
*   `$f_file_ecc_petsc=`: Filename for Fock matrix of ECC region in PETSc format.
*   `$s_file_ecc_petsc=`: Filename for Overlap matrix of ECC region in PETSc format.
*   `$d_file_ecc_petsc=`: Filename for Density matrix of ECC region in PETSc format.
*   `$f_file_elec_l_petsc=`: Filename for Fock matrix of left electrode in PETSc format.
*   `$s_file_elec_l_petsc=`: Filename for Overlap matrix of left electrode in PETSc format.
*   `$d_file_elec_l_petsc=`: Filename for Density matrix of left electrode in PETSc format.
*   `$f_file_elec_r_petsc=`: Filename for Fock matrix of right electrode in PETSc format.
*   `$s_file_elec_r_petsc=`: Filename for Overlap matrix of right electrode in PETSc format.
*   `$d_file_elec_r_petsc=`: Filename for Density matrix of right electrode in PETSc format.
*   `$scf_conv=`: SCF convergence criterion.
*   `$calculate_ti=`: Logical flag to calculate transmission integral.
*   `$ep_loe=`: Logical flag for electron-phonon Lowest Order Expansion.
*   `$ep_active_iat1=`: First active atom index for electron-phonon calculation.
*   `$ep_active_iat2=`: Second active atom index for electron-phonon calculation.
*   `$ep_bias_min=`: Minimum bias for electron-phonon calculation.
*   `$ep_bias_max=`: Maximum bias for electron-phonon calculation.
*   `$ep_bias_n=`: Number of bias points for electron-phonon calculation.
*   `$ep_reuse_lambda=`: Logical flag to reuse lambda in electron-phonon calculation.
*   `$eta_ph_cc=`: Broadening parameter (eta) for phonon-charge interaction.
*   `$eta_ph_elec=`: Broadening parameter (eta) for phonon-electron interaction.
*   `$temperature_ph=`: Temperature for phonon calculations.
*   `$solver_mode=`: Mode for the linear solver.
*   `$nsim_rhs=`: Number of simultaneous right-hand sides for solver.
*   `$ngroups=`: Number of groups for solver.
*   `$reuse_surface_gf=`: Logical flag to load/save surface Green's functions.
*   `$no_noneq=`: Logical flag to skip non-equilibrium part of calculation.
*   `$reaktor=`: Logical flag to enable Reaktor specific features/mode.
*   `$k_on_demand=`: Logical flag for k-point on demand generation.
*   `$dump_nzs=`: Logical flag to dump non-zero structure of matrices.
*   `$diag_fix_ef=`: Logical flag to fix Fermi energy during diagonalization.
*   `$diag_fixH=`: Logical flag to fix Hamiltonian during diagonalization.
*   `$diag_fixK=`: Logical flag to fix overlap matrix during diagonalization.
*   `$diag_dim=`: Dimension for diagonalization (e.g., number of states).
*   `$diag_skip_check=`: Logical flag to skip checks during diagonalization.
*   `$diag_nkz=`: Number of kz points for diagonalization.

### Other Input Files

Depending on the calculation mode and settings in `transp.ini`, TransOMat may require other input files. These can include:

*   **System Data Files:** Files containing information about the electronic structure of the central scattering region and the electrodes. These are often outputs from a preceding DFT calculation (e.g., from Conquest). Parameters like `$ecc_dir`, `$elec_left_dir`, `$elec_right_dir` in `transp.ini` specify directories where these files are located. Specific file names for Hamiltonian, overlap, and density matrices (e.g., `$f_file_ecc_petsc`, `$s_file_ecc_petsc`) are also set in `transp.ini`.
*   **Coordinate Files:** Files like `coord.in` (often used by Conquest) might be needed to define the atomic structure if TransOMat needs to process geometrical information directly, though often this information is already incorporated into the DFT output files.
*   **Conquest Input File:** If TransOMat directly invokes or reuses detailed settings from a Conquest run, a `Conquest_input` file might be read or referenced.

The exact nature and format of these auxiliary files depend on the workflow and how TransOMat interfaces with the DFT code that generated the primary electronic structure data.

## Usage

TransOMat is parallelized using MPI and potentially OpenMP. Therefore, it needs to be executed using an MPI launcher (e.g., `mpiexec`, `mpirun`, `srun`).

**General Command:**

A typical command to run TransOMat might look like this:

```bash
mpiexec -n <num_mpi_processes> <path_to_transomat_executable>/transomat < <input_control_file> > <output_log_file>
```

Where:
*   `<num_mpi_processes>`: Number of MPI processes to use.
*   `<path_to_transomat_executable>`: Path to the compiled `transomat` executable (e.g., `./bin/transomat` if compiled in the default location).
*   `<input_control_file>`: The main input control file (e.g., `transp.ini`). TransOMat reads this file from standard input.
*   `<output_log_file>`: A file to redirect standard output and error messages, which will contain logs and results from the run.

**Environment Variables:**

As seen in the example scripts (e.g., `example/al-hf-hfo2-al/run_sakura_dft_ecc.sh`), you might need to set environment variables:
*   `OMP_NUM_THREADS=1`: If using MPI for parallelism primarily, setting OpenMP threads to 1 can be beneficial. Adjust based on your intended MPI/OpenMP hybrid strategy.
*   `MKL_NUM_THREADS=1`: Similarly for the Intel MKL library.
*   Paths to other executables: If TransOMat interacts with other codes like Conquest, ensure the path to that executable (e.g., `conquest_bin`) is correctly set if TransOMat internally calls it.

**Example Script:**

The script `example/al-hf-hfo2-al/run_sakura_dft_ecc.sh` provides a concrete example of how to run TransOMat in a cluster environment (using SGE job scheduler commands, which might need adaptation for other schedulers like Slurm or PBS).

Key aspects from the example script:
*   **Module Loading:** It loads necessary modules for compilers and libraries (e.g., `gcc`, `intel`, `intelmpi`). This step is highly environment-specific.
*   **MPI Execution:** It uses `mpiexec.hydra` with specific options (`-ppn`, `-genvall`, etc.). The exact `mpiexec` command and its options can vary depending on your MPI implementation.
*   **Input/Output:** The example script for Conquest (`ecc/Conquest`) redirects output (`> conquest.out 2>conquest.err`). A similar approach should be used for TransOMat itself.
*   **Working Directory:** Ensure that the program is run from a directory where it can find all necessary input files, or provide full paths in `transp.ini`.

You will likely need to adapt the example run scripts to your specific cluster or computing environment, particularly the module loading and MPI execution commands.

## Output Files

TransOMat generates various output files depending on the calculations performed. The primary output is typically directed to standard output and should be redirected to a log file, as shown in the [Usage](#usage) section. This log file contains detailed information about the calculation progress, parameters used, and key results.

Specific output files that may be generated include:

*   **Transmission Data File:** The input parameter `$trans_file` in `transp.ini` suggests that a file is created to store the calculated transmission coefficients as a function of energy. This file is likely a plain text file with columns for energy and corresponding transmission values.
*   **Current Density Files:** If the `$currentdensity` parameter is enabled, TransOMat may produce files containing information about the calculated current density. The format and content of these files would be specific to this calculation mode.
*   **Green's Function Files:** If `$loadsave_gf` is enabled for saving, files containing Green's functions might be written to disk for reuse in subsequent calculations.
*   **Conquest Dump Files:** Parameters like `$dump_nzs` and module names like `write_conquest_dump.mod.f90` suggest that TransOMat might generate files in a format compatible with or similar to Conquest dump files, possibly for visualization or further analysis.

The naming and exact format of these files are determined by the program's internal logic and potentially influenced by input parameters. Users should check the standard output log for messages indicating the names of files being written.

## Example

The repository includes an example calculation in the `example/al-hf-hfo2-al/` directory. This example appears to simulate a system composed of Aluminum (Al), Hafnium (Hf), Hafnium Dioxide (HfO2), and Aluminum (Al), likely representing a metal-insulator-metal junction.

**Key Files in the Example Directory:**

*   `Conquest_input`: Input file for the Conquest DFT code, likely used to generate the initial electronic structure for the system.
*   `coord.in`: Atomic coordinates file, also likely for Conquest.
*   `pseudo/`: Directory containing pseudopotential files (`.ion`) for the different atomic species involved.
*   `transp.ini`: The main TransOMat input control file for this specific example. Examining this file will show which parameters are set for this particular Al-Hf-HfO2-Al system.
*   `run_sakura_dft_ecc.sh`, `run_sakura_dft_electrodes.sh`, `run_sakura_nscf.sh`: Shell scripts for running different stages of the calculation, possibly on a cluster using a job scheduler (e.g., SGE, as indicated by `#$` directives). These scripts handle setting up the environment, running Conquest for different parts of the system (electrodes, central region), and likely then running TransOMat.
*   `scripts/`: Contains helper scripts, possibly for pre- or post-processing or managing parts of the calculation workflow.

**What the Example Demonstrates:**

This example likely demonstrates a typical workflow for TransOMat:
1.  Performing DFT calculations (using Conquest) for the electrodes and the central scattering region.
2.  Running TransOMat to calculate transport properties (e.g., transmission) based on the DFT outputs, controlled by `transp.ini`.

Users can study the `transp.ini` file in this directory to see a practical application of the input parameters and use the run scripts as a template for their own calculations, adapting them to their specific system and computational environment.

## License

TransOMat is distributed under the terms of the license specified in the `LICENSE` file in the root directory of this repository. Please review this file for detailed licensing information.

## Contributing
Contributions are welcome! Please follow the standard GitHub flow:
1. Fork the repository.
2. Create a new branch for your feature or bug fix.
3. Commit your changes.
4. Push your branch to your fork.
5. Create a Pull Request.
