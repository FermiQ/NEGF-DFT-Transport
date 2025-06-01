# `transomat.f90` - Main Program for TransOMat

## Overview

This file contains the main program for the TransOMat simulation tool. It orchestrates the overall calculation flow, from initialization of parallel environments (MPI) and scientific libraries (PETSc/SLEPc) to reading control parameters, setting up the system, and invoking various computational modules. Depending on the input parameters, it can perform Self-Consistent Field (SCF) calculations, determine transmission properties, calculate Lowest Order Expansion (LOE) for electron-phonon coupling, or manage K-matrix operations. It also handles logging, timing, and finalization of the simulation.

## Key Components

*   **`program transomat`**: The main entry point of the TransOMat executable.
    *   **Initialization Phase:**
        *   `MPI_INIT_THREAD`: Initializes the Message Passing Interface (MPI) environment.
        *   `petsc_init()`: Initializes the Portable, Extensible Toolkit for Scientific Computation (PETSc).
        *   `init_control()`: Initializes and reads control parameters (likely from `transp.ini`). See `control.f90` and `init_control.f90`.
        *   `petsc_subcomm_init()`: Initializes PETSc sub-communicators.
        *   `petsc_solver_init()`: Initializes PETSc solvers.
        *   `init_timings()`: Sets up timers for performance tracking.
        *   `init_sys()`: Initializes the physical system based on input parameters. See `init_sys_mod.f90`.
    *   **Convergence Check & SCF Loop Control:**
        *   Reads `convergence.negf` file (if it exists) to check previous convergence status and resume calculations if possible.
        *   Manages `iter` (iteration count) and `conv` (convergence value).
        *   Broadcasts initialization and convergence status to all MPI processes.
    *   **Main Computational Logic (based on `k_mat_mode`):**
        *   **`k_mat_mode = 1` (Standard transport/SCF):**
            *   If calculation is `oneshot` or already converged:
                *   `transmission()` or `transmission_k_on_demand()`: Calculates transmission. See `transmission_mod.f90`.
                *   Optionally, `scf()`: Performs SCF if current density calculation is enabled. See `scf_mod.f90`.
            *   Else if Electron-Phonon (`l_ep`) is enabled:
                *   `loe_run()`: Runs LOE calculation. See `loe_mod.f90`.
            *   Else (not converged, not oneshot, not LOE):
                *   `scf()`: Performs SCF calculation. See `scf_mod.f90`.
        *   **`k_mat_mode = 2` (K-matrix generation/dumping):**
            *   If converged, writes `negf.converged`.
            *   Else, `get_kmat_from_diag()` and `dump_kmat()`: Generates K-matrix from diagonal part and dumps it. See `kmat_from_diag.f90` and `kmat_mod.f90`.
        *   **`k_mat_mode = 3` (Bond order):**
            *   `bond_order()`: Calls bond order calculation routine (likely in `analysis_mod.f90` or a dedicated module).
    *   **Cleanup and Finalization:**
        *   `petsc_cleanup()`: Cleans up PETSc objects.
        *   `slepcfinalize()`: Finalizes the Scalable Library for Eigenvalue Problem Computations (SLEPc).
        *   `MPI_FINALIZE()`: Finalizes the MPI environment.

## Important Variables/Constants

*   **`iter` (Integer):** Current iteration number, mainly for SCF loops.
*   **`conv` (Real(dp)):** Current convergence value (e.g., difference in density matrix between SCF iterations).
*   **`init_scf` (Logical):** Flag indicating if this is an initial SCF calculation (true) or a continuation (false).
*   **`converged` (Logical):** Flag indicating whether the main calculation (usually SCF) has converged.
*   **`k_mat_mode` (Integer):** Global parameter (from `globals.f90`, set via `init_control`) that dictates the main operational mode of the program.
*   **`oneshot` (Logical):** Global parameter; if true, performs a single calculation pass without an SCF loop.
*   **`l_ep` (Logical):** Global parameter; if true, enables electron-phonon calculations using LOE.
*   **`calc_current_density` (Logical):** Global parameter; if true, enables current density calculation.
*   **`pstr_out` (from `petsc_mod`):** PETSc string for standard output, ensuring output is typically handled by the master MPI process.

## Usage Examples

The program `transomat` is the main executable. It is typically run via an MPI executor:

```bash
mpiexec -n <num_processes> ./transomat < transp.ini > output.log
```
The behavior is controlled by parameters in `transp.ini`.

## Dependencies and Interactions

*   **Internal Dependencies:**
    *   `USE mpi`: Essential for parallel execution.
    *   `USE kinds`: For defining data types (e.g., `dp`).
    *   `USE misc`: For miscellaneous utility functions (e.g., `strln`).
    *   `USE globals`: For global parameters and variables (e.g., `k_mat_mode`, `l_ionode`, `inode`).
    *   `USE init_sys_mod`: For system setup (`init_sys`).
    *   `USE petsc_control`: For PETSc related control variables.
    *   `USE petsc_mod`: For PETSc utilities (`petsc_init`, `petsc_print_master`, `petsc_cleanup`).
    *   `USE scf_mod`: For SCF calculations (`scf`).
    *   `USE loe_mod`: For LOE calculations (`loe_run`).
    *   `USE phonon_mod`: (Currently appears unused directly in the main logic flow shown, but linked).
    *   `USE transmission_mod`: For transmission calculations (`transmission`, `transmission_k_on_demand`).
    *   `USE kmat_from_diag`: For `get_kmat_from_diag`.
    *   `USE kmat_mod`: For `dump_kmat`.
    *   `USE analysis`: (Linked, specific function `bond_order` likely comes from here or a related module).
    *   `USE timing`: For `init_timings`.
    *   Implicitly uses `init_control.f90` (via `call init_control()`) and `control.f90`.
*   **External Library Dependencies:**
    *   **MPI:** For parallel processing.
    *   **PETSc:** For linear algebra, solvers, and parallel utilities.
    *   **SLEPc:** For eigenvalue problems (used in finalization).
*   **Interactions:**
    *   Reads control parameters from a file (typically `transp.ini`) via `init_control`.
    *   Reads system information from files specified in `transp.ini` during `init_sys`.
    *   Writes output logs, convergence status (`convergence.negf`, `negf.converged`), and results of calculations (e.g., transmission data) to files.
    *   The flow of execution heavily depends on the `k_mat_mode` parameter.

---
*This documentation was auto-generated by Jules.*
