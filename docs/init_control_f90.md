# `init_control.f90` - Initialization of Control Parameters

## Overview

This subroutine is responsible for initializing the runtime configuration of the TransOMat simulation. It achieves this by reading key-value parameters from an external input file, typically named `transp.ini`. It uses the definitions and parsing routines from the `control` module to populate global variables (defined in `globals_mod`) that steer the behavior of the entire application. This includes setting up file paths, energy ranges, convergence criteria, calculation modes, and various other operational flags.

## Key Components

*   **`subroutine init_control()`**: The single public subroutine in this file.
    *   **Initialization:**
        *   Declares a local variable `cntrl` of type `cntrl_` (defined in the `control` module).
        *   Calls `init_cntrl(cntrl)` (from the `control` module) to populate the `kname` field (the string to search for in `transp.ini`, e.g., `$ecc_dir=`) for every parameter within the local `cntrl` structure.
    *   **File Handling:**
        *   Opens the input file `./transp.ini` for reading, associating it with the file unit `iunit_control` (a global variable).
    *   **Parameter Reading:**
        *   Sequentially calls the `get_cntrl_key` interface (which resolves to `get_skey`, `get_rkey`, `get_ikey`, `get_lkey` from the `control` module) for each expected parameter.
        *   The `kname` from the local `cntrl` structure is used to identify the parameter in `transp.ini`.
        *   The parsed value is stored directly into a corresponding global variable (e.g., `ecc_dir`, `estart`, `nkx_dez`, `oneshot`). Many of these global variables are defined in `globals.f90`.
    *   **Post-processing and Defaults:**
        *   Performs unit conversions for bias voltages (e.g., `vl = vl / eh`). The constant `eh` (Hartree energy) is expected to be globally available.
        *   Calculates `vb` (total bias) from `vr` and `vl`.
        *   Sets default values for some logical flags (e.g., `ldouble_contour = .false.`) before attempting to read them. This makes these parameters optional in `transp.ini`.
        *   Conditionally reads certain groups of parameters based on the value of other flags (e.g., DFT sigma parameters are only read if `dftsigma` is true; electron-phonon parameters if `l_ep_loe` is true).
    *   **File Closing:**
        *   Closes the `transp.ini` file.
    *   **Output:**
        *   If it's the master MPI process (`inode .eq. 0`), prints a separator line to standard output. The actual printing of each key-value pair as it's read is handled by the `get_key` subroutine in `control.f90`.

## Important Variables/Constants

This subroutine primarily populates global variables. Key global variables set here include:

*   **File/Directory Paths:** `ecc_dir`, `elec_l_dir`, `elec_r_dir`, `trans_file`.
*   **Energy & k-points:** `estart`, `eend`, `n_energy_steps`, `i_energy_start`, `i_energy_end`, `nkx_dez`, `nky_dez`.
*   **Voltages & Temperature:** `vl`, `vr`, `vb`, `temperature_el`, `temperature_ph`.
*   **Calculation Mode & Control Flags:** `k_mat_mode`, `oneshot`, `calc_current_density`, `dftsigma`, `l_ep_loe`, `l_loadsavegf`, `l_k_on_demand`, `solver_mode`.
*   **Convergence & Broadening:** `eta_cc`, `eta_elec`, `eps_geo`, `eps_real`, `conv_dez`, `scf_conv`, `delta_imag`, `epsfermi`.
*   **Algorithm Specifics:** `integrator`, `nint_order_el`, `maxsub`, `ldouble_contour`.
*   Many others as listed in `globals.f90` and defined in `type cntrl_` in `control.f90`.

Local variable:
*   **`cntrl` (type `cntrl_`):** A local instance used to hold the parameter key names. It is not used to store the values themselves after parsing; values go into global variables.

External Constants:
*   **`eh` (Real(dp)):** Expected to be a globally defined constant representing Hartree energy, used for voltage unit conversion.

## Usage Examples

This subroutine is called once at the beginning of a TransOMat run, typically from the main program `transomat.f90`, before any significant computations begin.

```fortran
program transomat
  use globals      ! Provides global variables
  use init_control_module ! Assuming init_control is in a module
  implicit none
  ! ...
  call init_control() ! Read and set all control parameters
  ! ... proceed with calculations using the now-populated global parameters ...
end program transomat
```
It requires `transp.ini` to be present in the execution directory.

## Dependencies and Interactions

*   **Internal Dependencies:**
    *   `USE globals`: Essential for accessing and populating the global variables that store the control parameters (e.g., `iunit_control`, `ecc_dir`, `estart`, `vl`, `eh`, `inode`).
    *   `USE control`: This is the primary dependency, used for:
        *   `type cntrl_`: Definition of the structure holding parameter key names.
        *   `subroutine init_cntrl()`: To initialize the key names in the local `cntrl` variable.
        *   `interface get_cntrl_key` (and its specific implementations like `get_skey`, `get_rkey`): To parse values from the input file.
*   **External Library Dependencies:**
    *   None directly within this file, but the `control` module it uses has a dependency on PETSc for `PetscPrintf`.
*   **Interactions:**
    *   **Reads from `./transp.ini`:** This is the primary input.
    *   **Writes to Global State:** Modifies numerous global variables defined in `globals.f90`.
    *   **Called by:** The main program `transomat.f90`.
    *   **Output:** Prints each parsed key-value pair to standard output (via `get_key` in `control.f90`) and a final separator line from the master process.

---
*This documentation was auto-generated by Jules.*
