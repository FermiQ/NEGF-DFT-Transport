# `globals.f90` - Global Variables and Constants

## Overview

This module serves as a central repository for global variables, constants, and PETSc matrix/solver type definitions used throughout the TransOMat application. It employs `USE kinds` for consistent data type definitions (e.g., `dp` for double precision) and `USE petscmat` along with an `#include <petsc/finclude/petscmat.h>` directive for PETSc functionalities. Many of these variables are populated by the `init_control` subroutine based on the `transp.ini` input file and then accessed by various computational modules.

## Key Components

*   **`module globals`**: The encompassing module for all global definitions.

## Important Variables/Constants

This module is primarily composed of variable declarations. They are grouped here by their apparent purpose:

**Physical Constants:**
*   `zione` (Complex(dp), Parameter): Complex number `(0.0, 1.0)`.
*   `pi` (Real(dp)): Mathematical constant &pi;.
*   `eh` (Real(dp)): Hartree energy in eV (27.211396132). Used for energy unit conversions.
*   `kb` (Real(dp)): Boltzmann constant in Hartree/K.

**System & MPI Flags:**
*   `ksym` (Logical): `.true.` if the system has k-space symmetry.
*   `l_ionode` (Logical): `.true.` if the current MPI process is designated for I/O operations.
*   `l_is_big_endian` (Logical): Flag for system endianness.
*   `inode` (Integer): MPI rank of the current process.
*   `nprocs` (Integer): Total number of MPI processes.
*   `psubcomm`, `inode_sub`, `inode_group`, `ngroups`, `group_range(2)`: Variables related to MPI communicators and process grouping.
*   `nodes_group(:)` (Integer, Allocatable): MPI nodes within a group.

**PETSc & Solver Settings:**
*   `solver_mode` (Integer): Selects the PETSc solver mode.
*   `nsim_rhs` (Integer): Number of simultaneous right-hand sides for PETSc solvers.
*   `mattype_surf`, `mattype_cc`, etc. (MatType): PETSc matrix type definitions.
*   `matsolvertype_surf`, `matsolvertype_cc` (MatSolverType): PETSc solver type definitions.
*   `cols_loc(:)`, `nzrow_loc(:)` (Integer, Allocatable): Local PETSc matrix column indexing and non-zero counts.
*   `nsim_inv` (Integer): Number of simultaneous inversions.

**Calculation Control & Flags (many set by `init_control`):**
*   `iigc` (Integer): Global integer, purpose not immediately clear from context (e.g. "Green's function calculation index").
*   `l_loadsavegf` (Logical): If true, load/save Green's functions.
*   `l_k_on_demand` (Logical): If true, k-points are generated as needed.
*   `integrator`, `nint_order_el`, `maxsub`, `contour_select`, `ineq`, `int_counter`: Parameters controlling numerical integration.
*   `epsfermi`, `delta_imag`, `eps_int`, `elow`, `eps_int_contour`: Tolerances and broadening factors for integration.
*   `ldouble_contour`, `lnoneq_int`, `l_output_progress`, `l_reaktor`, `l_no_noneq`: Logical flags for integration and operational modes.
*   `vb`, `mul`, `mur`, `temperature_el`, `ef_l`, `ef_c`, `ef_r`, `vl`, `vr`: Bias voltages, chemical potentials, temperature, Fermi levels.
*   `eta_ph_cc`, `eta_ph_elec`, `temperature_ph`: Parameters for phonon calculations.
*   `ep_active_atoms(2)`, `n_ep_active_atoms`, `ep_bias_n`, `ep_bias(2)`: Electron-phonon coupling settings.
*   `l_ep_loe`, `l_ep`, `l_ep_reuse_lambda`: Logical flags for electron-phonon calculations.
*   `eta_cc`, `eta_elec`, `conv_dez`, `eps_geo`, `eps_real`: Damping, convergence, and geometry tolerances.
*   `estart`, `eend`: Energy range for calculations.
*   `kk(3)` (Real(dp)): Current k-point coordinates.
*   `ecc_dir`, `elec_l_dir`, `elec_r_dir`, `trans_file` (Character(strln)): Key directory and file paths.
*   `species_ecc(:)`, `species_elec_l(:)`, `species_elec_r(:)` (Character(strln), Allocatable): Atomic species in different regions.
*   `nkx_dez`, `nky_dez`, `nkz_diag`, `n_energy_steps`, `i_energy_start`, `i_energy_end`: K-point mesh, energy steps, and current energy indices.
*   `iunit_control` (Integer): File unit for the main control input file (`transp.ini`).
*   `oldnegf`, `oneshot`, `calc_current_density`, `converged`, `init_scf`, `lget_ti`: Logical flags for calculation flow and status.
*   `k_mat_mode` (Integer): Main operational mode of TransOMat.
*   `l_diag_fix_ef`, `l_diag_fixH`, `l_diag_fixK`, `l_diag_skip_check`, `diag_dim`: Control flags and dimension for matrix diagonalization routines.
*   `dftsigma`, `d_occ`, `d_virt`, `imu_dftsigma`, `jmu_dftsigma`: DFT+Sigma parameters.
*   `scf_conv` (Real(dp)): SCF convergence criterion.
*   `l_dump_nzs` (Logical): Flag to dump non-zero structure of matrices.

**System Geometry & Basis Set:**
*   `dlat_l(3,3)`, `rlat_l(3,3)`, etc.: Direct and reciprocal lattice vectors for central, left, and right regions.
*   `xyz_ecc(:,:)`, `xyz_elec_l(:,:)`, etc. (Real(dp), Allocatable): Atomic coordinates.
*   `kp_l(:,:)`, `wkp_l(:)`, etc. (Real(dp)/Integer, Allocatable): K-points and their weights.
*   `kpoints_from_file` (Logical): If true, k-points are read from a file.
*   `nktot` (Integer): Total number of k-points.
*   `lcr_info(:)` (Character(1), Allocatable): Information about atom region (Left, Center, Right).
*   `nat_ecc`, `nmu_ecc`, etc.: Number of atoms and basis functions (orbitals) in different regions.
*   `n_electrons_ecc(2)`, etc.: Number of electrons.
*   `imu_ecc(:)`, `atoms_c(:)`, `tomat(:,:)`, etc.: Mapping arrays for basis functions, atoms, and Conquest format compatibility.
*   `ndim_total` (Integer): Total dimension of matrices.

**Temporary/Derived Variables:**
*   `iik` (Integer): Current k-point index (often loop variable).
*   `zdosl(2)`, `zdosr(2)` (Complex(dp)): Density of states for left/right leads.

## Usage Examples

Variables in this module are accessed via `USE globals` statement in other Fortran files.

```fortran
subroutine calculate_something
  use globals
  implicit none
  real(dp) :: energy_midpoint

  if (l_ionode) then
    print *, "Starting calculation on master node."
  end if

  energy_midpoint = (estart + eend) / 2.0_dp
  ! ... use other global variables like nkx_dez, vb, etc. ...
end subroutine calculate_something
```

## Dependencies and Interactions

*   **Internal Dependencies:**
    *   `USE kinds`: For `dp` (double precision) and `strln` (string length) type kinds.
*   **External Library Dependencies:**
    *   `USE petscmat`: For PETSc matrix type definitions (`MatType`, `MatSolverType`). This also brings in the need for the PETSc header file via `#include <petsc/finclude/petscmat.h>`.
*   **Interactions:**
    *   Populated primarily by `init_control.f90` which reads values from `transp.ini`.
    *   Read by almost all other modules in TransOMat to control their behavior and access system parameters.
    *   Some variables (like `converged`, `iter`) are modified during the computational flow by modules like `scf_mod.f90`.

---
*This documentation was auto-generated by Jules.*
