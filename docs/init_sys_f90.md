# `init_sys.f90` - System Initialization Module

## Overview

The `init_sys_mod` module, containing the `init_sys` subroutine, is responsible for the comprehensive setup of the physical system to be simulated by TransOMat. This involves reading detailed system information (atomic coordinates, species, lattice vectors, Fermi energies) for the central scattering region (ECC) and the left/right electrodes from their respective input files (typically `sys_info.dat` located in directories specified in `transp.ini`).

A significant part of this module's work is to load the real-space Hamiltonian (H), overlap (S), and Kohn-Sham/density (K) matrices from PETSc binary files for each region and for different unit cell offsets. If not operating in `k_on_demand` mode, it then performs Fourier transforms on these real-space matrices to obtain their k-space representations. It also handles the setup for electron-phonon (EP) coupling calculations if enabled, including reading phonon information and derivative of Hamiltonian (dH) matrices.

## Key Components

*   **`module init_sys_mod`**:
    *   **`subroutine init_sys()`**: The primary subroutine that orchestrates all system initialization tasks.
        *   **Read System Information:**
            *   Calls `read_sysinfo` (from `readsys.f90`) for the central region (ECC), left electrode (L), and right electrode (R) to populate global variables like atomic coordinates (`xyz_ecc`, `xyz_elec_l`, `xyz_elec_r`), species, basis function counts per atom (`imu_ecc`, etc.), lattice vectors (`dlat_c`, etc.), number of cells (`ncell_c`), and Fermi energies (`ef_c`, `ef_l`, `ef_r`).
            *   Calculates `nmu_c`, `nmu_l`, `nmu_r` (number of basis functions in the isolated central part of ECC, left interface, right interface).
            *   Creates `imu_to_at` mapping.
        *   **K-Point Setup:**
            *   Calls `get_rlat` (from `lattice_mod.f90`) to calculate reciprocal lattice vectors and generate k-points (`kp_c`, `kp_l`, `kp_r`) and weights (`wkp_c`, etc.) unless `kpoints_from_file` is true (in which case, `read_kpoint_file` is called).
        *   **Matrix Allocation:**
            *   Allocates numerous global PETSc matrix pointers for real-space matrices (e.g., `p_h00_cc(ix,iy)`, `p_s00_cc(ix,iy)`, `p_dmatxy(ix,iy,iz)`) representing Hamiltonian, overlap, and density/Kohn-Sham matrices for different cell offsets (ix, iy, iz).
            *   If `calc_current_density` is true, allocates `p_vnl00_cc`.
            *   If `k_mat_mode == 2`, allocates additional matrices for off-diagonal blocks in z-direction (`p_h01_cc`, `p_s01_cc`, etc.).
            *   If electron-phonon coupling (`l_ep`) is enabled, allocates `p_dh_cc` for Hamiltonian derivatives.
        *   **Phonon System Initialization (if `l_ep`):**
            *   Calls `init_phonons()` and `get_phonones()` (from `phonon_mod.f90`).
            *   Allocates k-dependent phonon matrices (`p_dhk_cc`, `p_Kphonon00k_cc`, etc.).
        *   **Real-Space Matrix Loading:**
            *   Determines PETSc parallel layout for matrices (`nzs_to_procs` or `rows_to_procs`).
            *   Iterates through unit cell offsets (`i1` for x, `i2` for y) and loads:
                *   `K_u_i1_i2_0_petsc.dat` into `p_dmatxy(i1,i2,0)` (and also `_petsc.dat` for `iz = -1, 1`).
                *   `H_u_i1_i2_0_petsc.dat` into `p_h00_cc(i1,i2)`.
                *   `S_i1_i2_0_petsc.dat` into `p_s00_cc(i1,i2)`.
                *   Similar loading for `p_h01_cc`, `p_s01_cc`, etc. if `k_mat_mode == 2`.
                *   If `l_ep`, loads `dH_iat_ixyz_i1_i2_0_petsc.dat` into `p_dh_cc(i1,i2,iat,ixyz)`.
                *   If `calc_current_density`, loads `Vnl_i1_i2_0_petsc.dat` into `p_vnl00_cc(i1,i2)`.
            *   Similar loading process for left (`p_k00_l`, `p_s00_l`, `p_h00_l`, `p_k01_l`, etc.) and right (`p_k00_r`, etc.) electrode matrices from their respective directories.
        *   **Conquest Information:**
            *   Calls `read_conquest_info` to read further structural details if needed.
        *   **K-Space Matrix Setup (if not `l_k_on_demand`):**
            *   Allocates k-space PETSc matrix pointers (e.g., `p_h00k_cc(ik)`, `p_s00k_cc(ik)`).
            *   Duplicates the structure of real-space matrices for k-space counterparts.
            *   Calls `fourier_trans` (from `ft_mod.f90`) to transform all loaded real-space H, S (and dH if `l_ep`) matrices into k-space for each k-point.
            *   **Electrode Embedding:** Adds the k-space electrode Hamiltonians and overlap matrices (`p_h00k_l`, `p_s00k_l`, `p_h00k_r`, `p_s00k_r`) into the corresponding blocks of the central region's k-space matrices (`p_h00k_cc`, `p_s00k_cc`).
            *   Deallocates real-space matrices after transformation.
        *   **K-Space Matrix Setup (if `l_k_on_demand`):**
            *   Initializes k-space matrix pointers (e.g., `p_h00_ik_cc(1)`, `p_s00_ik_cc(1)`) by performing a Fourier transform for a single k-point (usually the first one) to establish matrix structure and pre-allocate necessary structures like `p_invGr`.
            *   If `l_diag_fixH` is true, it also performs the electrode embedding for the real-space matrices `p_h00_cc(i1,i2)` and `p_s00_cc(i1,i2)`.
        *   **Electron-Phonon Lambda Calculation (if `l_ep`):**
            *   Allocates `p_ep_lambda_k`.
            *   Calls `get_lambda_ep` for each k-point to calculate electron-phonon coupling strengths.
        *   **Debugging Output:** Optionally dumps non-zero structure of key matrices if `l_dump_nzs` is true.

## Important Variables/Constants

This subroutine populates and uses many global variables from `globals.f90`. Key ones include:

*   **Input Data Holders:** `xyz_ecc`, `species_ecc`, `imu_ecc`, `dlat_c`, `ef_c`, (and similar for _l, _r).
*   **Matrix Pointers (Real Space):** `p_h00_cc`, `p_s00_cc`, `p_dmatxy`, `p_h01_cc`, `p_s01_cc`, `p_dh_cc`, `p_vnl00_cc`, and their electrode counterparts (`_l`, `_r`).
*   **Matrix Pointers (K-Space):** `p_h00k_cc`, `p_s00k_cc`, `p_k00k_cc` (density/KS in k-space for SCF), `p_dhk_cc` (dH in k-space), `p_vnl00k_cc`, and electrode/coupling counterparts. If `l_k_on_demand`: `p_h00_ik_cc`, `p_s00_ik_cc`, etc. (for a single k-point).
*   **Control Flags:** `k_mat_mode`, `l_diag_skip_check`, `calc_current_density`, `l_ep`, `l_ep_reuse_lambda`, `l_k_on_demand`, `l_dump_nzs`, `l_reaktor`, `l_diag_fixH`.
*   **System Dimensions:** `nmu_c`, `nmu_l`, `nmu_r`, `nmu_ecc`, `nat_ecc`, `ncell_c`.
*   **K-points:** `kp_c`, `kp_l`, `kp_r`, `wkp_l`, `nk_c`, `nk_l`, `nk_r`, `nkx_dez`, `nky_dez`, `nkz_diag`.
*   **PETSc variables:** `mattype_sparse`, `mattype_electrode`, `p_one`, `p_zero`.
*   **File Paths:** `ecc_dir`, `elec_l_dir`, `elec_r_dir`.

## Usage Examples

`init_sys` is called once by the main program `transomat.f90` after `init_control` has run and populated the necessary global control parameters.

```fortran
program transomat
  use init_sys_mod
  ! ... other use statements ...
  implicit none
  ! ...
  call init_control()
  call init_sys()      ! Initializes the system matrices and parameters
  ! ... proceed with calculations ...
end program transomat
```

## Dependencies and Interactions

*   **Internal Dependencies:**
    *   `USE slepceps`, `USE petscmat`, `USE petsc_mod`, `USE petsc_wrapper`: For PETSc/SLEPc matrix and vector operations, types, and utilities.
    *   `USE kinds`: For data type definitions.
    *   `USE globals`: Heavily used for accessing control parameters and storing system matrices/data.
    *   `USE readsys`: For `read_sysinfo` subroutine.
    *   `USE lattice_mod`: For `get_rlat` (reciprocal lattice and k-points).
    *   `USE ft_mod`: For `fourier_trans` and `fourier_trans_ik`, `fourier3d_trans_ik`.
    *   `USE conquest_mod`: For `read_conquest_info`.
    *   `USE kmat_mod`: (Potentially, e.g. `dump_kmat_lr` is commented out).
    *   `USE slepc_mod`: (Linked, specific usage not detailed in main flow).
    *   `USE init_ep_coupling`: For `init_ep_active_mat`, `check_lambdas_on_disk`.
    *   `USE phonon_mod`: For `init_phonons`, `get_phonones`, `get_lambda_ep`.
    *   `USE loe_mod`: (Linked, specific usage not detailed in main flow).
    *   `USE misc`: For utilities like `int2str`.
    *   `USE error_handler`: For `error()` routine.
    *   `USE k_on_demand_mod`: For `k_on_demand()` if that mode is active (though `init_sys` itself prepares for it rather than calling it).
    *   Calls `check_sys` (from `check_sys.f90`, not explicitly `USEd` but likely linked).
*   **External Library Dependencies:**
    *   **PETSc/SLEPc:** Core for all matrix and vector operations.
    *   **MPI:** Implicitly through PETSc and for `mpi_barrier`.
*   **Interactions:**
    *   Reads numerous data files (e.g., `sys_info.dat`, `K_u_..._petsc.dat`, `H_u_..._petsc.dat`, `S_..._petsc.dat`, `dH_..._petsc.dat`) from directories specified in `globals` (e.g., `ecc_dir`).
    *   Populates many global PETSc matrix variables (e.g., `p_h00k_cc`, `p_s00k_cc`) and other system parameters in `globals.f90`.
    *   Its execution flow is significantly affected by flags like `l_k_on_demand`, `l_ep`, `k_mat_mode`.

---
*This documentation was auto-generated by Jules.*
