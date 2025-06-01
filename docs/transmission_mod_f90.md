# `transmission_mod.f90` - Transmission Calculation Module

## Overview

The `transmission_mod` module is responsible for calculating the electron transmission coefficient `T(E)` through the system as a function of energy `E`. This is a key quantity in quantum transport and is determined using the Green's function formalism, specifically involving the formula `Tr[Gamma_L * G_r * Gamma_R * G_a]`, where `G_r` and `G_a` are the retarded and advanced Green's functions of the central scattering region, and `Gamma_L` and `Gamma_R` are broadening matrices due to coupling with the left and right electrodes, respectively. The module provides functionalities for standard transmission calculations over a range of energies and an on-demand version for specific k-points.

## Key Components

*   **`module transmission_mod`**:
    *   **`subroutine get_transmission(x, trans, t_i, dosc, dosl, dosr)`**:
        *   Calculates the transmission `trans` for a single energy `x`.
        *   `x` (Real(dp), in): Energy at which to calculate transmission.
        *   `trans` (Complex(dp), out): Calculated transmission coefficient.
        *   `t_i` (Vec, optional, out): Vector to store transmission eigenvalues/channels if `lget_ti` is true.
        *   `dosc`, `dosl`, `dosr` (Real(dp), optional, out): Placeholder for density of states (not fully implemented in the provided snippet).
        *   **Steps:**
            1.  Calls `init_gr(z, p_gr_inv)` (from `integrator_mod.f90`) to compute the retarded Green's function `p_gr_inv` (actually `(E*S - H - Sigma_L - Sigma_R)^-1`) at energy `z=x`.
            2.  Solves for specific blocks of the Green's function `p_g_block = p_gr_inv * p_gammar` using `petsc_call_solver`.
            3.  Extracts the `G_13` block (coupling left to right) into `p_g13`.
            4.  If `lget_ti` (calculate transmission eigenvalues/channels) is true:
                *   Computes `p_sqrt_gammar = Gamma_R^(1/2)`.
                *   Calculates `A = G_13 * Gamma_R^(1/2)`.
                *   Forms the transmission matrix `T_m = A^H * Gamma_L * A`.
                *   Calculates `trans` as `Trace(T_m)`.
                *   Diagonalizes `T_m` to get transmission eigenvalues into `t_i`.
            5.  Else (standard transmission):
                *   Computes `p_g13gammar = G_13 * Gamma_R`.
                *   Computes `p_gammalg13 = Gamma_L * G_13`.
                *   Calculates `trans = Trace(p_g13gammar * p_gammalg13^H)` using `petsc_matmatconj_restricted`.
    *   **`subroutine transmission()`**:
        *   Calculates transmission over a range of energies (`estart` to `eend` in `n_energy_steps`).
        *   Handles restarting calculations if a partial `trans_file` exists.
        *   **Loop over energies (`ie`):**
            *   Calculates current energy `energy = estart + de*real(ie, 8)`.
            *   Initializes `trans = 0.0`.
            *   **Loop over k-points (`ik` from 1 to `nk_c`):**
                *   Sets global `iik = ik`.
                *   Calls `get_transmission(energy, transk, p_tik)` to get transmission `transk` for the current energy and k-point.
                *   Accumulates `trans = trans + transk * real(wkp_r(ik), dp)`.
                *   If `lget_ti`, accumulates transmission channels `ti`.
            *   Averages `trans` (and `ti`) over k-points.
            *   Writes `energy`, `real(trans)`, and optionally `real(ti)` to `trans_file`.
            *   Also writes to DOS files (dos_r.dat, dos_l.dat, dos_cc.dat), though DOS calculation seems minimal here.
        *   Cleans up PETSc matrices.
    *   **`subroutine transmission_k_on_demand()`**:
        *   Similar to `transmission()` but designed to work with `l_k_on_demand = .true.`.
        *   If `trans_file` (final output) exists, it skips calculation.
        *   Manages temporary files (`trans.k.tmp`, `trans.tmp`) to accumulate results per k-point and then combine them.
        *   **Outer loop over k-points (`ik`):**
            *   Calls `k_on_demand(ik, 1)` (from `k_on_demand_mod.f90`) to set up matrices for the current k-point.
            *   **Inner loop over energies (`ie`):**
                *   Calls `get_transmission(energy, transk, p_tik)`.
                *   Writes `energy`, `real(transk)` (and `ti` if `lget_ti`) to `trans.k.tmp` for the current k-point.
            *   After each k-point's energy loop, reads the previous accumulated `trans.tmp`, adds the current k-point's contribution (weighted), and writes back to a new `trans.tmp`.
        *   After all k-points, renames the final `trans.tmp` to `trans_file`.
        *   Cleans up PETSc matrices.

## Important Variables/Constants

*   **Input (from `globals.f90`):**
    *   `estart`, `eend`, `n_energy_steps`, `i_energy_start`, `i_energy_end`: Energy range and steps.
    *   `nk_c`, `kp_l`, `wkp_r`: K-point information.
    *   `nmu_l`, `nmu_r`, `nmu_c`: Number of basis functions for left/right leads and central region.
    *   `p_h00k_cc`, `p_s00k_cc`, etc. (if not `l_k_on_demand`): K-space Hamiltonian and Overlap matrices.
    *   `p_h00_ik_cc`, `p_s00_ik_cc`, etc. (if `l_k_on_demand`): Real-space matrices for current k-point setup.
    *   `p_h00_l`, `p_s00_l`, `p_h01_l`, `p_s01_l`, `p_h10_l`, `p_s10_l` (and `_r` versions): Real-space electrode matrices (used by `k_on_demand` which then provides k-space versions).
    *   `trans_file` (Character(strln)): Name of the output transmission file.
    *   `lget_ti` (Logical): Flag to calculate transmission channels/eigenvalues.
    *   `l_k_on_demand` (Logical): Flag for on-demand k-point calculations.
    *   `l_ionode` (Logical): Flag for master MPI process.
    *   `solver_mode` (Integer): For PETSc solver selection in `get_transmission`.
*   **Key PETSc Matrices (Local or passed to `get_transmission`):**
    *   `p_gr_inv`: Retarded Green's function `(E*S - H_eff)^-1`.
    *   `p_gammal`, `p_gammar`: Broadening matrices for left/right leads.
    *   `p_g13`: `G_13` block of the Green's function.
    *   `p_trace`: Temporary matrix to hold product before taking trace.
    *   `p_tt`: Transmission matrix `A^H * Gamma_L * A` (if `lget_ti`).
*   **Output:**
    *   `trans` (Complex(dp)): Transmission coefficient `T(E)`.
    *   `t_i` (Vec): Transmission eigenvalues.

## Usage Examples

The `transmission()` or `transmission_k_on_demand()` subroutines are typically called from the main program `transomat.f90` after the system setup and potentially after an SCF cycle if the Hamiltonian needs to be self-consistent.

```fortran
! In transomat.f90 (conceptual)
if ( (oneshot .or. (converged .and. (.not. init_scf) .and. iter .ge. 2)) .and. (.not. l_ep) ) then
  write (pstr_out, fmt='(A)') "get transmission"; call petsc_print_master()
  if (l_k_on_demand) then
    call transmission_k_on_demand()
  else
    call transmission()
  end if
  ! ...
end if
```

## Dependencies and Interactions

*   **Internal Dependencies:**
    *   `USE petscmat`, `USE petsc_mod`, `USE petsc_wrapper`: For PETSc matrix/vector operations.
    *   `USE globals`: For accessing global parameters, k-space matrices, and control flags.
    *   `USE integrator_mod`: Crucially for `init_gr` which computes the Green's function.
    *   `USE slepc_mod`: For `diag_mat` if `lget_ti` is true (diagonalizing the transmission matrix).
    *   `USE k_on_demand_mod`: If `l_k_on_demand` is true, uses `k_on_demand` subroutine.
    *   `USE error_handler`: For `error()`.
*   **External Library Dependencies:**
    *   **PETSc/SLEPc:** Essential for all numerical computations.
    *   **MPI:** Implicitly via PETSc for parallel operations and `MPI_Barrier`, `MPI_Reduce`.
*   **Interactions:**
    *   Reads k-space Hamiltonian and overlap matrices from `globals` (prepared by `init_sys.f90` or `k_on_demand_mod.f90`).
    *   Calls `init_gr` from `integrator_mod.f90` to get the necessary Green's functions.
    *   Writes transmission data to the file specified by `trans_file` and potentially to DOS files.
    *   Manages temporary files if `transmission_k_on_demand()` is used.

---
*This documentation was auto-generated by Jules.*
