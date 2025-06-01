# `scf_mod.f90` - Self-Consistent Field (SCF) Module

## Overview

The `scf_mod` module contains the `scf` subroutine, which is responsible for performing Self-Consistent Field (SCF) calculations to determine the electron density matrix of the system. This is a crucial step in many electronic structure calculations, as it iteratively refines the density until a consistent solution is found between the electron density and the effective potential. The SCF procedure in TransOMat involves calculating the Green's functions and integrating them over energy to obtain the density matrix. It supports calculations with and without an applied bias voltage (`vb`) and can handle different k-point sampling strategies (`l_k_on_demand`).

## Key Components

*   **`module scf_mod`**:
    *   **`subroutine scf()`**: The main entry point for SCF calculations.
        *   **Initialization:**
            *   Sets chemical potentials: `mul = ef_l + vl`, `mur = ef_r + vr`.
            *   Initializes PETSc matrix pointers for Green's functions, self-energies, and broadening matrices (`p_grrr`, `p_gllr`, `p_sigmarr`, `p_gammar`, `p_sigmalr`, `p_gammal`).
            *   Initializes temporary matrices for k-space density (`p_tmpcc1`, `p_d_tmp1`) and, if `ldouble_contour` is active, additional matrices for the double contour part (`p_tmpcc2`, `p_d_tmp2`, `p_d_dc_k`).
            *   If `calc_current_density` and `l_k_on_demand` are true, allocates `p_pnl` for non-local potential contributions to current density.
            *   Initializes real-space density matrices `p_dmatxy(i1,i2,iz)` to zero.
        *   **Integration Setup:**
            *   Calls `init_eq_int` (from `integrator_mod.f90`) to set up parameters for equilibrium/non-equilibrium integration contours.
            *   Sets `ldouble_contour` flag based on whether `vb` is non-zero.
        *   **K-Point Loop:** Iterates over k-points (from `nk_c` down to 1).
            *   If `l_k_on_demand` is true, calls `k_on_demand(ik, 1)` to compute necessary k-dependent matrices.
            *   **Non-Equilibrium Contribution (if `mul != mur` and not `l_no_noneq`):**
                *   Calls `init_neq_int` to set up the non-equilibrium integration path.
                *   Calls `adaptive_int3` with `get_gr_cc_neq` as the integrand function. This calculates the non-equilibrium part of the density matrix (`p_d_tmp1(1)` and `p_d_tmp2(1)` if `ldouble_contour`).
                *   Transforms the k-space non-equilibrium density `p_d_tmp1` back to real space and adds it to `p_dmatxy(:,:,0)` using `fourier_back_add`.
                *   If `ldouble_contour`, calculates `dc_weights_k` using `get_alpha` and updates `p_d_dc_k`.
            *   **Equilibrium Contribution (Right Lead Reference):**
                *   Calls `init_eq_int(mur, mur, 2)` (sets Fermi level to right lead).
                *   Calculates contribution from the semi-circular contour path using `adaptive_int3` with `get_gr_cc_eq`.
                *   Transforms result back to real space (`p_dmatxy(:,:,0)`).
                *   Calculates contribution from the path along the real energy axis using `adaptive_int3` with `get_gr_cc_eq`.
                *   Transforms result back to real space (`p_dmatxy(:,:,0)`).
                *   Calculates contribution from Fermi poles using `get_fermipole_contribution(p_d_tmp1(1), mur)`.
                *   Transforms result back to real space (`p_dmatxy(:,:,0)`).
            *   **Double Contour Adjustments (if `ldouble_contour`):**
                *   Updates `p_d_dc_k` with the equilibrium contributions (negatively).
                *   Recalculates Fermi pole contribution with `mul` and adds to `p_d_dc_k`.
                *   Recalculates real-axis integral with average Fermi function and adds to `p_d_dc_k`.
                *   Scales `p_d_dc_k` using `scale_Dne_dc` with `dc_weights_k`.
                *   Transforms `p_d_dc_k` to real space and adds to `p_dmatxy(:,:,0)`.
        *   **Finalization & Output:**
            *   Scales the total real-space density matrix `p_dmatxy` by `1.0 / (pi * nktot)`.
            *   If `calc_current_density` is true:
                *   Calls `dump_j(p_dmatxy, "j_c", ecc_dir)` to save conventional current density.
                *   Calls `dump_j(p_pnl, "p_nl", ecc_dir)` to save non-local potential current density.
            *   Else (not calculating current density):
                *   Calls `dump_kmat` to save the density matrix (K_negf_matrix2...).
                *   If `l_reaktor` is true, calls `dump_kmat_lr` for left and right electrodes.
        *   **Cleanup:** Destroys all allocated PETSc matrices.

## Important Variables/Constants

This subroutine uses and modifies global variables, primarily:

*   **Input (from `globals.f90`, set by `init_control` and `init_sys`):**
    *   `ef_l`, `ef_r`, `vl`, `vr`, `vb`: Fermi levels and voltages.
    *   `nmu_c`, `ncell_c`, `nk_c`, `nktot`: System dimensions and k-point info.
    *   `l_k_on_demand`, `ldouble_contour`, `calc_current_density`, `l_no_noneq`, `l_reaktor`: Control flags.
    *   `p_h00k_cc`, `p_s00k_cc` (and `_ik_` versions): K-space Hamiltonian and Overlap matrices.
    *   `p_h00_l`, `p_s00_l`, etc. (and `_ik_` versions): Electrode matrices.
    *   `p_dmatxy`: Real-space density matrix blocks (updated by this routine).
    *   `p_vnl00_cc`, `p_vnl00k_cc`: Non-local potential matrices (if `calc_current_density`).
    *   `kp_l`, `wkp_l`, `dlat_l`: K-points, weights, and lattice vectors for Fourier transforms.
    *   `pi`: Mathematical constant.
    *   `ecc_dir`, `elec_l_dir`: Directory paths for output.
*   **Internal PETSc Matrices:** `p_grrr`, `p_gllr`, `p_sigmarr`, `p_gammar`, `p_sigmalr`, `p_gammal`, `p_tmpcc1`, `p_tmpcc2`, `p_d_tmp1`, `p_d_tmp2`, `p_d_dc_k`, `p_pnl`.
*   **Integration Control:** `nint_order_el`, `maxsub`, `eps_int`, `eps_int_contour`, `contour_select`.

## Usage Examples

The `scf` subroutine is typically called from the main program `transomat.f90` either iteratively until convergence is reached or once if `oneshot` is true.

```fortran
! In transomat.f90 (conceptual)
if ( (oneshot .or. (converged .and. (.not. init_scf) .and. iter .ge. 2)) .and. (.not. l_ep) ) then
  ! ... possibly call transmission ...
  if (calc_current_density .and. (vb .ne. 0d0) .and. (converged .or. oneshot)) then
    call scf() ! Called to get density for current calculation after transmission
  end if
else if (l_ep) then
  ! ... call loe_run ...
else
  calc_current_density = .false.
  call scf() ! Called in the main SCF loop
end if
```

## Dependencies and Interactions

*   **Internal Dependencies:**
    *   `USE petscmat`, `USE petsc_mod`, `USE petsc_wrapper`: For PETSc matrix operations.
    *   `USE kinds`: For data types.
    *   `USE globals`: For accessing global parameters and matrices.
    *   `USE integrator_mod`: For `init_eq_int`, `init_neq_int`, `adaptive_int3`, `get_gr_cc_neq`, `get_gr_cc_eq`, `get_fermipole_contribution`, `contour_x`, `phi_low`, `el`, `eu`, `get_alpha`, `scale_Dne_dc`.
    *   `USE ft_mod`: For `fourier_back_add`.
    *   `USE kmat_mod`: For `dump_kmat`, `dump_kmat_lr`.
    *   `USE timing`: For `timing_surf`, `timing_matsolve`, `timing_matmat_solve`.
    *   `USE k_on_demand_mod`: For `k_on_demand`.
    *   `USE j_dump_mod`: For `dump_j`.
    *   `USE error_handler`: For `error()`.
*   **External Library Dependencies:**
    *   **PETSc:** Core for all matrix operations.
*   **Interactions:**
    *   Reads k-dependent Hamiltonian and overlap matrices from `globals` (prepared by `init_sys` or `k_on_demand`).
    *   Computes and updates the real-space density matrix `p_dmatxy` in `globals`.
    *   If `calc_current_density`, it computes and stores current density contributions to `p_dmatxy` and `p_pnl`, then calls `dump_j`.
    *   Otherwise, it calls `dump_kmat` to save the converged density matrix.
    *   The Green's functions and self-energies (`p_grrr`, `p_sigmarr`, etc.) are intermediate results used in `get_gr_cc_neq` and `get_gr_cc_eq` (from `integrator_mod.f90`).

---
*This documentation was auto-generated by Jules.*
