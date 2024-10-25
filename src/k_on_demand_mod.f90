module k_on_demand_mod
  implicit none

contains

  subroutine k_on_demand(ik, iwhat)
    ! Purpose: This subroutine performs various calculations depending on the value of `iwhat`.  It primarily involves Fourier transforms of matrices related to Hamiltonian (H), overlap (S), and possibly current density (vnl) calculations.

    ! Parameters:
    !   ik: Integer, the k-point index.
    !   iwhat: Integer, determines the type of calculation to perform:
    !          1: H, S, and electrode calculations.
    !          2: D, S for bond order analysis.
    !          3: H, S calculations with Hermitian checks.

#include <petsc/finclude/petsc.h>
    use petscmat  ! PETSc matrix operations

    use globals   ! Module containing global variables
    use ft_mod    ! Module containing Fourier transform routines
    use petsc_mod ! Module containing PETSc-related variables and functions

    implicit none

    integer :: ik, iwhat

    integer :: ierr, i1, i2, i3  ! Error codes and temporary integer variables
    
    Mat :: A                ! PETSc matrix for Hermitian checks
    PetscScalar :: pp       ! PETSc scalar variable
    PetscReal :: norm       ! PETSc real variable for norm calculation

    if (iwhat .eq. 1) then ! H,S ecc and electrode
      ! Nullify pointers to avoid potential issues.  This is crucial for memory management, especially when dealing with dynamically allocated arrays.
      nullify (p_h00k_cc)
      nullify (p_s00k_cc)
      if (calc_current_density) nullify (p_vnl00k_cc)
      nullify (p_h00k_l)
      nullify (p_s00k_l)
      nullify (p_h10k_l)
      nullify (p_s10k_l)
      nullify (p_h01k_l)
      nullify (p_s01k_l)
      nullify (p_h00k_r)
      nullify (p_s00k_r)
      nullify (p_h10k_r)
      nullify (p_s10k_r)
      nullify (p_h01k_r)
      nullify (p_s01k_r)      

      ! Zero out the entries of the matrices.  This ensures that previous calculations don't interfere with the current ones.
      call MatZeroEntries(p_h00_ik_cc(1), ierr)
      call MatZeroEntries(p_s00_ik_cc(1), ierr)
      if (calc_current_density) call MatZeroEntries(p_vnl00_ik_cc(1), ierr)
      call MatZeroEntries(p_h00_ik_l(1), ierr)
      call MatZeroEntries(p_s00_ik_l(1), ierr)
      call MatZeroEntries(p_h10_ik_l(1), ierr)
      call MatZeroEntries(p_s10_ik_l(1), ierr)
      call MatZeroEntries(p_h01_ik_l(1), ierr)
      call MatZeroEntries(p_s01_ik_l(1), ierr)
      call MatZeroEntries(p_h00_ik_r(1), ierr)
      call MatZeroEntries(p_s00_ik_r(1), ierr)
      call MatZeroEntries(p_h10_ik_r(1), ierr)
      call MatZeroEntries(p_s10_ik_r(1), ierr)
      call MatZeroEntries(p_h01_ik_r(1), ierr)
      call MatZeroEntries(p_s01_ik_r(1), ierr)

      ! Perform Fourier transforms.  These transforms convert data from real space to reciprocal space (k-space).
      call fourier_trans_ik(p_h00_cc, p_h00_ik_cc(1), kp_r, ik, dlat_r, ncell_c(1), ncell_c(2), 0)
      call fourier_trans_ik(p_s00_cc, p_s00_ik_cc(1), kp_r, ik, dlat_r, ncell_c(1), ncell_c(2), 0)
      if (calc_current_density) &
        call fourier_trans_ik(p_vnl00_cc, p_vnl00_ik_cc(1), kp_r, ik, dlat_r, ncell_c(1), ncell_c(2), 0)
      call fourier_trans_ik(p_h00_l, p_h00_ik_l(1), kp_r, ik, dlat_r, ncell_c(1), ncell_c(2), 0)
      call fourier_trans_ik(p_s00_l, p_s00_ik_l(1), kp_r, ik, dlat_r, ncell_c(1), ncell_c(2), 0)
      call fourier_trans_ik(p_h10_l, p_h10_ik_l(1), kp_r, ik, dlat_r, ncell_c(1), ncell_c(2), 0)
      call fourier_trans_ik(p_s10_l, p_s10_ik_l(1), kp_r, ik, dlat_r, ncell_c(1), ncell_c(2), 0)
      call fourier_trans_ik(p_h01_l, p_h01_ik_l(1), kp_r, ik, dlat_r, ncell_c(1), ncell_c(2), 0)
      call fourier_trans_ik(p_s01_l, p_s01_ik_l(1), kp_r, ik, dlat_r, ncell_c(1), ncell_c(2), 0)
      call fourier_trans_ik(p_h00_r, p_h00_ik_r(1), kp_r, ik, dlat_r, ncell_c(1), ncell_c(2), 0)
      call fourier_trans_ik(p_s00_r, p_s00_ik_r(1), kp_r, ik, dlat_r, ncell_c(1), ncell_c(2), 0)
      call fourier_trans_ik(p_h10_r, p_h10_ik_r(1), kp_r, ik, dlat_r, ncell_c(1), ncell_c(2), 0)
      call fourier_trans_ik(p_s10_r, p_s10_ik_r(1), kp_r, ik, dlat_r, ncell_c(1), ncell_c(2), 0)
      call fourier_trans_ik(p_h01_r, p_h01_ik_r(1), kp_r, ik, dlat_r, ncell_c(1), ncell_c(2), 0)
      call fourier_trans_ik(p_s01_r, p_s01_ik_r(1), kp_r, ik, dlat_r, ncell_c(1), ncell_c(2), 0)

      ! Assign the results of the Fourier transforms to the appropriate variables.
      p_h00k_cc(ik:ik) => p_h00_ik_cc(1:1)
      p_s00k_cc(ik:ik) => p_s00_ik_cc(1:1)
      if (calc_current_density) p_vnl00k_cc(ik:ik) => p_vnl00_ik_cc(1:1)
      p_h00k_l(ik:ik) => p_h00_ik_l(1:1)
      p_s00k_l(ik:ik) => p_s00_ik_l(1:1)
      p_h10k_l(ik:ik) => p_h10_ik_l(1:1)
      p_s10k_l(ik:ik) => p_s10_ik_l(1:1)
      p_h01k_l(ik:ik) => p_h01_ik_l(1:1)
      p_s01k_l(ik:ik) => p_s01_ik_l(1:1)
      p_h00k_r(ik:ik) => p_h00_ik_r(1:1)
      p_s00k_r(ik:ik) => p_s00_ik_r(1:1)
      p_h10k_r(ik:ik) => p_h10_ik_r(1:1)
      p_s10k_r(ik:ik) => p_s10_ik_r(1:1)
      p_h01k_r(ik:ik) => p_h01_ik_r(1:1)
      p_s01k_r(ik:ik) => p_s01_ik_r(1:1)

    else if (iwhat .eq. 2) then ! D,S for bond order analysis
    
      nullify (p_k00k_cc)
      nullify (p_s00k_cc)


      call MatZeroEntries(p_k00_ik_cc(1), ierr)
      call MatZeroEntries(p_s00_ik_cc(1), ierr)


      call fourier_trans_ik(p_dmatxy(:, :, 0), p_k00_ik_cc(1), kp_r, ik, dlat_r, ncell_c(1), ncell_c(2), 0)
      call fourier_trans_ik(p_s00_cc, p_s00_ik_cc(1), kp_r, ik, dlat_r, ncell_c(1), ncell_c(2), 0)

      p_k00k_cc(ik:ik) => p_k00_ik_cc(1:1)
      p_s00k_cc(ik:ik) => p_s00_ik_cc(1:1)

    else if (iwhat .eq. 3) then ! H,S ecc
    
      nullify (p_h00k_cc)
      nullify (p_s00k_cc)


      call MatZeroEntries(p_h00_ik_cc(1), ierr)
      call MatZeroEntries(p_s00_ik_cc(1), ierr)

      ! 3D Fourier transform for case 3.
      call fourier3d_trans_ik(p_h00_diag, p_h00_ik_cc(1), kp_c, ik, dlat_c, ncell_c(1), &
        ncell_c(2), ncell_c(3), 0)
      call fourier3d_trans_ik(p_s00_diag, p_s00_ik_cc(1), kp_c, ik, dlat_c, ncell_c(1), &
        ncell_c(2), ncell_c(3), 0)

      p_h00k_cc(ik:ik) => p_h00_ik_cc(1:1)
      p_s00k_cc(ik:ik) => p_s00_ik_cc(1:1)
      
      pp = 0.5d0
      
      ! Hermitian checks and norm calculations for H and S matrices.
      call MatHermitianTranspose(p_h00_ik_cc(1), MAT_INITIAL_MATRIX, A, ierr)
      call MatAXPY(A,p_minus1,p_h00_ik_cc(1),DIFFERENT_NONZERO_PATTERN, ierr)
      call MatNorm(A, NORM_FROBENIUS, norm, ierr)
      call MatDestroy(A,ierr)
      write (pstr_out, fmt='(A,e24.12)') "norm(H(k)-H(k)')", norm      
      call petsc_print_master()
      
      call MatHermitianTranspose(p_h00_ik_cc(1), MAT_INITIAL_MATRIX, A, ierr)
      call MatAXPY(p_h00_ik_cc(1),p_one,A,DIFFERENT_NONZERO_PATTERN, ierr)
      call MatScale(p_h00_ik_cc(1), pp, ierr)
      call MatDestroy(A,ierr)

      call MatHermitianTranspose(p_h00_ik_cc(1), MAT_INITIAL_MATRIX, A, ierr)
      call MatAXPY(A,p_minus1,p_h00_ik_cc(1),DIFFERENT_NONZERO_PATTERN, ierr)
      call MatNorm(A, NORM_FROBENIUS, norm, ierr)
      call MatDestroy(A,ierr)
      write (pstr_out, fmt='(A,e24.12)') "norm(H(k)-H(k)')", norm      
      call petsc_print_master()



      call MatHermitianTranspose(p_s00_ik_cc(1), MAT_INITIAL_MATRIX, A, ierr)
      call MatAXPY(A,p_minus1,p_s00_ik_cc(1),DIFFERENT_NONZERO_PATTERN, ierr)
      call MatNorm(A, NORM_FROBENIUS, norm, ierr)
      call MatDestroy(A,ierr)
      write (pstr_out, fmt='(A,e24.12)') "norm(S(k)-S(k)')", norm
      call petsc_print_master()

      call MatHermitianTranspose(p_s00_ik_cc(1), MAT_INITIAL_MATRIX, A, ierr)
      call MatAXPY(p_s00_ik_cc(1),p_one,A,DIFFERENT_NONZERO_PATTERN, ierr)
      call MatScale(p_s00_ik_cc(1), pp, ierr)
      call MatDestroy(A,ierr)

      call MatHermitianTranspose(p_s00_ik_cc(1), MAT_INITIAL_MATRIX, A, ierr)
      call MatAXPY(A,p_minus1,p_s00_ik_cc(1),DIFFERENT_NONZERO_PATTERN, ierr)
      call MatNorm(A, NORM_FROBENIUS, norm, ierr)
      call MatDestroy(A,ierr)
      write (pstr_out, fmt='(A,e24.12)') "norm(S(k)-S(k)')", norm
      call petsc_print_master()

      

    else

    end if

  end subroutine k_on_demand

end module k_on_demand_mod