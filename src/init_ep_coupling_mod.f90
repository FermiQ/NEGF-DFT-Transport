module init_ep_coupling
  implicit none

contains

  subroutine check_lambdas_on_disk(l_load)
    ! Purpose: Checks if lambda files exist on disk for all k-points and modes.
    !
    ! Functionality: Iterates over all k-points and electron-phonon modes, 
    ! constructing filenames for lambda files. It then checks if these files exist. 
    ! If any file is missing, it sets l_load to .true., indicating that lambdas 
    ! need to be computed. The result is broadcast to all processes.
    !
    ! Parameters:
    !   l_load (logical, output):  A logical flag indicating whether lambda files need to be loaded (.true.) or not (.false.).

    use kinds
    use misc
    use globals
    implicit none

    logical :: l_load

    integer :: imode, ik, ierr
    character(strln) :: lambda_file, s_kx, s_ky
    logical :: l_exist

    if (l_ionode) then
      l_load = .false.
      ik_loop: do ik = 1, nk_c
        imode_loop: do imode = 1, n_ep_modes_k(ik)
          write (s_kx, fmt='(e16.9)') kp_r(1, ik)
          write (s_ky, fmt='(e16.9)') kp_r(2, ik)
          lambda_file = "ep_reaktor/lambda_"//trim(adjustl(s_kx))//"_"//trim(adjustl(s_ky))
          lambda_file = trim(lambda_file)//"_"//trim(int2str(imode))//".dat"
          inquire (file=trim(lambda_file), exist=l_exist)
!~             write(0,fmt='(A,l)') trim(lambda_file),l_exist
          if (.not. l_exist) then
            l_load = .true.
            exit ik_loop
          end if
        end do imode_loop
      end do ik_loop
    end if

    call MPI_bcast(l_load, 1, MPI_LOGICAL, 0, PETSC_COMM_WORLD, ierr)

  end subroutine check_lambdas_on_disk

  subroutine get_lambda_ep(ik)
    ! Purpose: Computes or loads the electron-phonon coupling matrix (lambda) for a given k-point.
    !
    ! Functionality: This subroutine either loads pre-computed lambda matrices from disk 
    ! or computes them if they don't exist.  It handles the loading and computation of 
    ! the electron-phonon coupling matrix for each mode at a given k-point.
    !
    ! Parameters:
    !   ik (integer, input): The index of the k-point.

#include <petsc/finclude/petsc.h>
    use petscmat
    use petsc_mod
    use petsc_wrapper
    use kinds
    use globals
    use misc

    implicit none

    integer, intent(in) :: ik

    integer :: iactive, mactive, icrt, imode, nl1, nl2, iat, ii, ierr
    PetscScalar :: prefac, mode_energy, mode_weight
    PetscReal :: norm
    logical :: l_exist
    character(strln) :: lambda_file, s_kx, s_ky

    l_exist = l_ep_reuse_lambda

    do imode = 1, n_ep_modes_k(ik)
      call petsc_vec_getvalue(imode - 1, mode_energy, p_phonon_EW(ik))
      prefac = zsqrt(0.5d0/mode_energy)
      ii = (ep_active_atoms(1) - 1)*3 !+1 is not need as we use C style boundarys i.e. 0 based

      if (l_ep_reuse_lambda) then
        write (s_kx, fmt='(e16.9)') kp_r(1, ik)
        write (s_ky, fmt='(e16.9)') kp_r(2, ik)
        lambda_file = "ep_reaktor/lambda_"//trim(adjustl(s_kx))//"_"//trim(adjustl(s_ky))
        lambda_file = trim(lambda_file)//"_"//trim(int2str(imode))//".dat"
        inquire (file=trim(lambda_file), exist=l_exist)
        if (l_exist) then
          call petsc_mat_direct_load(p_ep_lambda_k(imode, ik), lambda_file, ierr)
          l_exist = ierr .eq. 0
        end if
      end if

      if (.not. l_exist) then
        do iat = ep_active_atoms(1), ep_active_atoms(2)
          do icrt = 1, 3
            call petsc_mat_getvalue(ii, imode - 1, mode_weight, p_phonon_EV(ik), 1, PETSC_COMM_WORLD)
!~               write(pstr_out,fmt='(3i8,A,2e24.12)') iat,icrt,ii," mode_weight",mode_weight ; call petsc_print_master()
            call petsc_aXpY(p_ep_lambda_k(imode, ik), p_dhk_cc(ik, iat, icrt), mode_weight, PETSC_FALSE, .false.)
            ii = ii + 1
          end do
        end do
        call MatScale(p_ep_lambda_k(imode, ik), prefac, ierr)

        if (l_ep_reuse_lambda) call petsc_mat_direct_save(p_ep_lambda_k(imode, ik), lambda_file, ierr)

      end if

      call MatNorm(p_ep_lambda_k(imode, ik), NORM_FROBENIUS, norm, ierr)
      write (pstr_out, fmt='(i6,A,5e24.12,l)') imode, " mode energy, prefac, norm(lambda)", mode_energy, prefac, norm, l_exist; call petsc_print_master()
    end do

  end subroutine get_lambda_ep

  subroutine init_ep_active_mat(p_A, p_B, nat1, nat2, imu_to_at)
    ! Purpose: Creates a submatrix from a larger matrix, selecting only rows and columns corresponding to active atoms.
    !
    ! Functionality: This subroutine takes a matrix p_A and creates a submatrix p_B. 
    ! It selects only the rows and columns of p_A that correspond to atoms within the 
    ! range [nat1, nat2] using the mapping provided by the imu_to_at array.  It uses 
    ! PETSc's matrix manipulation functions for efficiency.
    !
    ! Parameters:
    !   p_A (Mat, input): The input matrix.
    !   p_B (Mat, output): The output submatrix containing only the active atoms.
    !   nat1 (integer, input): The starting index of the active atoms.
    !   nat2 (integer, input): The ending index of the active atoms.
    !   imu_to_at (integer array, input): Mapping between internal indices and atom indices.

#include <petsc/finclude/petsc.h>
    use petscmat
    use petsc_mod
    use petsc_wrapper
    use kinds
    use globals, only: nmu_c, mattype_sparse
    use error_handler

    implicit none

    Mat :: p_A, p_B
    integer :: nat1, nat2
    integer, allocatable :: imu_to_at(:)

    integer :: iat1, iat2, nl1, nl2, irow, ierr, nzcol, inz, ii(1)
    integer, allocatable :: icol(:)
    PetscScalar, allocatable :: p_vals(:)
    Mat :: p_tmp

    allocate (icol(nmu_c), p_vals(nmu_c), stat=ierr)
    if (ierr .ne. 0) then
      write (errormsg, fmt='(A,i8)') "allocation error icol ", ierr
      call error()
    end if
    p_vals = 0d0
    call petsc_get_a_with_b_c(p_tmp, p_A, p_A, mattype_sparse)
    call MatSetOption(p_tmp, MAT_NEW_NONZERO_LOCATIONS, PETSC_TRUE, ierr)

    call MatGetOwnershipRange(p_A, nl1, nl2, ierr)
    nl1 = nl1 + 1 ! to fortran numbering, nl2 is already correct as it is max_number of row+1

    do irow = nl1, nl2
      if ((imu_to_at(irow) .lt. nat1) .or. (imu_to_at(irow) .gt. nat2)) cycle
      call MatGetRow(p_A, irow - 1, nzcol, icol, PETSC_NULL_SCALAR, ierr)
      do inz = 1, nzcol
        if ((imu_to_at(icol(inz) + 1) .lt. nat1) .or. (imu_to_at(icol(inz) + 1) .gt. nat2)) icol(inz) = -1 ! negative coloumns are ignored by MatSetValues
      end do
      ii(1) = irow - 1 ! back to C numbering
      call MatSetValues(p_tmp, 1, ii, nzcol, icol(1:nzcol), p_vals(1:nzcol), INSERT_VALUES, ierr)
      call MatRestoreRow(p_A, irow - 1, nzcol, icol, PETSC_NULL_SCALAR, ierr)
    end do
    call MatAssemblyBegin(p_tmp, MAT_FINAL_ASSEMBLY, ierr)
    call MatAssemblyEnd(p_tmp, MAT_FINAL_ASSEMBLY, ierr)
    call MatDuplicate(p_tmp, MAT_COPY_VALUES, p_B, ierr)
    call MatDestroy(p_tmp, ierr)

  end subroutine init_ep_active_mat

  subroutine get_dS(p_tmp_dS, p_dS, i_d_at, imu_to_at, ij)
    ! Purpose: Extracts a submatrix from a larger matrix based on atom index and selection criteria.
    !
    ! Functionality: This subroutine extracts a submatrix p_dS from a temporary matrix p_tmp_dS. 
    ! The selection is based on the atom index i_d_at and a selection flag ij.  If ij=1, 
    ! it selects rows corresponding to i_d_at; if ij=2, it selects columns corresponding 
    ! to i_d_at.  It uses the imu_to_at array for mapping between indices.
    !
    ! Parameters:
    !   p_tmp_dS (Mat, input): The input temporary matrix.
    !   p_dS (Mat, output): The output submatrix.
    !   i_d_at (integer, input): The index of the atom to select.
    !   imu_to_at (integer array, input): Mapping between internal indices and atom indices.
    !   ij (integer, input): Selection flag (1 for rows, 2 for columns).

#include <petsc/finclude/petsc.h>
    use petscmat
    use petsc_mod
    use petsc_wrapper
    use kinds
    use globals, only: nmu_c, mattype_sparse
    use error_handler

    implicit none

    Mat :: p_tmp_dS, p_dS
    integer :: i_d_at, ij
    integer, allocatable :: imu_to_at(:)

    integer :: iat1, iat2, nl1, nl2, irow, ierr, nzcol, inz, ii(1)
    integer, allocatable :: icol(:)
    PetscScalar, allocatable :: p_vals(:)

    allocate (icol(nmu_c), p_vals(nmu_c), stat=ierr)
    if (ierr .ne. 0) then
      write (errormsg, fmt='(A,i8)') "allocation error icol ", ierr
      call error()
    end if

    call MatGetOwnershipRange(p_tmp_dS, nl1, nl2, ierr)
    nl1 = nl1 + 1 ! to fortran numbering, nl2 is already correct as it is max_number of row+1

    do irow = nl1, nl2
      if ((ij .eq. 1) .and. (imu_to_at(irow) .ne. i_d_at)) cycle
      call MatGetRow(p_tmp_dS, irow - 1, nzcol, icol, PETSC_NULL_SCALAR, ierr)
      do inz = 1, nzcol
        if ((ij .eq. 2) .and. (imu_to_at(icol(inz) + 1) .ne. i_d_at)) icol(inz) = -1 ! negative coloumns are ignored by MatSetValues
      end do
      ii(1) = irow - 1 ! back to C numbering
      call MatSetValues(p_dS, 1, ii, nzcol, icol(1:nzcol), p_vals(1:nzcol), INSERT_VALUES, ierr)
      call MatRestoreRow(p_tmp_dS, irow - 1, nzcol, icol, PETSC_NULL_SCALAR, ierr)
    end do
    call MatAssemblyBegin(p_dS, MAT_FINAL_ASSEMBLY, ierr)
    call MatAssemblyEnd(p_dS, MAT_FINAL_ASSEMBLY, ierr)

  end subroutine get_dS

end module init_ep_coupling