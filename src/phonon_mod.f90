module phonon_mod
  ! Module for phonon calculations.
  #include <petsc/finclude/petsc.h>
  use petscmat
  implicit none

  ! Integer variable for mode index.
  integer :: ii_mode
  ! PetscScalar variable for mode energy.
  PetscScalar :: mode_ii_energy

contains

  ! Subroutine to test the adaptive integrator with PETSc matrices.
  subroutine integrator_test()
    ! Tests the adaptive integrator using a 1x1 matrix.
    #include <petsc/finclude/petsc.h>
    use petscmat
    use petsc_mod
    use petsc_wrapper
    use kinds
    use globals
    use integrator_mod
    implicit none

    ! Real variables for integration limits, error tolerance, and integration error.
    real(dp) :: d_low, d_high, eps, dd_eps, err_eq_int
    ! PetscScalar variables for integration results and matrix values.
    PetscScalar :: p_zint1, p_zint2, pp(1), dx, x1, x2, x
    ! PETSc matrix variables.
    Mat :: p_d_tmp1(1), p_d_tmp2(1)
    ! Integer variables for error codes and indices.
    integer :: ierr, ii(1), jj(1), ik(2), i, n

    ! Create a 1x1 sparse matrix.
    call MatCreate(PETSC_COMM_WORLD, p_tmpcc1, ierr)
    call MatSetType(p_tmpcc1, mattype_sparse, ierr)
    call MatSetSizes(p_tmpcc1, PETSC_DECIDE, PETSC_DECIDE, 1, 1, ierr)
    ! Set matrix values based on processor rank.
    if (inode .eq. 0) then
      ik(1) = 0
      ik(2) = 1
      pp = 0d0
      call MatMPIAIJSetPreallocationCSR(p_tmpcc1, ik, ik, pp, ierr)
    else
      call MatMPIAIJSetPreallocationCSR(p_tmpcc1, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, pp, ierr)
    end if

    ! Assemble the matrix.
    call MatAssemblyBegin(p_tmpcc1, MAT_FINAL_ASSEMBLY, ierr)
    call MatAssemblyEnd(p_tmpcc1, MAT_FINAL_ASSEMBLY, ierr)
    ! Duplicate the matrix.
    call MatDuplicate(p_tmpcc1, MAT_DO_NOT_COPY_VALUES, p_d_tmp1(1), ierr)
    call MatDuplicate(p_tmpcc1, MAT_DO_NOT_COPY_VALUES, p_d_tmp2(1), ierr)

    ! Assemble the matrix again.
    call MatAssemblyBegin(p_tmpcc1, MAT_FINAL_ASSEMBLY, ierr)
    call MatAssemblyEnd(p_tmpcc1, MAT_FINAL_ASSEMBLY, ierr)
    ! Print matrix information.
    call petsc_mat_info(p_tmpcc1, "p_tmpcc1 ", ierr)

    ! Set integration parameters.
    eps = 1d-12
    dd_eps = 1d-8
    d_high = zsqrt(2d0*mode_ii_energy/dd_eps + mode_ii_energy*mode_ii_energy + eta_ph_cc*eta_ph_cc*0.25d0)
    d_low = -d_high
    nint_order = 64
    l_output_progress = .false.
    ! Perform adaptive integration.
    call adaptive_int3(dr_alpha_petsc, d_low, d_high, p_d_tmp1(1), p_d_tmp2(1), maxsub*1000, eps, err_eq_int)
    ! Get the integration result.
    call petsc_mat_getvalue(0, 0, p_zint1, p_d_tmp1(1), 1, PETSC_COMM_WORLD)

    ! Print the integration result.
    write (6, *) "p_zint1 ", p_zint1
    ! Destroy the matrix.
    call MatDestroy(p_tmpcc1, ierr)

  end subroutine integrator_test

  ! Function to calculate the equilibrium density of states.
  function rho_eq_alpha(x)
    ! Calculates the equilibrium density of states.
    use globals
    use kinds
    implicit none

    ! Real input variable for energy.
    real(dp) :: x
    ! Complex output variable for density of states.
    complex(dp) :: rho_eq_alpha

    ! Calculate the density of states.
    rho_eq_alpha = 1d0/pi*(eta_ph_cc*0.5d0/((x - mode_ii_energy)**2 + eta_ph_cc**2/4d0) - eta_ph_cc*0.5d0/((x + mode_ii_energy)**2 + eta_ph_cc**2/4d0))

  end function rho_eq_alpha

  ! Function to calculate the derivative of alpha.
  function dr_alpha(x)
    ! Calculates the derivative of alpha.
    use globals
    use kinds
    implicit none

    ! Real input variable for energy.
    real(dp) :: x
    ! Complex output variable for derivative of alpha.
    complex(dp) :: dr_alpha

    ! Calculate the derivative of alpha.
    dr_alpha = 2d0*mode_ii_energy/(x*x - mode_ii_energy*mode_ii_energy + zione*eta_ph_cc*x - eta_ph_cc*eta_ph_cc*0.25d0)

  end function dr_alpha

  ! Function to calculate the derivative of alpha using PETSc matrices.
  function dr_alpha_petsc(x)
    ! Calculates the derivative of alpha using PETSc matrices.
    use globals
    use kinds
    use petsc_wrapper, only: petsc_matassemble
    use petsc_mod
    implicit none

    ! Real input variable for energy.
    real(dp) :: x

    ! Real variable for temporary storage.
    real(dp) :: z
    ! Integer variables for error codes and indices.
    integer :: ii(1), ierr, dr_alpha_petsc
    ! PetscScalar variable for matrix value.
    PetscScalar :: pp(1)

    ! Set default return value.
    dr_alpha_petsc = -1

    ! Set matrix value on rank 0.
    if (inode .eq. 0) then
      z = x
      pp(1) = dr_alpha(z)
      ii = 0
      call MatSetValues(p_tmpcc1, 1, ii, 1, ii, pp, INSERT_VALUES, ierr)
    end if

    ! Assemble the matrix.
    call MatAssemblyBegin(p_tmpcc1, MAT_FINAL_ASSEMBLY, ierr)
    call MatAssemblyEnd(p_tmpcc1, MAT_FINAL_ASSEMBLY, ierr)

    ! Assemble the matrix using PETSc wrapper.
    call petsc_matassemble(p_tmpcc1)

    ! Set return value to 0.
    dr_alpha_petsc = 0

  end function dr_alpha_petsc

  ! Subroutine to initialize phonon data.
  subroutine init_phonons()
    ! Initializes phonon-related data by reading from binary files.
    #include <petsc/finclude/petsc.h>
    use petscmat
    use petsc_mod
    use petsc_wrapper
    use kinds
    use globals
    use misc

    implicit none

    ! Integer variables for error codes, indices, and file unit.
    integer :: ierr, i1, i2, iunit, ii(1)
    ! PETSc viewer for reading binary files.
    PetscViewer :: v_infile
    ! Character variable for file name.
    character(strln) :: infile, str_i1, str_i2
    ! PetscScalar variable for reading data.
    PetscScalar :: p_in(1)

    ! Allocate memory for phonon dynamical matrix.
    allocate (p_Kphonon00_cc(-ncell_c(1):ncell_c(1), -ncell_c(2):ncell_c(2)), &
              stat=ierr)

    ! Create PETSc viewer for binary files.
    call PetscViewerCreate(PETSC_COMM_WORLD, v_infile, ierr)
    call PetscViewerSetType(v_infile, PETSCVIEWERBINARY, ierr)
    call PetscViewerFileSetMode(v_infile, FILE_MODE_READ, ierr)

    ! Read dynamical matrices from binary files.
    write (pstr_out, fmt='(A)') "read dynamical matrix "; call petsc_print_master()
    do i1 = -ncell_c(1), ncell_c(1)
      do i2 = -ncell_c(2), ncell_c(2)
        str_i1 = int2str(i1)
        str_i2 = int2str(i2)
        infile = "FCmat_"//trim(str_i1)//"_"//trim(str_i2)//"_0_petsc.dat"
        call PetscViewerFileSetName(v_infile, trim(infile), ierr)
        call MatCreate(PETSC_COMM_WORLD, p_Kphonon00_cc(i1, i2), ierr)
        call MatSetType(p_Kphonon00_cc(i1, i2), mattype_sparse, ierr)
        call MatLoad(p_Kphonon00_cc(i1, i2), v_infile, ierr)
      end do
    end do

    ! Read inverse square root of mass matrix from binary file.
    write (pstr_out, fmt='(A)') "read MassMatrix^(-1/2) "; call petsc_print_master()
    infile = "sqrtMinv_petsc.dat"
    call PetscViewerFileSetName(v_infile, trim(infile), ierr)
    call MatCreate(PETSC_COMM_WORLD, p_invsqrt_mass, ierr)
    call MatSetType(p_invsqrt_mass, mattype_sparse, ierr)
    call MatLoad(p_invsqrt_mass, v_infile, ierr)

    ! Read mass matrix from binary file.
    write (pstr_out, fmt='(A)') "read MassMatrix "; call petsc_print_master()
    infile = "Massmat_petsc.dat"
    call PetscViewerFileSetName(v_infile, trim(infile), ierr)
    call MatCreate(PETSC_COMM_WORLD, p_mass_matrix, ierr)
    call MatSetType(p_mass_matrix, mattype_sparse, ierr)
    call MatLoad(p_mass_matrix, v_infile, ierr)

    ! Destroy PETSc viewer.
    call PetscViewerDestroy(v_infile, ierr)

  end subroutine init_phonons

  ! Subroutine to calculate phonon properties.
  subroutine get_phonones()
    ! Calculates phonon eigenenergies and eigenvectors.
    #include <petsc/finclude/petsc.h>
    use petscmat
    use slepceps
    use kinds
    use petsc_mod
    use ft_mod
    use globals
    use slepc_mod

    implicit none

    ! PETSc matrix variables.
    Mat :: p_tmp, p_tmp1, p_tmp2, p_tmp3
    ! Integer variables for error codes and indices.
    integer :: ia1, j, i, ii, ierr, i1, i2, jj, iroot
    ! Integer variables for matrix dimensions.
    PetscInt :: nl1, nl2
    ! Real variable for matrix norm.
    PetscReal :: norm
    ! PetscScalar array for phonon eigenenergies.
    PetscScalar, pointer :: p_x(:)

    ! Perform Fourier transform for each k-point.
    do i1 = 1, nk_c
      call MatDuplicate(p_Kphonon00_cc(0, 0), MAT_DO_NOT_COPY_VALUES, p_Kphonon00k_cc(i1), ierr)
      call MatCreateVecs(p_Kphonon00k_cc(i1), p_phonon_EW(i1), PETSC_NULL_vec, ierr)
      call fourier_trans(p_Kphonon00_cc, p_Kphonon00k_cc, kp_r, nk_r, dlat_r, ncell_c(1), ncell_c(2), 0)
    end do

    ! Set the number of phonon modes.
    n_ep_modes_k = nat_ecc*3
    ! Calculate phonon eigenenergies and eigenvectors for each k-point.
    do i1 = 1, nk_c
      call petsc_mat_info(p_Kphonon00k_cc(i1), "p_Kphonon00k_cc(i1) ", ierr)
      call MatPtAP(p_Kphonon00k_cc(i1), p_invsqrt_mass, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, p_tmp, ierr)
      call petsc_get_densemat(p_Kphonon00k_cc(i1), p_tmp2, mattype_dense)
      call diag_mat2(p_tmp, PETSC_NULL_MAT, p_phonon_EW(i1), p_tmp2, nat3, -1d0, EPS_LARGEST_REAL)
      call MatGetOwnershipRange(p_tmp, nl1, nl2, ierr)
      call VecGetArrayF90(p_phonon_EW(i1), p_x, ierr)
      p_x = zsqrt(p_x)
      write (pstr_out, fmt='(A8,A8,A48,A48)') "k-point", "Mode", "Energy(H)", "Energy(eV)"; call petsc_print_master()

      ! Print phonon eigenenergies.
      do j = 0, nprocs - 1
        if (j .eq. inode) then
          i = 0
          do i2 = nl1, nl2 - 1
            i = i + 1
            ii = i2
            jj = 0
            iroot = -1
            if ((real(p_x(i)) .le. 1d-8) .and. (aimag(p_x(i)) .le. 1d-12) .and.&
            &(n_ep_modes_k(i1) .eq. nat3)) n_ep_modes_k(i1) = i2
            write (6, fmt='(2i8,4e24.12E4)') i1, i2, p_x(i), p_x(i)*27.2114d0
            call flush (6)
          end do
        end if
        call MPI_BARRIER(PETSC_COMM_WORLD, ierr)
      end do

      ! Determine the number of non-zero modes.
      call MPI_BARRIER(PETSC_COMM_WORLD, ierr)
      call MPI_AllReduce(MPI_IN_PLACE, n_ep_modes_k(i1), 1, MPI_INTEGER, MPI_MIN, PETSC_COMM_WORLD, ierr)

      call MPI_BARRIER(PETSC_COMM_WORLD, ierr)
      write (pstr_out, fmt='(A,i8)') "number of nonzero modes=", n_ep_modes_k(i1); call petsc_print_master()
      call VecRestoreArrayF90(p_phonon_EW(i1), p_x, ierr)
      call MatDestroy(p_tmp, ierr)
      call MatMatMult(p_invsqrt_mass, p_tmp2, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, p_phonon_EV(i1), ierr)
      call MatDestroy(p_tmp2, ierr)
      call MatHermitianTranspose(p_phonon_EV(i1), MAT_INITIAL_MATRIX, p_tmp2, ierr)
      call MatMatMatMult(p_tmp2, p_mass_matrix, p_phonon_EV(i1), MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, p_tmp3, ierr)
      call MatShift(p_tmp3, p_minus1, ierr)
      call MatNorm(p_tmp3, NORM_FROBENIUS, norm, ierr)
      write (pstr_out, fmt='(A,e24.12)') "norm(A'*M*A-1)=", norm; call petsc_print_master()

      call MatDestroy(p_tmp3, ierr)
      call MatDestroy(p_tmp2, ierr)
      call MatDestroy(p_Kphonon00k_cc(i1), ierr)

    end do
    ! Determine the maximum number of non-zero modes.
    n_ep_modes_max = maxval(n_ep_modes_k)

    nullify (p_x)

  end subroutine get_phonones
end module phonon_mod