module kmat_from_diag
  ! This module contains subroutines for constructing the k-space matrix from diagonal elements.

  implicit none
  integer :: Ntot_diag ! Total number of diagonal elements considered.

  contains

    subroutine get_kmat_from_diag()
      ! This subroutine constructs the k-space matrix from the diagonal elements of the Hamiltonian.
      ! It calculates the Fermi energy, constructs the density matrix in k-space, and transforms it to real space.

#include <petsc/finclude/petsc.h>
      use slepceps
      use petsc
      use petsc_mod
      use petsc_wrapper
      use k_on_demand_mod
      use kinds
      use globals
      use ft_mod
      use slepc_mod
      use misc
      use error_handler

      implicit none

      ! PETSc objects
      Mat :: p_mat_coeff, p_tmp1, p_tmp2  ! Matrices for coefficients, temporary storage.
      Mat, allocatable :: p_mat_tmp(:)      ! Allocatable array of matrices for temporary storage.
      Vec, allocatable :: p_vec_ew(:)       ! Allocatable array of vectors for eigenvalues.
      Vec, pointer :: p_vec_stride(:)       ! Pointer to a vector for stride.
      Vec :: p_vec_ew_all                 ! Vector for all eigenvalues.
      ! Integer variables
      integer :: ierr, ik, nrow, ncol, n, nloc, nloc1, nloc2, i, j, N_ew_total, m, k, nz_proc, &
        ii(1), jj(1), i1, i2, iunit, ncols_nz, ncols_nz_max, i_max_ew, N_diag, nn, iamk, iew ! Various integer variables for indexing, error handling, and counting.
      integer, allocatable :: i_rows(:), i_rows_shift(:), idx(:), displs(:), counts(:), &
        cols_nz(:,:), i_max_ewk(:), i_row_nz(:) ! Allocatable integer arrays for indexing and storing data.
      ! Floating-point variables
      PetscScalar :: z_ew_n, c_jn, fac, p11, d11 ! Various complex and real variables for calculations.
      PetscScalar, pointer :: p_x1(:),p_x2(:), pp_coeff(:,:) ! Pointers to arrays for eigenvalues and matrix coefficients.
      PetscReal, allocatable :: ew(:) ! Allocatable array for eigenvalues.
      PetscScalar, allocatable :: coeff_batch(:), coeff_batch_j(:), cols_buf(:), local_Dmat(:) ! Allocatable arrays for coefficients and temporary storage.
      PetscReal :: N_elec, N_elec_target, Ef, E1, E2, N_left, fermi_fac, ew_n, norm, Ntrace ! Various real variables for electron count, Fermi energy, and other calculations.
      ! Integer(8) variables
      integer(8) :: counti, count_rate, countf, ipos, counter(3,256) ! Integer(8) variables for timing and counting.
      ! Character variable
      character(strln) :: file_ev, outstr ! Character variable for file name.

      counter = 0

      ! Get matrix size and ownership range
      call MatGetSize(p_h00_ik_cc(1), nrow, ncol, ierr)
      call MatGetOwnershipRange(p_h00_ik_cc(1), nloc1, nloc2, ierr)
      nloc = nloc2 - nloc1

      ! Allocate memory
      allocate(i_rows(nloc), i_rows_shift(nloc), coeff_batch(nloc), &
        coeff_batch_j(ncol), displs(nprocs), counts(nprocs), stat = ierr)

      ! Set counts and displacements for MPI_Allgather
      do i = 1, nloc
        i_rows(i) = nloc1 + i - 1
      end do

      counts(inode+1) = nloc
      displs(inode+1) = nloc1


      ! Gather counts and displacements from all processes
      call MPI_ALLGATHER(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, counts, 1, MPI_INTEGER, PETSC_COMM_WORLD, ierr)
      call MPI_ALLGATHER(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, displs, 1, MPI_INTEGER, PETSC_COMM_WORLD, ierr)

      ! Calculate total number of eigenvalues
      N_ew_total = nk_l * nrow
      allocate(p_vec_ew(nk_l), ew(N_ew_total), stat = ierr)

      ! Allocate temporary matrix
      allocate(p_mat_tmp(1), stat = ierr)

      ! Duplicate the Hamiltonian matrix
      call MatDuplicate(p_h00_ik_cc(1), MAT_DO_NOT_COPY_VALUES, p_mat_tmp(1), ierr)
      ! Set matrix options
      call MatSetOption(p_mat_tmp(1), MAT_NEW_NONZERO_LOCATION_ERR, PETSC_TRUE, ierr)
      call MatSetOption(p_mat_tmp(1), MAT_IGNORE_OFF_PROC_ENTRIES, PETSC_FALSE, ierr)
      call MatSetOption(p_mat_tmp(1), MAT_NO_OFF_PROC_ZERO_ROWS, PETSC_FALSE, ierr)

      ! Get dense matrix representation
      call petsc_get_densemat(p_mat_tmp(1), p_mat_coeff, mattype_dense)

      ! Set target number of electrons and diagonal elements
      N_elec_target = n_electrons_ecc(1)
      N_diag = ceiling(N_elec_target * 1.4d0)
      N_diag = min(N_diag, nrow)
      Ntot_diag = N_diag * nk_l

      ! Loop over k-points
      do ik = 1, nk_c
        ! Measure time
        call system_clock(counti, count_rate)
        ! Print k-point information
        write(pstr_out,fmt='(A,i8,3e24.12,i4)') "# k: ", nk_c-ik, kp_c(1:3,ik), wkp_c(ik)
        call petsc_print_master(.false.)

        ! If k-points are on demand, call k_on_demand subroutine
        if (l_k_on_demand) call k_on_demand(ik, 3)

        ! Create vector for eigenvalues
        call MatCreateVecs(p_h00k_cc(ik), p_vec_ew(ik), PETSC_NULL_VEC, ierr)
        ! Shift eigenvalues to avoid numerical issues
        call VecShift(p_vec_ew(ik), dcmplx(huge(1d0), 1d-64), ierr)
        ! Zero the coefficient matrix
        call MatZeroEntries(p_mat_coeff, ierr)
        ! Measure time for diagonalization
        counter(3,6) = 0
        call system_clock(counter(1,6), count_rate)
        ! Diagonalize the Hamiltonian matrix
        call diag_mat2(p_h00k_cc(ik), p_s00k_cc(ik), p_vec_ew(ik), p_mat_coeff, N_diag, eps_which_in=EPS_SMALLEST_REAL)
        call system_clock(counter(2,6))
        counter(3,6) = counter(3,6) + (counter(2,6) - counter(1,6))

        ! Save eigenvalues to file
        write(file_ev,fmt='(A,i0,A)') "ev_mat_", ik,".dat"
        call petsc_mat_direct_save(p_mat_coeff, file_ev, ierr)

        ! Print time information
        call system_clock(countf)
        write (pstr_out, fmt='(X,2e24.12)') real(countf - counti, 8)/real(count_rate, 8), &
          real(counter(3,6), 8)/real(count_rate, 8)
        call petsc_print_master()
      end do

      ! Print separator
      write (pstr_out, fmt='(A)') repeat("-",80)
      call petsc_print_master()

      ! Create vector for all eigenvalues
      n = 0
      if (inode .eq. 0) n = N_ew_total
      call VecCreateMPI(PETSC_COMM_WORLD, n, N_ew_total, p_vec_ew_all, ierr)

      ! Get array from vector
      call VecGetArrayF90(p_vec_ew_all, p_x2,ierr)

      ! Gather eigenvalues from all processes
      do ik = 1, nk_l

        call VecGetArrayF90(p_vec_ew(ik), p_x1, ierr)
        i_rows_shift = i_rows + nrow * (ik-1)
        call VecSetValues(p_vec_ew_all, nloc, i_rows_shift, p_x1, INSERT_VALUES, ierr)
        call VecRestoreArrayF90(p_vec_ew(ik), p_x1,ierr)

      end do
      call VecAssemblyBegin(p_vec_ew_all, ierr)
      call VecAssemblyEnd(p_vec_ew_all, ierr)

      ! Check for imaginary part of eigenvalues
      if (maxval(abs(aimag(p_x2))).ge.1d-14) then
        write (pstr_out, fmt='(A, 2e24.12)') "warning max(abs(imag(ew))) ", &
          maxval(abs(aimag(p_x2)))
        call petsc_print_master()
      end if

      ! Allocate array for indices
      allocate(idx(N_ew_total), stat = ierr)

      ! Sort eigenvalues and indices if on root process
      if (inode.eq.0) then

        ew=real(p_x2, dp)

        allocate(i_max_ewk(nk_l), stat = ierr)

        do i = 1, N_ew_total
          idx(i) = i - 1
        end do

        call PetscSortRealWithPermutation(N_ew_total, ew, idx, ierr)
        open (newunit = j, file = "eigall.dat", action="write", status = "replace")
        do i = 1, N_ew_total
          write(unit=j,fmt='(2i8,2e24.12,i8)') i, idx(i), p_x2(idx(i)+1),Ntot_diag
        end do
        close(j)
        ! Determine Fermi energy
        if (l_diag_fix_ef) then
          Ef = (ef_l + ef_r) * 0.5d0
        else
          E2 = ew(idx(Ntot_diag)+1)
          E1 = ew(1)
          Ef = ew(idx(Ntot_diag)*1) 
          k = 0

          do
            k=k+1
            N_elec = 0d0

            call Nelec_min_fct(Ef, N_left, N_elec_target, ew, wkp_l, idx, temperature_el)

            if (abs(N_left).le.1d-14) exit

            if (N_left.gt.0d0) then
              E1 = Ef
              Ef = (E2 + E1) * 0.5d0
            else
              E2 = Ef
              Ef = (E2 + E1) * 0.5d0
            end if
            if (k.ge.1000) exit
          end do
        end if

        ! Calculate number of electrons per k-point
        call Nelec_min_fct(Ef, N_left, 0d0, ew, wkp_l, idx, temperature_el, i_max_ewk)
        i_max_ew = maxval(i_max_ewk)
        i_max_ew = ceiling(i_max_ew * 1.1d0)
      end if

      ! Broadcast Fermi energy, electron count, and maximum number of eigenvalues
      call VecRestoreArrayF90(p_vec_ew_all, p_x2,ierr)

      call MPI_Bcast(Ef, 1, MPI_DOUBLE_PRECISION, 0, PETSC_COMM_WORLD, ierr)
      call MPI_Bcast(N_left, 1, MPI_DOUBLE_PRECISION, 0, PETSC_COMM_WORLD, ierr)
      call MPI_Bcast(i_max_ew, 1, MPI_INTEGER, 0, PETSC_COMM_WORLD, ierr)
      call MPI_Bcast(idx, N_ew_total, MPI_INTEGER, 0, PETSC_COMM_WORLD, ierr)


      ! Print Fermi energy and electron count
      write (pstr_out, fmt='(A, 2e24.12, i8)') "Fermi energy ", Ef, -N_left, i_max_ew
      call petsc_print_master()
      write (pstr_out, fmt='(A)') repeat("-",80)
      call petsc_print_master()

      ! Zero density matrices
      do i2 = -ncell_c(2), ncell_c(2)
        do i1 = -ncell_c(1), ncell_c(1)
          if (diag_dim.eq.3) then
            call MatZeroEntries(p_dmatxy(i1, i2, -1), ierr)
            call MatZeroEntries(p_dmatxy(i1, i2, 1), ierr)
          end if
          call MatZeroEntries(p_dmatxy(i1, i2, 0), ierr)
        end do
      end do

      ! Measure time for getting row information
      call system_clock(counti, count_rate)
      ncols_nz_max = 0
      do i = 1, nloc
        call MatGetRow(p_mat_tmp(1), i_rows(i), ncols_nz, PETSC_NULL_INTEGER, &
          PETSC_NULL_SCALAR, ierr)
        ncols_nz_max = max(ncols_nz_max, ncols_nz)
        call MatRestoreRow(p_mat_tmp(1), i_rows(i), ncols_nz, PETSC_NULL_INTEGER, &
          PETSC_NULL_SCALAR, ierr)
      end do

      ! Allocate memory for non-zero column indices
      allocate(cols_nz(ncols_nz_max, nloc), cols_buf(ncols_nz_max), &
        i_row_nz(nloc), stat = ierr)

      ! Get non-zero column indices for each row
      do i = 1, nloc
        call MatGetRow(p_mat_tmp(1), i_rows(i), i_row_nz(i), cols_nz(:,i), &
          PETSC_NULL_SCALAR, ierr)
        call MatRestoreRow(p_mat_tmp(1), i_rows(i), ncols_nz, PETSC_NULL_INTEGER, &
          PETSC_NULL_SCALAR, ierr)
      end do

      ! Print time information
      call system_clock(countf)
      write (pstr_out, fmt='(A,i8,e24.12)') "ncols_nz_max ", ncols_nz_max, &
        real(countf - counti, 8)/real(count_rate, 8)
      call petsc_print_master()
      write (pstr_out, fmt='(A)') repeat("-",80)
      call petsc_print_master()


      ! Initialize trace
      Ntrace = 0d0
      call MatDenseGetArrayF90(p_mat_coeff,pp_coeff, ierr)

      ! Open file for eigenvalues
      if (inode.eq.0) open(newunit=iunit,file="eigs.dat",action="write",status="replace")

      ! Set matrix options
      call MatSetOption(p_mat_tmp(1), MAT_NO_OFF_PROC_ENTRIES, PETSC_TRUE, ierr)
      call MatSetOption(p_mat_tmp(1), MAT_NEW_NONZERO_LOCATION_ERR , PETSC_TRUE, ierr)

      ! Calculate total number of non-zero elements
      nz_proc = sum(i_row_nz(1:nloc))
      allocate(local_Dmat(nz_proc))

      ! Loop over k-points
      do ik = 1, nk_c
        ! Initialize counter
        counter = 0d0
        ! Measure time
        call system_clock(counti, count_rate)
        ! Print k-point information
        write(pstr_out,fmt='(A,i8,3e24.12,i4)') "# k: ", nk_c-ik, kp_c(1:3,ik), wkp_c(ik)
        if (inode.eq.0) write(iunit,fmt='(A)') trim(pstr_out)
        call petsc_print_master(.false.)


        ! Load eigenvalues from file
        call MatZeroEntries(p_mat_coeff, ierr)
        write(file_ev,fmt='(A,i0,A)') "ev_mat_", ik,".dat"
        call system_clock(counter(1,5), count_rate)
        call petsc_mat_direct_load(p_mat_coeff, file_ev, ierr)
        call system_clock(counter(2,5))
        counter(3,5) = counter(3,5) + (counter(2,5) - counter(1,5))

        ! Initialize local density matrix
        local_Dmat = 0d0

        ! Zero temporary matrix
        call MatZeroEntries(p_mat_tmp(1), ierr)
        p11 = 0d0
        iew = 0

        ! Loop over eigenvalues
        do nn = 1, N_ew_total
          n = idx(nn) + 1

          iamk = ceiling(1d0/real(nrow,dp)*(n)-1d0/(real(nrow,dp)*2d0))

          if (ik.ne.iamk) cycle

          n = n - (nrow * (ik - 1))

          iew = iew + 1
          call petsc_vec_getvalue(n - 1, z_ew_n, p_vec_ew(ik), .true.)
          ew_n = z_ew_n
          fermi_fac = fermi_function(ew_n, Ef, temperature_el)
          if (inode.eq.0) write(iunit,fmt='(i8,2e24.12,E24.12E4)') iew, z_ew_n, fermi_fac
          if (fermi_fac .le. 1d-13) cycle
          call system_clock(counter(1,1), count_rate)
          coeff_batch_j = 0d0
          call MPI_Allgatherv(pp_coeff(:,n), nloc, MPI_DOUBLE_COMPLEX, coeff_batch_j, &
            counts, displs, MPI_DOUBLE_COMPLEX, PETSC_COMM_WORLD, ierr)
          call system_clock(counter(2,1))
          counter(3,1) = counter(3,1) + (counter(2,1) - counter(1,1))

          i1 = 0
          do i = 1, nloc
            call system_clock(counter(1,2), count_rate)
            do j = 1, i_row_nz(i)
              i1 = i1 + 1
              k = cols_nz(j,i) + 1
              c_jn = coeff_batch_j(k)
              cols_buf(j) = pp_coeff(i,n) * conjg(c_jn) * fermi_fac
              local_Dmat(i1) = local_Dmat(i1) + cols_buf(j)
            end do
            call system_clock(counter(2,2))
            counter(3,2) = counter(3,2) + (counter(2,2) - counter(1,2))
          end do
        end do

        ! Set values in temporary matrix
        call system_clock(counter(1,3), count_rate)
        i1 = 1
        do i = 1, nloc
          i2 = i1 + i_row_nz(i) - 1
          call MatSetValues(p_mat_tmp(1), 1, i_rows(i), i_row_nz(i), cols_nz(1:i_row_nz(i),i), &
            local_Dmat(i1:i2), INSERT_VALUES, ierr)
          i1 = i2 + 1
        end do
        call system_clock(counter(2,3))
        counter(3,3) = counter(3,3) +  (counter(2,3) - counter(1,3))

        ! Assemble temporary matrix
        call system_clock(counter(1,4), count_rate)
        call MatAssemblyBegin(p_mat_tmp(1), MAT_FINAL_ASSEMBLY, ierr)
        call MatAssemblyEnd(p_mat_tmp(1), MAT_FINAL_ASSEMBLY, ierr)
        call system_clock(counter(2,4))
        counter(3,4) = counter(3,4) + (counter(2,4) - counter(1,4))


        ! Perform inverse Fourier transform
        call fourier3d_back_add(p_mat_tmp, p_dmatxy(:, :, :), kp_c(1:3, ik:ik), &
          wkp_c(ik:ik), 1, dlat_c, ncell_c(1), ncell_c(2), ncell_c(3), 0, 0, PETSC_FALSE)

        ! Print time information
        call system_clock(countf)
        write (pstr_out, fmt='(X,e24.12)') real(countf - counti, 8)/real(count_rate, 8); call petsc_print_master()
        write (pstr_out, fmt='(X,e24.12)') real(counter(3,1), 8)/real(count_rate, 8); call petsc_print_master()
        write (pstr_out, fmt='(X,e24.12)') real(counter(3,2), 8)/real(count_rate, 8); call petsc_print_master()
        write (pstr_out, fmt='(X,e24.12)') real(counter(3,3), 8)/real(count_rate, 8); call petsc_print_master()
        write (pstr_out, fmt='(X,e24.12)') real(counter(3,4), 8)/real(count_rate, 8); call petsc_print_master()
        write (pstr_out, fmt='(X,e24.12)') real(counter(3,5), 8)/real(count_rate, 8); call petsc_print_master()

      end do
      ! Reset matrix option
      call MatSetOption(p_mat_tmp(1), MAT_NO_OFF_PROC_ENTRIES, PETSC_FALSE, ierr)
      if (inode.eq.0) close(iunit)
      call MatDenseRestoreArrayF90(p_mat_coeff,pp_coeff, ierr)

      ! Scale density matrices
      fac = 1d0/(nktot)
      Ntrace = Ntrace * fac
      do i2 = -ncell_c(2), ncell_c(2)
        do i1 = -ncell_c(1), ncell_c(1)
          if (diag_dim .eq. 3) then
            call MatScale(p_dmatxy(i1, i2, -1), fac, ierr)
            call MatScale(p_dmatxy(i1, i2, 1), fac, ierr)
          end if
          call MatScale(p_dmatxy(i1, i2, 0), fac, ierr)
        end do
      end do

      ! Destroy vectors
      do ik = 1, nk_l
        call VecDestroy(p_vec_ew(ik), ierr)
      end do

      deallocate(p_vec_ew, stat = ierr)

      call VecDestroy(p_vec_ew_all, ierr)

      ! Destroy matrices
      call MatDestroy(p_mat_coeff, ierr)
      call MatDestroy(p_mat_tmp(1), ierr)

      ! Save Fermi energy to file using MPI I/O
      call MPI_FILE_OPEN(MPI_COMM_WORLD, "sys_info.dat", MPI_MODE_WRONLY + MPI_MODE_CREATE, &
        MPI_INFO_NULL, iunit, ierr)

      if (l_ionode) then
        ipos = 3*3*24 + 3 + 3 * 22 + 1 + 22 + 1
        write(outstr,fmt='(e24.12, X, l, A)') Ef, .true., NEW_LINE('A')
        call MPI_FILE_WRITE_AT(iunit, ipos, trim(outstr), len(trim(outstr)), &
          MPI_CHARACTER, MPI_STATUS_IGNORE, ierr)
      end if


      call MPI_FILE_CLOSE(iunit, ierr)

    end subroutine get_kmat_from_diag

    function fermi_function(x,Ef,T)
      ! This function calculates the Fermi-Dirac distribution function.
      ! x: Energy
      ! Ef: Fermi energy
      ! T: Temperature

      use globals, only : kb
      use kinds
      implicit none

      real(dp) :: fermi_function
      real(dp), intent(in) :: x, Ef, T

      real(dp) :: beta, xx

      beta = 1d0/(T * kb)
      xx = (x-Ef)*beta
      fermi_function = 1d0 / ( exp(xx) + 1d0)

    end function fermi_function

    subroutine Nelec_min_fct(x, f, N, ew, wk, idx, T, i_max_ewk)
      ! This subroutine calculates the number of electrons for a given Fermi energy.
      ! x: Fermi energy
      ! f: Difference between target and calculated electron number
      ! N: Target number of electrons
      ! ew: Eigenvalues
      ! wk: k-point weights
      ! idx: Sorted indices of eigenvalues
      ! T: Temperature
      ! i_max_ewk: (Optional) Array to store the number of eigenvalues per k-point below the Fermi energy.

#include <petsc/finclude/petsc.h>
      use petsc
      use globals
      implicit none

      PetscReal :: x, f, N, T, ew(:)
      integer :: wk(:), idx(:)
      integer, optional :: i_max_ewk(:)

      integer :: N_ew_total, i, ik, nrow
      PetscReal :: N_elec, fac

      N_ew_total = size(ew)
      nrow  = N_ew_total / size(wk)

      N_elec = 0d0
      if (present(i_max_ewk)) i_max_ewk = 0
      do i = 1, Ntot_diag
        ik = ceiling(1d0/real(nrow,dp)*(idx(i)+1)-1d0/(real(nrow,dp)*2d0))
        fac = fermi_function(ew(idx(i)+1),x,T)
        if ((fac.ge.1d-13).and.(present(i_max_ewk))) &
          i_max_ewk(ik) = i_max_ewk(ik) + 1
        N_elec = N_elec + fac*wk(ik)
      end do
      N_elec = N_elec / nktot
      f = N - N_elec
    end subroutine Nelec_min_fct


end module kmat_from_diag