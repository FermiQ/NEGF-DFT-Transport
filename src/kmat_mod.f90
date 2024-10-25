module kmat_mod

  ! This module contains subroutines for dumping the K-matrix to files.

contains

  subroutine dump_kmat(lfixelectrodes, outfile, outdir)
    ! Purpose: Dumps the K-matrix to a Conquest-compatible format.
    !
    ! Functionality: This subroutine handles the dumping of the K-matrix, 
    !                considering whether electrodes are fixed or not.  If electrodes
    !                are fixed, it adds contributions from left and right electrodes
    !                to the central region's K-matrix before dumping.
    !
    ! Parameters:
    !   lfixelectrodes (logical, input): True if electrodes are fixed, False otherwise.
    !   outfile (character(*), input): The base name of the output file.
    !   outdir (character(*), input): The output directory.


    ! Include necessary PETSc header file.
    #include <petsc/finclude/petsc.h>
    ! Use necessary modules.
    use petscmat
    use petsc_mod
    use petsc_wrapper
    use kinds
    use globals
    use write_conquest_dump

    implicit none

    ! Declare variables.
    logical :: lfixelectrodes
    character(*) :: outfile, outdir

    real(dp) :: damp, dconv1, dconv2, deltal, deltar, b(3, 3), c(3), rhs(3)
    character(strln) :: k_negf_file
    integer :: iunitscf, ierr, iter, imu1, ijkmax(2), i1, i2
    PetscScalar :: p_scal
    PetscReal :: norm
    character(256) :: s1,s2
    Mat :: p_mat1



    ! Check if electrodes are fixed.
    if (lfixelectrodes) then

      ! Print message indicating addition of left electrode to central region.
      pstr_out = "add electrode L to C"; call petsc_print_master()
      ! Loop over the grid points of the left electrode.
      do i2 = -ncell_l(2), ncell_l(2)
        do i1 = -ncell_l(1), ncell_l(1)
          ! Zero out the entries of the submatrix.
          call MatZeroEntries(p_dmatxy(i1, i2, -1), ierr)
          ! Add the contribution of the left electrode to the central region's K-matrix.
          call petsc_add_sub_B_to_A(p_k10_l(i1, i2),p_dmatxy(i1, i2, -1), 0, nmu_ecc - nmu_l, p_one, INSERT_VALUES, PETSC_FALSE)
          ! Commented out lines likely related to adding k00_l contribution.
!~           call petsc_add_sub_B_to_A(p_k00_l(i1, i2),p_dmatxy(i1, i2, 0), 0, 0, p_one, INSERT_VALUES, PETSC_FALSE)
        end do
      end do

      ! Print message indicating addition of right electrode to central region.
      pstr_out = "add electrode R to C"; call petsc_print_master()
      ! Loop over the grid points of the right electrode.
      do i2 = -ncell_r(2), ncell_r(2)
        do i1 = -ncell_r(1), ncell_r(1)
          ! Zero out the entries of the submatrix.
          call MatZeroEntries(p_dmatxy(i1, i2, 1), ierr)
          ! Add the contribution of the right electrode to the central region's K-matrix.
          call petsc_add_sub_B_to_A(p_k01_r(i1, i2), p_dmatxy(i1, i2, 1), nmu_ecc - nmu_r, 0, p_one, INSERT_VALUES, PETSC_FALSE)
          ! Commented out lines likely related to adding k00_r contribution.
!~           call petsc_add_sub_B_to_A(p_k00_r(i1, i2), p_dmatxy(i1, i2, 0), nmu_ecc - nmu_r, nmu_ecc - nmu_r, p_one, INSERT_VALUES, PETSC_FALSE)
        end do
      end do

    end if

    ! Write the K-matrix to a Conquest-compatible file using write_conquest2.
    call write_conquest2(outfile, outdir, xyz_ecc, dlat_c, nat_ecc, tomat, tomatr, &
     & inao_c, neigh_c, nlist_c, ndim_c, atoms_c, p_dmatxy)


    ! Commented out alternative write_conquest call.
!~     call write_conquest(outfile, xyz_ecc, dlat_c, nat_ecc, tomat, inao_c, neigh_c,&
!~    &nlist_c, ndim_c, atoms_c, p_dmatxy)

  end subroutine dump_kmat

  subroutine dump_kmat_lr(outfile, outdir, lr)
    ! Purpose: Dumps the K-matrix for either the left ('l') or right ('r') electrode.
    !
    ! Functionality: Loads pre-calculated K-matrix parts for a specified electrode 
    !                (left or right), assembles them, and writes the result to a 
    !                Conquest-compatible file.
    !
    ! Parameters:
    !   outfile (character(*), input): The base name of the output file.
    !   outdir (character(*), input): The output directory containing pre-calculated K-matrix parts.
    !   lr (character(1), input): 'l' for left electrode, 'r' for right electrode.


    ! Include necessary PETSc header file.
    #include <petsc/finclude/petsc.h>
    ! Use necessary modules.
    use petscmat
    use petsc_mod
    use petsc_wrapper
    use kinds
    use globals
    use misc
    use write_conquest_dump

    implicit none

    character(1) :: lr
    character(*) :: outfile, outdir

    real(dp) :: damp, dconv1, dconv2, deltal, deltar, b(3, 3), c(3), rhs(3)
    character(strln) :: k_negf_file
    integer :: iunitscf, ierr, iter, imu1, ijkmax(2), i1, i2
    PetscScalar :: p_scal
    PetscReal :: norm
    character(256) :: s1,s2, filename
    Mat, allocatable :: p_mat1(:, :, :)
    Mat :: p_mat2

    ! Check if the electrode is left ('l').
    if (lr .eq. "l") then

      ! Print message indicating dumping of left electrode part.
      pstr_out = "dump part to fix L"; call petsc_print_master()
      ! Allocate memory for the K-matrix of the left electrode.
      allocate(p_mat1(-ncell_c(1):ncell_c(1), -ncell_c(2):ncell_c(2),-1:1))
      
      ! Loop over the grid points of the left electrode.
      do i2 = -ncell_l(2), ncell_l(2)
        do i1 = -ncell_l(1), ncell_l(1)

          ! Convert integer indices to strings.
          s1 = int2str(i1)
          s2 = int2str(i2)          

          ! Create and load pre-calculated K-matrix parts from files.
          call MatCreate(PETSC_COMM_WORLD, p_mat1(i1, i2, 0), ierr)
          call MatSetType(p_mat1(i1, i2, 0), mattype_sparse, ierr)
          filename = trim(outdir) // "/K_u_" // trim(s1) // "_" // trim(s2) //"_0_petsc.dat"
          call petsc_mat_direct_load(p_mat1(i1, i2, 0), filename, ierr)
          call MatZeroEntries(p_mat1(i1, i2, 0), ierr)
          
          call MatCreate(PETSC_COMM_WORLD, p_mat1(i1, i2, -1), ierr)
          call MatSetType(p_mat1(i1, i2, -1), mattype_sparse, ierr)
          filename = trim(outdir) // "/K_u_" // trim(s1) // "_" // trim(s2) //"_-1_petsc.dat"
          call petsc_mat_direct_load(p_mat1(i1, i2, -1), filename, ierr)
          call MatZeroEntries(p_mat1(i1, i2, -1), ierr)
          
          call MatCreate(PETSC_COMM_WORLD, p_mat1(i1, i2, 1), ierr)
          call MatSetType(p_mat1(i1, i2, 1), mattype_sparse, ierr)
          filename = trim(outdir) // "/K_u_" // trim(s1) // "_" // trim(s2) //"_-1_petsc.dat"
          call petsc_mat_direct_load(p_mat1(i1, i2, 1), filename, ierr)
          call MatZeroEntries(p_mat1(i1, i2, 1), ierr)
                              
          ! Assemble the K-matrix parts.
          call petsc_add_sub_B_to_A(p_k00_l(i1, i2), p_mat1(i1, i2, 0), 0, 0, p_one, INSERT_VALUES, PETSC_FALSE)
          call petsc_add_sub_B_to_A(p_k10_l(i1, i2), p_mat1(i1, i2, 0), nmu_l, 0, p_one, INSERT_VALUES, PETSC_FALSE)
          call petsc_add_sub_B_to_A(p_k01_l(i1, i2), p_mat1(i1, i2, 0), 0, nmu_l, p_one, INSERT_VALUES, PETSC_FALSE)
          call petsc_split_matrix(p_dmatxy(i1, i2, 0), p_mat2, &
         &  1, nmu_l * 2, 1, nmu_l * 2, MAT_INITIAL_MATRIX)
          call petsc_add_sub_B_to_A(p_mat2, p_mat1(i1, i2, 0), nmu_l , nmu_l, p_one, INSERT_VALUES, PETSC_FALSE)
          call MatDestroy(p_mat2, ierr)


        end do
      end do

      ! Write the assembled K-matrix to a Conquest-compatible file.
      call write_conquest2(outfile, outdir, xyz_elec_l3, dlat_l3, nat_elec_l3, tomat_l, tomatr_r, &
     &  inao_l, neigh_l, nlist_l, ndim_l, atoms_l, p_mat1)

      ! Destroy and deallocate matrices.
      do i2 = -ncell_l(2), ncell_l(2)
        do i1 = -ncell_l(1), ncell_l(1)    
          call MatDestroy(p_mat1(i1, i2, -1), ierr)
          call MatDestroy(p_mat1(i1, i2, 0), ierr)
          call MatDestroy(p_mat1(i1, i2, 1), ierr)
        end do
      end do
      
      deallocate(p_mat1)
      
    ! Check if the electrode is right ('r').
    else if (lr .eq. "r") then
    
      ! Print message indicating dumping of right electrode part.
      pstr_out = "dump part to fix R"; call petsc_print_master()
      ! Allocate memory for the K-matrix of the right electrode.
      allocate(p_mat1(-ncell_c(1):ncell_c(1), -ncell_c(2):ncell_c(2),-1:1))


      ! Loop over the grid points of the right electrode.
      do i2 = -ncell_r(2), ncell_r(2)
        do i1 = -ncell_r(1), ncell_r(1)

          ! Convert integer indices to strings.
          s1 = int2str(i1)
          s2 = int2str(i2)          

          ! Create and load pre-calculated K-matrix parts from files.
          call MatCreate(PETSC_COMM_WORLD, p_mat1(i1, i2, 0), ierr)
          call MatSetType(p_mat1(i1, i2, 0), mattype_sparse, ierr)
          filename = trim(outdir) // "/K_u_" // trim(s1) // "_" // trim(s2) //"_0_petsc.dat"
          call petsc_mat_direct_load(p_mat1(i1, i2, 0), filename, ierr)
          call MatZeroEntries(p_mat1(i1, i2, 0), ierr)
          
          call MatCreate(PETSC_COMM_WORLD, p_mat1(i1, i2, -1), ierr)
          call MatSetType(p_mat1(i1, i2, -1), mattype_sparse, ierr)
          filename = trim(outdir) // "/K_u_" // trim(s1) // "_" // trim(s2) //"_-1_petsc.dat"
          call petsc_mat_direct_load(p_mat1(i1, i2, -1), filename, ierr)
          call MatZeroEntries(p_mat1(i1, i2, -1), ierr)
          
          call MatCreate(PETSC_COMM_WORLD, p_mat1(i1, i2, 1), ierr)
          call MatSetType(p_mat1(i1, i2, 1), mattype_sparse, ierr)
          filename = trim(outdir) // "/K_u_" // trim(s1) // "_" // trim(s2) //"_-1_petsc.dat"
          call petsc_mat_direct_load(p_mat1(i1, i2, 1), filename, ierr)
          call MatZeroEntries(p_mat1(i1, i2, 1), ierr)
                              

          ! Assemble the K-matrix parts.
          call petsc_add_sub_B_to_A(p_k00_r(i1, i2), p_mat1(i1, i2, 0), 2 * nmu_r, 2 * nmu_r, p_one, INSERT_VALUES, PETSC_FALSE)
          call petsc_add_sub_B_to_A(p_k10_r(i1, i2), p_mat1(i1, i2, 0), 2 * nmu_r, nmu_r, p_one, INSERT_VALUES, PETSC_FALSE)
          call petsc_add_sub_B_to_A(p_k01_r(i1, i2), p_mat1(i1, i2, 0), nmu_r, 2 * nmu_r, p_one, INSERT_VALUES, PETSC_FALSE)
          call petsc_split_matrix(p_dmatxy(i1, i2, 0), p_mat2, &
         &  nmu_ecc - 2 * nmu_r + 1 , nmu_ecc, nmu_ecc - 2 * nmu_r + 1 , nmu_ecc, MAT_INITIAL_MATRIX)
          call petsc_add_sub_B_to_A(p_mat2, p_mat1(i1, i2, 0), 0 , 0, p_one, INSERT_VALUES, PETSC_FALSE)
          call MatDestroy(p_mat2, ierr)


        end do
      end do
      
      ! Write the assembled K-matrix to a Conquest-compatible file.
      call write_conquest2(outfile, outdir, xyz_elec_r3, dlat_r3, nat_elec_r3, tomat_r, tomatr_r, &
     &  inao_r, neigh_r, nlist_r, ndim_r, atoms_r, p_mat1)

      ! Destroy and deallocate matrices.
      do i2 = -ncell_l(2), ncell_l(2)
        do i1 = -ncell_l(1), ncell_l(1)    
          call MatDestroy(p_mat1(i1, i2, -1), ierr)
          call MatDestroy(p_mat1(i1, i2, 0), ierr)
          call MatDestroy(p_mat1(i1, i2, 1), ierr)
        end do
      end do
       
      deallocate(p_mat1)

    end if



  end subroutine dump_kmat_lr

end module kmat_mod