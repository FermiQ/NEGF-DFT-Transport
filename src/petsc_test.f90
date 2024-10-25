module petsc_test
  ! This module contains a subroutine for testing PETSc functionalities.
  #include <slepc/finclude/slepceps.h>
  use petscmat  ! Use the PETSc matrix module.
  use slepceps  ! Use the SLEPc eigenvalue solver module.
  use slepc_mod ! Use a custom SLEPc module (likely containing wrappers or additional functionalities).
  use petsc_mod ! Use a custom PETSc module (likely containing wrappers or additional functionalities).
  use petsc_wrapper ! Use a custom PETSc wrapper module.
  implicit none

contains

  subroutine test_petsc()
    ! This subroutine demonstrates a simple linear system solve using PETSc.
    use globals  ! Use a module named 'globals', likely containing global variables.
    use petsc_mod ! Use the custom PETSc module.
    implicit none

    ! Declare PETSc objects.
    Mat :: a, c, f, y, e, d, g, h  ! Declare various matrix variables of type Mat.
    Vec :: b, x                    ! Declare vector variables of type Vec.
    KSP :: ksp                     ! Declare a Krylov subspace solver object.
    PC :: pc                       ! Declare a preconditioner object.
    PetscScalar, pointer ::  xx_v(:), bb_v(:) ! Declare pointers to arrays for vector data.
    PetscScalar, pointer :: mat_p(:,:) ! Declare a pointer to a 2D array for matrix data.

    ! Declare PETSc data types.
    MatSolverType :: matsolvertype1 ! Type of matrix solver to be used.
    MatType ::  mattype            ! Type of matrix.
    PetscReal :: norm, tol         ! Variables for norm and tolerance.

    ! Declare integer variables.
    integer :: ierr, nl1, nl2, irow ! Error code and indices.

    ! Commenting out alternative matrix solver types.
!~     matsolvertype1=MATSOLVERMUMPS
!~     matsolvertype1=MATSOLVERPASTIX
!~     matsolvertype1=MATSOLVERSUPERLU_DIST
    matsolvertype1 = MATSOLVERMKL_CPARDISO ! Selecting MKL Pardiso as the matrix solver.

    ! Get information about a matrix (likely from the 'globals' module).
    call petsc_mat_info(p_h00_cc(0, 0), "p_h00_cc(0,0) ", ierr)

    ! Duplicate the matrix p_h00_cc(0,0) to matrix 'a'.
    call MatDuplicate(p_h00_cc(0, 0), MAT_COPY_VALUES, a, ierr)

    ! Get the ownership range of matrix 'a'.
    call MatGetOwnershipRange(a, nl1, nl2, ierr); CHKERRQ(ierr)

    ! Create vectors 'b' and 'x' compatible with matrix 'a'.
    call MatCreateVecs(a, b, PETSC_NULL_VEC, ierr)
    call MatCreateVecs(a, x, PETSC_NULL_VEC, ierr)

    ! Set irow to 0.  This seems to be an index for a specific row operation.
    irow = 0

    ! Set a specific element of vector 'b' if irow is within the ownership range.
    if ((irow .ge. nl1) .and. (irow .le. nl2 - 1)) then
      call VecGetArrayF90(b, xx_v, ierr)
      xx_v(irow - nl1 + 1) = 1d0
      call VecRestoreArrayF90(b, xx_v, ierr)
    end if

    ! Create and setup the Krylov subspace solver.
    call KSPCreate(PETSC_COMM_WORLD, ksp, ierr)
    call KSPSetOperators(ksp, a, a, ierr)  ! Set the matrix for the linear system.
    call KSPSetType(ksp, KSPPREONLY, ierr) ! Set the KSP type to pre-only.
    call KSPGetPC(ksp, pc, ierr)           ! Get the preconditioner from the KSP.
    call PCSetType(pc, PCLU, ierr)         ! Set the preconditioner type to LU.
    call PCFactorSetMatSolverType(pc, matsolvertype1, ierr) ! Set the matrix solver type for the preconditioner.
    call KSPSetUp(ksp, ierr)               ! Set up the KSP solver.
    call PCFactorGetMatrix(pc, F, ierr)     ! Get the factored matrix from the preconditioner.
    call KSPGetPC(ksp, pc, ierr)           ! Get the preconditioner again (redundant?).
    call KSPSolve(ksp, b, x, ierr)         ! Solve the linear system.

    ! Set the element back to 0 (cleanup).
    if ((irow .ge. nl1) .and. (irow .le. nl2 - 1)) then
      call VecGetArrayF90(b, xx_v, ierr)
      xx_v(irow - nl1 + 1) = 0d0
      call VecRestoreArrayF90(b, xx_v, ierr)
    end if

  end subroutine test_petsc

end module petsc_test