module timing
  use kinds  ! Use the kinds module (presumably defining data types like dp)
  implicit none

  ! Declare variables to store timing information.  'dp' is assumed to be a double precision data type.
  real(dp) :: timing_surf       ! Timing for surface calculations
  real(dp) :: timing_matsolve   ! Timing for matrix solves
  real(dp) :: timing_petsc_load ! Timing for PETSc load operations
  real(dp) :: timing_petsc_save ! Timing for PETSc save operations
  real(dp) :: timing_matmat_solve(2) ! Timing for matrix-matrix solves (2 elements, possibly for different methods)

contains

  subroutine init_timings()
    implicit none

    ! Initialize all timing variables to zero.
    timing_surf = 0d0
    timing_matsolve = 0d0
    timing_petsc_load = 0d0
    timing_petsc_save = 0d0
    timing_matmat_solve = 0d0

  end subroutine init_timings

end module timing