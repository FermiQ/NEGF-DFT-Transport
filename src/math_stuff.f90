module mathstuff

  implicit none

  ! Interface block for matrix inversion routines.  Provides a consistent way to call either real or complex inversion.
  interface invert
    ! Subroutine to invert a complex double precision matrix.
    subroutine cinvert(A)
      implicit none
      ! The complex double precision matrix to be inverted.  Modified in place.
      double complex :: A(:, :)
    end subroutine cinvert

    ! Subroutine to invert a real double precision matrix.
    subroutine rinvert(A)
      implicit none
      ! The real double precision matrix to be inverted. Modified in place.
      double precision :: A(:, :)
    end subroutine rinvert

  end interface invert

  ! Interface block for matrix norm calculation. Allows calling either real or complex version.
  interface matnorm_my
    ! Function to calculate the Frobenius norm of a real double precision matrix.
    function rmatnorm_my(a)
      implicit none
      ! The real double precision input matrix.
      double precision :: a(:, :)
      ! The Frobenius norm of the input matrix.
      double precision :: rmatnorm_my
    end function rmatnorm_my

    ! Function to calculate the Frobenius norm of a complex double precision matrix.
    function cmatnorm_my(a)
      implicit none
      ! The complex double precision input matrix.
      double complex :: a(:, :)
      ! The Frobenius norm of the input matrix.
      double precision :: cmatnorm_my
    end function cmatnorm_my
  end interface matnorm_my

contains

  ! Function to check if an interval [j1, j2] overlaps with interval [i1, i2].
  function ininterval(i1, i2, j1, j2)
    integer :: i1, i2, j1, j2
    logical ::ininterval

    ! Initialize the result to false.
    ininterval = .false.

    ! Check for overlap.  The conditions are not entirely comprehensive and might need review.
    if ((j1 .ge. i1) .and. (j1 .le. i2)) ininterval = .true.

    if ((j2 .ge. i1) .and. (j1 .le. i1)) ininterval = .true.

  end function

  ! Function to compute a modified Frobenius norm between two complex matrices.
  function matnorm_myab(a, b)
    implicit none

    ! The two complex double precision input matrices.
    double complex :: a(:, :), b(:, :)
    ! The computed modified Frobenius norm.
    double precision :: matnorm_myab
    ! Variable to accumulate the trace.
    double precision :: trace
    ! Temporary variable for intermediate calculations.
    double complex :: d
    ! Loop counters.
    integer :: n, i, j

    ! Get the size of the matrices. Assumes both are square and of the same size.
    n = size(a, 1)

    ! Initialize the trace.
    trace = 0d0
    ! Compute the trace of (a-b)*(a-b)^H.
    do j = 1, n
      do i = 1, n
        d = a(i, j) - b(i, j)
        trace = trace + d*conjg(d)
      end do
    end do
    ! Return the square root of the trace.
    matnorm_myab = sqrt(trace)

  end function

  ! Function to compute the trace of the product of two real double precision matrices.
  function rmatdotab(a, b)
    implicit none

    ! The two real double precision input matrices.
    double precision :: a(:, :), b(:, :)
    ! The computed trace.
    double precision :: rmatdotab
    ! Variable to accumulate the trace.
    double precision :: trace
    ! Loop counters.
    integer :: n, i, j

    ! Get the size of the matrices. Assumes both are square and of the same size.
    n = size(a, 1)

    ! Initialize the trace.
    trace = 0d0
    ! Compute the trace of a*b.
    do j = 1, n
      do i = 1, n
        trace = trace + a(i, j)*b(i, j)
      end do
    end do
    ! Return the trace.
    rmatdotab = trace

  end function

  ! Function to compute the cross product of two 3D vectors.
  function vector_product(a, b)
    use kinds
    implicit none
    ! The two input 3D vectors.
    real(dp), intent(in) :: a(3), b(3)
    ! The resulting cross product vector.
    real(dp) :: vector_product(3)

    ! Compute the cross product using the standard formula.
    vector_product(1) = a(2)*b(3) - a(3)*b(2)
    vector_product(2) = a(3)*b(1) - a(1)*b(3)
    vector_product(3) = a(1)*b(2) - a(2)*b(1)

  end function vector_product

end module