subroutine line_ncc_rule(n, a, b, x, w)

  !*****************************************************************************80
  !
  !! LINE_NCC_RULE computes a Newton-Cotes Closed (NCC) quadrature rule.
  !
  !  Discussion:
  !
  !    The integral:
  !
  !      Integral ( A <= X <= B ) F(X) dx
  !
  !    The quadrature rule:
  !
  !      Sum ( 1 <= I <= N ) W(I) * F ( X(I) ).
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    09 April 2014
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer ( kind = 4 ) N, the order.
  !
  !    Input, real ( kind = 8 ) A, B, the endpoints of the interval.
  !
  !    Input, real ( kind = 8 ) X(N), the abscissas.
  !
  !    Output, real ( kind = 8 ) W(N), the weights.
  !
  implicit none

  integer(kind=4) n

  real(kind=8) a
  real(kind=8) b
  real(kind=8) d(n)
  integer(kind=4) i
  integer(kind=4) j
  integer(kind=4) k
  real(kind=8) w(n)
  real(kind=8) x(n)
  real(kind=8) y_a
  real(kind=8) y_b
  !
  !  Define the points X.
  !
  call r8vec_linspace(n, a, b, x)
  !
  !  Compute the Lagrange basis polynomial which is 1 at X(I),
  !  and zero at the other nodes.
  !
  do i = 1, n

    d(1:n) = 0.0D+00
    d(i) = 1.0D+00

    do j = 2, n
      do k = j, n
        d(n + j - k) = (d(n + j - k - 1) - d(n + j - k))/(x(n + 1 - k) - x(n + j - k))
      end do
    end do

    do j = 1, n - 1
      do k = 1, n - j
        d(n - k) = d(n - k) - x(n - k - j + 1)*d(n - k + 1)
      end do
    end do
    !
    !  Evaluate the antiderivative of the polynomial at the endpoints.
    !
    y_a = d(n)/real(n, kind=8)
    do j = n - 1, 1, -1
      y_a = y_a*a + d(j)/real(j, kind=8)
    end do
    y_a = y_a*a

    y_b = d(n)/real(n, kind=8)
    do j = n - 1, 1, -1
      y_b = y_b*b + d(j)/real(j, kind=8)
    end do
    y_b = y_b*b

    w(i) = y_b - y_a

  end do

  return
end subroutine line_ncc_rule

subroutine r8vec_linspace(n, a, b, x)

  !*****************************************************************************80
  !
  !! R8VEC_LINSPACE creates a vector of linearly spaced values.
  !
  !  Discussion:
  !
  !    An R8VEC is a vector of R8's.
  !
  !    4 points evenly spaced between 0 and 12 will yield 0, 4, 8, 12.
  !
  !    In other words, the interval is divided into N-1 even subintervals,
  !    and the endpoints of intervals are used as the points.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    14 March 2011
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer ( kind = 4 ) N, the number of entries in the vector.
  !
  !    Input, real ( kind = 8 ) A, B, the first and last entries.
  !
  !    Output, real ( kind = 8 ) X(N), a vector of linearly spaced data.
  !
  implicit none

  integer(kind=4) n

  real(kind=8) a
  real(kind=8) b
  integer(kind=4) i
  real(kind=8) x(n)

  if (n == 1) then

    x(1) = (a + b)/2.0D+00

  else

    do i = 1, n
      x(i) = (real(n - i, kind=8)*a &
              + real(i - 1, kind=8)*b) &
             /real(n - 1, kind=8)
    end do

  end if

  return
end subroutine r8vec_linspace