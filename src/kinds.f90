module kinds

  ! Module for defining kind parameters for various data types.

  ! sp: Kind parameter for single precision real numbers.
  !     Provides at least 6 decimal digits of precision and a range of at least 10**37.
  integer, parameter :: sp = selected_real_kind(6, 37)

  ! dp: Kind parameter for double precision real numbers.
  !     Provides at least 15 decimal digits of precision and a range of at least 10**307.
  integer, parameter :: dp = selected_real_kind(15, 307)

  ! qp: Kind parameter for quadruple precision real numbers (if available).
  !     Provides at least 33 decimal digits of precision and a range of at least 10**4931.
  integer, parameter :: qp = selected_real_kind(33, 4931)

  ! strln: Kind parameter for string length.  Defines a maximum string length of 256 characters.
  integer, parameter :: strln = 256

end module kinds