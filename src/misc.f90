module misc
  implicit none

  ! Interface block for string-to-various-type conversion subroutines.
  interface str2

    ! Converts a string to a 4-byte integer.
    !
    ! Args:
    !   instr: Input string containing the integer value.
    !   outint: Output 4-byte integer.
    subroutine str2int4(instr, outint)
      implicit none
      character(*) :: instr
      integer(4) :: outint
    end subroutine str2int4

    ! Converts a string to an 8-byte integer.
    !
    ! Args:
    !   instr: Input string containing the integer value.
    !   outint: Output 8-byte integer.
    subroutine str2int8(instr, outint)
      implicit none
      character(*) :: instr
      integer(8) :: outint
    end subroutine str2int8

    ! Converts a string to a double-precision real number.
    !
    ! Args:
    !   instr: Input string containing the real number value.
    !   outreal: Output double-precision real number.  Uses dp kind from kinds module.
    subroutine str2real(instr, outreal)
      use kinds
      implicit none
      character(*) :: instr
      real(dp) :: outreal
    end subroutine str2real

    ! Copies a string to another string variable.
    !
    ! Args:
    !   instr: Input string.
    !   outstr: Output string (copy of instr).
    subroutine str2str(instr, outstr)
      use kinds
      implicit none
      character(*) :: instr
      character(*) :: outstr
    end subroutine str2str

    ! Converts a string to a logical value.  '1' is considered true, anything else is false.
    !
    ! Args:
    !   instr: Input string containing '1' or any other value.
    !   outlogic: Output logical value.
    subroutine str2logic(instr, outlogic)
      use kinds
      implicit none
      character(*) :: instr
      logical :: outlogic
    end subroutine str2logic

  end interface str2

  ! Interface block for integer-to-string conversion function.
  interface int2str
    ! Converts an integer to a string.
    !
    ! Args:
    !   i: Input integer.
    ! Returns:
    !   i2s: Output string representation of the integer.
    function i2s(i)
      character(256) :: i2s
      integer :: i
    end function i2s
  end interface int2str

end module misc

! Converts a string to a double-precision real number.
subroutine str2real(str, outreal)
  use kinds
  implicit none

  character(*) :: str
  real(dp) :: outreal

  read (str, *) outreal

end subroutine str2real

! Converts a string to a 4-byte integer.
subroutine str2int4(str, outint)

  implicit none

  character(*) :: str
  integer(4) :: outint

  read (str, *) outint

end subroutine str2int4

! Converts a string to an 8-byte integer.
subroutine str2int8(str, outint)

  implicit none

  character(*) :: str
  integer(8) :: outint

  read (str, *) outint

end subroutine str2int8

! Copies a string to another string variable.
subroutine str2str(str, outstr)

  implicit none

  character(*) :: str
  character(*) :: outstr

  outstr = str

end subroutine str2str

! Converts a string to a logical value. '1' is considered true, anything else is false.
subroutine str2logic(str, outlogic)

  implicit none

  character(*) :: str
  logical :: outlogic
  integer :: istr

  read (str, *) istr

  outlogic = 1 .eq. istr

end subroutine str2logic

! Converts an integer to a string.
function i2s(i)
  character(256) :: i2s
  integer :: i

  write (i2s, fmt='(i8)') i
  i2s = adjustl(i2s)

end function i2s