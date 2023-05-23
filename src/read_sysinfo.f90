module readsys
contains

  subroutine read_sysinfo(sys_xyz, sys_species, sys_nmu, natoms, dlat, ncell, nmat, ef, ne, sysfile)

    use kinds
    use lattice_mod, only: read_cell
    use globals, only: inode
    use petsc_mod, only: pstr_out, petsc_print_master
    use error_handler

    implicit none

    real(dp), allocatable, intent(out) :: sys_xyz(:, :)
    real(dp) :: dlat(3, 3), ef, ne(2)
    integer :: ncell(3), nmat
    integer, allocatable, intent(out) :: sys_nmu(:)
    integer, intent(out) :: natoms
    character(strln), allocatable, intent(out)  :: sys_species(:)

    character(*), intent(in) :: sysfile

    integer :: iunit, io, ierr, i, iat, ityp

    open (newunit=iunit, file=trim(sysfile), action="read", status="old")
    call read_cell(dlat, iunit)
    read (iunit, *) ncell(1:3)
    read (iunit, *) nmat
    read (iunit, *) ef
    read (iunit, *) ne(1)
    read (iunit, *) ne(2)
    read (iunit, *) natoms
    allocate (sys_xyz(3, natoms), sys_nmu(natoms), sys_species(natoms), stat=ierr)
    if (ierr .ne. 0) then
      write (errormsg, fmt='(A,i8)') "allocation error ", ierr
      call error()
    end if

    do i = 1, natoms
      read (iunit, *) iat, ityp, sys_nmu(i), sys_xyz(1:3, i), sys_species(i)
    end do

    close (iunit)

    ncell(3) = 1
  end subroutine read_sysinfo

end module readsys
