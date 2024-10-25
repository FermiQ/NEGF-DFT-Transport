module write_conquest_dump
  implicit none

contains

  subroutine write_conquest(outfile, xyz, cell, nat, tomat, inao, neigh, nlist, ndim, atoms, p_mat)
    ! Purpose: Writes Conquest-format dump files containing atomic coordinates, connectivity, and density matrix information.
    !
    ! Functionality: This subroutine writes data to a file in a format compatible with the Conquest DFT code.  It handles parallel I/O using MPI and PETSc for efficient processing of large datasets.  The subroutine writes atomic positions, connectivity information (neighbors, atom types), and the density matrix.
    !
    ! Parameters:
    !   outfile (character*): Name of the output file.
    !   xyz (real(dp), allocatable(:, :)): Array of atomic coordinates (x, y, z for each atom).
    !   cell (real(dp)(3, 3)): Unit cell vectors.
    !   nat (integer): Total number of atoms.
    !   tomat (integer, allocatable(:, :)): Connectivity matrix defining the relationship between atoms.
    !   inao (integer, allocatable(:)): Array containing the number of atomic orbitals for each atom.
    !   neigh (integer, allocatable(:)): Array specifying the number of neighbors for each atom.
    !   nlist (integer, allocatable(:)): Neighbor list array.
    !   ndim (integer, allocatable(:)): Array containing dimensionality information for each atom.
    !   atoms (integer, allocatable(:)): Array mapping global atom indices to local indices.
    !   p_mat (Mat, allocatable(:, :, :)): 3D array representing the density matrix.

#include <petsc/finclude/petsc.h>
    use petscmat
    use mpi
    use kinds
    use globals, only: imu_ecc, inode
    use petsc_mod, only: pstr_out

    Mat, allocatable :: p_mat(:, :, :)
    real(dp), allocatable :: xyz(:, :)
    real(dp) :: cell(3, 3)
    integer, allocatable :: atoms(:)
    integer, allocatable :: tomat(:, :), inao(:), neigh(:), ndim(:), nlist(:)
    integer :: nat
    character(*) :: outfile
    character(256) :: outstr

    integer :: iunit, ierr, ik1, ik2, ik3, imu1, imu2, iat1, iat2, inu1, inu2, junit, iat3, jmu, &
               nl1, nl2, p_inu1(1), p_inu2(1), iipos, jjpos, buflen_i, buflen_j, ikind
    integer(KIND=MPI_OFFSET_KIND) :: ipos, jpos, ioff, ipos2, jpos2, nwrites, ndwrites
    real(dp) :: dist(3), r_out
    PetscScalar :: p_scal(1)
    logical :: ldebug_dump

    ikind = kind(ierr)
    ldebug_dump = .true.

    call MatGetOwnershipRange(p_mat(0, 0, 0), nl1, nl2, ierr)

    nl1 = nl1 + 1
    call MPI_File_delete(trim(outfile), MPI_INFO_NULL, ierr)
    if (ldebug_dump) call MPI_File_delete(trim(outfile)//".dat", MPI_INFO_NULL, ierr)
    call MPI_FILE_OPEN(PETSC_COMM_WORLD, trim(outfile), MPI_MODE_WRONLY + MPI_MODE_CREATE, MPI_INFO_NULL, iunit, ierr)
    if (ldebug_dump) call MPI_FILE_OPEN(PETSC_COMM_WORLD, trim(outfile)//".dat", MPI_MODE_WRONLY + MPI_MODE_CREATE, MPI_INFO_NULL, junit, ierr)

    if (inode .eq. 0) then
      ipos = 0
      buflen_i = 1
      call MPI_File_write_at(iunit, ipos, 1, buflen_i, MPI_INTEGER, MPI_STATUS_IGNORE, ierr)
      ipos = ipos + buflen_i*ikind
      buflen_i = 1
      call MPI_File_write_at(iunit, ipos, nat, buflen_i, MPI_INTEGER, MPI_STATUS_IGNORE, ierr)
      ipos = ipos + buflen_i*ikind
      buflen_i = size(inao)
      call MPI_File_write_at(iunit, ipos, inao, buflen_i, MPI_INTEGER, MPI_STATUS_IGNORE, ierr)
      ipos = ipos + buflen_i*ikind
      buflen_i = size(neigh)
      call MPI_File_write_at(iunit, ipos, neigh, buflen_i, MPI_INTEGER, MPI_STATUS_IGNORE, ierr)
      ipos = ipos + buflen_i*ikind
      buflen_i = size(ndim)
      call MPI_File_write_at(iunit, ipos, ndim, buflen_i, MPI_INTEGER, MPI_STATUS_IGNORE, ierr)
      ipos = ipos + buflen_i*ikind
      buflen_i = 1
      call MPI_File_write_at(iunit, ipos, 1, buflen_i, MPI_INTEGER, MPI_STATUS_IGNORE, ierr) ! spin identification number
      ipos = ipos + buflen_i*ikind
    end if
    call MPI_Bcast(ipos, 1, MPI_INTEGER8, 0, PETSC_COMM_WORLD, ierr)

    imu1 = 0
    imu2 = 0
    write (pstr_out, fmt='(A,2i8)') "dumped data for atom ", 0, nat; call PetscPrintf(PETSC_COMM_WORLD, trim(pstr_out), ierr)

    if (ldebug_dump) jpos = 0

    do iat1 = 1, nat

      call MPI_Barrier(MPI_COMM_WORLD, ierr)

      if (inode .eq. 0) then
        buflen_i = 1
        call MPI_File_write_at(iunit, ipos, atoms(iat1), buflen_i, MPI_INTEGER, MPI_STATUS_IGNORE, ierr)
        ipos = ipos + buflen_i*ikind
        do iat2 = 1, neigh(iat1)
          buflen_i = 1
          call MPI_File_write_at(iunit, ipos, imu_ecc(nlist(sum(neigh(1:iat1 - 1)) + iat2)), buflen_i, MPI_INTEGER, MPI_STATUS_IGNORE, ierr)
          ipos = ipos + buflen_i*ikind
        end do
        buflen_i = size(nlist(1 + sum(neigh(1:iat1 - 1)):sum(neigh(1:iat1))))
        call MPI_File_write_at(iunit, ipos, nlist(1 + sum(neigh(1:iat1 - 1)):sum(neigh(1:iat1))), buflen_i, MPI_INTEGER, MPI_STATUS_IGNORE, ierr)
        ipos = ipos + buflen_i*ikind
        imu2 = imu1
        do iat2 = 1, neigh(iat1)
          imu1 = imu1 + imu_ecc(nlist(sum(neigh(1:iat1 - 1)) + iat2))*inao(iat1)
          ik1 = -tomat(3, imu1)
          ik2 = -tomat(4, imu1)
          ik3 = -tomat(5, imu1)
          dist = ik1*cell(1:3, 1) + ik2*cell(1:3, 2) + ik3*cell(1:3, 3) + xyz(1:3, atoms(iat1)) - xyz(1:3, nlist(sum(neigh(1:iat1 - 1)) + iat2))
          dist(1) = dist(1)/cell(1, 1)
          dist(2) = dist(2)/cell(2, 2)
          dist(3) = dist(3)/cell(3, 3)
          buflen_i = size(dist)
          call MPI_File_write_at(iunit, ipos, dist, buflen_i, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, ierr)
          ipos = ipos + buflen_i*8
        end do
      end if

      call MPI_Bcast(imu2, 1, MPI_INTEGER, 0, PETSC_COMM_WORLD, ierr)
      call MPI_Bcast(imu1, 1, MPI_INTEGER, 0, PETSC_COMM_WORLD, ierr)
      call MPI_Bcast(ipos, 1, MPI_INTEGER8, 0, PETSC_COMM_WORLD, ierr)

      buflen_i = 0
      ioff = 0

      call MPI_Barrier(MPI_COMM_WORLD, ierr)

      nwrites = 0

      if (ldebug_dump) ndwrites = 0
      do inu1 = imu2 + 1, imu1

        buflen_i = 0
        ioff = 0

        if (ldebug_dump) buflen_j = 0

        if ((tomat(1, inu1) .ge. nl1) .and. (tomat(1, inu1) .le. nl2)) then
          p_inu1 = tomat(1, inu1) - 1
          p_inu2 = tomat(2, inu1) - 1
          call MatGetValues(p_mat(tomat(3, inu1), tomat(4, inu1), tomat(5, inu1)), 1, p_inu1, 1, p_inu2, p_scal, ierr)
          r_out = real(p_scal(1), 8)
          buflen_i = 1
          ioff = inu1 - imu2 - 1
          nwrites = nwrites + 1
          ipos2 = ipos + ioff*8
          call MPI_FILE_WRITE_AT(iunit, ipos2, r_out, buflen_i, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, ierr)
          if (ldebug_dump) then
            write (outstr, fmt='(5i7,es45.24e5,A)') tomat(3, inu1), tomat(4, inu1), tomat(5, inu1), p_inu1 + 1, p_inu2 + 1, r_out, new_line('A')
            buflen_j = len(trim(outstr))
          end if
        end if


        if (ldebug_dump) then
          jpos2 = jpos + buflen_j*ioff
          ndwrites = ndwrites + buflen_j
          call MPI_FILE_WRITE_AT(junit, jpos2, outstr, buflen_j, MPI_CHARACTER, MPI_STATUS_IGNORE, ierr)
        end if

      end do

      call MPI_ALLREDUCE(MPI_IN_PLACE, nwrites, 1, MPI_INTEGER8, MPI_SUM, MPI_COMM_WORLD, ierr)
      ipos = ipos + nwrites*8

      if (ldebug_dump) then
        call MPI_ALLREDUCE(MPI_IN_PLACE, ndwrites, 1, MPI_INTEGER8, MPI_SUM, MPI_COMM_WORLD, ierr)
        jpos = jpos + ndwrites
      end if
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
      write (pstr_out, fmt='(A,2i8)') repeat(ACHAR(8), 16), iat1, nat; call PetscPrintf(PETSC_COMM_WORLD, trim(pstr_out), ierr)
    end do

    call MPI_File_Close(iunit, ierr)
    if (ldebug_dump) call MPI_File_Close(junit, ierr)
    call PetscPrintf(PETSC_COMM_WORLD, new_line('A'), ierr)

    if (inode .eq. 0) then
      open (newunit=iunit, file="InfoGlobal_negf.i00", action="write", status="replace")
      write (iunit, fmt='(A)') "F F   = flag_velocity, flag_MDstep"
      write (iunit, fmt='(2i22)') nat, 1
      write (iunit, fmt='(3i22)') 0, 0, 0
      write (iunit, fmt='(3e24.12)') cell(1, 1), cell(2, 2), cell(3, 3)
      do iat1 = 1, nat
        if (mod(iat1, 6) .eq. 0 .or. iat1 .eq. nat) then
          write (iunit, fmt='(i22)') 1
        else
          write (iunit, fmt='(i22)', advance="no") 1
        end if
      end do
      write (iunit, fmt='(A)') "0   MD step"
      do iat1 = 1, nat
        write (iunit, fmt='(5x,i8,3x,3e22.13)') iat1, xyz(1:3, iat1)
      end do
      write (iunit, *) "0  = index_local in dump_InfoMatGlobal"
      close (iunit)
    end if

    call MPI_Barrier(MPI_COMM_WORLD, ierr)

  end subroutine write_conquest

  subroutine write_conquest2(outfile, outdir, xyz, cell, nat, tomat, tomatr, inao, neigh, &
 &  nlist, ndim, atoms, p_mat)
    ! Purpose: Writes Conquest-format dump files with improved efficiency and handling of large density matrices.
    !
    ! Functionality: This subroutine is a more efficient version of `write_conquest`, designed to handle larger density matrices and improve I/O performance. It uses buffered writing to reduce the number of MPI calls.  The output filename is constructed by combining `outdir` and `outfile`.
    !
    ! Parameters:
    !   outfile (character*): Base name of the output file.
    !   outdir (character*): Directory where the output file will be written.
    !   xyz (real(dp), allocatable(:, :)): Array of atomic coordinates.
    !   cell (real(dp)(3, 3)): Unit cell vectors.
    !   nat (integer): Total number of atoms.
    !   tomat (integer, allocatable(:, :)): Connectivity matrix.
    !   tomatr (real(dp), allocatable(:, :)): Matrix containing additional real-valued data related to connectivity.
    !   inao (integer, allocatable(:)): Number of atomic orbitals per atom.
    !   neigh (integer, allocatable(:)): Number of neighbors for each atom.
    !   nlist (integer, allocatable(:)): Neighbor list.
    !   ndim (integer, allocatable(:)): Dimensionality information.
    !   atoms (integer, allocatable(:)): Mapping of global to local atom indices.
    !   p_mat (Mat, allocatable(:, :, :)): Density matrix.

#include <petsc/finclude/petsc.h>
    use petscmat
    use mpi
    use kinds
    use globals, only: imu_ecc, inode
    use error_handler
    use petsc_mod, only: pstr_out, petsc_mat_getvalue, petsc_print_master

    Mat, allocatable :: p_mat(:, :, :)
    real(dp), allocatable :: xyz(:, :), tomatr(:, :)
    real(dp) :: cell(3, 3)
    integer, allocatable :: atoms(:)
    integer, allocatable :: tomat(:, :), inao(:), neigh(:), ndim(:), nlist(:)
    integer :: nat
    character(*) :: outfile, outdir
    character(256) :: outstr

    integer :: iunit, ierr, ik1, ik2, ik3, imu1, imu2, iat1, iat2, inu1, inu2, junit, &
      iat3, jmu, nl1, nl2, p_inu1(1), p_inu2(1), iipos, jjpos, buflen_i, buflen_j, &
      ikind, maxblock, j, i
    integer, allocatable :: imu1_at(:), imu2_at(:)
    integer(KIND=MPI_OFFSET_KIND) :: ipos, jpos, ioff, ipos2, jpos2, nwrites, ndwrites
    integer(KIND=MPI_OFFSET_KIND), allocatable :: ipos_at(:), jpos_at(:), ipos2_at(:)
    real(dp) :: dist(3), r_out
    real(dp), allocatable :: buf_out(:)
    PetscScalar :: p_scal(1)
    logical :: ldebug_dump
    character(256) :: filename

    filename = trim(outdir)//"/"//trim(outfile)
    write (pstr_out, fmt='(A,A,X)') "filename: ",trim(filename) 
    call petsc_print_master()
    
    
    ikind = kind(ierr)
    ldebug_dump = .false.

    call MatGetOwnershipRange(p_mat(0, 0, 0), nl1, nl2, ierr)

    nl1 = nl1 + 1
    call MPI_File_delete(trim(filename), MPI_INFO_NULL, ierr)
    if (ldebug_dump) call MPI_File_delete(trim(filename)//".dat", MPI_INFO_NULL, ierr)
    call MPI_FILE_OPEN(PETSC_COMM_WORLD, trim(filename), MPI_MODE_WRONLY + MPI_MODE_CREATE, MPI_INFO_NULL, iunit, ierr)
    if (ldebug_dump) call MPI_FILE_OPEN(PETSC_COMM_WORLD, trim(filename)//".dat", MPI_MODE_WRONLY + MPI_MODE_CREATE, MPI_INFO_NULL, junit, ierr)

    if (inode .eq. 0) then
      ipos = 0
      buflen_i = 1
      call MPI_File_write_at(iunit, ipos, 1, buflen_i, MPI_INTEGER, MPI_STATUS_IGNORE, ierr)
      ipos = ipos + buflen_i*ikind
      buflen_i = 1
      call MPI_File_write_at(iunit, ipos, nat, buflen_i, MPI_INTEGER, MPI_STATUS_IGNORE, ierr)
      ipos = ipos + buflen_i*ikind
      buflen_i = size(inao)
      call MPI_File_write_at(iunit, ipos, inao, buflen_i, MPI_INTEGER, MPI_STATUS_IGNORE, ierr)
      ipos = ipos + buflen_i*ikind
      buflen_i = size(neigh)
      call MPI_File_write_at(iunit, ipos, neigh, buflen_i, MPI_INTEGER, MPI_STATUS_IGNORE, ierr)
      ipos = ipos + buflen_i*ikind
      buflen_i = size(ndim)
      call MPI_File_write_at(iunit, ipos, ndim, buflen_i, MPI_INTEGER, MPI_STATUS_IGNORE, ierr)
      ipos = ipos + buflen_i*ikind
      buflen_i = 1
      call MPI_File_write_at(iunit, ipos, 1, buflen_i, MPI_INTEGER, MPI_STATUS_IGNORE, ierr) ! spin identification number
      ipos = ipos + buflen_i*ikind
    end if
    call MPI_Bcast(ipos, 1, MPI_INTEGER8, 0, PETSC_COMM_WORLD, ierr)

    imu1 = 0
    imu2 = 0
    write (pstr_out, fmt='(A,2i8)') "dumped data for atom ", 0, nat; call PetscPrintf(PETSC_COMM_WORLD, trim(pstr_out), ierr)

    if (ldebug_dump) jpos = 0

  
    allocate(ipos_at(nat), jpos_at(nat), imu1_at(nat), imu2_at(nat), stat = ierr)
    if (ierr .ne. 0) then
      write (errormsg, *) "allocation error write_conquest2 ", ierr
      call error()
    end if

    maxblock = 0

    do iat1 = 1, nat

      call MPI_Barrier(MPI_COMM_WORLD, ierr)

      if (inode .eq. 0) then
        buflen_i = 1
        call MPI_File_write_at(iunit, ipos, atoms(iat1), buflen_i, MPI_INTEGER, MPI_STATUS_IGNORE, ierr)
        ipos = ipos + buflen_i*ikind
        do iat2 = 1, neigh(iat1)
          buflen_i = 1
          call MPI_File_write_at(iunit, ipos, imu_ecc(nlist(sum(neigh(1:iat1 - 1)) + iat2)), buflen_i, MPI_INTEGER, MPI_STATUS_IGNORE, ierr)
          ipos = ipos + buflen_i*ikind
        end do
        buflen_i = size(nlist(1 + sum(neigh(1:iat1 - 1)):sum(neigh(1:iat1))))
        call MPI_File_write_at(iunit, ipos, nlist(1 + sum(neigh(1:iat1 - 1)):sum(neigh(1:iat1))), buflen_i, MPI_INTEGER, MPI_STATUS_IGNORE, ierr)
        ipos = ipos + buflen_i*ikind
        imu2 = imu1
        do iat2 = 1, neigh(iat1)
          imu1 = imu1 + imu_ecc(nlist(sum(neigh(1:iat1 - 1)) + iat2))*inao(iat1)
          ik1 = tomat(3, imu1)
          ik2 = tomat(4, imu1)
          ik3 = tomat(5, imu1)
          dist = ik1*cell(1:3, 1) + ik2*cell(1:3, 2) + ik3*cell(1:3, 3) + xyz(1:3, atoms(iat1)) - xyz(1:3, nlist(sum(neigh(1:iat1 - 1)) + iat2))
          dist(1) = dist(1)/cell(1, 1)
          dist(2) = dist(2)/cell(2, 2)
          dist(3) = dist(3)/cell(3, 3)
!~           write(0, fmt = '(2i4,i8,X,A,X,3e12.4,A,3e12.4,A,3i4,A,3i4)') iat1,iat2,imu1,trim(outdir),dist(1:3), " , ", tomatr(1:3, imu1), " , ", tomat(3:5, imu1), " , ", tomat(6:8, imu1)
          dist(1:3) = tomatr(1:3, imu1)
          buflen_i = size(dist)
          call MPI_File_write_at(iunit, ipos, dist, buflen_i, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, ierr)
          ipos = ipos + buflen_i*8
        end do
      end if

      call MPI_Bcast(imu2, 1, MPI_INTEGER, 0, PETSC_COMM_WORLD, ierr)
      call MPI_Bcast(imu1, 1, MPI_INTEGER, 0, PETSC_COMM_WORLD, ierr)
      call MPI_Bcast(ipos, 1, MPI_INTEGER8, 0, PETSC_COMM_WORLD, ierr)

      buflen_i = 0
      ioff = 0

      call MPI_Barrier(MPI_COMM_WORLD, ierr)

      nwrites = 0

      if (ldebug_dump) ndwrites = 0
      
      ipos_at(iat1) = ipos
      if (ldebug_dump) jpos_at(iat1) = jpos
      imu1_at(iat1) = imu1
      imu2_at(iat1) = imu2
      
      do inu1 = imu2 + 1, imu1

        if (ldebug_dump) buflen_j = 0

        if ((tomat(1, inu1) .ge. nl1) .and. (tomat(1, inu1) .le. nl2)) then
          nwrites = nwrites + 1
          if (ldebug_dump) then
!~             write (outstr, fmt='(5i7,es45.24e5,A)') tomat(3, inu1), tomat(4, inu1), &
            write (outstr, fmt=*) tomat(3, inu1), tomat(4, inu1), &
              tomat(5, inu1), p_inu1 + 1, p_inu2 + 1, "666d0", new_line('A')
            buflen_j = len(trim(outstr))
          end if
        end if

        if (ldebug_dump) then
          ndwrites = ndwrites + buflen_j
        end if

      end do

      maxblock = max(maxblock, nwrites)
      call MPI_ALLREDUCE(MPI_IN_PLACE, nwrites, 1, MPI_INTEGER8, MPI_SUM, MPI_COMM_WORLD, ierr)
      ipos = ipos + nwrites*8

      if (ldebug_dump) then
        call MPI_ALLREDUCE(MPI_IN_PLACE, ndwrites, 1, MPI_INTEGER8, MPI_SUM, MPI_COMM_WORLD, ierr)
        jpos = jpos + ndwrites
      end if
      
      call MPI_Barrier(MPI_COMM_WORLD, ierr)      
    end do


    allocate(ipos2_at(maxblock), buf_out(maxblock), stat = ierr)
    if (ierr .ne. 0) then
      write (errormsg, *) "allocation error write_conquest2 ", ierr
      call error()
    end if


    do iat1 = 1, nat
      

      buflen_i = 0
      ioff = 0

      nwrites = 0
      if (ldebug_dump) ndwrites = 0
      
      ipos = ipos_at(iat1)
      jpos = jpos_at(iat1)
      imu1 = imu1_at(iat1)
      imu2 = imu2_at(iat1)
      
      do inu1 = imu2 + 1, imu1

        buflen_i = 0
        ioff = 0

        if (ldebug_dump) buflen_j = 0

        if ((tomat(1, inu1) .ge. nl1) .and. (tomat(1, inu1) .le. nl2)) then
          p_inu1 = tomat(1, inu1) - 1
          p_inu2 = tomat(2, inu1) - 1
!~           call MatGetValues(p_mat(tomat(3, inu1), tomat(4, inu1), tomat(5, inu1)), 1, &
!~             p_inu1, 1, p_inu2, p_scal, ierr)
!~           call petsc_mat_getvalue(p_inu1(1), p_inu2(1), p_scal(1), &
!~          & p_mat(tomat(3, inu1), tomat(4, inu1), tomat(5, inu1)), &
!~          & 1, PETSC_COMM_WORLD, lcast_in = .false.)
          call petsc_mat_getvalue(p_inu1(1), p_inu2(1), p_scal(1), &
         & p_mat(tomat(6, inu1), tomat(7, inu1), tomat(8, inu1)), &
         & 1, PETSC_COMM_WORLD, lcast_in = .false.)
          r_out = real(p_scal(1), 8)
          buflen_i = 1
          ioff = inu1 - imu2 - 1          
          ipos2 = ipos + ioff*8
          nwrites = nwrites + 1
          ipos2_at(nwrites) = ipos2
          buf_out(nwrites) = r_out
          if (ldebug_dump) then
            write (outstr, fmt='(5i7,es45.24e5,A)') tomat(3, inu1), tomat(4, inu1), &
              tomat(5, inu1), p_inu1 + 1, p_inu2 + 1, r_out, new_line('A')
            buflen_j = len(trim(outstr))
          end if
        end if

!~         ipos2 = ipos + ioff*8
!~         call MPI_FILE_WRITE_AT_ALL(iunit, ipos2, r_out, buflen_i, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, ierr)

        if (ldebug_dump) then
          jpos2 = jpos + buflen_j*ioff
          ndwrites = ndwrites + buflen_j
          call MPI_FILE_WRITE_AT_ALL(junit, jpos2, outstr, buflen_j, MPI_CHARACTER, MPI_STATUS_IGNORE, ierr)
        end if

      end do

      j = -1
      do i = 1, nwrites
        if (j .eq. -1) j = i        
        if (i .eq. nwrites) then          
          call MPI_FILE_WRITE_AT(iunit, ipos2_at(j:i), buf_out(j:i), (i-j)+1, &
          MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, ierr)
        else        
          if ((ipos2_at(i+1) - ipos2_at(i)) .ne. 8) then
            call MPI_FILE_WRITE_AT(iunit, ipos2_at(j:i), buf_out(j:i), (i-j)+1, &
            MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, ierr)
            j = -1
          end if
        end if
        call MPI_File_sync(iunit, ierr)
      end do

      write (pstr_out, fmt='(A,2i8)') repeat(ACHAR(8), 16), iat1, nat; call PetscPrintf(PETSC_COMM_WORLD, trim(pstr_out), ierr)
      
    end do

    call MPI_File_sync(iunit, ierr)
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    
    call MPI_File_Close(iunit, ierr)
    if (ldebug_dump) call MPI_File_Close(junit, ierr)
    call PetscPrintf(PETSC_COMM_WORLD, new_line('A'), ierr)

    if (inode .eq. 0) then
      open (newunit=iunit, file=trim(outdir)//"/InfoGlobal_negf.i00", action="write", status="replace")
      write (iunit, fmt='(A)') "F F   = flag_velocity, flag_MDstep"
      write (iunit, fmt='(2i22)') nat, 1
      write (iunit, fmt='(3i22)') 0, 0, 0
      write (iunit, fmt='(3e24.12)') cell(1, 1), cell(2, 2), cell(3, 3)
      do iat1 = 1, nat
        if (mod(iat1, 6) .eq. 0 .or. iat1 .eq. nat) then
          write (iunit, fmt='(i22)') 1
        else
          write (iunit, fmt='(i22)', advance="no") 1
        end if
      end do
      write (iunit, fmt='(A)') "0   MD step"
      do iat1 = 1, nat
        write (iunit, fmt='(5x,i8,3x,3e22.13)') iat1, xyz(1:3, iat1)
      end do
      write (iunit, *) "0  = index_local in dump_InfoMatGlobal"
      close (iunit)
    end if

    call MPI_Barrier(MPI_COMM_WORLD, ierr)

  end subroutine write_conquest2

end module write_conquest_dump