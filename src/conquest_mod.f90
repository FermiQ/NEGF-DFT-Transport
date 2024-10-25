module conquest_mod
  implicit none

contains

  ! Purpose: Reads Conquest output file to get information about the atomic structure and connectivity.
  ! Functionality: This subroutine reads data from a Conquest output file ("mat_info.dat") located in the specified directory. 
  !   It extracts information about atom connectivity, neighbor lists, and atomic basis function indices. This information is then 
  !   broadcast to all MPI processes.
  !
  ! Parameters:
  !   data_dir (input, character(*) ): Directory containing the Conquest output files.
  !   xyz (input, real(dp), dimension(:,:)): Cartesian coordinates of atoms.
  !   cell (input, real(dp), dimension(3,3)): Unit cell vectors.
  !   nat (input, integer): Total number of atoms.
  !   tomat (output, integer, dimension(:,:)): Connectivity information between basis functions.
  !   tomatr (output, real(dp), dimension(:,:)): Reduced coordinates of neighboring atoms.
  !   inao (output, integer, dimension(:)): Number of atomic orbitals per atom.
  !   neigh (output, integer, dimension(:)): Number of neighbors for each atom.
  !   nlist (output, integer, dimension(:)): Neighbor list for each atom.
  !   ndim (output, integer, dimension(:)): Dimensionality of each atom.
  !   atoms (output, integer, dimension(:)): Global atom indices.

  subroutine read_conquest_info(data_dir, xyz, cell, nat, tomat, tomatr, inao, neigh, nlist, ndim, atoms)
    use MPI
    use kinds
    use error_handler

    use globals, only: imu_ecc, l_ionode, ecc_dir
    implicit none

    character(*) :: data_dir
    real(dp) :: xyz(:, :), tomatr(:, :)
    real(dp) :: cell(3, 3)
    integer, allocatable :: atoms(:)
    integer, allocatable :: tomat(:, :), inao(:), neigh(:), ndim(:), nlist(:), tmp_nlist(:)
    integer :: nat, nfiles
    integer :: i1, nneigh, iat1, iat2, idum, natfile, &
      ifile, iunit, ioff, ierr, imat, imu1, imu2, &
      inu1, inu2, iat3, intvec(3)
    integer, allocatable :: natms(:), nnao(:)
    integer, allocatable :: niks(:, :), miks(:, :)
    real(dp), allocatable :: riks(:,:)
    real(dp) :: vec(3), dist(3), vec2(3)
    
    ! Only the ionode (rank 0) reads the data from the file.
    if (l_ionode) then
      ! Allocate a dummy array to avoid issues with unallocated arrays.
      allocate (nlist(0))
      ioff = 0; inao = 0; neigh = 0; ndim = 0
      ! Construct the file path.
      write(0,*) trim(data_dir)//"/mat_info.dat"
      ! Open the Conquest output file.
      open (newunit=iunit, file=trim(data_dir)//"/mat_info.dat", status="old")
      ! Read the number of files.
      read (iunit, *) nfiles

      iat1 = 0; ioff = 0; imat = 0; nneigh = 0
      ! Loop over each file.
      do ifile = 1, nfiles

        ! Read the file information.
        read (iunit, *) idum, natfile
        ioff = ioff + natfile

        ! Read the number of atomic orbitals, neighbors, and dimensionality for each atom.
        read (iunit, *) inao(ioff - natfile + 1:ioff)
        read (iunit, *) neigh(ioff - natfile + 1:ioff)
        read (iunit, *) ndim(ioff - natfile + 1:ioff)

        ! Calculate the total number of neighbors.
        nneigh = nneigh + sum(neigh)
        ! Allocate temporary array for neighbor list.
        allocate (tmp_nlist(nneigh), stat=ierr)
        ! Check for allocation errors.
        if (ierr .ne. 0) then
          write (errormsg, *) "allocation error ", ierr
          call error()
        end if
        ! Copy existing neighbor list to temporary array.
        tmp_nlist(1:size(nlist)) = nlist(1:size(nlist))
        ! Move allocated memory from temporary to main array.
        call move_alloc(tmp_nlist, nlist)

        ! Loop over each atom in the file.
        do i1 = 1, natfile
          iat1 = iat1 + 1
          ! Read the global atom index.
          read (iunit, *) iat2
          atoms(iat1) = iat2
          ! Allocate arrays for neighbor information.
          allocate (nnao(neigh(iat1)), natms(neigh(iat1)), niks(3, neigh(iat1)), &
         &  riks(3, neigh(iat1)),  miks(3, neigh(iat1)), stat=ierr)
          ! Check for allocation errors.
          if (ierr .ne. 0) then
            write (errormsg, *) "allocation error ", ierr
            call error()
          end if
          ! Read neighbor information.
          read (iunit, *) nnao(1:neigh(iat1))
          read (iunit, *) natms(1:neigh(iat1))
          ! Update the neighbor list.
          nlist(1 + sum(neigh(1:iat1 - 1)):sum(neigh(1:iat1))) = natms(1:neigh(iat1))
          ! Loop over each neighbor.
          do iat2 = 1, neigh(iat1)
            ! Read the reduced coordinates and integer vectors.
            read (iunit, *) vec2(1:3), intvec(1:3)
            ! Convert reduced coordinates to Cartesian coordinates.
            vec = vec2(1)*cell(1:3, 1) + vec2(2)*cell(1:3, 2) + vec2(3)*cell(1:3, 3)
            ! Calculate the distance vector.
            dist = vec - xyz(1:3, atoms(iat1)) + xyz(1:3, natms(iat2))
            ! Calculate the integer components of the distance vector.
            niks(1, iat2) = nint(dist(1)/cell(1, 1))
            niks(2, iat2) = nint(dist(2)/cell(2, 2))
            niks(3, iat2) = nint(dist(3)/cell(3, 3))
            ! Store the integer vectors and reduced coordinates.
            miks(1:3, iat2) = intvec(1:3)            
            riks(1:3, iat2) = vec2(1:3)
          end do
          ! Loop over basis functions of the atom and its neighbors to construct the connectivity matrix.
          do iat2 = 1, neigh(iat1)
            do imu2 = 1, nnao(iat2)
              do imu1 = 1, inao(iat1)
                ! Calculate the global indices of the basis functions.
                inu1 = 0
                do iat3 = 1, atoms(iat1) - 1
                  inu1 = inu1 + imu_ecc(iat3)
                end do
                inu1 = inu1 + imu1
                inu2 = 0
                do iat3 = 1, natms(iat2) - 1
                  inu2 = inu2 + imu_ecc(iat3)
                end do
                inu2 = inu2 + imu2
                ! Update the connectivity matrix.
                imat = imat + 1
                tomat(1, imat) = inu1
                tomat(2, imat) = inu2
                tomat(3:5, imat) = niks(1:3, iat2)
                tomat(6:8, imat) = miks(1:3, iat2)
                tomatr(1:3, imat) = riks(1:3, iat2)
              end do
            end do
          end do

          ! Deallocate temporary arrays.
          deallocate (nnao, natms, niks, miks, riks)

        end do
      end do
      ! Close the file.
      close (iunit)
    end if

    ! Broadcast the total number of neighbors to all processes.
    call MPI_Bcast(nneigh, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    ! Allocate the neighbor list on non-ionode processes.
    if (.not. l_ionode) then
      allocate (nlist(nneigh), stat=ierr)
      ! Check for allocation errors.
      if (ierr .ne. 0) then
        write (errormsg, *) "allocation error ", ierr
        call error()
      end if
    end if
    ! Synchronize all processes.
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    ! Broadcast the neighbor list, number of atomic orbitals, dimensionality, and connectivity matrix to all processes.
    call MPI_Bcast(nlist, size(nlist), MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(inao, size(inao), MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(ndim, size(ndim), MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(tomat, size(tomat), MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

  end subroutine read_conquest_info

end module conquest_mod