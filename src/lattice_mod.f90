module lattice_mod
  implicit none
contains

  subroutine read_cell(brav, iunit)
    ! Purpose: Reads the lattice vectors from a file.
    !
    ! Functionality: Reads three lattice vectors from the specified file unit.
    !                It then checks if the vectors are orthogonal and raises an error if not.
    !
    ! Parameters:
    !   brav (out): A 3x3 array containing the lattice vectors. Each column represents a vector.
    !   iunit (in): The Fortran unit number of the file containing the lattice vectors.

    use petscmat ! Provides PETSc matrix operations.
    use kinds ! Provides definitions for kind types (e.g., dp for double precision).
    use globals, only: eps_real ! Provides access to the global variable eps_real (a small tolerance value).
    use error_handler ! Provides error handling routines.
    implicit none

    real(dp), intent(out) :: brav(3, 3) ! Output: 3x3 array of lattice vectors (double precision).
    integer, intent(in) :: iunit ! Input: Fortran unit number.

    integer :: i, j, ierr ! Loop counters and error code.
    character(256) :: str ! String variable for output.

    ! Read lattice vectors from file.
    do i = 1, 3
      read (iunit, *) brav(1:3, i) ! Read the i-th lattice vector.
      write (str, fmt='(3f24.12)') brav(1:3, i) ! Format the vector for output.
      call PetscPrintf(PETSC_COMM_WORLD, trim(str)//New_Line('A'), ierr) ! Print the vector.
    end do

    ! Check orthogonality of lattice vectors.
    do i = 1, 2
      do j = i + 1, 3
        if (.not. abs(dot_product(brav(1:3, i), brav(1:3, j))) .lt. eps_real) then
          ! If the dot product is not close to zero, the vectors are not orthogonal.
          write (errormsg, *) "super cell vector ", i, " and ", j, " are not orthogonal"
          call error() ! Raise an error.
        end if
      end do
    end do

  end subroutine read_cell

  subroutine read_kpoint_file(kps, wkps, nkp)
    ! Purpose: Reads k-points and their weights from a file.
    !
    ! Functionality: Reads k-points and their corresponding weights from the file "k_points.dat".
    !                It handles both fractional and Cartesian coordinates.
    !
    ! Parameters:
    !   kps (out): A 2xnkp array containing the k-points.
    !   wkps (out): A nkp array containing the weights of the k-points.
    !   nkp (out): The number of k-points.

    use petscmat
    use kinds
    use globals, only: rlat_l ! Provides access to the reciprocal lattice vectors.
    use error_handler
    implicit none

    real(dp), allocatable :: kps(:, :) ! Output: Array of k-points (double precision).
    integer, allocatable :: wkps(:) ! Output: Array of k-point weights.
    integer :: nkp ! Output: Number of k-points.

    integer :: iunit, ik, ierr ! File unit, loop counter, and error code.
    character(256) :: instr, frac, str ! Variables for reading file and formatting output.
    logical :: lfrac ! Flag indicating whether k-points are in fractional coordinates.

    ! Open the k-points file.
    open (newunit=iunit, file="k_points.dat", status="old", action="read")

    ! Read the header line and number of k-points.
    read (iunit, fmt='(A)') instr
    call PetscPrintf(PETSC_COMM_WORLD, trim(instr)//New_Line('A'), ierr)
    read (instr, *) nkp, frac
    frac = adjustl(frac)
    lfrac = .false.
    if (index(frac, "frac") .ge. 1) lfrac = .true.

    ! Allocate memory for k-points and weights.
    allocate (kps(2, nkp), wkps(nkp), stat=ierr)
    if (ierr .ne. 0) then
      write (errormsg, *) "allocation error ", ierr
      call error()
    end if

    ! Print header for k-points and weights.
    write (str, fmt='(A)') "k-points and weights from file:"
    call PetscPrintf(PETSC_COMM_WORLD, trim(str)//New_Line('A'), ierr)

    ! Read k-points and weights.
    do ik = 1, nkp
      read (iunit, *) kps(1:2, ik), wkps(ik)
      if (lfrac) kps(1:2, ik) = kps(1, ik)*rlat_l(1:2, 1) + kps(2, ik)*rlat_l(1:2, 2)
      write (str, fmt='(2e24.12,i8)') kps(1:2, ik), wkps(ik)
      call PetscPrintf(PETSC_COMM_WORLD, trim(str)//New_Line('A'), ierr)
    end do

    ! Close the k-points file.
    close (iunit)

  end subroutine read_kpoint_file

  subroutine get_rlat(dlat, rlat, kps, wkps, nkp, nkx, nky, nkz, makeks)
    ! Purpose: Computes the reciprocal lattice vectors and generates k-points.
    !
    ! Functionality: Calculates the reciprocal lattice vectors from the direct lattice vectors.
    !                If makeks is true, it generates a set of k-points based on nkx, nky, nkz,
    !                considering symmetry if ksym is true.
    !
    ! Parameters:
    !   dlat (in): A 3x3 array containing the direct lattice vectors.
    !   rlat (out): A 3x3 array containing the reciprocal lattice vectors.
    !   kps (out): A 3xnkp array containing the generated k-points.
    !   wkps (out): A nkp array containing the weights of the generated k-points.
    !   nkp (out): The number of generated k-points.
    !   nkx (in): Number of k-points in the x-direction.
    !   nky (in): Number of k-points in the y-direction.
    !   nkz (in): Number of k-points in the z-direction.
    !   makeks (in): A logical flag indicating whether to generate k-points.

    use petscmat
    use kinds
    use globals, only: pi, ksym ! Provides access to pi and ksym (symmetry flag).
    use mathstuff, only: vector_product ! Provides vector product function.
    use error_handler
    implicit none

    real(dp), intent(in) :: dlat(3, 3) ! Input: 3x3 array of direct lattice vectors.
    integer, intent(in) :: nkx, nky, nkz ! Input: Number of k-points in each direction.
    logical, intent(in) :: makeks ! Input: Flag to generate k-points.
    real(dp), allocatable :: kps(:, :) ! Output: Array of k-points.
    integer, allocatable :: wkps(:) ! Output: Array of k-point weights.
    real(dp), intent(out) :: rlat(3, 3) ! Output: 3x3 array of reciprocal lattice vectors.
    integer :: nkp ! Output: Number of k-points.
    integer :: i, j, nxk, nyk, nk1, nk2, nk3, ierr, n1, &
      ik1, ik2, ik3, i1, in1, in2, in3 ! Loop counters and variables.
    integer, allocatable :: ikps(:, :) ! Array to store k-point indices.
    real(dp) :: dk(3), deltak(3, 3) ! Variables for k-point spacing and shift.
    character(256) :: str ! String variable for output.
    logical :: l_withgamma ! Flag indicating whether the Gamma point is included.

    ! Calculate reciprocal lattice vectors.
    rlat(1:3, 1) = vector_product(dlat(1:3, 2), dlat(1:3, 3))
    rlat(1:3, 1) = rlat(1:3, 1)/(dot_product(dlat(1:3, 1), rlat(1:3, 1)))

    rlat(1:3, 2) = vector_product(dlat(1:3, 3), dlat(1:3, 1))
    rlat(1:3, 2) = rlat(1:3, 2)/(dot_product(dlat(1:3, 2), rlat(1:3, 2)))

    rlat(1:3, 3) = vector_product(dlat(1:3, 1), dlat(1:3, 2))
    rlat(1:3, 3) = rlat(1:3, 3)/(dot_product(dlat(1:3, 3), rlat(1:3, 3)))

    rlat = rlat*2d0*pi
    call PetscPrintf(PETSC_COMM_WORLD, "reciprocal lattice:"//New_Line('A'), ierr)
    do i = 1, 3
      write (str, fmt='(3f24.12)') rlat(1:3, i)
      call PetscPrintf(PETSC_COMM_WORLD, trim(str)//New_Line('A'), ierr)
    end do
    deltak(1:3, 1) = rlat(1:3, 1)/real(max(1, nkx), dp)
    deltak(1:3, 2) = rlat(1:3, 2)/real(max(1, nky), dp)
    deltak(1:3, 3) = rlat(1:3, 3)/real(max(1, nkz), dp)

    ! If makeks is false, return without generating k-points.
    if (.not. makeks) return

    call PetscPrintf(PETSC_COMM_WORLD, "steps reciprocal lattice:"//New_Line('A'), ierr)
    do i = 1, 3
      write (str, fmt='(3f24.12)') deltak(1:3, i)
      call PetscPrintf(PETSC_COMM_WORLD, trim(str)//New_Line('A'), ierr)
    end do

    ! Determine the number of k-points considering symmetry.
    l_withgamma = .true.
    if (ksym) then
      nkp = floor(nkx*nky*nkz*0.5d0) + 1
      if (mod(nkx*nky*nkz, 2) .eq. 0) then
        nkp = nkp - 1
        l_withgamma = .false.
      end if
    else
      nkp = nkx*nky*nkz
    end if

    ! Calculate k-point shift for symmetry.
    dk = 0d0
    if (mod(nkx, 2) .eq. 0) dk = dk + deltak(1:3, 1)*0.5d0
    if (mod(nky, 2) .eq. 0) dk = dk + deltak(1:3, 2)*0.5d0
    if (mod(nkz, 2) .eq. 0) dk = dk + deltak(1:3, 3)*0.5d0

    ! Allocate memory for k-points, weights, and indices.
    allocate (kps(3, nkp), wkps(nkp), ikps(3, nkp), stat=ierr)
    if (ierr .ne. 0) then
      write (errormsg, *) "allocation error ", ierr
      call error()
    end if

    ! Generate k-points and weights.
    kps = 0d0
    wkps = 0d0
    n1 = 0d0

    in1 = -1
    ik1_loop: do ik1 = floor(nkx*0.5d0), -floor(nkx*0.5d0), -1
      if ((mod(nkx, 2) .eq. 0) .and. (ik1 .eq. 0)) cycle ik1_loop
      in1 = in1 + 1
      in2 = -1
      ik2_loop: do ik2 = floor(nky*0.5d0), -floor(nky*0.5d0), -1
        if ((mod(nky, 2) .eq. 0) .and. (ik2 .eq. 0)) cycle ik2_loop
        in2 = in2 + 1
        in3 = -1        
        ik3_loop: do ik3 = floor(nkz*0.5d0), -floor(nkz*0.5d0), -1
          if ((mod(nkz, 2) .eq. 0) .and. (ik3 .eq. 0)) cycle ik3_loop
  
          ! Handle symmetry.
          if (ksym) then
            do i1 = 1, n1
              if (all((/ik1, ik2, ik3/) + ikps(1:3, i1) .eq. 0)) then
                wkps(i1) = wkps(i1) + 1
                cycle ik3_loop
              end if
            end do
          end if
          
          in3 = in3 + 1
          n1 = n1 + 1
          kps(1:3, n1) = -dk(1:3) + (floor(nkx*0.5d0) - in1)*deltak(1:3, 1) +&
            (floor(nky*0.5d0) - in2)*deltak(1:3, 2) + (floor(nkz*0.5d0) - in3)*deltak(1:3, 3)
          ikps(1:3, n1) = (/ik1, ik2, ik3/)
          wkps(n1) = 1
        end do ik3_loop
      end do ik2_loop
    end do ik1_loop


    ! Print k-point shift and k-points with weights.
    call PetscPrintf(PETSC_COMM_WORLD, "k-point shift"//New_Line('A'), ierr)
    write (str, fmt='(3e24.12)') dk(1:3)
    call PetscPrintf(PETSC_COMM_WORLD, trim(str)//New_Line('A'), ierr)
    call PetscPrintf(PETSC_COMM_WORLD, "k-points and weights:"//New_Line('A'), ierr)
    do ik1 = 1, nkp
      write (str, fmt='(4i8,3e24.12,i4)') ik1, ikps(1:3, ik1), kps(1:3, ik1), wkps(ik1)
      call PetscPrintf(PETSC_COMM_WORLD, trim(str)//New_Line('A'), ierr)
    end do

  end subroutine get_rlat

end module lattice_mod