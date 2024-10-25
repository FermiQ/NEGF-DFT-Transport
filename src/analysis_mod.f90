module analysis
  implicit none
  
  contains
  
  subroutine bond_order()
    ! Purpose: Calculate and output the bond order matrix.
    ! Functionality: This subroutine calculates the bond order matrix using a k-point summation. 
    !                It utilizes PETSc for matrix operations and MPI for parallel computation.
    !                The resulting bond order matrix is written to a file named "bond_order.dat".

    ! Include necessary PETSc header file
#include <petsc/finclude/petsc.h>
    use petsc  ! Use PETSc module
    use petsc_mod ! Use PETSc module (likely a custom extension)
    use globals ! Use module containing global variables
    use k_on_demand_mod ! Use module for on-demand k-point calculations

    implicit none
    
    ! p_DS: PETSc matrix to store intermediate results.
    Mat :: p_DS
    ! ik: k-point index.
    ! iat1, iat2: Atom indices.
    ! imu1, imu2: Orbital indices within atoms.
    ! ipos, jpos: Linear indices for matrix elements.
    ! ierr: PETSc error code.
    ! iunit: Unit number for output file.
    ! iam: MPI rank.
    ! nao_max: Maximum number of orbitals per atom.
    ! nn: Total number of orbital pairs.
    ! n1, n2: Number of orbitals in atoms iat1 and iat2, respectively.
    integer :: ik, iat1, iat2, imu1, imu2, ipos, jpos, ierr, iunit, iam, &
      nao_max, nn, n1, n2
    ! p_Bmunu, p_Bnumu: Temporary variables for matrix elements.
    ! Bk_ij: Bond order contribution from a single k-point.
    PetscScalar :: p_Bmunu, p_Bnumu, Bk_ij
    ! B_ij: Bond order matrix.
    ! p_BBmunu, p_BBnumu: Temporary matrices to store intermediate results.
    PetscScalar, allocatable :: B_ij(:,:), p_BBmunu(:,:), p_BBnumu(:,:)
    ! counti, countf: Counters for timing.
    ! count_rate: Clock rate.
    integer(8) :: counti, count_rate, countf
    
    ! Determine the maximum number of orbitals per atom.
    nao_max = maxval(imu_ecc)
    ! Allocate memory for matrices.
    allocate(B_ij(nat_ecc,nat_ecc), p_BBmunu(nao_max, nao_max), &
      p_BBnumu(nao_max, nao_max), stat = ierr)
    ! Initialize the bond order matrix to zero.
    B_ij = 0d0
    ! Loop over k-points in reverse order.
    do ik = nk_c, 1, -1
      ! Perform on-demand k-point calculation.
      call k_on_demand(ik, 2)
      ! Perform matrix multiplication: p_k00k_cc(ik) * p_s00k_cc(ik) -> p_DS
      call MatMatMult(p_k00k_cc(ik), p_s00k_cc(ik),  MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, p_DS, ierr)      
      ! Print k-point and atom information.
      write (pstr_out, fmt='(A, 2i8)') "k-point, atom ", ik, nat_ecc ; call petsc_print_master(.false.)
      ! Start timer.
      call system_clock(counti, count_rate)
      ! Loop over atom pairs.
      do iat1 = 1, nat_ecc
        ! Print progress indicator.
        write (pstr_out, fmt='(A)') repeat(ACHAR(8), 8); call petsc_print_master(.false.)            
        write (pstr_out, fmt='(i8)') nat_ecc - iat1; call petsc_print_master(.false.)
        do iat2 = 1, nat_ecc
          ! Initialize variables.
          Bk_ij = 0d0
          p_BBmunu = 0d0
          p_BBnumu = 0d0
          ! Loop over orbitals within atom pairs.
          do imu1 = 1, imu_ecc(iat1)
            ! Calculate linear index for orbital imu1.
            ipos = sum(imu_ecc(1:iat1-1)) + imu1 - 1
            do imu2 = 1, imu_ecc(iat2)
              ! Calculate linear index for orbital imu2.
              jpos = sum(imu_ecc(1:iat2-1)) + imu2 - 1
              ! Get matrix element from p_DS.
              call petsc_mat_getvalue(ipos, jpos, p_BBmunu(imu1, imu2), p_DS, 1,&
                PETSC_COMM_WORLD, iam, .false.)
              ! Get matrix element from p_DS.
              call petsc_mat_getvalue(jpos, ipos, p_BBnumu(imu2, imu1), p_DS, 1,&
                PETSC_COMM_WORLD, iam, .false.)
            end do           
          end do
          ! Set sizes for MPI reduction.
          n1 = imu_ecc(iat1)
          n2 = imu_ecc(iat2)
          nn = n1 * n2
          ! Perform MPI AllReduce to sum matrix elements across processors.
          call MPI_AllReduce(MPI_IN_PLACE, p_BBmunu(1:n1,1:n2), nn, MPI_DOUBLE_COMPLEX, &
            MPI_SUM, PETSC_COMM_WORLD, ierr)  
          call MPI_AllReduce(MPI_IN_PLACE, p_BBnumu(1:n2,1:n1), nn, MPI_DOUBLE_COMPLEX, &
            MPI_SUM, PETSC_COMM_WORLD, ierr)  
          ! Calculate k-point contribution to bond order.
          do imu1 = 1, imu_ecc(iat1)
            do imu2 = 1, imu_ecc(iat2)
              Bk_ij = Bk_ij + p_BBmunu(imu1, imu2)*p_BBnumu(imu2, imu1)
            end do           
          end do
          
          
          ! Accumulate bond order contribution from current k-point.
          B_ij(iat1, iat2) = B_ij(iat1, iat2) + Bk_ij
          ! Add complex conjugate if k-point is of type 2.
          if (wkp_l(ik).eq.2) B_ij(iat1, iat2) = B_ij(iat1, iat2) + conjg(Bk_ij) 
        end do         
      end do
      ! Destroy the intermediate matrix.
      call MatDestroy(p_DS, ierr)
      ! Stop timer and print elapsed time.
      call system_clock(countf)
      write (pstr_out, fmt='(X,e24.12)') real(countf - counti, 8)/real(count_rate, 8); call petsc_print_master()
    end do
    
    ! Normalize the bond order matrix.
    B_ij = B_ij / (nktot)
    
    ! Write the bond order matrix to a file if on the ionode.
    if (l_ionode) then
      open(newunit = iunit, file = "bond_order.dat", action = "write", status = "replace")
      do iat1 = 1, nat_ecc
        do iat2 = 1, nat_ecc
          write(unit = iunit, fmt = '(E24.12E4)', advance="no") real(B_ij(iat1, iat2),dp)
        end do
        write(iunit,fmt=*)
      end do      
      close(iunit)
    end if
    
    ! Ensure all processes have finished before exiting.
    call MPI_BARRIER(PETSC_COMM_WORLD, ierr)
    
  end subroutine bond_order

end module analysis