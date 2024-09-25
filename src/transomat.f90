!> The main program for the Trans'o'mat simulation tool. This program initializes the MPI and PETSc libraries, sets up the system, and performs various calculations such as self-consistent field (SCF) calculations, linear optical excitation (LOE) calculations, and transmission calculations. The program also handles convergence and output management.
program transomat
  ! Main program for Trans'o'mat, a simulation tool.

  use mpi ! Use the MPI library for parallel computing.
  use kinds ! Use the kinds module for defining data types.
  use misc ! Use the misc module for miscellaneous utility functions.
  use globals ! Use the globals module for global variables.
  use init_sys_mod ! Use the init_sys_mod module for system initialization.
  use petsc_control ! Use the petsc_control module for PETSc control parameters.
  use petsc_mod, only: pstr_out, petsc_print_master, petsc_cleanup ! Use selected components from petsc_mod.
  use scf_mod ! Use the scf_mod module for self-consistent field calculations.
  use loe_mod ! Use the loe_mod module for linear optical excitation calculations.
  use phonon_mod ! Use the phonon_mod module for phonon calculations.  (Likely unused in this specific code flow)
  use transmission_mod ! Use the transmission_mod module for transmission calculations.
  use kmat_from_diag ! Use the kmat_from_diag module for constructing the K matrix from diagonal elements.
  use kmat_mod ! Use the kmat_mod module for K matrix operations.
  use analysis ! Use the analysis module for analysis functions. (Likely unused in this specific code flow)
  use timing ! Use the timing module for performance measurements.

  implicit none
  integer :: dateandtime(8) ! Array to store date and time information.
  character(strln) :: hostname ! String variable to store the hostname.
  integer :: iunitscf, iter, ierr, provided, iin ! Integer variables for unit numbers, iteration count, error codes, and input.
  real(dp) :: conv, din ! Double-precision real variables for convergence criteria and input data.
  logical :: lflag ! Logical variable for flags. (Unused in this specific code flow)

  call MPI_INIT_THREAD(MPI_THREAD_MULTIPLE, provided, ierr) ! Initialize MPI with multiple threads.
!~   call MPI_INIT_THREAD(MPI_THREAD_SERIALIZED, provided, ierr) ! Alternative initialization with serialized threads.
!~  call MPI_INIT(ierr)  ! call MPI_INIT here because for some reason initializing MPI in petsc breaks mpiP profiler
  call petsc_init() ! Initialize PETSc library.

  call date_and_time(values=dateandtime) ! Get the current date and time.
  CALL HOSTNM(hostname) ! Get the hostname.
  write (pstr_out, fmt='(A)'); call petsc_print_master() ! Print a blank line.
  write (pstr_out, fmt='(A)') "Trans'o'mat v1.1 petsc beta"; call petsc_print_master() ! Print the program version.
  write (pstr_out, fmt='(A)') repeat("-", 80); call petsc_print_master() ! Print a separator line.
  write (pstr_out, fmt='(i2,"/",i2,"/",i4,X,i2,":",i2," @ ",A)') & ! Print the date, time, and hostname.
    dateandtime(2), dateandtime(3), dateandtime(1), dateandtime(5), &
    dateandtime(6), trim(hostname); call petsc_print_master()
  write (pstr_out, fmt='(A)') repeat("-", 80); call petsc_print_master() ! Print a separator line.

  call init_control() ! Initialize control parameters.
  call petsc_subcomm_init() ! Initialize PETSc subcommunicators.
  call petsc_solver_init() ! Initialize PETSc solvers.
  call init_timings() ! Initialize timing variables.

  iigc=0 ! Initialize a variable (likely related to Green's function calculations).

  write (pstr_out, fmt='(A)') repeat("-", 80); call petsc_print_master() ! Print a separator line.

  if (l_ionode .and. l_loadsavegf) then
!~ this is not portable at all. probaly using some C_BINDING would be better.
    call system("mkdir reaktor >/dev/null 2>/dev/null") ! Create a directory (non-portable).
  end if

  if (l_ionode .and. l_ep_reuse_lambda) then
    call system("mkdir ep_reaktor >/dev/null 2>/dev/null") ! Create another directory (non-portable).
  end if

  iter = 1 ! Initialize iteration counter.
  conv = 666d0 ! Initialize convergence criterion.
  init_scf = .true. ! Initialize flag for SCF calculation.

  if (inode .eq. 0) then
!    inquire(file="K_negf_matrix2.i00.p000000",exist=oldnegf) ! Check for existence of a file (commented out).
    inquire (file="convergence.negf", exist=init_scf) ! Check for existence of convergence file.
    init_scf = .not. init_scf ! Invert the flag based on file existence.
    converged = .false. ! Initialize convergence flag.
    if (init_scf) then
      iter = 0 ! Reset iteration counter if no convergence file exists.
      conv = 666d0 ! Reset convergence criterion.
      init_scf = .true. ! Set the init_scf flag.
    else
      open (newunit=iunitscf, file="convergence.negf", action="read", status="old") ! Open convergence file for reading.
      do
        read (unit=iunitscf, fmt=*, iostat=ierr) iin, din ! Read iteration number and convergence value.
        if (ierr .ne. 0) exit ! Exit loop if read error occurs.
        iter = iin ! Update iteration counter.
        conv = din ! Update convergence value.
      end do
      close (iunitscf) ! Close convergence file.
      converged = conv .le. scf_conv ! Check if convergence criterion is met.
    end if
    write (pstr_out, fmt='(A,i8,e24.12,X,e24.12,X,l)') "iter,conv,scf_conv", iter, conv, scf_conv, converged; call petsc_print_master() ! Print iteration, convergence, and convergence flag.

  end if

  call MPI_bcast(init_scf, 1, MPI_logical, 0, PETSC_COMM_WORLD, ierr) ! Broadcast init_scf flag.
  call MPI_bcast(oldnegf, 1, MPI_logical, 0, PETSC_COMM_WORLD, ierr) ! Broadcast oldnegf flag. (Likely unused in this specific code flow)
  call MPI_bcast(iter, 1, MPI_INTEGER, 0, PETSC_COMM_WORLD, ierr) ! Broadcast iteration counter.
  call MPI_bcast(conv, 1, MPI_DOUBLE_PRECISION, 0, PETSC_COMM_WORLD, ierr) ! Broadcast convergence value.
  call MPI_bcast(converged, 1, MPI_LOGICAL, 0, PETSC_COMM_WORLD, ierr) ! Broadcast convergence flag.
  write (pstr_out, fmt='(A,l)') "continue NEGF SCF: ", oldnegf; call petsc_print_master() ! Print NEGF SCF continuation flag.
  write (pstr_out, fmt='(A,l)') "init_scf: ", init_scf; call petsc_print_master() ! Print init_scf flag.

  call init_sys() ! Initialize the system.

  if (k_mat_mode .eq. 1) then
    if ((oneshot .or. (converged .and. (.not. init_scf) .and. iter .ge. 2)) .and. (.not. l_ep)) then
      write (pstr_out, fmt='(A)') "get transmission"; call petsc_print_master() ! Print message indicating transmission calculation.
      if (l_k_on_demand) then
        call transmission_k_on_demand() ! Perform on-demand transmission calculation.
      else
        call transmission() ! Perform transmission calculation.
      end if
      if (.not. oneshot) then
        open (newunit=iunitscf, file="negf.converged", action="write", status="replace") ! Create a file to indicate convergence.
        close (iunitscf) ! Close the file.
      end if
      if (calc_current_density .and. (vb .ne. 0d0) .and. (converged .or. oneshot)) then
        write (pstr_out, fmt='(A)') "start scf for current density"; call petsc_print_master() ! Print message indicating SCF calculation for current density.
        call scf() ! Perform SCF calculation.
      end if
    else if (l_ep) then
      calc_current_density = .false. ! Set current density calculation flag to false.
      call loe_run() ! Run linear optical excitation calculation.
    else
      calc_current_density = .false. ! Set current density calculation flag to false.
      call scf() ! Perform SCF calculation.
    end if
  else if (k_mat_mode .eq. 2) then
    if (converged) then
      open (newunit=iunitscf, file="negf.converged", action="write", status="replace") ! Create a file to indicate convergence.
      close (iunitscf) ! Close the file.
    else
      call get_kmat_from_diag() ! Get K matrix from diagonal elements.
      write (pstr_out, fmt='(A)') "...dump matrix..."; call petsc_print_master() ! Print message indicating matrix dump.
      call dump_kmat(l_diag_fixK, "K_negf_matrix2.i00.p000000" , ecc_dir) ! Dump the K matrix.
    end if
  else if (k_mat_mode .eq. 3) then
    call bond_order() ! Calculate bond order.
  end if


  write (pstr_out, fmt='(A)') "...matrix clean up..."; call petsc_print_master() ! Print message indicating matrix cleanup.
  if (.not. (k_mat_mode .eq. 3)) call petsc_cleanup() ! Clean up PETSc objects.

  write (pstr_out, fmt='(A)') "...end"; call petsc_print_master() ! Print end message.
!~   call MPI_FINALIZE(ierr) ! Finalize MPI (commented out).
  call slepcfinalize(ierr) ! Finalize SLEPc library.
  call MPI_FINALIZE(ierr) ! Finalize MPI.
end program transomat