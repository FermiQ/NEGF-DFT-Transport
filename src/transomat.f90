program transomat

  use mpi
  use kinds
  use misc
  use globals
  use init_sys_mod
  use petsc_control
  use petsc_mod, only: pstr_out, petsc_print_master, petsc_cleanup
  use scf_mod
  use loe_mod
  use phonon_mod
  use transmission_mod
  use kmat_from_diag
  use kmat_mod
  use analysis
  use timing

  implicit none
  integer :: dateandtime(8)
  character(strln) :: hostname
  integer :: iunitscf, iter, ierr, provided, iin
  real(dp) :: conv, din
  logical :: lflag

  call MPI_INIT_THREAD(MPI_THREAD_MULTIPLE, provided, ierr)
!~   call MPI_INIT_THREAD(MPI_THREAD_SERIALIZED, provided, ierr)
!~  call MPI_INIT(ierr)  ! call MPI_INIT here because for some reason initializing MPI in petsc breaks mpiP profiler
  call petsc_init()

  call date_and_time(values=dateandtime)
  CALL HOSTNM(hostname)
  write (pstr_out, fmt='(A)'); call petsc_print_master()
  write (pstr_out, fmt='(A)') "Trans'o'mat v1.1 petsc beta"; call petsc_print_master()
  write (pstr_out, fmt='(A)') repeat("-", 80); call petsc_print_master()
  write (pstr_out, fmt='(i2,"/",i2,"/",i4,X,i2,":",i2," @ ",A)') &
    dateandtime(2), dateandtime(3), dateandtime(1), dateandtime(5), &
    dateandtime(6), trim(hostname); call petsc_print_master()
  write (pstr_out, fmt='(A)') repeat("-", 80); call petsc_print_master()

  call init_control()
  call petsc_subcomm_init()
  call petsc_solver_init()
  call init_timings()

  iigc=0

  write (pstr_out, fmt='(A)') repeat("-", 80); call petsc_print_master()

  if (l_ionode .and. l_loadsavegf) then
!~ this is not portable at all. probaly using some C_BINDING would be better.
    call system("mkdir reaktor >/dev/null 2>/dev/null")
  end if

  if (l_ionode .and. l_ep_reuse_lambda) then
    call system("mkdir ep_reaktor >/dev/null 2>/dev/null")
  end if

  iter = 1
  conv = 666d0
  init_scf = .true.

  if (inode .eq. 0) then
!    inquire(file="K_negf_matrix2.i00.p000000",exist=oldnegf)
    inquire (file="convergence.negf", exist=init_scf)
    init_scf = .not. init_scf
    converged = .false.
    if (init_scf) then
      iter = 0
      conv = 666d0
      init_scf = .true.
    else
      open (newunit=iunitscf, file="convergence.negf", action="read", status="old")
      do
        read (unit=iunitscf, fmt=*, iostat=ierr) iin, din
        if (ierr .ne. 0) exit
        iter = iin
        conv = din
      end do
      close (iunitscf)
      converged = conv .le. scf_conv
    end if
    write (pstr_out, fmt='(A,i8,e24.12,X,e24.12,X,l)') "iter,conv,scf_conv", iter, conv, scf_conv, converged; call petsc_print_master()

  end if

  call MPI_bcast(init_scf, 1, MPI_logical, 0, PETSC_COMM_WORLD, ierr)
  call MPI_bcast(oldnegf, 1, MPI_logical, 0, PETSC_COMM_WORLD, ierr)
  call MPI_bcast(iter, 1, MPI_INTEGER, 0, PETSC_COMM_WORLD, ierr)
  call MPI_bcast(conv, 1, MPI_DOUBLE_PRECISION, 0, PETSC_COMM_WORLD, ierr)
  call MPI_bcast(converged, 1, MPI_LOGICAL, 0, PETSC_COMM_WORLD, ierr)
  write (pstr_out, fmt='(A,l)') "continue NEGF SCF: ", oldnegf; call petsc_print_master()
  write (pstr_out, fmt='(A,l)') "init_scf: ", init_scf; call petsc_print_master()

  call init_sys()

  if (k_mat_mode .eq. 1) then
    if ((oneshot .or. (converged .and. (.not. init_scf) .and. iter .ge. 2)) .and. (.not. l_ep)) then
      write (pstr_out, fmt='(A)') "get transmission"; call petsc_print_master()
      if (l_k_on_demand) then
        call transmission_k_on_demand()
      else
        call transmission()
      end if
      if (.not. oneshot) then
        open (newunit=iunitscf, file="negf.converged", action="write", status="replace")
        close (iunitscf)
      end if
      if (calc_current_density .and. (vb .ne. 0d0) .and. (converged .or. oneshot)) then
        write (pstr_out, fmt='(A)') "start scf for current density"; call petsc_print_master()
        call scf()
      end if
    else if (l_ep) then
      calc_current_density = .false.
      call loe_run()
    else
      calc_current_density = .false.
      call scf()
    end if
  else if (k_mat_mode .eq. 2) then
    if (converged) then
      open (newunit=iunitscf, file="negf.converged", action="write", status="replace")
      close (iunitscf)
    else
      call get_kmat_from_diag()
      write (pstr_out, fmt='(A)') "...dump matrix..."; call petsc_print_master()
      call dump_kmat(l_diag_fixK, "K_negf_matrix2.i00.p000000" , ecc_dir)
    end if
  else if (k_mat_mode .eq. 3) then
    call bond_order()
  end if


  write (pstr_out, fmt='(A)') "...matrix clean up..."; call petsc_print_master()
  if (.not. (k_mat_mode .eq. 3)) call petsc_cleanup()

  write (pstr_out, fmt='(A)') "...end"; call petsc_print_master()
!~   call MPI_FINALIZE(ierr)
  call slepcfinalize(ierr)
  call MPI_FINALIZE(ierr)
end program transomat
