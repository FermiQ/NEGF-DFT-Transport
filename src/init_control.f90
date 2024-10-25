! subroutine init_control
!
! Purpose:
!   Initializes the control parameters for the simulation by reading them from the "transp.ini" file.
!   This subroutine reads various parameters related to the simulation setup, including directory paths,
!   energy range, convergence criteria, solver settings, and other options.  It also performs some
!   unit conversions and sets default values for certain parameters.
!
! Functionality:
!   1. Initializes a control structure (cntrl) using init_cntrl.
!   2. Opens the "transp.ini" file for reading.
!   3. Reads various control parameters from the "transp.ini" file using get_cntrl_key subroutine.
!   4. Performs unit conversions (e.g., converting bias voltages).
!   5. Sets default values for some boolean flags.
!   6. Handles conditional parameter reading based on other parameters (e.g., DFT sigma related parameters).
!   7. Closes the "transp.ini" file.
!   8. Prints a separator line if it's the main process (inode == 0).
!
! Parameters:
!   None.  It uses common blocks (implied through "use globals" and "use control") to access and modify global variables.

subroutine init_control()

  ! Modules used for global variables and control parameters.
  use globals
  use control

  implicit none

  ! Variable to hold the control parameters.
  type(cntrl_) :: cntrl
  ! Variable for file I/O (not explicitly used, possibly a leftover).
  integer :: io

  ! Initialize the control structure.
  call init_cntrl(cntrl)

  ! Open the control file.  newunit=iunit_control is assumed to be defined elsewhere and represents a unit number.
  open (newunit=iunit_control, file="./transp.ini", action="read", &
    position="rewind", FORM='FORMATTED', status="old")

  ! Read control parameters from the file using get_cntrl_key subroutine.  Each call reads a specific parameter.
  call get_cntrl_key(cntrl%ecc_dir%kname, ecc_dir)
  call get_cntrl_key(cntrl%elec_l_dir%kname, elec_l_dir)
  call get_cntrl_key(cntrl%elec_r_dir%kname, elec_r_dir)

  call get_cntrl_key(cntrl%trans_file%kname, trans_file)
  call get_cntrl_key(cntrl%eta_cc%kname, eta_cc)
  call get_cntrl_key(cntrl%eta_elec%kname, eta_elec)
  call get_cntrl_key(cntrl%eps_geo%kname, eps_geo)
  call get_cntrl_key(cntrl%eps_real%kname, eps_real)
  call get_cntrl_key(cntrl%conv_dez%kname, conv_dez)
  call get_cntrl_key(cntrl%kx_dez%kname, nkx_dez)
  call get_cntrl_key(cntrl%ky_dez%kname, nky_dez)
  call get_cntrl_key(cntrl%estart%kname, estart)
  call get_cntrl_key(cntrl%eend%kname, eend)
  call get_cntrl_key(cntrl%n_energy_steps%kname, n_energy_steps)
  call get_cntrl_key(cntrl%i_energy_start%kname, i_energy_start)
  call get_cntrl_key(cntrl%i_energy_end%kname, i_energy_end)
!~   call get_cntrl_key(cntrl%bias%kname, vb) ! Commented out line
  call get_cntrl_key(cntrl%bias_l%kname, vl)
  call get_cntrl_key(cntrl%bias_r%kname, vr)

!~   vb = vb/eh ! Commented out line
  vl = vl / eh ! Convert left bias voltage.  eh is assumed to be defined elsewhere.
  vr = vr / eh ! Convert right bias voltage.
  vb = vr - vl ! Calculate total bias voltage.

  call get_cntrl_key(cntrl%integrator%kname, integrator)
  call get_cntrl_key(cntrl%nint_order%kname, nint_order_el)
  call get_cntrl_key(cntrl%maxsub%kname, maxsub)
  
  ! Initialize to false, then potentially overwritten by get_cntrl_key
  ldouble_contour = .false. 
  call get_cntrl_key(cntrl%double_contour%kname, ldouble_contour, l_iscritical = .false.)
  
  ! Initialize to false, then potentially overwritten by get_cntrl_key
  l_no_noneq = .false. 
  call get_cntrl_key(cntrl%no_noneq%kname, l_no_noneq, l_iscritical = .false.)

  call get_cntrl_key(cntrl%epsfermi%kname, epsfermi)
  call get_cntrl_key(cntrl%delta_imag%kname, delta_imag)
  call get_cntrl_key(cntrl%eps_int%kname, eps_int)
  call get_cntrl_key(cntrl%eps_int_contour%kname, eps_int_contour)
  call get_cntrl_key(cntrl%elow%kname, elow)
  call get_cntrl_key(cntrl%temperature%kname, temperature_el)
  call get_cntrl_key(cntrl%oneshot%kname, oneshot)
  call get_cntrl_key(cntrl%currentdensity%kname, calc_current_density)
  call get_cntrl_key(cntrl%dftsigma%kname, dftsigma)

  ! Conditional reading of DFT sigma related parameters.
  if (dftsigma) then
    call get_cntrl_key(cntrl%d_occ%kname, d_occ)
    call get_cntrl_key(cntrl%d_virt%kname, d_virt)
    call get_cntrl_key(cntrl%imu_dftsigma%kname, imu_dftsigma)
    call get_cntrl_key(cntrl%jmu_dftsigma%kname, jmu_dftsigma)
  end if

  ! Initialize to false, then potentially overwritten by get_cntrl_key
  l_reaktor = .false. 
  call get_cntrl_key(cntrl%reaktor%kname, l_reaktor, l_iscritical = .false. )

  call get_cntrl_key(cntrl%scf_conv%kname, scf_conv)  
  
  ! Default values, potentially overwritten later.
  k_mat_mode = 1
  l_diag_skip_check = .false.
  nkz_diag = 1 
  l_diag_fixH = .true.
  call get_cntrl_key(cntrl%diag_fixH%kname, l_diag_fixH, l_iscritical = .false.)
  call get_cntrl_key(cntrl%k_mat_mode%kname, k_mat_mode, l_iscritical = .false.)
  ! Conditional reading based on k_mat_mode.
  if (k_mat_mode .eq. 2) then
    call get_cntrl_key(cntrl%diag_skip_check%kname, l_diag_skip_check)
    call get_cntrl_key(cntrl%diag_fix_ef%kname, l_diag_fix_ef)
    call get_cntrl_key(cntrl%diag_fixH%kname, l_diag_fixH)
    call get_cntrl_key(cntrl%diag_fixK%kname, l_diag_fixK)
    call get_cntrl_key(cntrl%diag_dim%kname, diag_dim)
    call get_cntrl_key(cntrl%diag_nkz%kname, nkz_diag)
  end if

  call get_cntrl_key(cntrl%calculate_ti%kname, lget_ti)

  
  call get_cntrl_key(cntrl%solver_mode%kname, solver_mode)

  ! Default values, potentially overwritten later.
  nsim_rhs = 1
  ngroups = 1
  call get_cntrl_key(cntrl%nsim_rhs%kname, nsim_rhs, l_iscritical=.false.)
  call get_cntrl_key(cntrl%ngroups%kname, ngroups, l_iscritical=.false.)
  nsim_rhs = max(1, nsim_rhs)
  ngroups = max(1, ngroups)
  call get_cntrl_key(cntrl%ep_loe%kname, l_ep_loe)
  ! Conditional reading based on l_ep_loe.
  if (l_ep_loe) then
    l_ep = .true.
    call get_cntrl_key(cntrl%ep_active_iat1%kname, ep_active_atoms(1))
    call get_cntrl_key(cntrl%ep_active_iat2%kname, ep_active_atoms(2))
    call get_cntrl_key(cntrl%ep_bias_min%kname, ep_bias(1))
    call get_cntrl_key(cntrl%ep_bias_max%kname, ep_bias(2))
    call get_cntrl_key(cntrl%ep_bias_n%kname, ep_bias_n)
    call get_cntrl_key(cntrl%ep_reuse_lambda%kname, l_ep_reuse_lambda)
    call get_cntrl_key(cntrl%eta_ph_cc%kname, eta_ph_cc)
    call get_cntrl_key(cntrl%temperature_ph%kname, temperature_ph)
!~     call get_cntrl_key(cntrl%eta_ph_elec%kname,eta_ph_elec) ! Commented out line
  end if

  call get_cntrl_key(cntrl%loadsave_gf%kname, l_loadsavegf)

  call get_cntrl_key(cntrl%dump_nzs%kname, l_dump_nzs, l_iscritical=.false.)

  ! Initialize to true, then potentially overwritten by get_cntrl_key
  l_k_on_demand = .true. 
  call get_cntrl_key(cntrl%k_on_demand%kname, l_k_on_demand, l_iscritical=.false.)

  ! Close the control file.
  close (iunit_control)

  ! Print a separator line if it's the main process. inode is assumed to be defined elsewhere.
  if (inode .eq. 0) write (6, fmt='(A)') repeat("-", 80)

end subroutine init_control