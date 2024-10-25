module control

  ! Module for defining and managing control parameters.

  use kinds  ! Use the kinds module for defining data types.

  implicit none

  ! Structure for integer keys.
  type ikey
    character(strln) :: kname  ! Key name.
    integer :: def  ! Default value.
  end type ikey

  ! Structure for real keys.
  type rkey
    character(strln) :: kname  ! Key name.
    real(dp) :: def  ! Default value.
  end type rkey

  ! Structure for logical keys.
  type lkey
    character(strln) :: kname  ! Key name.
    logical :: def  ! Default value.
  end type lkey

  ! Structure for string keys.
  type skey
    character(strln) :: kname  ! Key name.
    character(strln) :: def  ! Default value.
  end type skey

  ! Structure containing all control parameters.
  type cntrl_

    type(skey) :: ecc_dir  ! Directory for ECC data.
    type(skey) :: elec_l_dir  ! Directory for left electrode data.
    type(skey) :: elec_r_dir  ! Directory for right electrode data.

    type(skey) :: trans_file  ! Transmission file.

    type(skey) :: cell_l_file  ! Left cell file.
    type(skey) :: cell_r_file  ! Right cell file.

    type(rkey) :: eta_cc  ! Eta for charge-charge interaction.
    type(rkey) :: eta_elec  ! Eta for electron-electron interaction.
    type(rkey) :: eps_geo  ! Geometric epsilon.
    type(rkey) :: eps_real  ! Real epsilon.
    type(rkey) :: conv_dez  ! Convergence criterion.
    type(rkey) :: estart  ! Start energy.
    type(rkey) :: eend  ! End energy.
    
    type(ikey) :: k_mat_mode  ! k-point matrix mode.

    type(ikey) :: kx_dez  ! Number of kx points.
    type(ikey) :: ky_dez  ! Number of ky points.
    type(ikey) :: n_energy_steps  ! Number of energy steps.
    type(ikey) :: i_energy_start  ! Starting energy index.
    type(ikey) :: i_energy_end  ! Ending energy index.

    type(rkey) :: ef_c  ! Fermi energy of the central region.
    type(rkey) :: ef_l  ! Fermi energy of the left electrode.
    type(rkey) :: ef_r  ! Fermi energy of the right electrode.
!~     type(rkey) :: bias  ! Bias voltage.
    type(rkey) :: bias_l  ! Bias voltage on the left electrode.
    type(rkey) :: bias_r  ! Bias voltage on the right electrode.

    type(ikey) :: integrator  ! Integrator type.
    type(ikey) :: nint_order  ! Integration order.
    type(ikey) :: maxsub  ! Maximum number of subdivisions.
    type(lkey) :: double_contour  ! Use double contour integration.
    type(rkey) :: epsfermi  ! Fermi energy broadening.
    type(rkey) :: delta_imag  ! Imaginary part of energy.
    type(rkey) :: eps_int  ! Integration epsilon.
    type(rkey) :: eps_int_contour  ! Contour integration epsilon.
    type(rkey) :: elow  ! Lower energy bound.
    type(rkey) :: temperature  ! Temperature.

    type(lkey) :: oneshot  ! One-shot calculation.

    type(lkey) :: currentdensity  ! Calculate current density.

    type(lkey) :: dftsigma  ! Use DFT sigma.
    type(rkey) :: d_occ  ! Occupation broadening for DFT sigma.
    type(rkey) :: d_virt  ! Virtual broadening for DFT sigma.
    type(ikey) :: imu_dftsigma  ! Mu index for DFT sigma.
    type(ikey) :: jmu_dftsigma  ! Mu index for DFT sigma.

    type(skey) :: f_file_ecc_petsc  ! File for ECC data in PETSc format.
    type(skey) :: s_file_ecc_petsc  ! File for ECC data in PETSc format.
    type(skey) :: d_file_ecc_petsc  ! File for ECC data in PETSc format.

    type(skey) :: f_file_elec_l_petsc  ! File for left electrode data in PETSc format.
    type(skey) :: s_file_elec_l_petsc  ! File for left electrode data in PETSc format.
    type(skey) :: d_file_elec_l_petsc  ! File for left electrode data in PETSc format.

    type(skey) :: f_file_elec_r_petsc  ! File for right electrode data in PETSc format.
    type(skey) :: s_file_elec_r_petsc  ! File for right electrode data in PETSc format.
    type(skey) :: d_file_elec_r_petsc  ! File for right electrode data in PETSc format.

    type(rkey) :: scf_conv  ! SCF convergence criterion.

    type(lkey) :: calculate_ti  ! Calculate transmission integral.

    type(lkey) :: ep_loe  ! Use linear order expansion.
    type(ikey) :: ep_active_iat1  ! Active index 1 for EP calculation.
    type(ikey) :: ep_active_iat2  ! Active index 2 for EP calculation.

    type(rkey) :: ep_bias_min  ! Minimum bias for EP calculation.
    type(rkey) :: ep_bias_max  ! Maximum bias for EP calculation.
    type(ikey) :: ep_bias_n  ! Number of bias points for EP calculation.

    type(lkey) :: ep_reuse_lambda  ! Reuse lambda for EP calculation.

    type(rkey) :: eta_ph_cc  ! Eta for phonon-charge interaction.
    type(rkey) :: eta_ph_elec  ! Eta for phonon-electron interaction.
    type(rkey) :: temperature_ph  ! Temperature for phonon calculation.

    type(ikey) :: solver_mode  ! Solver mode.
    type(ikey) :: nsim_rhs  ! Number of simultaneous RHS.
    type(ikey) :: ngroups  ! Number of groups.

    type(lkey) :: loadsave_gf  ! Load/save Green's function.
    type(lkey) :: no_noneq  ! No non-equilibrium calculation.
    
    type(lkey) :: reaktor  ! Reaktor flag.

    type(lkey) :: k_on_demand  ! k-point on demand.

    type(lkey) :: dump_nzs  ! Dump non-zero elements.
    
    type(lkey) :: diag_fix_ef  ! Fix Fermi energy during diagonalization.
    type(lkey) :: diag_fixH  ! Fix Hamiltonian during diagonalization.
    type(lkey) :: diag_fixK  ! Fix overlap matrix during diagonalization.
    type(lkey) :: diag_skip_check  ! Skip check during diagonalization.
    type(ikey) :: diag_dim  ! Dimension for diagonalization.
    type(ikey) :: diag_nkz  ! Number of kz points for diagonalization.

  end type cntrl_

  ! Interface for getting control parameters.
  interface get_cntrl_key

    ! Subroutine to get a string key.
    subroutine get_skey(keyname, keystr, l_iscritical)
      use kinds
      character(strln) :: keyname  ! Name of the key.
      character(strln) :: keystr  ! Variable to store the string value.
      logical, optional :: l_iscritical  ! Flag indicating if the key is critical.
    end subroutine get_skey

    ! Subroutine to get a real key.
    subroutine get_rkey(keyname, kvalue, l_iscritical)
      use kinds
      character(strln) :: keyname  ! Name of the key.
      real(dp) :: kvalue  ! Variable to store the real value.
      logical, optional :: l_iscritical  ! Flag indicating if the key is critical.
    end subroutine get_rkey

    ! Subroutine to get an integer key.
    subroutine get_ikey(keyname, kvalue, l_iscritical)
      use kinds
      character(strln) :: keyname  ! Name of the key.
      integer :: kvalue  ! Variable to store the integer value.
      logical, optional :: l_iscritical  ! Flag indicating if the key is critical.
    end subroutine get_ikey

    ! Subroutine to get a logical key.
    subroutine get_lkey(keyname, kvalue, l_iscritical)
      use kinds
      character(strln) :: keyname  ! Name of the key.
      logical :: kvalue  ! Variable to store the logical value.
      logical, optional :: l_iscritical  ! Flag indicating if the key is critical.
    end subroutine get_lkey

  end interface get_cntrl_key

contains

  ! Subroutine to initialize the control parameters.
  subroutine init_cntrl(cntrl)

    implicit none

    type(cntrl_) :: cntrl  ! Control parameter structure.

    ! Initialize key names for each control parameter.
    cntrl%ecc_dir%kname = "$ecc_dir="
    cntrl%elec_l_dir%kname = "$elec_left_dir="
    cntrl%elec_r_dir%kname = "$elec_right_dir="

    cntrl%trans_file%kname = "$trans_file="

    cntrl%eta_cc%kname = "$eta_cc="

    cntrl%eps_geo%kname = "$eps_geo="
    cntrl%eps_real%kname = "$eps_real="

    cntrl%eta_elec%kname = "$eta_elec="
    cntrl%conv_dez%kname = "$conv_dez="
        
    cntrl%k_mat_mode%kname = "$k_mat_mode="

    cntrl%estart%kname = "$e_start="
    cntrl%eend%kname = "$e_end="
    cntrl%n_energy_steps%kname = "$n_steps="
    cntrl%i_energy_start%kname = "$i_start="
    cntrl%i_energy_end%kname = "$i_end="

    cntrl%kx_dez%kname = "$nkx_dez="
    cntrl%ky_dez%kname = "$nky_dez="

!~     cntrl%bias%kname = "$bias="
    cntrl%bias_l%kname = "$bias_l="
    cntrl%bias_r%kname = "$bias_r="
    cntrl%integrator%kname = "$integrator="
    cntrl%nint_order%kname = "$nint_order="
    cntrl%maxsub%kname = "$maxsub="
    cntrl%double_contour%kname = "$double_contour="
    cntrl%epsfermi%kname = "$epsfermi="
    cntrl%delta_imag%kname = "$delta_imag="
    cntrl%eps_int%kname = "$eps_int="
    cntrl%eps_int_contour%kname = "$eps_int_contour="
    cntrl%elow%kname = "$elow="
    cntrl%temperature%kname = "$temperature="

    cntrl%oneshot%kname = "$oneshot="

    cntrl%currentdensity%kname = "$current_density="

    cntrl%dftsigma%kname = "$dftsigma="
    cntrl%d_occ%kname = "$dftsigma_occ="
    cntrl%d_virt%kname = "$dftsigma_virt="
    cntrl%imu_dftsigma%kname = "$dftsigma_imu="
    cntrl%jmu_dftsigma%kname = "$dftsigma_jmu="

    cntrl%f_file_ecc_petsc%kname = "$f_file_ecc_petsc="
    cntrl%s_file_ecc_petsc%kname = "$s_file_ecc_petsc="
    cntrl%d_file_ecc_petsc%kname = "$d_file_ecc_petsc="

    cntrl%f_file_elec_l_petsc%kname = "$f_file_elec_l_petsc="
    cntrl%s_file_elec_l_petsc%kname = "$s_file_elec_l_petsc="
    cntrl%d_file_elec_l_petsc%kname = "$d_file_elec_l_petsc="

    cntrl%f_file_elec_r_petsc%kname = "$f_file_elec_r_petsc="
    cntrl%s_file_elec_r_petsc%kname = "$s_file_elec_r_petsc="
    cntrl%d_file_elec_r_petsc%kname = "$d_file_elec_r_petsc="

    cntrl%scf_conv%kname = "$scf_conv="

    cntrl%calculate_ti%kname = "$calculate_ti="

    cntrl%ep_loe%kname = "$ep_loe="
    cntrl%ep_active_iat1%kname = "$ep_active_iat1="
    cntrl%ep_active_iat2%kname = "$ep_active_iat2="

    cntrl%ep_bias_min%kname = "$ep_bias_min="
    cntrl%ep_bias_max%kname = "$ep_bias_max="
    cntrl%ep_bias_n%kname = "$ep_bias_n="

    cntrl%ep_reuse_lambda%kname = "$ep_reuse_lambda="

    cntrl%eta_ph_cc%kname = "$eta_ph_cc="
    cntrl%eta_ph_elec%kname = "$eta_ph_elec="
    cntrl%temperature_ph%kname = "$temperature_ph="

    cntrl%solver_mode%kname = "$solver_mode="
    cntrl%nsim_rhs%kname = "$nsim_rhs="
    cntrl%ngroups%kname = "$ngroups="

    cntrl%loadsave_gf%kname = "$reuse_surface_gf="
    
    cntrl%no_noneq%kname = "$no_noneq="
    
    cntrl%reaktor%kname = "$reaktor="

    cntrl%k_on_demand%kname = "$k_on_demand="

    cntrl%dump_nzs%kname = "$dump_nzs="
    
    cntrl%diag_fix_ef%kname = "$diag_fix_ef="
    cntrl%diag_fixH%kname = "$diag_fixH="
    cntrl%diag_fixK%kname = "$diag_fixK="
    cntrl%diag_dim%kname = "$diag_dim="
    cntrl%diag_skip_check%kname = "$diag_skip_check="
    cntrl%diag_nkz%kname = "$diag_nkz="

  end subroutine init_cntrl

  ! Subroutine to get a control parameter from a file.
  subroutine get_key(keyname, keystr, keyint, keyreal, keylogic, l_iscritical)
    use petsc, only : PETSC_COMM_WORLD  ! Use PETSc for parallel communication.
    use kinds  ! Use kinds module for data types.
    use misc  ! Use misc module for utility functions.
    use globals, only: inode, iunit_control  ! Use globals module for global variables.
    implicit none

    character(strln) :: keyname  ! Name of the key.
    character(strln), optional :: keystr  ! String value (optional).
    integer, optional :: keyint  ! Integer value (optional).
    real(dp), optional :: keyreal  ! Real value (optional).
    logical, optional :: keylogic  ! Logical value (optional).
    logical :: l_iscritical  ! Flag indicating if the key is critical.

    integer :: io, ierr  ! I/O status and error code.
    character(strln) :: instr  ! Input string.
    logical ::found  ! Flag indicating if the key is found.
  
    character(256) :: outstr  ! Output string.

    rewind (iunit_control)  ! Rewind the control file.

    io = 0  ! Initialize I/O status.
    found = .false.  ! Initialize found flag.

    ! Read the control file line by line.
    do while (io .eq. 0)

      read (iunit_control, fmt='(A256)', IOSTAT=io) instr  ! Read a line from the file.
      if (io .ne. 0) exit  ! Exit if end of file is reached.
      instr = adjustl(instr)  ! Adjust the string to left justify.
      ! Check if the key is found in the line.
      if (index(instr, trim(keyname)) .ge. 1) then
        found = .true.  ! Set found flag to true.
        instr = instr(1:index(instr//"#", "#") - 1)  ! Remove comments from the line.
        instr = instr(index(instr, trim(keyname)) + len(trim(keyname)):strln)  ! Extract the value.
        ! Assign the value to the appropriate variable.
        if (present(keystr)) call str2(instr, keystr)
        if (present(keyint)) call str2(instr, keyint)
        if (present(keyreal)) call str2(instr, keyreal)
        if (present(keylogic)) call str2(instr, keylogic)
        write (outstr, fmt='(A)') trim(keyname)//" "//trim(instr)        ! Create output string.
        call PetscPrintf(PETSC_COMM_WORLD, trim(outstr)//New_line('A'), ierr)  ! Print the key and value.
        return  ! Return from the subroutine.
      end if

    end do

    ! If the key is not found and it is critical, stop the program.
    if ((.not. found) .and. (l_iscritical)) then
      close (iunit_control)  ! Close the control file.
      write (0, *) "non optional keyword ", trim(keyname), " not found in transp.ini"  ! Print an error message.
      stop  ! Stop the program.
    else
      ! If the key is not found and it is not critical, print a message.
      if (present(keystr)) write(instr, fmt='(A)') keystr
      if (present(keyint)) write(instr, fmt='(i16)') keyint
      if (present(keyreal)) write(instr, fmt='(e24.12)') keyreal
      if (present(keylogic)) write(instr, fmt='(l)') keylogic
      write (outstr, fmt='(A)') trim(keyname)//" "//trim(instr)        ! Create output string.
      call PetscPrintf(PETSC_COMM_WORLD, trim(outstr)//New_line('A'), ierr)  ! Print the key and value.
    end if

  end subroutine get_key

end module control

! Subroutine to get a string key.
subroutine get_skey(keyname, kvalue, l_iscritical)
  use kinds  ! Use kinds module for data types.
  use control, only: get_key  ! Use get_key subroutine from control module.
  implicit none
  character(strln) :: keyname  ! Name of the key.
  character(strln) :: kvalue  ! Variable to store the string value.
  logical, optional :: l_iscritical  ! Flag indicating if the key is critical.

  logical :: l_critical  ! Local variable for critical flag.

  l_critical = .true.  ! Set critical flag to true by default.
  if (present(l_iscritical)) l_critical = l_iscritical  ! Override if l_iscritical is present.

  call get_key(keyname, keystr=kvalue, l_iscritical=l_critical)  ! Call get_key subroutine.

end subroutine get_skey

! Subroutine to get a real key.
subroutine get_rkey(keyname, kvalue, l_iscritical)
  use kinds  ! Use kinds module for data types.
  use control, only: get_key  ! Use get_key subroutine from control module.
  implicit none
  character(strln) :: keyname  ! Name of the key.
  real(dp) :: kvalue  ! Variable to store the real value.
  logical, optional :: l_iscritical  ! Flag indicating if the key is critical.

  logical :: l_critical  ! Local variable for critical flag.

  l_critical = .true.  ! Set critical flag to true by default.
  if (present(l_iscritical)) l_critical = l_iscritical  ! Override if l_iscritical is present.

  call get_key(keyname, keyreal=kvalue, l_iscritical=l_critical)  ! Call get_key subroutine.

end subroutine get_rkey

! Subroutine to get an integer key.
subroutine get_ikey(keyname, kvalue, l_iscritical)
  use kinds  ! Use kinds module for data types.
  use control, only: get_key  ! Use get_key subroutine from control module.
  implicit none
  character(strln) :: keyname  ! Name of the key.
  integer :: kvalue  ! Variable to store the integer value.
  logical, optional :: l_iscritical  ! Flag indicating if the key is critical.

  logical :: l_critical  ! Local variable for critical flag.

  l_critical = .true.  ! Set critical flag to true by default.
  if (present(l_iscritical)) l_critical = l_iscritical  ! Override if l_iscritical is present.

  call get_key(keyname, keyint=kvalue, l_iscritical=l_critical)  ! Call get_key subroutine.

end subroutine get_ikey

! Subroutine to get a logical key.
subroutine get_lkey(keyname, kvalue, l_iscritical)
  use kinds  ! Use kinds module for data types.
  use control, only: get_key  ! Use get_key subroutine from control module.
  implicit none
  character(strln) :: keyname  ! Name of the key.
  logical :: kvalue  ! Variable to store the logical value.
  logical, optional :: l_iscritical  ! Flag indicating if the key is critical.

  logical :: l_critical  ! Local variable for critical flag.

  l_critical = .true.  ! Set critical flag to true by default.
  if (present(l_iscritical)) l_critical = l_iscritical  ! Override if l_iscritical is present.

  call get_key(keyname, keylogic=kvalue, l_iscritical=l_critical)  ! Call get_key subroutine.

end subroutine get_lkey