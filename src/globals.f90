module globals
  ! Include the PETSc matrix module for Fortran
  #include <petsc/finclude/petscmat.h>
  use petscmat
  use kinds
  implicit none

  ! Define a complex constant with value (0, 1)
  complex(dp), parameter :: zione = cmplx(0d0, 1d0, 8)

  ! Integer variable for global use.  Purpose not explicitly defined.
  integer :: iigc

  ! Define physical constants: pi, Hartree energy (eh), and Boltzmann constant (kb)
  real(dp) :: pi = atan(1d0)*4d0, eh = 27.211396132d0
  real(dp) :: kb = (8.6173303d-5)/27.211396132d0

  ! Flag to indicate if the system is k-symmetric (true if k-symmetric)
  logical :: ksym = .true.

  ! MPI related flags and variables
  logical :: l_ionode  ! Flag indicating if the current process is the I/O node
  logical :: l_is_big_endian  ! Flag indicating if the system uses big-endian byte order

  ! PETSc related variables
  integer :: solver_mode, nsim_rhs  ! Solver mode and number of simultaneous right-hand sides
  integer :: inode, nprocs, psubcomm, inode_sub, inode_group, ngroups, group_range(2)  ! MPI process and communication variables
  integer, allocatable :: nrow_pproc(:, :), nrow_pproc_elec_l(:, :), nrow_pproc_elec_r(:, :), nodes_group(:)  ! Allocatable arrays for distributing rows among processes
  logical :: l_loadsavegf  ! Flag to control loading/saving of Green's functions

  ! Matrix types and solver types for different components of the system
  MatType :: mattype_surf, mattype_cc, mattype_surf2, mattype_dense, mattype_sparse, mattype_electrode  ! PETSc matrix types
  MatSolverType :: matsolvertype_surf, matsolvertype_cc  ! PETSc matrix solver types
  Integer, allocatable :: cols_loc(:), nzrow_loc(:)  ! Allocatable arrays for local matrix column indices and number of non-zero elements per row

  integer :: nsim_inv  ! Number of simultaneous inversions

  ! K-point related variables
  integer :: iik  ! Index of the current k-point
  logical :: l_k_on_demand  ! Flag to indicate whether k-points are calculated on demand

  ! Integration control and parameters
  integer :: nint_order, maxsub, integrator, contour_select, ineq, int_counter, nint_order_el  ! Integration parameters: order, maximum subintervals, integrator type, contour selection, inequality type, counter, electron integration order
  real(dp) :: epsfermi, delta_imag, eps_int, elow, eps_int_contour  ! Integration parameters: Fermi energy, imaginary part of energy, integration tolerance, lower energy bound, contour integration tolerance
  logical :: ldouble_contour, lnoneq_int, l_output_progress, l_reaktor, l_no_noneq  ! Flags for integration control: double contour, non-equilibrium integration, output progress, reaktor mode, no non-equilibrium calculation

  ! Transport parameters: bias, temperature, chemical potentials
  real(dp) :: vb, mul, mur, temperature_el, ef_l, ef_c, ef_r, vl, vr  ! Transport parameters: bias voltage, chemical potential left/right, electron temperature, Fermi energy left/center/right, voltage left/right

  ! Transport parameters for phonons
  real(dp) :: eta_ph_cc, eta_ph_elec, temperature_ph  ! Phonon transport parameters: damping central/electrode, phonon temperature

  ! Transport parameters for electron-phonon coupling
  integer :: ep_active_atoms(2), n_ep_active_atoms, nat3, n_ep_modes_max, ep_bias_n  ! Electron-phonon coupling parameters: active atoms, number of active atoms, total number of atoms, maximum number of electron-phonon modes, number of bias points
  real(dp) :: ep_bias(2)  ! Electron-phonon bias
  integer, allocatable :: n_ep_modes_k(:)  ! Allocatable array for number of electron-phonon modes per k-point
  logical :: l_ep_loe, l_ep, l_ep_reuse_lambda  ! Flags for electron-phonon coupling: linear order expansion, electron-phonon coupling active, reuse lambda

  ! Transport control parameters
  real(dp) :: eta_cc, eta_elec, conv_dez, eps_geo, eps_real, estart, eend, kk(3)  ! Transport control parameters: damping central/electrode, convergence tolerance, geometric tolerance, real part tolerance, energy start/end, k-point
  character(strln) :: ecc_dir, elec_l_dir, elec_r_dir, trans_file  ! Directories and filenames for different components
  character(strln), allocatable :: species_ecc(:), species_elec_l(:), species_elec_r(:)  ! Allocatable arrays for species in different regions
  integer :: nkx_dez, nky_dez, nkz_diag, n_energy_steps, i_energy_start, i_energy_end  ! Transport control parameters: number of k-points, energy steps, energy start/end indices
  integer :: nfiles_c, nfiles_l, nfiles_r  ! Number of files for central region, left lead, right lead
  integer :: iunit_control  ! Unit number for control file

  logical :: oldnegf, oneshot, calc_current_density, converged, init_scf, lget_ti  ! Flags for transport control: old NEGF, oneshot calculation, current density calculation, convergence, SCF initialization, get transmission time

  ! Lattice parameters
  real(dp) :: dlat_l(3, 3), rlat_l(3, 3), dlat_r(3, 3), rlat_r(3, 3), dlat_c(3, 3), rlat_c(3, 3), dlat_l3(3,3), dlat_r3(3, 3)  ! Direct and reciprocal lattice vectors for left/right leads and central region
  real(dp), allocatable :: xyz_ecc(:, :), xyz_elec_l(:, :), xyz_elec_r(:, :)  ! Allocatable arrays for atomic coordinates in different regions
  real(dp), allocatable :: xyz_elec_l3(:, :), xyz_elec_r3(:, :)  ! Allocatable arrays for atomic coordinates (additional)
  real(dp), allocatable :: kp_l(:, :), kp_r(:, :), kp_c(:,:)  ! Allocatable arrays for k-points in different regions
  integer, allocatable :: wkp_r(:), wkp_l(:), wkp_c(:)  ! Allocatable arrays for k-point weights in different regions
  logical :: kpoints_from_file  ! Flag indicating whether k-points are read from a file
  integer :: nktot  ! Total number of k-points
  integer :: nred_l(2), nred_r(2)  ! Reduction factors for left and right leads

  ! System information
  character(1), allocatable :: lcr_info(:)  ! Allocatable array for system information (L, C, R)
  integer :: nat_ecc, nat_elec_l, nat_elec_r, nmu_ecc, nmu_elec_l, nmu_elec_r, nmu_c, nmu_l, nmu_r  ! Number of atoms and basis functions in different regions
  integer :: nmat_l3, nmat_r3, nat_elec_l3, nat_elec_r3  ! Additional system information
  integer :: nx_l_max, ny_l_max, nx_r_max, ny_r_max, nk_l, nk_r, nk_c, nx_max_c, ny_max_c  ! System dimensions in different regions
  integer :: ncell_c(3), ncell_l(3), ncell_r(3), nmat_c, nmat_l, nmat_r  ! Number of cells and matrices in different regions

  real(dp) :: n_electrons_ecc(2), n_electrons_l(2), n_electrons_r(2)  ! Number of electrons in different regions

  ! Mapping to Conquest format
  integer, allocatable :: imu_ecc(:), imu_elec_l(:), imu_elec_r(:), atoms_c(:), imu_to_at(:)  ! Allocatable arrays for mapping basis functions to atoms
  integer, allocatable :: tomat(:, :), inao_c(:), neigh_c(:), ndim_c(:), nlist_c(:)  ! Allocatable arrays for mapping and connectivity
  real(dp), allocatable :: tomatr(:, :), tomatr_l(:, :), tomatr_r(:, :)  ! Allocatable arrays for mapping (real part)
  integer, allocatable :: tomat_l(:, :), inao_l(:), neigh_l(:), ndim_l(:), nlist_l(:), atoms_l(:)  ! Allocatable arrays for mapping left lead
  integer, allocatable :: tomat_r(:, :), inao_r(:), neigh_r(:), ndim_r(:), nlist_r(:), atoms_r(:)  ! Allocatable arrays for mapping right lead

  integer :: ndim_total  ! Total number of dimensions

  complex(dp) :: zdosl(2), zdosr(2)  ! Complex variables for density of states left/right

  ! DFT+Sigma related variables
  real(dp) ::  d_occ, d_virt  ! Occupation and virtual parameters for DFT+Sigma
  integer :: imu_dftsigma, jmu_dftsigma  ! Basis function indices for DFT+Sigma
  logical :: dftsigma  ! Flag for DFT+Sigma calculation

  ! SCF related variables
  real(dp) :: scf_conv  ! SCF convergence parameter

  ! Debug output flag
  logical :: l_dump_nzs  ! Flag to control dumping of non-zero matrix elements

  ! K-matrix mode selection
  integer :: k_mat_mode  ! K-matrix mode

  ! Diagonalization control flags and parameters
  logical :: l_diag_fix_ef, l_diag_fixH, l_diag_fixK, l_diag_skip_check  ! Flags for diagonalization control: fix Fermi energy, Hamiltonian, K-matrix, skip check
  integer :: diag_dim  ! Dimension for diagonalization

end module globals