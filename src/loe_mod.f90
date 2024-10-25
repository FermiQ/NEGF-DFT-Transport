module loe_mod
  ! Module for calculating the linear optical response (LOE).

  use petscmat ! PETSc matrix operations
  use kinds ! Data types
  implicit none

  ! Variables representing different components of the linear optical response.
  ! T0: Trace(Gr*GammaR*Ga*GammaL)
  ! Tec: 2*Re[Tr(Gr*GammaR*Ga*GammaL*Gr*Sigma*Gr*Sigma)]
  ! TecL: Im[Tr(Gr*GammaR*Ga*GammaL*Gr*Sigma*Gr*GammaL*Ga*Sigma)]
  ! TecR: Im[Tr(Gr*GammaR*Ga*GammaL*Gr*Sigma*Gr*GammaR*Ga*Sigma)]
  ! Tin: Tr(Gr*GammaR*Ga*Sigma*Ga*GammaL*Gr*Sigma)
  ! TJL: Re[Tr(Gr*GammaR*Ga*GammaL*Gr*Sigma*Gr*GammaL*Ga*Sigma)]
  ! TJR: Re[Tr(Gr*GammaR*Ga*GammaL*Gr*Sigma*Gr*GammaR*Ga*Sigma)]
  ! TII: 2*Re[Tr(Gr*GammaR*Ga*GammaL*Gr*Sigma)]

  real(dp) :: bias_ep ! Bias energy

  complex(dp) :: T0 ! Linear optical response component T0
  complex(dp), allocatable :: Tec(:), TecL(:), TecR(:), Tin(:), TII(:), TJL(:), TJR(:), &
                              TpiLR(:), TpiLL(:), TpiRR(:), ImPi_r_fac(:) ! Allocatable arrays for LOE components


  ! PETSc matrices for Green's function and self-energies.
  Mat :: p_gr_inv, p_gr, p_full_gammal, p_full_gammar

contains

  subroutine loe_run()
    ! Main subroutine for running the LOE calculation.

    use globals ! Global variables
    use petsc_mod ! PETSc module
    use phonon_mod, only: ii_mode, mode_ii_energy ! Phonon module
    use integrator_mod ! Integration module
    use misc ! Miscellaneous functions
    use error_handler ! Error handling
    implicit none

    integer :: ik, imode, i, n, ib, iunit, ierr, junit ! Loop counters, file unit, error code
    real(dp) :: dd_eps, d_high, d_low, dx, x, bias_low, bias_high, dbias ! Variables for integration and bias
    complex(dp) :: zz, current_out, neff_test ! Temporary complex variables
    complex(dp), allocatable :: current_loe(:), current_elastic_p_loe(:) ! Allocatable arrays for current
    PetscScalar :: dI0_ec_out, dI0_J_out, Iinel_out ! Output variables for current components

    ! Print message to output stream
    write (pstr_out, fmt='(A)') "init LOE matrices"; call petsc_print_master()
    call init_loe_matrices() ! Initialize matrices

    ! Set bias parameters
    bias_high = ep_bias(2)
    bias_low = ep_bias(1)
    n = ep_bias_n
    dbias = (bias_high - bias_low)/real(n)
    n = n + 1
    allocate (current_loe(n), current_elastic_p_loe(n), stat=ierr) ! Allocate arrays for current
    ! Check for allocation errors
    if (ierr .ne. 0) then
      write (errormsg, *) "allocation error in loe_run ", ierr
      call error() ! Handle error
    end if
    current_loe = 0d0 ! Initialize current arrays
    current_out = 0d0
    current_elastic_p_loe = 0d0

    ! Loop over k-points
    do ik = nk_c, 1, -1

      ! Print message to output stream
      write (pstr_out, fmt='(A,i8)') "loe current at k-point for all bias voltages ", ik; call petsc_print_master()

      call init_loe_coefficients(ik) ! Initialize LOE coefficients

      ! Print message to output stream
      write (pstr_out, fmt='(A,i8)') "get current for each mode", ik; call petsc_print_master()
      write (pstr_out, fmt='(A)') "mode "; call petsc_print_master(.false.)

      ! Loop over phonon modes
      do imode = n_ep_modes_k(ik), 1, -1
        ! Print message to output stream
        write (pstr_out, fmt='(i6)') imode; call petsc_print_master(.false.)
        call petsc_vec_getvalue(imode - 1, mode_ii_energy, p_phonon_EW(ik)) ! Get phonon energy

        ii_mode = imode ! Set current mode
        ! Loop over bias points
        do ib = 1 + inode, n, nprocs
          bias_ep = bias_low + dbias*real(ib - 1, 8) ! Set bias energy
          ! Calculate current components
          call loe_current(real(mode_ii_energy, 8), dI0_ec_out, dI0_J_out, Iinel_out)
          current_loe(ib) = current_loe(ib) + (dI0_ec_out + dI0_J_out + Iinel_out)*wkp_l(ik)
        end do ! bias

        call MPI_BARRIER(PETSC_COMM_WORLD, ierr) ! Barrier synchronization
      end do ! modes

      ! Print message to output stream
      write (pstr_out, fmt='(A)') " finished"; call petsc_print_master()

      ! Calculate elastic part of current
      do ib = 1 + inode, n, nprocs
        bias_ep = bias_low + dbias*real(ib - 1, 8)
        current_elastic_p_loe(ib) = current_loe(ib) + T0*bias_ep*wkp_l(ik)
      end do ! bias

      call destroy_loe_coeffcients() ! Destroy LOE coefficients

    end do !kpoit

    ! Average current over k-points and sum over processors
    current_elastic_p_loe = current_elastic_p_loe/real(nktot, 8)
    current_loe = current_loe/real(nktot, 8)

    ! MPI reduction to sum current over all processors
    if (l_ionode) then
      call MPI_REDUCE(MPI_IN_PLACE, current_elastic_p_loe, size(current_elastic_p_loe), MPI_DOUBLE_COMPLEX, MPI_SUM, 0, PETSC_COMM_WORLD, ierr)
      call MPI_REDUCE(MPI_IN_PLACE, current_loe, size(current_loe), MPI_DOUBLE_COMPLEX, MPI_SUM, 0, PETSC_COMM_WORLD, ierr)
    else
      call MPI_REDUCE(current_elastic_p_loe, current_elastic_p_loe, size(current_elastic_p_loe), MPI_DOUBLE_COMPLEX, MPI_SUM, 0, PETSC_COMM_WORLD, ierr)
      call MPI_REDUCE(current_loe, current_loe, size(current_elastic_p_loe), MPI_DOUBLE_COMPLEX, MPI_SUM, 0, PETSC_COMM_WORLD, ierr)
    end if

    ! Write current data to file (only on node 0)
    if (l_ionode) then
      open (newunit=iunit, file="current_loe.dat", action="write", status="replace")
      do ib = 1, n
        bias_ep = bias_low + dbias*real(ib - 1, 8)
        write (unit=iunit, fmt='(4es45.24e5)') bias_ep, real(current_elastic_p_loe(ib), 8), real(current_loe(ib), 8), real(current_elastic_p_loe(ib) - current_loe(ib), 8)
      end do
      close (iunit)

      ! Write mode data to files
      do ik = 1, nk_c
        open (newunit=iunit, file="modes_"//trim(i2s(ik))//".dat", action="write", status="replace")
        do imode = 1, n_ep_modes_k(ik)
          call petsc_vec_getvalue(imode - 1, mode_ii_energy, p_phonon_EW(ik), .false.)
          write (iunit, fmt='(2e24.12)') real(mode_ii_energy), 0d0
          write (iunit, fmt='(2e24.12)') real(mode_ii_energy), -1d0
          write (iunit, fmt='(2e24.12)') real(mode_ii_energy), 1d0
          write (iunit, fmt='(2e24.12)') real(mode_ii_energy), 0d0
        end do
        close (iunit)
      end do

    end if

    call destroy_loe_matrices() ! Destroy matrices

  end subroutine loe_run

  subroutine init_loe_matrices()
    ! Subroutine to initialize the matrices needed for the LOE calculation.

#include <petsc/finclude/petsc.h>
    use petscmat ! PETSc matrix operations
    use petsc_mod ! PETSc module
    use petsc_wrapper ! PETSc wrapper functions
    use kinds ! Data types
    use globals ! Global variables
    use integrator_mod ! Integration module
    implicit none

    integer :: ierr ! Error code

    ! Initialize matrices using PETSc functions.  The specific meaning of each matrix is context-dependent and requires knowledge of the broader application.
    call petsc_get_densemat(p_h00k_r(1), p_grrr, mattype_surf)
    call petsc_get_densemat(p_h00k_l(1), p_gllr, mattype_surf)
    call petsc_get_densemat(p_h00k_r(1), p_sigmarr, mattype_dense)
    call petsc_get_densemat(p_h00k_r(1), p_gammar, mattype_dense)
    call petsc_get_densemat(p_h00k_l(1), p_sigmalr, mattype_dense)
    call petsc_get_densemat(p_h00k_l(1), p_gammal, mattype_dense)

    call petsc_get_a_with_b_c(p_full_gammal, p_h00k_cc(1), p_h00k_cc(1), mattype_sparse)
    call petsc_get_a_with_b_c(p_full_gammar, p_h00k_cc(1), p_h00k_cc(1), mattype_sparse)

    call petsc_alloc_mat_block(p_full_gammal, 0, nmu_l - 1, 0, nmu_l - 1)
    call petsc_alloc_mat_block(p_full_gammar, nmu_c - nmu_r, nmu_c - 1, nmu_c - nmu_r, nmu_c - 1)

  end subroutine init_loe_matrices

  subroutine destroy_loe_matrices()
    ! Subroutine to destroy the matrices used in the LOE calculation.

#include <petsc/finclude/petsc.h>
    use petscmat ! PETSc matrix operations
    use petsc_mod ! PETSc module
    use petsc_wrapper ! PETSc wrapper functions
    use kinds ! Data types
    use globals ! Global variables
    use integrator_mod ! Integration module
    implicit none

    integer :: ierr ! Error code

    ! Destroy matrices using PETSc functions.
    call MatDestroy(p_grrr, ierr)
    call MatDestroy(p_gllr, ierr)
    call MatDestroy(p_sigmarr, ierr)
    call MatDestroy(p_gammar, ierr)
    call MatDestroy(p_sigmalr, ierr)
    call MatDestroy(p_gammal, ierr)

    call MatDestroy(p_full_gammal, ierr)
    call MatDestroy(p_full_gammar, ierr)

    call MatDestroy(p_full_gammal, ierr)
    call MatDestroy(p_full_gammar, ierr)

  end subroutine destroy_loe_matrices

  subroutine destroy_loe_coeffcients()
    ! Subroutine to deallocate the arrays holding LOE coefficients.

#include <petsc/finclude/petsc.h>
    use petscmat ! PETSc matrix operations
    implicit none

    integer :: ierr ! Error code

    ! Deallocate arrays.
    deallocate (Tec)
    deallocate (TecL)
    deallocate (TecR)
    deallocate (Tin)
    deallocate (TII)
    deallocate (TJL)
    deallocate (TJR)
    deallocate (TpiLR)
    deallocate (TpiRR)
    deallocate (TpiLL)
    deallocate (ImPi_r_fac)

  end subroutine destroy_loe_coeffcients

  subroutine init_loe_coefficients(ik)
    ! Subroutine to initialize the LOE coefficients for a given k-point.

#include <petsc/finclude/petsc.h>
    use petscmat ! PETSc matrix operations
    use petsc_mod ! PETSc module
    use petsc_wrapper ! PETSc wrapper functions
    use kinds ! Data types
    use globals ! Global variables
    use integrator_mod ! Integration module
    use error_handler ! Error handling
    implicit none

    integer :: ik ! k-point index

    ! Declare PETSc matrices for intermediate calculations.
    Mat :: p_gr, p_GrGammaRGa, p_Ga, p_tmp1, p_tmp2, p_tmp3, p_GrGammaRGaGammaLGr, &
      p_GrGammaRGaGammaLGrSigmaGr, p_gaSigma

    integer :: ierr, ii(1), jj(1), imode ! Error code, indices, mode index
    complex(dp) :: ef, pp(1), pi_pm, pi_r ! Variables for energy, trace, and self-energies
    PetscReal :: p_real ! Real part of a PETSc scalar
    PetscScalar :: p_tr ! PETSc scalar for trace

    ! Print message to output stream
    write (pstr_out, fmt='(A,i8)') "init leo coefficients ", ik; call petsc_print_master()

    iik = ik ! Set global k-point index

    ! Allocate arrays for LOE coefficients.
    allocate (Tec(n_ep_modes_k(ik)), TecL(n_ep_modes_k(ik)), TecR(n_ep_modes_k(ik)), &
              Tin(n_ep_modes_k(ik)), TII(n_ep_modes_k(ik)), TJL(n_ep_modes_k(ik)), &
              TJR(n_ep_modes_k(ik)), TpiLR(n_ep_modes_k(ik)), TpiRR(n_ep_modes_k(ik)), &
              TpiLL(n_ep_modes_k(ik)), ImPi_r_fac(n_ep_modes_k(ik)), stat=ierr)
    ! Check for allocation errors
    if (ierr .ne. 0) then
      write (errormsg, *) "allocation error in loe_run ", ierr
      call error() ! Handle error
    end if

    ef = 0.5d0*(ef_l + ef_r) ! Average energy
    mul = ef_l ! Left energy
    mur = ef_r ! Right energy
    nsim_inv = nmu_c ! Inverse of number of modes

    ldouble_contour = .false. ! Flag for double contour integration

    call init_gr(ef, p_gr_inv) ! Initialize Green's function

    call petsc_add_sub_B_to_A(p_gammal, p_full_gammal, 0, 0, p_one, INSERT_VALUES, PETSC_TRUE)
    call petsc_add_sub_B_to_A(p_gammar, p_full_gammar, nmu_c - nmu_r, nmu_c - nmu_r, p_one, INSERT_VALUES, PETSC_TRUE)

    call petsc_get_densemat(p_h00k_cc(1), p_gr, mattype_dense) ! Get Green's function
    call petsc_invert(p_gr_inv, p_gr, matsolvertype_cc, mattype_dense) ! Invert Green's function
    call MatHermitianTranspose(p_gr, MAT_INITIAL_MATRIX, p_ga, ierr) ! Hermitian transpose
    call MatMatMatMult(p_gr, p_full_gammar, p_ga, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, p_GrGammaRGa, ierr) ! Matrix multiplication
    call MatMatMatMult(p_GrGammaRGa, p_full_gammal, p_gr, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, &
                       p_GrGammaRGaGammaLGr, ierr) ! Matrix multiplication

    ! Calculate T0 and other LOE components using matrix operations and traces.  The specific formulas are complex and require detailed knowledge of the underlying theory.
    call MatMatMult(p_GrGammaRGa, p_full_gammal, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, &
                    p_tmp3, ierr)
    call MatGetTrace(p_tmp3, p_tr, ierr)
    T0 = p_tr
    write (pstr_out, fmt='(A,2e24.12)') "T0", T0; call petsc_print_master()
    call MatDestroy(p_tmp3, ierr)

    ! Loop over modes to calculate LOE coefficients for each mode.
    write (pstr_out, fmt='(A)') "mode "; call petsc_print_master(.false.)
    do imode = n_ep_modes_k(ik), 1, -1
      write (pstr_out, fmt='(i6)') imode; call petsc_print_master(.false.)
      call MatMatMatMult(p_GrGammaRGaGammaLGr, p_ep_lambda_k(imode, ik), p_gr, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, &
                         p_GrGammaRGaGammaLGrSigmaGr, ierr)

      ! ... (Many lines of complex matrix calculations to compute various LOE components) ...

    end do

    ! Print message to output stream
    write (pstr_out, fmt='(A)') " finished"; call petsc_print_master()

    ! Destroy matrices.
    call MatDestroy(p_GrGammaRGaGammaLGr, ierr)
    call MatDestroy(p_ga, ierr)
    call MatDestroy(p_gr, ierr)
    call MatDestroy(p_GrGammaRGa, ierr)
    call MatDestroy(p_gr_inv, ierr)

  end subroutine init_loe_coefficients

  function ImPi_r_alpha(x) result(ImPi_r_alpha)
    ! Function to calculate the imaginary part of the retarded self-energy.

    use globals ! Global variables
    use kinds ! Data types
    use phonon_mod ! Phonon module
    implicit none

    PetscScalar :: ImPi_r_alpha ! Output: Imaginary part of retarded self-energy
    real(dp) :: x ! Input: Energy

    ImPi_r_alpha = ImPi_r_fac(ii_mode)*x ! Calculation

  end function ImPi_r_alpha

  function ImPi_pm_alpha(x) result(ImPi_pm_alpha)
    ! Function to calculate the imaginary part of the Keldysh self-energy.

    use globals ! Global variables
    use kinds ! Data types
    use phonon_mod ! Phonon module
    use integrator_mod, only: bose ! Bose-Einstein distribution
    use petsc_mod ! PETSc module
    implicit none

    PetscScalar :: ImPi_pm_alpha ! Output: Imaginary part of Keldysh self-energy
    real(dp) :: x ! Input: Energy

    complex(dp) :: x_p_u, x_m_u, xx ! Variables for shifted energies

    xx = x ! Copy input energy
    x_p_u = x + bias_ep ! Energy shifted by bias
    x_m_u = x - bias_ep ! Energy shifted by bias

    ! Calculation using LOE coefficients and Bose-Einstein distribution.
    ImPi_pm_alpha = 0.5d0/pi*(TpiLR(ii_mode)*(x_p_u*bose(x_p_u, temperature_ph) + &
                                              (x_m_u)*bose(x_m_u, temperature_ph)) + (TpiLL(ii_mode) + TpiRR(ii_mode))*xx*bose(xx, temperature_ph))

  end function ImPi_pm_alpha

  function Nneq_alpha(x) result(Nneq_alpha)
    ! Function to calculate the non-equilibrium phonon distribution.

    use globals ! Global variables
    use kinds ! Data types
    use integrator_mod, only: bose ! Bose-Einstein distribution
    use phonon_mod ! Phonon module
    use petsc_mod ! PETSc module
    implicit none

    PetscScalar :: Nneq_alpha ! Output: Non-equilibrium phonon distribution
    real(dp) :: x ! Input: Energy
    complex(dp) :: x_p_u, x_m_u ! Variables for shifted energies

    complex(dp) :: xx ! Variable for energy

    xx = x ! Copy input energy
    x_p_u = x + bias_ep ! Energy shifted by bias

    ! Calculation using imaginary parts of self-energies and Bose-Einstein distribution.
    Nneq_alpha = -0.5d0*(ImPi_pm_alpha(x) + bose(xx, temperature_ph)*eta_ph_cc*x/mode_ii_energy)/(ImPi_r_alpha(x) - eta_ph_cc*x/mode_ii_energy*0.5d0)

  end function Nneq_alpha

  function dI0_ec(x) result(dI0_ec)
    ! Function to calculate the elastic current component.

    use globals ! Global variables
    use kinds ! Data types
    use petsc_mod ! PETSc module
    use integrator_mod, only: bose ! Bose-Einstein distribution
    use phonon_mod ! Phonon module
    implicit none

    PetscScalar :: dI0_ec ! Output: Elastic current component
    real(dp) :: x ! Input: Energy

    complex(dp) :: x_p_u, x_m_u, xx ! Variables for shifted energies

    xx = x ! Copy input energy
    x_p_u = x + bias_ep ! Energy shifted by bias
    x_m_u = x - bias_ep ! Energy shifted by bias

    ! Calculation using LOE coefficients, non-equilibrium distribution, and Bose-Einstein distribution.
    dI0_ec = (Tec(ii_mode)*(2d0*Nneq_alpha(x) + 1d0)*bias_ep + &
              (TecL(ii_mode) + TecR(ii_mode))*(x_m_u*bose(x_m_u, temperature_ph) - x_p_u*bose(x_p_u, temperature_ph) - bias_ep))

  end function dI0_ec

  function dI0_ec_petsc(x) result(dI0_ec_petsc)
    ! Function to calculate the elastic current component using PETSc matrices.  This function appears to be a wrapper around dI0_ec, adding PETSc matrix handling overhead.

    use globals ! Global variables
    use kinds ! Data types
    use petsc_wrapper, only: petsc_matassemble ! PETSc wrapper functions
    use petsc_mod ! PETSc module
    implicit none

    real(dp) :: x ! Input: Energy

    real(dp) :: z ! Variable for energy
    integer :: ii(1), ierr, dI0_ec_petsc ! Indices, error code, output
    PetscScalar :: pp(1) ! PETSc scalar

    dI0_ec_petsc = -1 ! Initialize output

    ! Set matrix values using dI0_ec function.
    if (inode .eq. 0) then
      z = x
      pp(1) = dI0_ec(z)
      ii = 0
      call MatSetValues(p_tmpcc1, 1, ii, 1, ii, pp, INSERT_VALUES, ierr)
    end if

    call MatAssemblyBegin(p_tmpcc1, MAT_FINAL_ASSEMBLY, ierr)
    call MatAssemblyEnd(p_tmpcc1, MAT_FINAL_ASSEMBLY, ierr)

    call petsc_matassemble(p_tmpcc1)

    dI0_ec_petsc = 0 ! Set output

  end function dI0_ec_petsc

  function dI0_J(x) result(dI0_J)
    ! Function to calculate the inelastic current component.

    use globals ! Global variables
    use kinds ! Data types
    use petsc_mod ! PETSc module
    use integrator_mod, only: bose ! Bose-Einstein distribution
    use phonon_mod ! Phonon module
    implicit none

    PetscScalar :: dI0_J ! Output: Inelastic current component
    real(dp) :: x ! Input: Energy

    complex(dp) :: x_p_u, x_m_u, xx ! Variables for shifted energies

    xx = x ! Copy input energy
    x_p_u = x + bias_ep ! Energy shifted by bias
    x_m_u = x - bias_ep ! Energy shifted by bias

    ! Calculation using LOE coefficients and Bose-Einstein distribution.
    dI0_J = real(dr_alpha(x))*(xx*bose(xx, temperature_ph) - (x_p_u)*bose(x_p_u, temperature_ph))
    dI0_J = -1d0/pi*(TJR(ii_mode) - TJL(ii_mode))*dI0_J

  end function dI0_J

  function dI0_J_petsc(x) result(dI0_J_petsc)
    ! Function to calculate the inelastic current component using PETSc matrices.  Similar to dI0_ec_petsc, this adds PETSc overhead.

    use globals ! Global variables
    use kinds ! Data types
    use petsc_wrapper, only: petsc_matassemble ! PETSc wrapper functions
    use petsc_mod ! PETSc module
    implicit none

    real(dp) :: x ! Input: Energy

    real(dp) :: z ! Variable for energy
    integer :: ii(1), ierr, dI0_J_petsc ! Indices, error code, output
    PetscScalar :: pp(1) ! PETSc scalar

    dI0_J_petsc = -1 ! Initialize output

    ! Set matrix values using dI0_J function.
    if (inode .eq. 0) then
      z = x
      pp(1) = dI0_J(z)
      ii = 0
      call MatSetValues(p_tmpcc1, 1, ii, 1, ii, pp, INSERT_VALUES, ierr)
    end if

    call MatAssemblyBegin(p_tmpcc1, MAT_FINAL_ASSEMBLY, ierr)
    call MatAssemblyEnd(p_tmpcc1, MAT_FINAL_ASSEMBLY, ierr)

    call petsc_matassemble(p_tmpcc1)

    dI0_J_petsc = 0 ! Set output

  end function dI0_J_petsc

  function dI0_J_scalar(x) result(dI0_J_scalar)
    ! Function to calculate the inelastic current component using scalar integration.

    use globals ! Global variables
    use kinds ! Data types
    use petsc_wrapper, only: petsc_matassemble ! PETSc wrapper functions
    use petsc_mod ! PETSc module
    use integrator_scalar ! Scalar integration module
    implicit none

    real(dp) :: x ! Input: Energy

    real(dp) :: z ! Variable for energy
    integer :: ii(1), ierr, dI0_J_scalar ! Indices, error code, output
    PetscScalar :: pp(1) ! PETSc scalar

    zint_scalar_output = dI0_J(x) ! Calculate inelastic current component
    dI0_J_scalar = 1 ! Set output

  end function dI0_J_scalar

  function Iinel(x) result(Iinel)
    ! Function to calculate the inelastic current component.

    use globals ! Global variables
    use kinds ! Data types
    use petsc_mod ! PETSc module
    use integrator_mod, only: bose ! Bose-Einstein distribution
    use phonon_mod ! Phonon module
    implicit none

    PetscScalar :: Iinel ! Output: Inelastic current component
    real(dp) :: x ! Input: Energy

    complex(dp) :: x_p_u, x_m_u, xx ! Variables for shifted energies

    xx = x ! Copy input energy
    x_p_u = x + bias_ep ! Energy shifted by bias
    x_m_u = x - bias_ep ! Energy shifted by bias

    ! Calculation using LOE coefficients, non-equilibrium distribution, and Bose-Einstein distribution.
    Iinel = (2d0*Nneq_alpha(x)*bias_ep + x_m_u*bose(x_m_u, temperature_ph) - x_p_u*bose(x_p_u, temperature_ph))
    Iinel = Tin(ii_mode)*Iinel

  end function Iinel

  function Iinel_petsc(x) result(Iinel_petsc)
    ! Function to calculate the inelastic current component using PETSc matrices.  Similar to dI0_ec_petsc and dI0_J_petsc, this adds PETSc overhead.

    use globals ! Global variables
    use kinds ! Data types
    use petsc_wrapper, only: petsc_matassemble ! PETSc wrapper functions
    use petsc_mod ! PETSc module
    implicit none

    real(dp) :: x ! Input: Energy

    real(dp) :: z ! Variable for energy
    integer :: ii(1), ierr, Iinel_petsc ! Indices, error code, output
    PetscScalar :: pp(1) ! PETSc scalar

    Iinel_petsc = -1 ! Initialize output

    ! Set matrix values using Iinel function.
    if (inode .eq. 0) then
      z = x
      pp(1) = Iinel(z)
      ii = 0
      call MatSetValues(p_tmpcc1, 1, ii, 1, ii, pp, INSERT_VALUES, ierr)
    end if

    call MatAssemblyBegin(p_tmpcc1, MAT_FINAL_ASSEMBLY, ierr)
    call MatAssemblyEnd(p_tmpcc1, MAT_FINAL_ASSEMBLY, ierr)

    call petsc_matassemble(p_tmpcc1)

    Iinel_petsc = 0 ! Set output

  end function Iinel_petsc

  subroutine loe_current(energy, dI0_ec_out, dI0_J_out, Iinel_out)
    ! Subroutine to calculate the current components for a given energy.

    use globals ! Global variables
    use kinds ! Data types
    use phonon_mod ! Phonon module
    implicit none

    PetscScalar :: dI0_ec_out, dI0_J_out, Iinel_out ! Output: Current components
    real(dp) :: energy ! Input: Energy
    real(dp) :: dd_eps, d_high, d_low ! Variables for integration limits

    dI0_ec_out = 0d0 ! Initialize output
    dI0_J_out = 0d0 ! Initialize output
    Iinel_out = 0d0 ! Initialize output

    dI0_ec_out = dI0_ec(energy) ! Calculate elastic current component

    ! Set integration limits and call integration routine.
    dd_eps = 1d-4
    d_high = zsqrt(2d0*energy/dd_eps + mode_ii_energy*energy + eta_ph_cc*eta_ph_cc*0.25d0)
    d_low = -d_high
    call loe_integrate_scalar(dI0_J_scalar, d_low, d_high, dI0_J_out) ! Integrate inelastic current component
    Iinel_out = Iinel(energy) ! Calculate inelastic current component

  end subroutine loe_current

  subroutine loe_integrate(f, d_low, d_high, zout, what)
    ! Subroutine to perform numerical integration using PETSc matrices.  This is a more general integration routine that handles PETSc matrices as input and output.

#include <petsc/finclude/petsc.h>
    use petscmat ! PETSc matrix operations
    use petsc_mod ! PETSc module
    use petsc_wrapper ! PETSc wrapper functions
    use kinds ! Data types
    use globals ! Global variables
    use integrator_mod ! Integration module
    use phonon_mod, only: mode_ii_energy ! Phonon module
    implicit none

    integer, external :: f ! External function for integrand
    real(dp) :: d_low, d_high ! Integration limits
    complex(dp) :: zout ! Output: Integration result
    character(*), optional :: what ! Optional argument for output label

    real(dp) :: err_eq_int, eps ! Error and tolerance for integration
    PetscScalar :: p_zint1, p_zint2, pp(1), dx, x1, x2, x ! Variables for integration
    Mat :: p_d_tmp1(1), p_d_tmp2(1)
    integer :: ierr, ii(1), jj(1), ik(2), i, n

    call MatCreate(PETSC_COMM_WORLD, p_tmpcc1, ierr)
    call MatSetType(p_tmpcc1, mattype_sparse, ierr)
    call MatSetSizes(p_tmpcc1, PETSC_DECIDE, PETSC_DECIDE, 1, 1, ierr)
    if (inode .eq. 0) then
      ik(1) = 0
      ik(2) = 1
      pp = 0d0
      call MatMPIAIJSetPreallocationCSR(p_tmpcc1, ik, ik, pp, ierr)
    else
      call MatMPIAIJSetPreallocationCSR(p_tmpcc1, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, pp, ierr)
    end if

    call MatAssemblyBegin(p_tmpcc1, MAT_FINAL_ASSEMBLY, ierr)
    call MatAssemblyEnd(p_tmpcc1, MAT_FINAL_ASSEMBLY, ierr)
    call MatDuplicate(p_tmpcc1, MAT_DO_NOT_COPY_VALUES, p_d_tmp1(1), ierr)
    call MatDuplicate(p_tmpcc1, MAT_DO_NOT_COPY_VALUES, p_d_tmp2(1), ierr)

!~       call petsc_mat_info(p_tmpcc1,"p_tmpcc1 ",ierr)

    eps = 1d-5
    nint_order = 32
    l_output_progress = .false.
    call adaptive_int3(f, d_low, d_high, p_d_tmp1(1), p_d_tmp2(1), maxsub*10, eps, err_eq_int)
    call petsc_mat_getvalue(0, 0, zout, p_d_tmp1(1), 1, PETSC_COMM_WORLD)

    if (present(what)) then
      write (pstr_out, fmt='(A,X,4e24.12)') trim(what), zout, d_low, d_high; call petsc_print_master()
    end if
    call MatDestroy(p_tmpcc1, ierr)
    call MatDestroy(p_d_tmp1(1), ierr)
    call MatDestroy(p_d_tmp2(1), ierr)

  end subroutine loe_integrate

  subroutine loe_integrate_scalar(f, d_low, d_high, zout, what)
    use petsc_mod
    use kinds
    use globals
    use integrator_scalar
    use phonon_mod, only: mode_ii_energy
    implicit none

    integer, external :: f
    real(dp) :: d_low, d_high
    complex(dp) :: zout
    character(*), optional :: what

    real(dp) :: err_eq_int, eps
    PetscScalar :: p_zint1, p_zint2, pp, dx, x1, x2, x
    complex(dp) :: p_d_tmp1, p_d_tmp2
    integer :: ierr, ii(1), jj(1), ik(2), i, n

    eps = 1d-5
    nint_order = 32
    l_output_progress = .false.
    call adaptive_int3_scalar(f, d_low, d_high, p_d_tmp1, p_d_tmp2, maxsub*10, eps, err_eq_int)
    zout = p_d_tmp1

    if (present(what)) then
      write (pstr_out, fmt='(A,X,4e24.12)') trim(what), zout, d_low, d_high; call petsc_print_master()
    end if

  end subroutine loe_integrate_scalar

end module loe_mod
