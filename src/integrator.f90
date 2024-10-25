module integrator_mod
  use kinds
  use misc
  use globals
!~   use blas95
!~   use lapack95
  use mathstuff
  use surf_gf_mod
  use petsc_mod, only: pstr_out, petsc_print_master
  implicit none

  ! Global variables for the integrator module
  integer :: intlr, n_fermipols  ! intlr: Integer variable, likely representing left/right electrode. n_fermipols: Number of Fermi poles.
  real(dp) :: r_eq, e0, rpfix, ipfix, phi_low, eu, el ! r_eq: Equilibrium distance. e0: Energy offset. rpfix, ipfix: Real and imaginary fixes. phi_low: Lower phi value. eu, el: Upper and lower energies.
  complex(dp) :: y, x_low, y_low, dx_low, dy_low ! y: Complex variable. x_low, y_low: Lower x and y values. dx_low, dy_low: Differences in x and y.


contains

  subroutine get_alpha(p_Dl, p_Dr, dc_weights)
#include <petsc/finclude/petscmat.h>
    use petsc_mod
    implicit none

    ! Input: PETSc matrices representing left and right contributions
    Mat :: p_Dl, p_Dr 
    ! Output: Allocatable array of PETSc scalars representing weights
    PetscScalar, allocatable :: dc_weights(:)

    PetscScalar :: a_l, a_r, p_val ! a_l, a_r: Sum of diagonal elements for left and right matrices. p_val: Temporary variable to store matrix values.
    integer :: imu1, imu2, iat1, ierr, imu3 ! Loop indices and error code.
    Mat :: p_tmp_l, p_tmp_r ! Temporary PETSc matrices (not used in the current implementation).
    
    imu2 = 0
    imu3 = 0
    ! Loop over atoms
    do iat1 = 1, nat_ecc
      a_l = 0d0
      a_r = 0d0
      ! Loop over modes for each atom
      do imu1 = 1, imu_ecc(iat1)
        imu2 = imu2 + 1
        ! Get diagonal element from p_Dl
        call petsc_mat_getvalue(imu2 - 1, imu2 - 1, p_val, p_Dl, 1, PETSC_COMM_WORLD)
        a_l = a_l + p_val
        ! Get diagonal element from p_Dr
        call petsc_mat_getvalue(imu2 - 1 , imu2 - 1 , p_val, p_Dr, 1, PETSC_COMM_WORLD)
        a_r = a_r + p_val
      end do
      ! Calculate weight
      p_val = a_l / (a_l + a_r)
      ! Assign weight to dc_weights
      do imu1 = 1, imu_ecc(iat1)
        imu3 = imu3 + 1       
        dc_weights(imu3) = p_val
      end do
    end do
    
!~     do imu2 = 1, nmu_c
!~       write(0, fmt='(i8,2e24.12)') imu2, dc_tmp(imu2)
!~       do imu1 = 1, nmu_c
!~         dc_weights(imu1, imu2) = dc_tmp(imu1) * dc_tmp(imu2)
!~       end do
!~     end do
!~     write(0, fmt='(4e24.12)') &
!~    &  maxval(real(dc_weights)),minval(real(dc_weights)),&
!~    &  maxval(aimag(dc_weights)),minval(aimag(dc_weights))

  end subroutine get_alpha

  subroutine get_weight(wl, wr, alphal, alphar)

    implicit none

    ! Input/Output: Complex matrices and vectors representing weights
    complex(dp), allocatable :: wl(:, :), wr(:, :), alphal(:), alphar(:)
    complex(dp) :: znorm ! Temporary complex variable for normalization
    integer :: iat1, iat2 ! Loop indices

    ! Calculate weights
    do iat2 = 1, nat_ecc
      do iat1 = 1, nat_ecc
        znorm = zsqrt(alphal(iat1)*alphal(iat2)) + zsqrt(alphar(iat1)*alphar(iat2))
        wl(iat1, iat2) = zsqrt(alphal(iat1))*zsqrt(alphal(iat2))/znorm
        wr(iat1, iat2) = zsqrt(alphar(iat1))*zsqrt(alphar(iat2))/znorm
      end do
    end do

  end subroutine get_weight

  subroutine init_eq_int(mu1, mu2, lr)

    implicit none

    ! Input: Chemical potentials and left/right indicator
    real(dp) :: mu1, mu2 ! mu1, mu2: Chemical potentials.
    integer :: lr ! lr: Integer representing left/right electrode.
    real(dp) :: eoff_l, eoff_u, mu_l, mu_u ! eoff_l, eoff_u: Lower and upper energy offsets. mu_l, mu_u: Lower and upper chemical potentials.

    ! Initialize equilibrium integration parameters
    mu_l = min(mu1, mu2)
    mu_u = max(mu1, mu2)
    intlr = lr
    eoff_l = log(epsfermi/(1d0 - epsfermi))*kb*temperature_el + mu_l
    eoff_u = log(1d0/epsfermi - 1d0)*kb*temperature_el + mu_u
    r_eq = abs(elow - eoff_l)*0.5d0
    e0 = eoff_l - r_eq
    phi_low = pi - asin(delta_imag/r_eq)
    x_low = contour(phi_low)
    dx_low = real(eoff_u - x_low)
    n_fermipols = (delta_imag/(pi*kb*temperature_el) - 1)*0.5d0

    lnoneq_int = .false.

  end subroutine init_eq_int

  subroutine init_neq_int(lr)
#include <petsc/finclude/petscmat.h>
    use globals, only: l_ionode
    implicit none
    ! Input: Left/right indicator
    integer :: lr ! lr: Integer representing left/right electrode.
    complex(dp) :: z1, z2 ! Temporary complex variables.

    ! Initialize non-equilibrium integration parameters
    intlr = lr
    el = log(epsfermi/(1d0 - epsfermi))*kb*temperature_el + min(mul, mur)
    eu = log(1d0/epsfermi - 1d0)*kb*temperature_el + max(mul, mur)
    z1 = el - mul
    z2 = el - mur
    write (pstr_out, fmt='(A,4e24.12)') "el, eu", el, eu; call petsc_print_master()
    write (pstr_out, fmt='(A,2e24.12)') "low :", fermi(z1, temperature_el) - fermi(z2, temperature_el); call petsc_print_master()
    z1 = eu - mul
    z2 = eu - mur
    write (pstr_out, fmt='(A,2e24.12)') "up :", fermi(z1, temperature_el) - fermi(z2, temperature_el); call petsc_print_master()

    lnoneq_int = .true.

  end subroutine init_neq_int

  subroutine scale_Dne_dc(p_Dneq_dc, dc_weights)
#include <petsc/finclude/petscmat.h>
    use globals
      
    implicit none
    
    ! Input: PETSc matrix and weights
    Mat :: p_Dneq_dc ! p_Dneq_dc: PETSc matrix.
    PetscScalar :: dc_weights(:) ! dc_weights: Array of PETSc scalars representing weights.
    
    ! Local variables
    integer :: ierr, imu1, imu2,  nl1, nl2, i1, nlc1, nlc2, &
      nrow, ncol, nzrow, ii(1), nz ! Loop indices, matrix dimensions, number of non-zero elements, and error code.
    integer, allocatable :: cols(:) ! cols: Array to store column indices.
    PetscScalar, allocatable :: row_vals(:) ! row_vals: Array to store row values.
    PetscScalar :: p_val ! Temporary variable to store matrix values.
    PetscReal :: weight ! Temporary variable to store weight.
      
    ! Get matrix dimensions
    call MatGetSize(p_Dneq_dc, nrow, ncol, ierr)
    call MatGetOwnershipRange(p_Dneq_dc, nl1, nl2, ierr)
    
    ! Allocate arrays
    allocate(cols(ncol), row_vals(ncol), stat = ierr)
    
    ! Loop over rows
    do imu1 = 0, nmu_c - 1
      if ((imu1 .ge. nl1) .and. (imu1 .le. nl2 - 1)) then
        ! Get row
        call MatGetRow(p_Dneq_dc, imu1, nzrow, cols, row_vals, ierr)
        ! Loop over non-zero elements in the row
        do i1 = 1, nzrow
          imu2 = cols(i1)
!~           write(0, fmt='(2i8, 2e24.12)') imu1, imu2, dc_weights(imu2 + 1, imu1 + 1)
!~           row_vals(i1) = row_vals(i1) * dc_weights(imu2 + 1, (imu1 + 1))
          ! Scale row values by weights
          if (((real(dc_weights(imu1), 8).lt.0d0).or.(real(dc_weights(imu2), 8).lt.0d0)).and.&
         &  ((abs((aimag(dc_weights(imu1)))).ge.1d-10).or.&
         &   (abs((aimag(dc_weights(imu2)))).ge.1d-10))) then
            write(0, fmt='(A,2i8,4e24.12)') "warning: ", imu1, imu2, dc_weights(imu1), dc_weights(imu2)
          end if
          weight = dsqrt(real(dc_weights(imu1), 8) * real(dc_weights(imu2), 8))
          if (weight.lt.0d0) weight = 0d0
          row_vals(i1) = row_vals(i1) * weight
        end do
        ii(1) = imu1
        nz = nzrow
        ! Restore row and set values
        call MatRestoreRow(p_Dneq_dc, imu1, nz, PETSC_NULL_INTEGER, PETSC_NULL_SCALAR, ierr)        
        call MatSetValues(p_Dneq_dc, 1, ii, nzrow, cols, row_vals, INSERT_VALUES, ierr)
      end if
      call MatAssemblyBegin(p_Dneq_dc, MAT_FINAL_ASSEMBLY, ierr)
      call MatAssemblyEnd(p_Dneq_dc, MAT_FINAL_ASSEMBLY, ierr)
    end do
    
  end subroutine scale_Dne_dc

  subroutine blockscale_mat(a, w, nat, imuecc)
    implicit none

    ! Input/Output: Complex matrix a, input complex matrix w, number of atoms, and number of modes per atom
    complex(dp), intent(inout) :: a(:, :) ! a: Matrix to be scaled.
    complex(dp), intent(in) :: w(:, :) ! w: Weight matrix.
    integer, intent(in) :: nat, imuecc(:) ! nat: Number of atoms. imuecc: Number of modes per atom.
    integer :: imu1, imu2, iat1, iat2, i1, i2, j1, j2 ! Loop indices.

    ! Scale matrix a by blocks using matrix w
    i2 = 1
    j2 = 1
    do iat2 = 1, nat
      j2 = i2 + imuecc(iat2) - 1
      i1 = 1
      j1 = 1
      do iat1 = 1, nat
        j1 = i1 + imuecc(iat1) - 1
        a(i1:j1, i2:j2) = w(iat1, iat2)*a(i1:j1, i2:j2)
        i1 = j1 + 1
      end do
      i2 = j2 + 1
    end do

  end subroutine blockscale_mat

  subroutine hpsort_aa(ra, rb, reverse)

    implicit none

! input,output:
    real(dp) :: ra(:) ! ra: Array to be sorted.
    real(dp) :: rb(:, :) ! rb: Matrix whose rows are sorted along with ra.
! local :
    real(dp), allocatable :: rrb(:) ! rrb: Temporary array.
    real(dp) :: rra ! Temporary variable.
    integer :: l, ir, n, j, i, m, minus ! Loop indices and other variables.
    logical, optional :: reverse ! Optional logical variable to specify reverse sorting.

    ! Heap sort algorithm for ra and rb
    minus = 1
    if (present(reverse)) then
      if (reverse) minus = -1
    end if

    n = size(ra, 1)
    m = size(rb, 1)
    allocate (rrb(m))
    l = n/2 + 1
    ir = n
    ra = ra*minus
    do
      if (l > 1) then
        l = l - 1
        rra = ra(l)
        rrb(:) = rb(:, l)
      else
        rra = ra(ir)
        rrb(:) = rb(:, ir)
        ra(ir) = ra(1)
        rb(:, ir) = rb(:, 1)
        ir = ir - 1
        if (ir .eq. 1) then
          ra(1) = rra
          rb(:, 1) = rrb(:)
          ra = ra*minus
          return
        end if
      end if
      i = l
      j = l + l
      do while (j .le. ir)
        if (j .lt. ir) then
          if (ra(j) .lt. ra(j + 1)) j = j + 1
        end if
        if (rra < ra(j)) then
          ra(i) = ra(j)
          rb(:, i) = rb(:, j)
          i = j
          j = j + j
        else
          j = ir + 1
        end if

      end do
      ra(i) = rra
      rb(:, i) = rrb(:)
    end do
    ra = ra*minus
  end subroutine hpsort_aa

  subroutine get_quadrule(integrator, xi, v1, v2, nint_order, nlow)
    use kinds
    use integrations_weights
    use error_handler

    implicit none

    ! Input: Integrator type and order, Output: quadrature rule
    integer :: integrator ! integrator: Integer specifying the type of integrator.
    integer :: nint_order, nlow ! nint_order: Order of integration. nlow: Number of lower order points.
    real(dp), allocatable :: xi(:), v1(:), v2(:) ! xi: Quadrature points. v1, v2: Weights.

    ! Local variables
    integer :: off, i, ierr ! Loop index and error code.
    real(dp), allocatable :: x2(:), w1(:), w2(:), vdummy(:, :) ! Temporary arrays.

    ! Select quadrature rule based on integrator type
    if (integrator .eq. 1) then
      ! Gauss-Kronrod quadrature
      if (mod(nint_order, 2) .eq. 0) nint_order = nint_order + 1
      nlow = (nint_order - 1)/2
      if (l_output_progress) then
        write (pstr_out, fmt='(A,2i8)') " Adaptive Gauss Kronrod Integrator ", nlow, nint_order; call petsc_print_master()
      end if
      allocate (xi(nint_order), v1(nint_order), v2(nint_order), &
                x2(nlow + 1), w1(nlow + 1), w2(nlow + 1), vdummy(1, nint_order), stat=ierr)
      if (ierr .ne. 0) then
        write (errormsg, fmt='(A,i8)') "allocation error ", ierr
        call error()
      end if
      call kronrod(nlow, 1d-12, x2, w1, w2)
      off = 0

      do i = 1, nlow
        xi(i) = x2(i)
        v1(i) = w1(i)
        v2(i) = w2(i)
        xi(nint_order - i + 1) = -x2(i)
        v1(nint_order - i + 1) = w1(i)
        v2(nint_order - i + 1) = w2(i)
      end do
      xi(nlow + 1) = x2(nlow + 1)
      v1(nlow + 1) = w1(nlow + 1)
      v2(nlow + 1) = w2(nlow + 1)
    else if (integrator .eq. 2) then
      ! Clenshaw-Curtis quadrature
      if ((mod(nint_order, 2)) .eq. 0) nint_order = nint_order + 1
      nlow = (nint_order - 1)/2 + 1
      if (l_output_progress) then
        write (pstr_out, fmt='(2i8,A)') nlow, nint_order, " Adaptive Clenshaw-Curtis integrator"; call petsc_print_master()
      end if
      allocate (xi(nint_order), v1(nint_order), v2(nint_order), &
                x2(nint_order), w1(nint_order), w2(nint_order), vdummy(1, nint_order), stat=ierr)
      x2 = 0d0
      w1 = 0d0
      xi = 0d0
      v1 = 0d0
      v2 = 0d0
      call clenshaw_curtis_compute(nlow, x2(1:nlow), w2(1:nlow))
      call clenshaw_curtis_compute(nint_order, x2, w1)
      do i = 1, nlow
        v1(i) = w1(2*(i - 1) + 1)
        xi(i) = x2(2*(i - 1) + 1)
        v2(i) = w2(i)
        if (i .eq. nlow) exit
        v1(nlow + i) = w1(2*i)
        xi(nlow + i) = x2(2*i)
      end do

    else if (integrator .eq. 3) then
      ! Newton-Cotes quadrature
      if (mod(nint_order, 2) .eq. 0) nint_order = nint_order + 1
      nlow = (nint_order)/2 + 1
      if (l_output_progress) then
        write (pstr_out, fmt='(2i8,A)') nlow, nint_order, " Adaptive Newton-Cotes integrator"; call petsc_print_master()
      end if
      allocate (xi(nint_order), v1(nint_order), v2(nint_order), &
                x2(nint_order), w1(nint_order), w2(nint_order), vdummy(1, nint_order), stat=ierr)
      if (ierr .ne. 0) then
        write (errormsg, fmt='(A,i8)') "allocation error ", ierr
        call error()
      end if
      x2 = 0d0
      w2 = 0d0
      call line_ncc_rule(nlow, -1d0, 1d0, x2(1:nlow), w2(1:nlow))
      call line_ncc_rule(nint_order, -1d0, 1d0, x2, w1)
      do i = 1, nlow
        v1(i) = w1(2*(i - 1) + 1)
        xi(i) = x2(2*(i - 1) + 1)
        v2(i) = w2(i)
        if (i .eq. nlow) exit
        v1(nlow + i) = w1(2*i)
        xi(nlow + i) = x2(2*i)
      end do

    end if

    deallocate (x2, w1, w2)
    allocate (x2(nint_order))
    x2 = xi
    vdummy(1, :) = v1(:)
    call hpsort_aa(xi, vdummy)
    v1(:) = vdummy(1, :)
    vdummy(1, :) = v2(:)
    call hpsort_aa(x2, vdummy)
    v2(:) = vdummy(1, :)

!~  if (inode.eq.0) then
!~         do i=1,nint_order
!~           write(0,fmt='(i8,3e24.12)') i,xi(i),v1(i),v2(i)
!~         end do
!~  end if

    if (l_output_progress) then
      write (pstr_out, fmt='(A)') repeat("-", 4); call petsc_print_master()
    end if

  end subroutine get_quadrule

  function fermi(e, t)

    use kinds

    implicit none

    complex(dp) :: ii = dcmplx(0d0, 1d0) ! Imaginary unit
    complex(dp) :: fermi, e, x, ze ! Fermi function and temporary variables
    real(dp) :: t, kt, infinity ! Temperature, kt, and infinity

    ! Fermi-Dirac distribution function
    kt = t*kb
    x = e/kt

    ze = zexp(x)
    infinity = HUGE(0_dp)
    if (real(ze, 8) .gt. infinity) then
      fermi = 1d-32
    else
      fermi = 1d0/(1d0 + zexp(x))
    end if

  end function fermi

  function bose(e, t)

    use kinds

    implicit none

    complex(dp) :: ii = dcmplx(0d0, 1d0) ! Imaginary unit
    complex(dp) :: bose, e, x, ze, b ! Bose function and temporary variables
    real(dp) :: t, kt, infinity ! Temperature, kt, and infinity

    ! Bose-Einstein distribution function
    kt = t*kb
    x = e/kt

    ze = zexp(x)
    b = 1d0/(zexp(x) - 1d0)
    if (isnan(real(b))) b = 0d0
    if (isnan(aimag(b))) b = 0d0
    bose = b

  end function bose

  function contour(x)

    use kinds

    implicit none

    ! Input: x value, Output: complex contour value
    real(dp) :: x ! x: Real variable.
    complex(dp) :: contour ! contour: Complex variable representing the contour.

    ! Contour function
    contour = cmplx(-cos(x), sin(x), 8)
    contour = r_eq*contour + e0

  end function contour

  function dcontour(x)

    use kinds

    implicit none

    ! Input: x value, Output: complex derivative of contour
    real(dp) :: x ! x: Real variable.
    complex(dp) :: dcontour ! dcontour: Complex variable representing the derivative of the contour.

    ! Derivative of the contour function
    dcontour = r_eq*cmplx(-sin(x), -cos(x), 8)

  end function dcontour

  function contour_x(x)

    use kinds

    implicit none

    ! Input: x value, Output: complex contour value
    real(dp) :: x ! x: Real variable.
    complex(dp) :: contour_x ! contour_x: Complex variable representing the contour.

    ! Contour function (alternative form)
    contour_x = x_low + dx_low*(1d0 - x)

  end function contour_x

  function dcontour_x(x)

    use kinds

    implicit none

    ! Input: x value, Output: complex derivative of contour
    real(dp) :: x ! x: Real variable.
    complex(dp) :: dcontour_x ! dcontour_x: Complex variable representing the derivative of the contour.

    ! Derivative of the contour function (alternative form)
    dcontour_x = -dx_low

  end function dcontour_x

  function contour_y(x)

    use kinds

    implicit none

    ! Input: x value, Output: complex contour value
    real(dp) :: x ! x: Real variable.
    complex(dp) :: contour_y ! contour_y: Complex variable representing the contour.

    ! Contour function (alternative form)
    contour_y = y_low + dy_low*x

  end function contour_y

  function dcontour_y(x)

    use kinds

    implicit none

    ! Input: x value, Output: complex derivative of contour
    real(dp) :: x ! x: Real variable.
    complex(dp) :: dcontour_y ! dcontour_y: Complex variable representing the derivative of the contour.

    ! Derivative of the contour function (alternative form)
    dcontour_y = dy_low

  end function dcontour_y

  subroutine get_fermipole_contribution(p_dmat, mu)
#include <petsc/finclude/petsc.h>
    use petscmat
    use petsc_mod
    use petsc_wrapper
    implicit none

    ! Input: PETSc matrix and chemical potential, Output: updated matrix
    Mat :: p_dmat ! p_dmat: PETSc matrix.
    real(dp) :: mu ! mu: Chemical potential.

    complex(dp) :: z, fac ! z: Complex variable. fac: Scaling factor.
    integer :: ifp, ierr ! Loop index and error code.

    ! Calculate Fermi pole contributions and add to the matrix
    call MatZeroEntries(p_dmat, ierr)

    do ifp = 0, n_fermipols
      z = cmplx(mu, (2d0*ifp + 1)*pi*kb*temperature_el)
      write (pstr_out, fmt='(i8,2e24.12)') ifp, z; call petsc_print_master()
      ierr = eval_gr(z)
      fac = zione*2d0*pi*kb*temperature_el
      call MatAXPY(p_dmat, fac, p_tmpcc1, SAME_NONZERO_PATTERN, ierr)
!~         dmat=dmat+zione*2d0*pi*kb*temp*tmpcc
    end do

  end subroutine get_fermipole_contribution

  subroutine evalqr(f, a, b, x, w1, w2, p_zint, p_zint2, nhigh, nlow, errout, errout2, d1, d2)
#include <petsc/finclude/petsc.h>
    use petscmat
    use petsc_mod
    use petsc_wrapper
    use kinds
    use mathstuff
    use globals, only: ldouble_contour, lnoneq_int

    implicit none

    ! Input: Function, integration limits, quadrature points and weights, matrixes, number of high and low order points, Output: results and error
    real(dp) ::a, b, errout, errout2 ! a, b: Integration limits. errout, errout2: Errors.
    real(dp), optional :: d1, d2 ! d1, d2: Optional subdivision points.
    real(dp), allocatable :: x(:), w1(:), w2(:) ! x: Quadrature points. w1, w2: Weights.
    Mat :: p_zint, p_zint2 ! p_zint, p_zint2: PETSc matrices to store integration results.
    integer :: nhigh, nlow, ierr ! nhigh, nlow: Number of high and low order points. ierr: Error code.
    integer, external :: f ! f: External function to be integrated.

    ! Local variables
    Mat :: p_zintg, p_zintg2 ! p_zintg, p_zintg2: Temporary PETSc matrices.
    integer :: i, n, j, i_low, i_high, i_max, k1, k2 ! Loop indices.
    real(dp) :: xx, norm, m1, m2, timing_local ! Temporary variables.
    integer(8) :: counti, count_rate, countf ! Variables for timing.
    PetscScalar, allocatable :: ws1(:), ws2(:) ! Temporary arrays.

    ! Perform quadrature integration
    allocate (ws1(nhigh), ws2(nhigh))

    ws1 = w1*(b - a)*0.5d0
    ws2 = w2*(b - a)*0.5d0

    call MatZeroEntries(p_zint, ierr)

    call MatDuplicate(p_zint, MAT_SHARE_NONZERO_PATTERN, p_zintg, ierr)

    if ((ldouble_contour) .and. (lnoneq_int)) then
      call MatDuplicate(p_zint2, MAT_SHARE_NONZERO_PATTERN, p_zintg2, ierr)
      call MatZeroEntries(p_zint2, ierr)
    end if
    timing_local = 0d0
    if (l_output_progress) then
      write (pstr_out, fmt='(i3,e16.6)') 0, timing_local; call PetscPrintf(PETSC_COMM_WORLD, trim(pstr_out), ierr)
    end if
    do i = 1, nhigh
      call system_clock(counti, count_rate)
      xx = x(i)*(b - a)*0.5d0 + (a + b)*0.5d0

      ierr = f(xx)

      call MatAXPY(p_zint, ws1(i), p_tmpcc1, SAME_NONZERO_PATTERN, ierr)

      if (ws2(i) .ne. 0d0) then
        call MatAXPY(p_zintg, ws2(i), p_tmpcc1, SAME_NONZERO_PATTERN, ierr)
      end if

      if ((ldouble_contour) .and. (lnoneq_int)) then
        call MatAXPY(p_zint2, ws1(i), p_tmpcc2, SAME_NONZERO_PATTERN, ierr)
        if (ws2(i) .ne. 0d0) call MatAXPY(p_zintg2, ws2(i), p_tmpcc2, SAME_NONZERO_PATTERN, ierr)
      end if

      int_counter = int_counter + 1
      call system_clock(countf)
      timing_local = timing_local + real(countf - counti, 8)/real(count_rate, 8)
      if (l_output_progress) then
        write (pstr_out, fmt='(A,i3,e16.6)') repeat(ACHAR(8), 19), i, timing_local; call PetscPrintf(PETSC_COMM_WORLD, trim(pstr_out), ierr)
      end if
    end do
    call PetscPrintf(PETSC_COMM_WORLD, " int eps", ierr)

    call MatAXPY(p_zintg, p_minus1, p_zint, SAME_NONZERO_PATTERN, ierr)
    call MatNorm(p_zintg, NORM_FROBENIUS, errout, ierr)
    if ((ldouble_contour) .and. (lnoneq_int)) then
      call MatAXPY(p_zintg2, p_minus1, p_zint2, SAME_NONZERO_PATTERN, ierr)
      call MatNorm(p_zintg2, NORM_FROBENIUS, errout2, ierr)
!~       write(0,*) errout, errout2
!~       errout = (errout + errout2)*0.5d0
      errout = max(errout, errout2)
      errout2 = errout
    else
      errout2 = errout
    end if

    call MatDestroy(p_zintg, ierr)

    if ((ldouble_contour) .and. (lnoneq_int)) then
      call MatDestroy(p_zintg2, ierr)
    end if

    if (.not. present(d2)) return

    d1 = a + (b - a)/3d0
    d2 = a + (b - a)*2d0/3d0

  end subroutine evalqr

  function get_gr_cc_neq(x)
#include <petsc/finclude/petscmat.h>
    use kinds
    use petsc_mod
    use petsc_wrapper
    use globals, only: eta_elec, iik, vb, nmu_l, nmu_r, nmu_c, wkp_r, eta_cc

    implicit none

    real(dp) :: x

    integer:: get_gr_cc_neq

    Mat :: p_gammal_block, p_gammar_block, p_x, b, p_g_block, p_rhs
    Mat, pointer :: p_gr_inv

    complex(dp) :: z

    integer :: ierr
    real(dp) :: vbb
    logical :: l_solvemode_5
    
    PetscScalar :: p_scal
    PetscReal :: norm

    get_gr_cc_neq = 0

    z = x

    p_gr_inv => p_invGr ! see eval_gr
    
    call init_gr(z, p_gr_inv)

    l_solvemode_5 = solver_mode .eq. 5
    
      if (l_solvemode_5) solver_mode = 1
      call petsc_get_a_with_b_c(p_g_block, p_gr_inv, p_gammal, mattype_dense)    
      if (solver_mode .eq. 2) then
        call petsc_get_a_with_b_c(p_rhs, p_gr_inv, p_gammal, mattype_dense)      
        call petsc_one_mat(p_rhs, 0, nmu_l - 1) ! prepare RHS nmu_c,nmu_l (1 0 0 .. 0,0 1 0 .. 0, ..) C style 1=0
        call petsc_call_solver(p_gr_inv, p_rhs, p_g_block, matsolvertype_cc, solver_mode)
      else if ((solver_mode .eq. 1).or.(solver_mode .eq. 3).or.(solver_mode .eq. 4)) then
        call petsc_get_a_with_b_c(p_rhs, p_gr_inv, p_gammal, mattype_sparse)
        call petsc_one_mat(p_rhs, 0, nmu_l - 1,2) ! prepare RHS nmu_c,nmu_l (1 0 0 .. 0,0 1 0 .. 0, ..) C style 1=0
        call petsc_call_solver(p_gr_inv, p_rhs, p_g_block, matsolvertype_cc, solver_mode)
      end if
      if (l_solvemode_5) solver_mode = 5


      call MatDestroy(p_rhs, ierr)

      call MatMatMult(p_g_block, p_gammal, MAT_INITIAL_MATRIX, 1d0, p_x, ierr) ! X_block=(G11*Gamma_L,G12*Gamma_L,G13*Gamma_L)
      call petsc_matmatconj_restricted(p_x, p_g_block, p_tmpcc1)      

!~ ! debug start --
!~         call MatNorm(p_gammal,NORM_FROBENIUS,norm,ierr)
!~         if (inode.eq.0)  write(0,fmt='(2e16.8,A,e16.8)') z, " Gammal ",norm
!~         call MatNorm(p_tmpcc2,NORM_FROBENIUS,norm,ierr)
!~         if (inode.eq.0)  write(0,fmt='(2e16.8,A,e16.8)') z, " G*Gammal*G' ",norm
!~ ! debug end --

      p_scal = 0.5d0*(fermi(z - mul, temperature_el) - fermi(z - mur, temperature_el)) ! add weight functions, fermi window and absorb factor 0.5d0
      call MatScale(p_tmpcc1, p_scal, ierr)
      
      call MatDestroy(p_x, ierr)
      call MatDestroy(p_rhs, ierr)
      call MatDestroy(p_g_block, ierr)    
    
     

    if (intlr .eq. 2) then  
    
      if (l_solvemode_5) solver_mode = 1
      call petsc_get_a_with_b_c(p_g_block, p_gr_inv, p_gammar, mattype_dense)    
      if (solver_mode .eq. 2) then
        call petsc_get_a_with_b_c(p_rhs, p_gr_inv, p_gammar, mattype_dense)      
        call petsc_one_mat(p_rhs, nmu_c - nmu_r + 1 - 1, nmu_c - 1) ! prepare RHS nmu_c,nmu_r ( .., 0 0 .. 1 0, 0 0 .. 0 1 )   C style 1=0
        call petsc_call_solver(p_gr_inv, p_rhs, p_g_block, matsolvertype_cc, solver_mode)
      else if ((solver_mode .eq. 1).or.(solver_mode .eq. 3).or.(solver_mode .eq. 4)) then
        call petsc_get_a_with_b_c(p_rhs, p_gr_inv, p_gammar, mattype_sparse)
        call petsc_one_mat(p_rhs, nmu_c - nmu_r + 1 - 1, nmu_c - 1, 2) ! prepare RHS nmu_c,nmu_r ( .., 0 0 .. 1 0, 0 0 .. 0 1 )   C style 1=0
        call petsc_call_solver(p_gr_inv, p_rhs, p_g_block, matsolvertype_cc, solver_mode)
      end if
      if (l_solvemode_5) solver_mode = 5
  
      call MatDestroy(p_rhs, ierr)
  
      call MatMatMult(p_g_block, p_gammar, MAT_INITIAL_MATRIX, 1d0, p_x, ierr) ! X_block=(G13*Gamma_R,G23*Gamma_R,G33*Gamma_R)    
      call petsc_matmatconj_restricted(p_x, p_g_block, p_tmpcc2)
!~      ! debug start --
!~        call MatNorm(p_gammar,NORM_FROBENIUS,norm,ierr)
!~        if (inode.eq.0)  write(0,fmt='(2e16.8,A,e16.8)') z," Gammar ",norm
!~        call MatNorm(p_tmpcc1,NORM_FROBENIUS,norm,ierr)
!~        if (inode.eq.0)  write(0,fmt='(2e16.8,A,e16.8)') z," G*Gammar*G' ",norm
!~  ! debug end --
  
      p_scal = 0.5d0*(fermi(z - mul, temperature_el) - fermi(z - mur, temperature_el)) ! add weight functions, fermi window and absorb factor 0.5d0
      call MatScale(p_tmpcc2, p_scal, ierr)
  
      call MatDestroy(p_x, ierr)
      call MatDestroy(p_rhs, ierr)
      call MatDestroy(p_g_block, ierr)    

    end if

!~     call MatDestroy(p_gr_inv, ierr)

  end function get_gr_cc_neq

  function get_gr_cc_eq(x)

    use kinds
    use petsc_mod
    use petsc_wrapper
    use globals, only: eta_elec, iik, vb, nmu_l, nmu_r, nmu_c, wkp_r, eta_cc, ngroups
    use error_handler
    use pexsi_wrapper

    implicit none

    real(dp) :: x

    integer:: get_gr_cc_eq

    Mat :: p_mat_one
    Mat, pointer :: p_gr_inv

    complex(dp) :: z, dc
    integer :: ierr
    PetscScalar :: p_scal
    PetscReal :: norm

    get_gr_cc_eq = 0

    if (contour_select .eq. 1) then
      z = contour(x)
      dc = dcontour(x)
    else if (contour_select .eq. 2) then
      z = contour_x(x)
      dc = dcontour_x(x)
    else
      write (errormsg, fmt='(A,i8)') "contour_select invalid", contour_select
      call error()
    end if

    p_gr_inv => p_invGr ! see eval_gr

    call init_gr(z, p_gr_inv)

    if (solver_mode .eq. 2) then
      call petsc_get_densemat(p_gr_inv, p_mat_one, mattype_dense)
      call petsc_one_mat(p_mat_one, 0, nmu_c - 1)
      call petsc_call_solver(p_gr_inv, p_mat_one, p_tmpcc1, matsolvertype_cc, solver_mode)
      call MatDestroy(p_mat_one, ierr)
    else if ((solver_mode .eq. 1).or.(solver_mode .eq. 3).or.(solver_mode .eq. 4)) then
      call petsc_get_a_with_b_c(p_mat_one, p_gr_inv, p_gr_inv, mattype_sparse)
      call petsc_one_mat(p_mat_one, 0, nmu_c - 1, 2)
      call petsc_call_solver(p_gr_inv, p_mat_one, p_tmpcc1, matsolvertype_cc, solver_mode)
      call MatDestroy(p_mat_one, ierr)
    else if (solver_mode .eq. 5) then
      call inv_sel(p_gr_inv, p_tmpcc1)
    end if


    if (intlr .eq. 1) then
      p_scal = dc*fermi(z - mul, temperature_el)
    else if (intlr .eq. 2) then
      p_scal = dc*fermi(z - mur, temperature_el)
    else if (intlr .eq. 3) then
      p_scal = dc*(fermi(z - mul, temperature_el) + fermi(z - mur, temperature_el))*0.5d0
    else if (intlr .eq. 4) then
      p_scal = dc*(fermi(z - mul, temperature_el) - fermi(z - mur, temperature_el))
    else
      write (errormsg, fmt='(A,i8)') "only left=1 and right=2 electrodes", intlr
      call error()
    end if

    call MatScale(p_tmpcc1, p_scal, ierr)

!~     call MatDestroy(p_gr_inv, ierr)

  end function get_gr_cc_eq

!~     subroutine init_gr(z,p_tmp1)
!~ initialises G^r ,i.e. p_tmp1=E*S-H-Sigma_L-Sigma_R
  subroutine init_gr(z, p_tmp1)
#include <petsc/finclude/petscmat.h>
    use petsc_mod
    use petsc_wrapper
    use globals, only: eta_elec, iik
    use kinds
    use error_handler

    implicit none

    complex(dp) :: z
    Mat :: p_tmp1

    integer :: ierr
    real(dp) :: vbb
    complex(dp) :: zz
    PetscReal :: norm

    ierr = eval_sigma(z)
    call MatZeroEntries(p_tmp1, ierr)
    call MatAXPY(p_tmp1, p_minus1, p_h00k_cc(iik), DIFFERENT_NONZERO_PATTERN, ierr)
    zz = (z + eta_cc*zione)
    call MatAXPY(p_tmp1, zz, p_s00k_cc(iik), DIFFERENT_NONZERO_PATTERN, ierr)

    call petsc_add_sub_B_to_A(p_sigmalr, p_tmp1, 0, 0, p_minus1,&
   &ADD_VALUES, PETSC_FALSE)
    call petsc_add_sub_B_to_A(p_sigmarr, p_tmp1, nmu_c - nmu_r,&
   &nmu_c - nmu_r, p_minus1, ADD_VALUES, PETSC_FALSE)

!bias shift
     if (((vl.ne.0d0).or.(vr.ne.0d0))) then
       zz = -vr
!~        call petsc_add_sub_B_to_A(p_s10k_r(iik), p_tmp1, nmu_c - nmu_r, nmu_c - nmu_r*2, zz, ADD_VALUES, PETSC_FALSE)       
!~        call petsc_add_sub_B_to_A(p_s01k_r(iik), p_tmp1, nmu_c - nmu_r*2, nmu_c - nmu_r, zz, ADD_VALUES, PETSC_FALSE)      
       call petsc_add_sub_B_to_A(p_s00k_r(iik), p_tmp1, nmu_c - nmu_r, nmu_c - nmu_r, zz, ADD_VALUES, PETSC_FALSE)            
!~        call petsc_add_sub_B_to_A(p_s00k_r(iik), p_tmp1, nmu_c - nmu_r*2, nmu_c - nmu_r*2, zz, ADD_VALUES, PETSC_FALSE)
       
     
       zz = -vl ! this shall be change to genral bias voltage in left and right if vl and vr       
       call petsc_add_sub_B_to_A(p_s00k_l(iik), p_tmp1, 0, 0, zz, ADD_VALUES, PETSC_FALSE)
!~        call petsc_add_sub_B_to_A(p_s00k_l(iik), p_tmp1, nmu_l, nmu_l, zz, ADD_VALUES, PETSC_FALSE)
!~        call petsc_add_sub_B_to_A(p_s10k_l(iik), p_tmp1, nmu_l, 0, zz, ADD_VALUES, PETSC_FALSE)
!~        call petsc_add_sub_B_to_A(p_s01k_l(iik), p_tmp1, 0, nmu_l, zz, ADD_VALUES, PETSC_FALSE)

     end if

  end subroutine init_gr

  function eval_gr(z)
#include <petsc/finclude/petscmat.h>
    use petsc_mod
    use petsc_wrapper
    use globals, only: nsim_inv
    use kinds
    use pexsi_wrapper

    implicit none

    complex(dp) :: z

    integer :: eval_gr

    Mat :: p_tmp1, p_tmp2

    Mat, pointer :: p_gr_inv
    Mat :: p_mat_one
    integer :: ierr

    eval_gr = 0

!~ Bug in PETSC 3.16 causes preallocator matrix to be messed up after first call.
!~ So we initialize p_invGr in init_sys.f90 and reuse it in each call to init_gr,
!~ where it's entries are set to 0 at the beginning. To keep the structure 
!~ of the code and in case we want to go back use pointer to p_invG and leave most
!~ of the code unchanged, save the call to preallocated in init_gr.
    p_gr_inv => p_invGr
    
    call init_gr(z, p_gr_inv)

    if (solver_mode .eq. 2) then
      call petsc_get_densemat(p_gr_inv, p_mat_one, mattype_dense)
      call petsc_one_mat(p_mat_one, 0, nmu_c - 1)
      call petsc_call_solver(p_gr_inv, p_mat_one, p_tmpcc1, matsolvertype_cc, solver_mode)
      call MatDestroy(p_mat_one, ierr)
    else if ((solver_mode .eq. 1).or.(solver_mode .eq. 3).or.(solver_mode .eq. 4)) then
      call petsc_get_a_with_b_c(p_mat_one, p_gr_inv, p_gr_inv, mattype_sparse)
      call petsc_one_mat(p_mat_one, 0, nmu_c - 1, 2)
      call petsc_call_solver(p_gr_inv, p_mat_one, p_tmpcc1, matsolvertype_cc, solver_mode)
      call MatDestroy(p_mat_one, ierr)
    else if (solver_mode .eq. 5) then
      call inv_sel(p_gr_inv, p_tmpcc1)
    end if

!~     call MatDestroy(p_gr_inv, ierr)

  end function eval_gr

  function eval_sigma(z)
#include <petsc/finclude/petscmat.h>
    use petsc_mod
    use petsc_wrapper
    use kinds
!~       use blas95
    use globals, only: eta_elec, iik

    implicit none

    complex(dp) :: z

    integer:: eval_sigma

    complex(dp) :: zz
    integer :: ierr
    PetscReal :: norml, normr, normgrr, normgll
    PetscScalar :: pl11,plnn,pr11,prnn

    Mat :: p_tmp1, p_tmp2, p_tmp3

    eval_sigma = 0

! sigma_l

!~       call petsc_get_densemat(p_h00k_l(iik),p_tmp1,mattype_surf)
!~       call petsc_get_densemat(p_h00k_l(iik),p_tmp2,mattype_surf)
!~       call petsc_get_densemat(p_h00k_l(iik),p_tmp3,mattype_surf)

    call surf_gf3_petsc(z - vl + cmplx(0d0, eta_elec), iik, "l")
!~     call MatNorm(p_gllr, NORM_FROBENIUS, normgll, ierr)
    
    ! pull out -1 factor to avoid calling MatScale twice !!! (z*S_CL-H_CL)*g_LL*(z*S_LC-H_LC)
    zz = -(z - vl + cmplx(0d0, eta_elec))

    call MatDuplicate(p_s10k_l(iik), MAT_COPY_VALUES, p_tmp1, ierr)
    call MatScale(p_tmp1, zz, ierr)
    call MatAXPY(p_tmp1, p_one, p_h10k_l(iik), DIFFERENT_NONZERO_PATTERN, ierr)

    call MatDuplicate(p_s01k_l(iik), MAT_COPY_VALUES, p_tmp2, ierr)
    call MatScale(p_tmp2, zz, ierr)
    call MatAXPY(p_tmp2, p_one, p_h01k_l(iik), DIFFERENT_NONZERO_PATTERN, ierr)

    call MatConvert(p_tmp1, mattype_surf, MAT_INPLACE_MATRIX, p_tmp1, ierr)
    call MatConvert(p_tmp2, mattype_surf, MAT_INPLACE_MATRIX, p_tmp2, ierr)
    call MatMatMult(p_tmp1, p_gllr, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, p_tmp3, ierr)
    call MatDestroy(p_tmp1, ierr)
    call MatMatMult(p_tmp3, p_tmp2, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, p_tmp1, ierr) ! there is no ELEMENTAL*ELEMENTAL=DENSE so we need MatConvert
    call MatConvert(p_tmp1, mattype_dense, MAT_REUSE_MATRIX, p_sigmalr, ierr)

! gammaL=i*(sigmaLr-sigmaLa)
!~     call MatNorm(p_sigmalr, NORM_FROBENIUS, norml, ierr)
    call MatTranspose(p_sigmalr, MAT_REUSE_MATRIX, p_gammal, ierr)
    call MatConjugate(p_gammal, ierr)
    call MatAYPX(p_gammal, p_minus1, p_sigmalr, SAME_NONZERO_PATTERN, ierr) !p_gammal=(sigmaLr-sigmaLa)
    call MatScale(p_gammal, p_zione, ierr)

    call MatDestroy(p_tmp1, ierr)
    call MatDestroy(p_tmp2, ierr)
    call MatDestroy(p_tmp3, ierr)

! sigma_r

!~       call petsc_get_densemat(p_h00k_r(iik),p_tmp1,mattype_surf)
!~       call petsc_get_densemat(p_h00k_r(iik),p_tmp2,mattype_surf)
!~       call petsc_get_densemat(p_h00k_r(iik),p_tmp3,mattype_surf)

    call surf_gf3_petsc(z - vr + cmplx(0d0, eta_elec), iik, "r")
!~     call MatNorm(p_grrr, NORM_FROBENIUS, normgrr, ierr)
    ! pull out -1 factor to avoid calling MatScale twice !!! (z*S_CR-H_CR)*g_RR*(z*S_RC-H_RC)
    zz = -(z - vr + cmplx(0d0, eta_elec))

    call MatDuplicate(p_s01k_r(iik), MAT_COPY_VALUES, p_tmp1, ierr)
    call MatScale(p_tmp1, zz, ierr)
    call MatAXPY(p_tmp1, p_one, p_h01k_r(iik), DIFFERENT_NONZERO_PATTERN, ierr)

    call MatDuplicate(p_s10k_r(iik), MAT_COPY_VALUES, p_tmp2, ierr)
    call MatScale(p_tmp2, zz, ierr)
    call MatAXPY(p_tmp2, p_one, p_h10k_r(iik), DIFFERENT_NONZERO_PATTERN, ierr)

    call MatConvert(p_tmp1, mattype_surf, MAT_INPLACE_MATRIX, p_tmp1, ierr)
    call MatConvert(p_tmp2, mattype_surf, MAT_INPLACE_MATRIX, p_tmp2, ierr)
    call MatMatMult(p_tmp1, p_grrr, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, p_tmp3, ierr)
    call MatDestroy(p_tmp1, ierr)
    call MatMatMult(p_tmp3, p_tmp2, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, p_tmp1, ierr) ! there is no ELEMENTAL*ELEMENTAL=DENSE so we need MatConvert
    call MatConvert(p_tmp1, mattype_dense, MAT_REUSE_MATRIX, p_sigmarr, ierr)

! gammaR=i*(sigmaRr-sigmaRa)
!~     call MatNorm(p_sigmarr, NORM_FROBENIUS, normr, ierr)
    call MatTranspose(p_sigmarr, MAT_REUSE_MATRIX, p_gammar, ierr)
    call MatConjugate(p_gammar, ierr)
    call MatAYPX(p_gammar, p_minus1, p_sigmarr, SAME_NONZERO_PATTERN, ierr) !p_gammar=(sigmaRr-sigmaRa)
    call MatScale(p_gammar, p_zione, ierr)

    call MatDestroy(p_tmp1, ierr)
    call MatDestroy(p_tmp2, ierr)
    call MatDestroy(p_tmp3, ierr)
!~     call petsc_mat_getvalue(0, 0, pl11, p_sigmalr, 1, PETSC_COMM_WORLD, lcast_in=.true.)
!~     call petsc_mat_getvalue(0, 0, pr11, p_sigmarr, 1, PETSC_COMM_WORLD, lcast_in=.true.)
!~     call petsc_mat_getvalue(nmu_l - 1 - 9 + 1, nmu_l - 1 - 9 + 1, plnn, p_sigmalr, 1, PETSC_COMM_WORLD, lcast_in=.true.)
!~     call petsc_mat_getvalue(nmu_r - 1 - 9 + 1, nmu_r - 1 - 9 + 1, prnn, p_sigmarr, 1, PETSC_COMM_WORLD, lcast_in=.true.)
!~     if (inode.eq.0) then
!~       write(0,fmt='(A,14e15.6)') "check sigma ",zz, normgll, normgrr, norml, normr, pl11, pr11, plnn, prnn
!~     end if

  end function eval_sigma

!    subroutine adaptive_int(f,a,b,zint1,zint2,nsub,eps,abserr)
!
!      use kinds
!      use mklfi_include
!
!      implicit none
!
!      integer, external :: f
!      integer :: nsub,i1,i2,norder
!      real(dp) :: abserr,a,b,eps
!      complex(dp), pointer :: zint1(:,:),zint2(:,:)
!
!      integer :: ierr,isub,n,nlow,iunit,i
!
!      complex(dp), allocatable, target :: tmp1_int1(:,:),tmp1_int2(:,:)
!      complex(dp), allocatable, target :: tmp2_int1(:,:),tmp2_int2(:,:)
!
!      complex(dp), pointer :: tmp1_int1_ptr(:,:),tmp1_int2_ptr(:,:)
!      complex(dp), pointer :: tmp2_int1_ptr(:,:),tmp2_int2_ptr(:,:)
!
!      real(dp), allocatable :: nodes(:),w1(:),w2(:),subint(:,:),suberr(:),&
!                               subtmp1(:,:),subtmp2(:,:),work(:)
!      real(dp) :: m,subcurrent(2,2),globalnorm,errorg,error2,error_1,error_2
!
!      int_counter=0
!      call get_quadrule(integrator,nodes,w1,w2,nint_order,nlow)
!
!
!
!      n=size(zint1,1)
!
!      allocate(tmp1_int1(n,n),tmp1_int2(n,n),work(1),stat=ierr)
!      if (ierr.ne.0) then
!        write(0,*) "allocation error ",ierr
!      end if
!
!      tmp1_int1_ptr=>tmp1_int1
!      tmp1_int2_ptr=>tmp1_int2
!      write(0,*) ldouble_contour,lnoneq_int
!      if (ldouble_contour.and.lnoneq_int) then
!        allocate(tmp2_int1(n,n),tmp2_int2(n,n),stat=ierr)
!        if (ierr.ne.0) then
!          write(0,*) "allocation error " ,ierr
!        end if
!        tmp2_int1_ptr=>tmp2_int1
!        tmp2_int2_ptr=>tmp2_int2
!        zint1=0d0
!        zint2=0d0
!      else
!        tmp2_int1_ptr=>null()
!        tmp2_int2_ptr=>null()
!        zint1=0d0
!      end if
!
!      allocate(subint(2,nsub),suberr(nsub))
!      isub=1
!      subint(1,1)=a
!      subint(2,1)=b
!      suberr=0d0
!
!      i=0
!      call evalqr(f,subint(1,1),subint(2,1),nodes,w1,w2,tmp1_int1_ptr,tmp2_int1_ptr,nint_order,nlow,error_1,error_2)
!      errorg=0d0
!      suberr(1)=max(error_1,error_2) !/globalnorm
!!~       error2=suberr(1)
!      subcurrent(1:2,1)=subint(1:2,1)
!      zint1=tmp1_int1_ptr
!      if (ldouble_contour.and.lnoneq_int) zint2=tmp2_int1_ptr
!
!      globalnorm=zlange('F',n,n,zint1,n,work)
!      if (ldouble_contour.and.lnoneq_int) globalnorm=(globalnorm+zlange('F',n,n,zint2,n,work))*0.5d0
!      abserr=suberr(1)
!!~       suberr(1)=0d0
!
!
!
!
!      do
!
!        abserr=maxval(suberr)/globalnorm
!        write(0,fmt='(i8,6e24.12,e10.3,i8)') isub,suberr(1:3)/globalnorm,subint(1:2,1),globalnorm,eps,int_counter
!        if ((abserr.gt.eps).and.(isub+3.le.nsub)) then
!!~           error2=error2-suberr(1)
!          suberr(1)=0d0
!          isub=isub+1
!          m=subint(2,1)-subint(1,1)
!          m=m*0.5d0
!          i1=isub+1
!          i2=isub+2
!          subint(1,i1)=subint(1,1)
!          subint(2,i1)=subint(1,1)+m
!          subint(1,i2)=subint(1,1)+m
!          subint(2,i2)=subint(2,1)
!!~           write(0,fmt='(A,3e24.12)') "new devision ",subint(1:2,i1),subint(2,i1)-subint(1,i1)
!!~           write(0,fmt='(A,3e24.12)') "new devision ",subint(1:2,i2),subint(2,i2)-subint(1,i2)
!          isub=isub+2
!        else
!!~           write(0,*) "abserr",abserr,isub
!          exit
!        end if
!
!
!        if (all(subcurrent(1:2,1)-subint(1:2,1).eq.0)) then
!          zint1=zint1-tmp1_int1_ptr
!          if (ldouble_contour.and.lnoneq_int) zint2=zint2-tmp2_int1_ptr
!        else if (all(subcurrent(1:2,2)-subint(1:2,1).eq.0)) then
!          zint1=zint1-tmp1_int2_ptr
!          if (ldouble_contour.and.lnoneq_int) zint2=zint2-tmp2_int2_ptr
!        else
!!~           write(0,*) "recalc",subint(1,1),subint(2,1)
!          call evalqr(f,subint(1,1),subint(2,1),nodes,w1,w2,tmp1_int1_ptr,tmp2_int1_ptr,nint_order,nlow,error_1,error_2)
!          zint1=zint1-tmp1_int1_ptr
!          if (ldouble_contour.and.lnoneq_int) zint2=zint2-tmp2_int1_ptr
!        end if
!
!        call evalqr(f,subint(1,i1),subint(2,i1),nodes,w1,w2,tmp1_int1_ptr,tmp2_int1_ptr,nint_order,nlow,error_1,error_2)
!        suberr(i1)=max(error_1,error_2)
!        zint1=zint1+tmp1_int1_ptr
!        if (ldouble_contour.and.lnoneq_int) zint2=zint2+tmp2_int1_ptr
!        call evalqr(f,subint(1,i2),subint(2,i2),nodes,w1,w2,tmp1_int2_ptr,tmp2_int2_ptr,nint_order,nlow,error_1,error_2)
!        suberr(i2)=max(error_1,error_2)
!        zint1=zint1+tmp1_int2_ptr
!        if (ldouble_contour.and.lnoneq_int) zint2=zint2+tmp2_int2_ptr
!
!        globalnorm=zlange('F',n,n,zint1,n,work)
!        if (ldouble_contour.and.lnoneq_int) globalnorm=(globalnorm+zlange('F',n,n,zint2,n,work))*0.5d0
!!~         globalnorm=dlange('F',n,n,aimag(zint),n,work)
!        suberr(i1)=suberr(i1)!/globalnorm
!        suberr(i2)=suberr(i2)!/globalnorm
!!~         error2=error2+suberr(i1)+suberr(i2)
!!~         write(6,*) "sub ",subint(1,isub+1),subint(2,isub+1),suberr(isub+1)
!
!!~         write(6,*) "sub ",subint(1,isub+2),subint(2,isub+2),suberr(isub+2)
!
!!~         write(0,fmt='(i8,4e24.12)') i1,subint(1:2,i1),suberr(i1),zlange('F',n,n,tmp1_int1,n,work)
!!~         write(0,fmt='(i8,4e24.12)') i2,subint(1:2,i2),suberr(i2),zlange('F',n,n,tmp1_int2,n,work)
!
!        subcurrent(1:2,1)=subint(1:2,i1)
!        subcurrent(1:2,2)=subint(1:2,i2)
!
!        call hpsort_aa (suberr,subint,.true.)
!
!
!
!
!      end do
!
!
!      nullify(tmp1_int1_ptr)
!      nullify(tmp1_int2_ptr)
!      nullify(tmp2_int1_ptr)
!      nullify(tmp2_int2_ptr)
!
!    end subroutine adaptive_int

  subroutine adaptive_int3(f, a, b, p_zint1, p_zint2, nsub, eps, abserr, subfile)
#include <petsc/finclude/petsc.h>
    use petscmat
    use petsc_mod
    use petsc_wrapper
    use kinds

    implicit none

    integer, external :: f
    integer :: nsub, i1, i2, i3, norder
    real(dp) :: abserr, a, b, eps, abserr_tot
    character(*), optional :: subfile
    logical :: l_subfile
    Mat :: p_zint1, p_zint2

    integer :: ierr, isub, n, nlow, iunit, i, ii(1)

    Mat :: p_tmp1_int1, p_tmp1_int2, p_tmp1_int3
    Mat :: p_tmp2_int1, p_tmp2_int2, p_tmp2_int3

    real(dp), allocatable :: nodes(:), w1(:), w2(:), subint(:, :), suberr(:), suberr2(:), &
                             subdiv(:, :)
    real(dp) :: m, subcurrent(2, 3), globalnorm, errorg, error2, error_1, error_2, norm, sub(2)

    int_counter = 0

    call get_quadrule(integrator, nodes, w1, w2, nint_order, nlow)

!~       n=size(zint1,1)

! setup necessary temporary matrices
    call MatDuplicate(p_zint1, MAT_SHARE_NONZERO_PATTERN, p_tmp1_int1, ierr)
    call MatDuplicate(p_zint1, MAT_SHARE_NONZERO_PATTERN, p_tmp1_int2, ierr)
    call MatDuplicate(p_zint1, MAT_SHARE_NONZERO_PATTERN, p_tmp1_int3, ierr)

! initialize matrices
    call MatZeroEntries(p_zint1, ierr)
    if (ldouble_contour .and. lnoneq_int) call MatZeroEntries(p_zint2, ierr)

    if (ldouble_contour .and. lnoneq_int) then
      call MatDuplicate(p_zint1, MAT_SHARE_NONZERO_PATTERN, p_tmp2_int1, ierr)
      call MatDuplicate(p_zint1, MAT_SHARE_NONZERO_PATTERN, p_tmp2_int2, ierr)
      call MatDuplicate(p_zint1, MAT_SHARE_NONZERO_PATTERN, p_tmp2_int3, ierr)
    end if

    allocate (subint(2, nsub), suberr(nsub), suberr2(nsub), subdiv(2, nsub))
    subint = 0d0
    suberr = 0d0
    suberr2 = 0d0
    subdiv = 0d0
    l_subfile = .false.
    if (present(subfile)) inquire (file=trim(subfile), exist=l_subfile)
    if (l_output_progress) then
      write (pstr_out, fmt='(A,l)') "read subdivision from file "//trim(subfile)//" ", l_subfile; call petsc_print_master()
    end if

    if (l_subfile) then
      suberr = 0d0
      if (l_ionode) then
        open (newunit=iunit, file=trim(subfile), action="read", status="old")
        read (iunit, *) isub
        do i1 = 1, isub
          read (iunit, *) i3, subint(1:2, i1), suberr(i1), subdiv(1:2, i1)
        end do
        close (iunit)
      end if
      call MPI_bcast(isub, 1, MPI_INTEGER, 0, PETSC_COMM_WORLD, ierr)
      call MPI_bcast(subint(1:2, 1:isub), 2*isub, MPI_DOUBLE_PRECISION, 0, PETSC_COMM_WORLD, ierr)
      call MPI_bcast(suberr(1:isub), isub, MPI_DOUBLE_PRECISION, 0, PETSC_COMM_WORLD, ierr)
      call MPI_bcast(subdiv(1:2, 1:isub), 2*isub, MPI_DOUBLE_PRECISION, 0, PETSC_COMM_WORLD, ierr)
      do i2 = isub, 1, -1
!~           i1=isub-i2+1
        if (l_output_progress) then
          write (pstr_out, fmt='(A,i8,3e24.12)') "doing sub ", i2, subint(1:2, i2), suberr(i2); call petsc_print_master(.false.)
        end if
        call evalqr(f, subint(1, i2), subint(2, i2), nodes, w1, w2, p_tmp1_int1, p_tmp2_int1, nint_order, nlow, error_1, error_2, subdiv(1, i2), subdiv(2, i2))
        suberr(i2) = max(error_1, error_2)
        if (l_output_progress) then
          write (pstr_out, fmt='(2X,e24.12)') suberr(i2); call petsc_print_master()
        end if
        subcurrent(1:2, 1) = subint(1:2, i2)
        call MatAXPY(p_zint1, p_one, p_tmp1_int1, SAME_NONZERO_PATTERN, ierr)
        if (ldouble_contour .and. lnoneq_int) call MatAXPY(p_zint2, p_one, p_tmp2_int1, SAME_NONZERO_PATTERN, ierr)       !     p_zint2=p_tmp2_int1
      end do
      call MatNorm(p_zint1, NORM_FROBENIUS, norm, ierr)
      globalnorm = norm
      if (ldouble_contour .and. lnoneq_int) then
        call MatNorm(p_zint2, NORM_FROBENIUS, norm, ierr)
        globalnorm = (globalnorm + norm)*0.5d0
      end if
      subcurrent(1:2, 1) = subint(1:2, isub)
      subcurrent(1:2, 2) = 0d0
      subcurrent(1:2, 3) = 0d0
      suberr2 = suberr
      call hpsort_aa(suberr, subint, .true.)
      call hpsort_aa(suberr2, subdiv, .true.)
      isub = isub + 3
    else
      isub = 1
      subint(1, isub) = a
      subint(2, isub) = b
      suberr = 0d0

      i = 0
      call evalqr(f, subint(1, 1), subint(2, 1), nodes, w1, w2, p_tmp1_int1, p_tmp2_int1, nint_order, nlow, error_1, error_2, subdiv(1, 1), subdiv(2, 1))
      errorg = 0d0
      suberr(1) = max(error_1, error_2)
      subcurrent(1:2, 1) = subint(1:2, 1)
      call MatAXPY(p_zint1, p_one, p_tmp1_int1, SAME_NONZERO_PATTERN, ierr)
      if (ldouble_contour .and. lnoneq_int) call MatAXPY(p_zint2, p_one, p_tmp2_int1, SAME_NONZERO_PATTERN, ierr)       !     p_zint2=p_tmp2_int1

      call MatNorm(p_zint1, NORM_FROBENIUS, norm, ierr)
      globalnorm = norm
      if (ldouble_contour .and. lnoneq_int) then
        call MatNorm(p_zint2, NORM_FROBENIUS, norm, ierr)
        globalnorm = (globalnorm + norm)*0.5d0
      end if
      abserr = suberr(1)
    end if
!~       suberr(1)=0d0

    if (l_output_progress) then
      call PetscPrintf(PETSC_COMM_WORLD, NEW_LINE('A'), ierr)
    end if
    do

      abserr = maxval(suberr)/globalnorm
      if (l_output_progress) then
        write (pstr_out, fmt='(i4,7e16.8,e10.3,i8)') &
          isub, suberr(1:3)/globalnorm, subint(1:2, 1), globalnorm, abserr, eps, int_counter; call petsc_print_master()
      end if
      if ((abserr .gt. eps) .and. (isub*3 + 3 .le. nsub)) then
!~           error2=error2-suberr(1)
        suberr(1) = 0d0
        m = subint(2, 1) - subint(1, 1)
        m = m/3d0
        i1 = isub + 1
        i2 = isub + 2
        i3 = isub + 3
        subint(1, i1) = subint(1, 1)
        subint(2, i1) = subdiv(1, 1) !subint(1,1)+1d0*m

        subint(1, i2) = subdiv(1, 1) !subint(1,1)+1d0*m
        subint(2, i2) = subdiv(2, 1) !subint(1,1)+2d0*m

        subint(1, i3) = subdiv(2, 1) !subint(1,1)+2d0*m
        subint(2, i3) = subint(2, 1)
        suberr(1) = 0d0
        isub = isub + 3
!~           write(0,fmt='(A,3e24.12)') "new devision ",subint(1:2,i1),subint(2,i1)-subint(1,i1)
!~           write(0,fmt='(A,3e24.12)') "new devision ",subint(1:2,i2),subint(2,i2)-subint(1,i2)
!~           write(0,fmt='(A,3e24.12)') "new devision ",subint(1:2,i3),subint(2,i3)-subint(1,i3)
      else
        exit
      end if

      if (all(subcurrent(1:2, 1) - subint(1:2, 1) .eq. 0)) then
        call MatAXPY(p_zint1, p_minus1, p_tmp1_int1, SAME_NONZERO_PATTERN, ierr)
        if (ldouble_contour .and. lnoneq_int) call MatAXPY(p_zint2, p_minus1, p_tmp2_int1, SAME_NONZERO_PATTERN, ierr)
      else if (all(subcurrent(1:2, 2) - subint(1:2, 1) .eq. 0)) then
        call MatAXPY(p_zint1, p_minus1, p_tmp1_int2, SAME_NONZERO_PATTERN, ierr)
        if (ldouble_contour .and. lnoneq_int) call MatAXPY(p_zint2, p_minus1, p_tmp2_int2, SAME_NONZERO_PATTERN, ierr)
      else if (all(subcurrent(1:2, 3) - subint(1:2, 1) .eq. 0)) then
        call MatAXPY(p_zint1, p_minus1, p_tmp1_int3, SAME_NONZERO_PATTERN, ierr)
        if (ldouble_contour .and. lnoneq_int) call MatAXPY(p_zint2, p_minus1, p_tmp2_int3, SAME_NONZERO_PATTERN, ierr)
      else
        call evalqr(f, subint(1, 1), subint(2, 1), nodes, w1, w2, p_tmp1_int1, p_tmp2_int1, nint_order, nlow, error_1, error_2, subdiv(1, 1), subdiv(2, 1))
        call MatAXPY(p_zint1, p_minus1, p_tmp1_int1, SAME_NONZERO_PATTERN, ierr)
        if (ldouble_contour .and. lnoneq_int) call MatAXPY(p_zint2, p_minus1, p_tmp2_int1, SAME_NONZERO_PATTERN, ierr)
      end if

! a..a+1/3
      call evalqr(f, subint(1, i1), subint(2, i1), nodes, w1, w2, p_tmp1_int1, p_tmp2_int1, nint_order, nlow, error_1, error_2, subdiv(1, i1), subdiv(2, i1))
      suberr(i1) = max(error_1, error_2)
      call MatAXPY(p_zint1, p_one, p_tmp1_int1, SAME_NONZERO_PATTERN, ierr)
      if (ldouble_contour .and. lnoneq_int) call MatAXPY(p_zint2, p_one, p_tmp2_int1, SAME_NONZERO_PATTERN, ierr)

! a+1/3*d..a+2/3*d
      call evalqr(f, subint(1, i2), subint(2, i2), nodes, w1, w2, p_tmp1_int2, p_tmp2_int2, nint_order, nlow, error_1, error_2, subdiv(1, i2), subdiv(2, i2))
      suberr(i2) = max(error_1, error_2)
      call MatAXPY(p_zint1, p_one, p_tmp1_int2, SAME_NONZERO_PATTERN, ierr)
      if (ldouble_contour .and. lnoneq_int) call MatAXPY(p_zint2, p_one, p_tmp2_int2, SAME_NONZERO_PATTERN, ierr)

! a+2/3d..a+3/3*d=b
      call evalqr(f, subint(1, i3), subint(2, i3), nodes, w1, w2, p_tmp1_int3, p_tmp2_int3, nint_order, nlow, error_1, error_2, subdiv(1, i3), subdiv(2, i3))
      suberr(i3) = max(error_1, error_2)
      call MatAXPY(p_zint1, p_one, p_tmp1_int3, SAME_NONZERO_PATTERN, ierr)
      if (ldouble_contour .and. lnoneq_int) call MatAXPY(p_zint2, p_one, p_tmp2_int3, SAME_NONZERO_PATTERN, ierr)

      call MatNorm(p_zint1, NORM_FROBENIUS, norm, ierr)
      globalnorm = norm
      if (ldouble_contour .and. lnoneq_int) then
        call MatNorm(p_zint2, NORM_FROBENIUS, norm, ierr)
        globalnorm = (globalnorm + norm)*0.5d0
      end if

      if (l_output_progress) call PetscPrintf(PETSC_COMM_WORLD, NEW_LINE('A'), ierr)
      if (l_output_progress) then
        write (pstr_out, fmt='(i8,6e24.12)') i1, subint(1:2, i1),&
        &subint(2, i1) - subint(1, i1), suberr(i1), subdiv(1:2, i1)
        call petsc_print_master()
      end if
      if (l_output_progress) then
        write (pstr_out, fmt='(i8,6e24.12)') i2, subint(1:2, i2),&
        &subint(2, i2) - subint(1, i2), suberr(i2), subdiv(1:2, i2)
        call petsc_print_master()
      end if
      if (l_output_progress) then
        write (pstr_out, fmt='(i8,6e24.12)') i3, subint(1:2, i3),&
        &subint(2, i3) - subint(1, i3), suberr(i3), subdiv(1:2, i3)
        call petsc_print_master()
      end if

      if (abs((subint(2, i1) - subint(1, i1))) .le. 1d-3) then
        suberr(i1) = 0d0
      end if
      if (abs((subint(2, i2) - subint(1, i2))) .le. 1d-3) then
        suberr(i2) = 0d0
      end if
      if (abs((subint(2, i3) - subint(1, i3))) .le. 1d-3) then
        suberr(i3) = 0d0
      end if

      subint(1:2, 1) = -1d0
      suberr(1) = -1d0
      subcurrent(1:2, 1) = subint(1:2, i1)
      subcurrent(1:2, 2) = subint(1:2, i2)
      subcurrent(1:2, 3) = subint(1:2, i3)
      suberr2 = suberr
      call hpsort_aa(suberr, subint, .true.)
      call hpsort_aa(suberr2, subdiv, .true.)

    end do

    call MatDestroy(p_tmp1_int1, ierr)
    call MatDestroy(p_tmp1_int2, ierr)
    call MatDestroy(p_tmp1_int3, ierr)
    if (ldouble_contour .and. lnoneq_int) then
      call MatDestroy(p_tmp2_int1, ierr)
      call MatDestroy(p_tmp2_int2, ierr)
      call MatDestroy(p_tmp2_int3, ierr)
    end if

    if (present(subfile)) then
      if (l_ionode) then
!~           suberr=-suberr

!~           call hpsort_aa (suberr(1:i3),subint_s(1:2,1:i3),.true.)
!~           call hpsort_aa (suberr2(1:i3),subdiv(1:2,1:i3),.true.)
!~           suberr=-suberr

        ii = minloc(abs(suberr))
        i3 = ii(1) - 1
        if (i3 .eq. 1) return
        do i1 = 1, i3 - 1
          do i2 = i1 + 1, i3
            if (subint(1, i1) .ge. subint(1, i2)) then
              sub = subint(1:2, i1)
              subint(1:2, i1) = subint(1:2, i2)
              subint(1:2, i2) = sub
              sub = subdiv(1:2, i1)
              subdiv(1:2, i1) = subdiv(1:2, i2)
              subdiv(1:2, i2) = sub
              error2 = suberr(i1)
              suberr(i1) = suberr(i2)
              suberr(i2) = error2
            end if
          end do
        end do

!~           abserr=suberr(i3)
!~           i2=i3
!~           do i=i3-1,1,-1
!~             abserr=abserr+suberr(i)
!~             write(pstr_out,fmt='(2i4,2e24.12)') i,i2,abserr,eps*1d-1 ; call petsc_print_master()
!~             if (abserr.gt.eps*1d-1) then
!~               suberr(i2)=sum(suberr(i+1:i2))
!~               subint(1,i2)=subint(2,i)
!~               suberr(i+1:i2-1)=0d0
!~               subint(1:2,i+1:i2-1)=0d0
!~               abserr=suberr(i)
!~               i2=i
!~             end if
!~           end do
        suberr2 = suberr
        call hpsort_aa(suberr(1:i3), subint(1:2, 1:i3), .true.)
        call hpsort_aa(suberr2(1:i3), subdiv(1:2, 1:i3), .true.)
        ii = minloc(abs(suberr))
        i3 = ii(1) - 1
        open (newunit=iunit, file=trim(subfile), action="write", status="replace")
        write (iunit, fmt='(i8)') i3
        do i1 = i3, 1, -1
          write (iunit, fmt='(i8,5es45.24e5)') i1, subint(1:2, i1), suberr(i1), subdiv(1:2, i1)
        end do
        close (iunit)
      end if
    end if

!~       nullify(p_tmp1_int1_ptr)
!~       nullify(p_tmp1_int2_ptr)
!~       nullify(p_tmp1_int3_ptr)
!~       nullify(p_tmp2_int1_ptr)
!~       nullify(p_tmp2_int2_ptr)
!~       nullify(p_tmp2_int3_ptr)

  end subroutine adaptive_int3

end module integrator_mod
